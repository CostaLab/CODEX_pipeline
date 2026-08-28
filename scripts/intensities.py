"""
Step 1 of the CODEX pipeline.

Loads cropped multiplex images + one or two segmentation masks per region
(merging the two, if a second one is configured), quantifies marker
intensity per cell, and writes one zarr store per region.

All parameters come from config.yaml - see py_config.py. Set
paths.seg2_masks_subdir to null/blank for single-mask mode (no merge step).
"""
import os
from glob import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tifffile
import dask.array as da
import spatialproteomics

from py_config import load_config, resolve_paths, parse_args, get_marker_channels


def fix_for_zarr(obj):
    # 1. Fix Global Attributes
    obj.attrs = {str(k): (v.item() if hasattr(v, 'item') else v)
                 for k, v in obj.attrs.items()}
    # 2. Fix Variable Attributes AND Encoding (the "Dask/Tiff" fix)
    for var in obj.variables:
        obj[var].encoding = {}
        obj[var].attrs = {str(k): (v.item() if hasattr(v, 'item') else v)
                           for k, v in obj[var].attrs.items()}
    return obj


### different quantification methods of marker intensity per cell
def intensity_sum(regionmask, intensity):
    return np.sum(intensity[regionmask])


def intensity_median(regionmask, intensity):
    return np.median(intensity[regionmask])


def intensity_sum_arcsinh_scaled(regionmask, intensity):
    x = np.sum(intensity[regionmask])
    return np.arcsinh(1 + x / 5)


CUSTOM_INTENSITY_FUNCS = {
    "intensity_sum": intensity_sum,
    "intensity_median": intensity_median,
    "intensity_sum_arcsinh_scaled": intensity_sum_arcsinh_scaled,
}


def resolve_intensity_func(name):
    """'intensity_mean' (and other skimage/spatialproteomics builtins)
    are passed through as a string; anything in CUSTOM_INTENSITY_FUNCS
    is passed as the actual function object."""
    return CUSTOM_INTENSITY_FUNCS.get(name, name)


def main():
    args = parse_args("Step 1: build per-region zarr stores with quantified intensities")
    cfg = load_config(args.config)
    paths = resolve_paths(cfg)
    seg_cfg = cfg["segmentation"]

    channels = get_marker_channels(paths["marker_list_file"])

    seg_suffix = seg_cfg["seg_files_suffix"]
    seg2_suffix = seg_cfg["seg2_files_suffix"]
    img_suffix = seg_cfg["img_files_suffix"]
    suffix_replace = seg_cfg["files_suffix_replace"]
    segm_mks = seg_cfg["segmentation_markers"]
    seg_colors = seg_cfg["segmentation_colors"]
    merge_labels_order = seg_cfg["merge_labels_order"]
    intensity_func = resolve_intensity_func(seg_cfg["intensity_func"])
    merge_threshold = seg_cfg["merge_threshold"]
    plot_dpi = seg_cfg["plot_dpi"]

    run_ids = {str(r["id"]) for r in cfg["runs"]}

    # --- build {region_key: filepath} dicts for images and both masks ---
    data = {f.split('/')[-1].replace('.tif', ''): f for f in sorted(glob(os.path.join(paths["input_images"], '*.tif')))}
    for key in list(data.keys()):
        data[key.replace(img_suffix, suffix_replace)] = data.pop(key)

    use_second_mask = paths["seg2_masks"] is not None

    seg = {f.split('/')[-1].replace('.tif', ''): f for f in sorted(glob(os.path.join(paths["seg1_masks"], '*.tif')))}
    for key in list(seg.keys()):
        seg[key.replace(seg_suffix, suffix_replace)] = seg.pop(key)

    if use_second_mask:
        seg2 = {f.split('/')[-1].replace('.tif', ''): f for f in sorted(glob(os.path.join(paths["seg2_masks"], '*.tif')))}
        for key in list(seg2.keys()):
            seg2[key.replace(seg2_suffix, suffix_replace)] = seg2.pop(key)
    else:
        seg2 = {}
        print("No second segmentation mask configured (seg2_masks_subdir is blank) - "
              "running in single-mask mode, no merge_segmentation step.")

    # only process regions that belong to a run listed in config.yaml
    def run_of(region_key):
        # region_key looks like 'run21_...' after the suffix replace above
        for rid in run_ids:
            if region_key.startswith(f"{suffix_replace}{rid}"):
                return rid
        return None

    regions = [r for r in data.keys() if run_of(r) is not None]
    skipped = [r for r in data.keys() if run_of(r) is None]
    if skipped:
        print(f"Skipping {len(skipped)} region(s) not listed under `runs` in config.yaml: {skipped}")

    dat = {}

    for reg in regions:
        print(f"Processing {reg} ...")
        segm = tifffile.imread(seg[reg])

        with tifffile.TiffFile(data[reg], mode='r') as tif:
            mmap = tif.asarray(out='memmap')
            d_array = da.from_array(mmap, chunks=(1, mmap.shape[1], mmap.shape[2]))

        print(f"Dask array chunks: {d_array.chunks}")

        sdata = spatialproteomics.load_image_data(d_array, channel_coords=channels, segmentation=segm)

        if use_second_mask:
            segm2 = tifffile.imread(seg2[reg])
            # two segmentation masks
            sdata['_seg1'] = sdata._segmentation
            sdata = sdata.pp.drop_layers('_segmentation')
            sdata = sdata.pp.add_segmentation(segm2)
            sdata['_seg2'] = sdata._segmentation
            sdata = sdata.pp.drop_layers('_segmentation')

            ### Plot DAPI & second-marker segmentation masks side by side
            fig, ax = plt.subplots(1, 2, figsize=(10, 5))
            _ = (
                sdata.pp[segm_mks]
                .pl.colorize(seg_colors)
                .pl.show(render_segmentation=True, segmentation_kwargs={"layer_key": "_seg1"}, ax=ax[0])
            )
            _ = (
                sdata.pp[segm_mks]
                .pl.colorize(seg_colors)
                .pl.show(render_segmentation=True, segmentation_kwargs={"layer_key": "_seg2"}, ax=ax[1])
            )
            plot_path = os.path.join(paths["visual_output_dir"], f"segmentation_plot_pro_ind_{reg}.png")
            plt.savefig(plot_path, dpi=plot_dpi, bbox_inches='tight')
            plt.close()

            ### merge segmentation masks
            sdata = sdata.pp.merge_segmentation(
                layer_key=["_seg2", "_seg1"],  # big-cells mask goes first
                labels=merge_labels_order,
                threshold=merge_threshold,
            ).pp.add_segmentation("_merged_segmentation")
        else:
            ### single mask: nothing to merge - sdata._segmentation (set by
            ### load_image_data above) is already the mask to quantify against.
            fig, ax = plt.subplots(1, 1, figsize=(6, 5))
            _ = (
                sdata.pp[segm_mks]
                .pl.colorize(seg_colors)
                .pl.show(render_segmentation=True, ax=ax)
            )
            plot_path = os.path.join(paths["visual_output_dir"], f"segmentation_plot_pro_ind_{reg}.png")
            plt.savefig(plot_path, dpi=plot_dpi, bbox_inches='tight')
            plt.close()

        sdata = sdata.pp.add_quantification(func=intensity_func).pp.transform_expression_matrix(method="arcsinh")
        sdata = sdata.pp.add_observations('area')

        sdata = fix_for_zarr(sdata)
        dat[reg] = sdata.copy()

        ### Plot final segmentation (merged, or the single mask)
        plt.figure(figsize=(12, 12))
        if use_second_mask:
            # set_label_colors needs the per-cell origin-mask labels that
            # merge_segmentation's `labels=` param creates as a side effect -
            # there's no such label layer in single-mask mode.
            sdata = sdata.la.set_label_colors(segm_mks, seg_colors)
        _ = sdata.pp[segm_mks].pl.colorize(colors=seg_colors).pl.show(render_segmentation=True)
        plt.axis('off')
        plot_path = os.path.join(paths["visual_output_dir"], f"segmentation_plot_pro_{reg}.png")
        plt.savefig(plot_path, dpi=plot_dpi, bbox_inches='tight')
        plt.close()

        del d_array, segm, sdata

    for k, v in dat.items():
        out_path = os.path.join(paths["zarr_dir"], f"{paths['zarr_prefix']}{k}.zarr")
        v.to_zarr(out_path)
        print(f"Wrote {out_path}")

    print("intensities.py: DONE")


if __name__ == "__main__":
    main()
