"""
Step 3 of the CODEX pipeline.

For every region zarr, subtracts each marker's per-cell median from itself
(clipped at 0) and plots the resulting "background-filtered" spatial
marker maps. Scanpy writes these into a `figures/` folder next to wherever
this script is run from.
"""
import os
from glob import glob

import pandas as pd
import xarray as xr
import scanpy as sc
import anndata as ad

from py_config import load_config, resolve_paths, parse_args, get_marker_channels, get_run_dicts


def main():
    args = parse_args("Step 3: background-filtered spatial marker plots")
    cfg = load_config(args.config)
    paths = resolve_paths(cfg)
    vis_cfg = cfg["visualization"]

    sc.set_figure_params(dpi=vis_cfg["scanpy_dpi"], dpi_save=vis_cfg["scanpy_dpi_save"])

    channels = get_marker_channels(paths["marker_list_file"])
    condition_vecs, _ = get_run_dicts(cfg)
    spot_size = vis_cfg["spatial_spot_size"]

    for run in condition_vecs.keys():
        pattern = os.path.join(paths["zarr_dir"], f"{paths['zarr_prefix']}run{run}*.zarr")
        dat = {f.split('/')[-1].split('.')[0]: xr.open_zarr(f) for f in sorted(glob(pattern))}
        if not dat:
            print(f"No zarr stores found for run {run} ({pattern}), skipping.")
            continue

        condition_vec = condition_vecs[run]

        j = 0
        for i in dat.keys():
            print(i)

            anndata = ad.AnnData(pd.DataFrame(dat[i]._intensity, columns=channels).iloc[:, 1:len(channels)])
            sp = pd.DataFrame(dat[i]._obs).iloc[:, -2:]  # last 2 columns
            sp.columns = ['0', '1']
            sp.index = anndata.obs.index
            anndata.obsm['spatial'] = sp.to_numpy()

            print(anndata.obs.shape)
            anndata.obs_names = [f"{n}_{condition_vec[j]}_{j}" for n in anndata.obs_names]

            ### Median filtering: subtract each marker's median, clip at 0
            thr = anndata.to_df().median()
            dat1 = anndata.to_df().copy()
            for c in dat1.columns:
                dat1[c] = dat1[c] - thr[c]
                dat1.loc[dat1[c] < 0, c] = 0
            anndata.X = dat1.to_numpy()

            sc.pl.spatial(
                anndata, color=anndata.var_names, spot_size=spot_size, use_raw=False,
                save=f"run{run}_{condition_vec[j]}_{j}_mks_bgFilt.png"
            )
            j += 1

    print("codex_vis_bgFilt.py: DONE")


if __name__ == "__main__":
    main()
