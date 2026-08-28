"""
Step 2 of the CODEX pipeline.

For every region zarr produced by intensities.py:
  - exports per-cell marker intensities to CSV
  - plots a marker histogram grid
  - plots pairwise scatter QC plots against a configurable set of markers
    (useful to check markers that shouldn't co-occur in the same cell)
"""
import os
from glob import glob

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import math

from py_config import load_config, resolve_paths, parse_args, get_marker_channels, get_run_dicts

sc.set_figure_params(dpi=100, dpi_save=900)


def main():
    args = parse_args("Step 2: export intensities + QC plots")
    cfg = load_config(args.config)
    paths = resolve_paths(cfg)
    vis_cfg = cfg["visualization"]

    sc.set_figure_params(dpi=vis_cfg["scanpy_dpi"], dpi_save=vis_cfg["scanpy_dpi_save"])

    channels = get_marker_channels(paths["marker_list_file"])
    condition_vecs, _ = get_run_dicts(cfg)
    intensities_file_suffix = vis_cfg["intensities_file_suffix"]
    scatter_plots_mkrs = vis_cfg["scatter_markers"]

    ### new
    ## Grids for figures
    n_channels = len(channels)

    # 1. Define desired number of columns and calculate required rows
    ncols = 6
    nrows = math.ceil(n_channels / ncols)
    
    # 2. Scale figure size dynamically (width_per_col, height_per_row)
    cell_width = 2.6  # inches per subplot
    cell_height = 3.0
    figsize = (ncols * cell_width, nrows * cell_height)
    ###

    for run in condition_vecs.keys():
        pattern = os.path.join(paths["zarr_dir"], f"{paths['zarr_prefix']}run{run}*.zarr")
        dat = {f.split('/')[-1].split('.')[0]: xr.open_zarr(f) for f in sorted(glob(pattern))}
        if not dat:
            print(f"No zarr stores found for run {run} ({pattern}), skipping.")
            continue

        condition_vec = condition_vecs[run]

        # store intensity values per cell
        for i in dat.keys():
            out_csv = os.path.join(paths["quant_dir"], f"{i}{intensities_file_suffix}")
            pd.DataFrame(dat[i]._intensity.data, columns=channels).to_csv(out_csv)

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

            ### Histogram per protein
            df = pd.DataFrame(dat[i]._intensity.data, columns=channels)
            #fig, axes = plt.subplots(8, 6, figsize=(16, 24), sharex=True, sharey=True)
            fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
            ax = axes.flatten()
            for x in range(len(channels)):
                _ = df.iloc[:, x].hist(bins=100, ax=ax[x], legend=True)
            hist_path = os.path.join(paths["visual_output_dir"], f"allmks_hist_run{run}_{condition_vec[j]}_{j}.png")
            plt.savefig(hist_path)
            plt.close(fig)

            ### Scatter plots for markers that shouldn't co-occur in the same cell
            for mkr in scatter_plots_mkrs:
                grid_coords = [[i2, j2] for i2 in range(8) for j2 in range(6)]
                #fig, axes = plt.subplots(8, 6, figsize=(16, 24))
                fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
                for idx in range(len(anndata.var_names)):
                    r, c = grid_coords[idx]
                    axes[r][c].scatter(anndata.X[:, anndata.var_names.values == mkr], anndata.X[:, idx])
                    axes[r][c].set_title(anndata.var_names[idx])
                scatter_path = os.path.join(
                    paths["visual_output_dir"], f"run{run}_{condition_vec[j]}_{j}{mkr}_scatter.png"
                )
                fig.savefig(scatter_path)
                plt.close(fig)

            j += 1

    print("codex_vis.py: DONE")


if __name__ == "__main__":
    main()
