### CODEX data processing and visualization helper module
import numpy as np
import os
import xarray as xr
import pandas as pd
from glob import glob
import re

def get_cells_mask_and_coords_old(adata, condition_vecs_dict, paths_dict, working_samples_dict, sample_col_name="sample"):

    n_cells = adata.n_obs
    spatial_coords = np.full((n_cells, 2), np.nan)
    mask_origin = np.full(n_cells, np.nan)  # numeric, not dtype=object - object-dtype floats break h5ad's string inference on write
    
    obs_names = adata.obs_names.to_numpy()
    sample_col = adata.obs[sample_col_name].to_numpy()
    
    matched = 0
    unmatched_samples = []
    
    for run_id, conditions in condition_vecs_dict.items():
        run_key = f"run{run_id}"
        pattern = os.path.join(paths_dict["zarr_dir"], f"{paths_dict['zarr_prefix']}{run_key}*.zarr")
        # sorted() so this matches R's list.files() (alphabetical by default) -
        # same fix applied throughout the rest of the pipeline's Python scripts
        zarr_files = sorted(glob(pattern))
    
        for i in working_samples_dict[run_id]:  # 1-based, matches R's `i`
            sample_label = f"{run_key}_{conditions[i - 1]}_{i}"
            zarr_path = zarr_files[i - 1]  # i-th alphabetically-sorted file for this run (1-based)
    
            dat = xr.open_zarr(zarr_path)
            obs_df = pd.DataFrame(dat._obs)
            coords = obs_df.iloc[:, -2:].to_numpy()   # last 2 columns = spatial x, y
            mask_ids = obs_df.iloc[:, 0].to_numpy()   # first column = which seg. mask this cell came from
    
            sample_mask = sample_col == sample_label
            sample_cell_names = obs_names[sample_mask]
    
            if len(sample_cell_names) == 0:
                unmatched_samples.append(sample_label)
                continue
    
            for name in sample_cell_names:
                orig_idx = int(name.split("_")[0])
                row = np.where(obs_names == name)[0][0]
                spatial_coords[row] = coords[orig_idx]
                mask_origin[row] = mask_ids[orig_idx]
                matched += 1
    
    adata.obsm["spatial"] = spatial_coords
    adata.obs["mask"] = mask_origin

    print(f"Matched {matched} / {n_cells} cells to spatial coordinates.")
    if unmatched_samples:
        print("WARNING: no h5ad cells found for these samples (check config.yaml runs/conditions):")
        for s in unmatched_samples:
            print(" -", s)
    n_missing = np.isnan(spatial_coords).any(axis=1).sum()
    if n_missing > 0:
        print(f"WARNING: {n_missing} cells still have no spatial coordinates - "
              f"double check config.yaml's runs/conditions/working_samples match what P1.r used.")
    else:
        print("All cells matched successfully.")
    
    return adata

def get_cells_mask_and_coords(
    adata,
    condition_vecs_dict,
    paths_dict,
    working_samples_dict,
    sample_col_name="sample",
):

  n_cells = adata.n_obs
  spatial_coords = np.full((n_cells, 2), np.nan)
  mask_origin = np.full(n_cells, np.nan)

  obs_names = adata.obs_names.to_numpy()
  sample_col = adata.obs[sample_col_name].to_numpy()

  matched = 0
  unmatched_samples = []

  for run_id, conditions in condition_vecs_dict.items():
    run_key = f"run{run_id}"
    pattern = os.path.join(
        paths_dict["zarr_dir"], f"{paths_dict['zarr_prefix']}{run_key}*.zarr"
    )
    zarr_files = sorted(glob(pattern))

    for i in working_samples_dict[run_id]:  # 1-based index
      sample_label = f"{run_key}_{conditions[i - 1]}_{i}"
      zarr_path = zarr_files[i - 1]

      dat = xr.open_zarr(zarr_path)
      obs_df = pd.DataFrame(dat._obs)
      coords = obs_df.iloc[:, -2:].to_numpy()  # last 2 columns = spatial x, y
      mask_ids = (
          obs_df.iloc[:, 0].to_numpy()
      )  # first column = segmentation mask ID

      # 1. Get array row indices directly from the boolean mask in O(N) once per sample
      sample_mask = sample_col == sample_label
      sample_indices = np.where(sample_mask)[0]

      if len(sample_indices) == 0:
        unmatched_samples.append(sample_label)
        continue

      # 2. Extract original index integers in a single list comprehension
      orig_indices = [
          int(name.split("_", 1)[0]) for name in obs_names[sample_indices]
      ]

      # 3. Vectorized assignment directly into target array slices
      spatial_coords[sample_indices] = coords[orig_indices]
      mask_origin[sample_indices] = mask_ids[orig_indices]
      matched += len(sample_indices)

  adata.obsm["spatial"] = spatial_coords
  adata.obs["mask"] = mask_origin

  print(f"Matched {matched} / {n_cells} cells to spatial coordinates.")
  if unmatched_samples:
    print(
        "WARNING: no h5ad cells found for these samples (check config.yaml"
        " runs/conditions):"
    )
    for s in unmatched_samples:
      print(" -", s)

  n_missing = np.isnan(spatial_coords[:, 0]).sum()
  if n_missing > 0:
    print(
        f"WARNING: {n_missing} cells still have no spatial coordinates - "
        "double check config.yaml's runs/conditions/working_samples match what"
        " P1.r used."
    )
  else:
    print("All cells matched successfully.")

  return adata


def to_r_safe_name(raw_name):
    """Mirrors R's make.names() default transformation, e.g.
    '127 CD61' -> 'X127.CD61', 'Osteopontin 166' -> 'Osteopontin.166'.
    Double check against anndata.var_names if a name doesn't match -
    make.names() has a few more edge cases this doesn't cover."""
    safe = re.sub(r"[^0-9A-Za-z.]", ".", raw_name)
    if re.match(r"^[0-9]", safe):
        safe = "X" + safe
    return safe


def get_spatialproteomics_obj(adata, condition_vecs_dict, paths_dict, working_samples_dict, cell_colors, celltype_col, sample_col="sample"):
    # label -> color, from build_final_anndata.ipynb's uns, + gray for filtered-out cells
    clist = cell_colors + ["lightgray"]
    clist_series = pd.Series(
        clist, index=list(adata.obs[celltype_col].cat.categories) + ["Removed"]
    )
    obs_names = adata.obs_names.to_numpy()
    
    dat = {}          # sample_label -> spatialproteomics object with labels attached
    CTs_by_sample = {} # sample_label -> the cell-type dataframe (used by the ROI section below)
    
    for run_id, conditions in condition_vecs_dict.items():
        run_key = f"run{run_id}"
        pattern = os.path.join(paths_dict["zarr_dir"], f"{paths_dict['zarr_prefix']}{run_key}*.zarr")
        zarr_files = sorted(glob(pattern))  # sorted() to match R's list.files() ordering
    
        for i in working_samples_dict[run_id]:  # 1-based, matches R's `i`
            sample_label = f"{run_key}_{conditions[i - 1]}_{i}"
            zarr_path = zarr_files[i - 1]
    
            d = xr.open_zarr(zarr_path)
            n_cells_total = len(pd.DataFrame(d._obs))
    
            h5ad_names = obs_names[adata.obs[sample_col] == sample_label]
            if len(h5ad_names) == 0:
                print(f"WARNING: no h5ad cells found for {sample_label}, skipping")
                continue
    
            # cell types: full-size, "Removed"-filled, real values for surviving cells
            CTs = pd.DataFrame({"cell": np.arange(1, n_cells_total + 1), "label": "Removed"})
    
            for name in h5ad_names:
                orig_idx = int(name.split("_")[0])
                CTs.loc[orig_idx, "label"] = adata.obs.loc[name, celltype_col]
    
            CTs["label"] = CTs["label"].astype(str)
            CTs_by_sample[sample_label] = CTs
    
            if "_la_properties" in d.data_vars:
                d = d.pp.drop_layers("_la_properties")
            d = d.la.add_labels_from_dataframe(CTs)
    
            # recolor to match build_final_anndata.ipynb's palette
            lab_props = pd.DataFrame(d._la_properties, columns=["color", "label"])
            for lp_i in range(len(lab_props)):
                lbl = lab_props["label"][lp_i]
                mask = d._la_properties[:, 1] == lbl
                d._la_properties[mask, 0] = clist_series.get(lbl, "lightgray")
    
            dat[sample_label] = d
            print(f"{sample_label}: {n_cells_total} total cells, {len(h5ad_names)} labeled, "
                  f"{n_cells_total - len(h5ad_names)} marked Removed")

    return dat, CTs_by_sample
    print("\nSamples ready for plotting:", list(dat.keys()))


def filter_high(ds, q=0.99):
    """Clip each channel's values above quantile q to 0."""
    cutoff = ds.quantile(q, dim=["x", "y"])
    return ds.where(ds <= cutoff, 0)


def filter_high_multi(ds, quantiles=[0.99, 0.95]):
    """Same as filter_high, but with one quantile per channel."""
    img_var = "_image" if "_image" in ds.data_vars else "image"
    for i, q in enumerate(quantiles):
        channel_slice = ds[img_var].isel(channels=i)
        cutoff = channel_slice.quantile(q).compute()
        ds[img_var].loc[dict(channels=ds.channels.values[i])] = channel_slice.where(channel_slice <= cutoff, 0)
    return ds


def filter_bone(ds, filts=[0, 1], cutoffs_q=[0.5, 0.9]):
    """Two-marker-gated filter (e.g. keep only pixels that are simultaneously
    above cutoffs_q[0] for channel filts[0] and below cutoffs_q[1] for
    channel filts[1]) - example use case was a collagen/DAPI 'bone' mask."""
    img_var = "_image" if "_image" in ds.data_vars else "image"
    f1_slice = ds[img_var].isel(channels=filts[0])
    f2_slice = ds[img_var].isel(channels=filts[1])
    f1_cut = f1_slice.quantile(cutoffs_q[0]).compute()
    f2_cut = f2_slice.quantile(cutoffs_q[1]).compute()
    for i, _ in enumerate(list(ds.channels.to_dataframe().index)):
        channel_slice = ds[img_var].isel(channels=i)
        ds[img_var].loc[dict(channels=ds.channels.values[i])] = channel_slice.where(f1_slice >= f1_cut, 0)
        ds[img_var].loc[dict(channels=ds.channels.values[i])] = channel_slice.where(f2_slice <= f2_cut, 0)
    return ds