"""
Generates a tiny synthetic CODEX dataset (1 run, 2 samples, 3 markers,
40x40 px images with 4 toy "cells" each) that runs end-to-end through the
whole pipeline in a few seconds. Useful as:

  - a smoke test after editing the pipeline scripts
  - a quick way to sanity-check your config.yaml structure and both conda
    environments before submitting real SLURM jobs
  - a teaching example for how the input files need to be named/shaped

Usage:
    python generate_toy_data.py [output_dir]

    output_dir defaults to ./toy_data next to this script. This becomes
    the `paths.wdir` used by example/toy_config.yaml.
"""
import os
import sys

import numpy as np
import tifffile

RNG = np.random.default_rng(42)

# marker order matters: index 0 must be the nuclear marker used for
# segmentation (DAPI), matching segmentation.segmentation_markers[0]
# in the config.

#MARKERS = ["DAPI", "127 CD61", "010 CD14"]

# Generate 30 markers (Index 0 = nuclear segmentation marker)
BASE_MARKERS = ["DAPI", "127 CD61", "010 CD14", "CD3", "CD4", "CD8a", "CD20", "CD11b", "CD45", "FoxP3"]
EXTRA_MARKERS = [f"Marker_{i:02d}" for i in range(1, 21)]
MARKERS = BASE_MARKERS + EXTRA_MARKERS  # Total: 30 markers

RUN_ID = "99"
SAMPLES = ["sample1", "sample2"]
#IMG_H, IMG_W = 40, 40
IMG_H, IMG_W = 200, 200  # Enlarged grid canvas
N_CELLS = 50             # 50 cells per sample -> 100 cells total


def make_quadrant_mask(h, w, n_labels):
    """Simple non-overlapping label mask: split the image into n_labels
    horizontal bands, each with a unique integer label (1..n_labels)."""
    mask = np.zeros((h, w), dtype=np.int32)
    band = h // n_labels
    for i in range(n_labels):
        r0 = i * band
        r1 = (i + 1) * band if i < n_labels - 1 else h
        mask[r0:r1, :] = i + 1
    return mask

### New
def make_grid_mask(h, w, n_labels):
    """Creates a grid of n_labels cell masks across the canvas."""
    mask = np.zeros((h, w), dtype=np.int32)
    rows_cols = int(np.ceil(np.sqrt(n_labels)))
    dh, dw = h // rows_cols, w // rows_cols
    
    label = 1
    for i in range(rows_cols):
        for j in range(rows_cols):
            if label > n_labels:
                break
            r0, r1 = i * dh, min((i + 1) * dh, h)
            c0, c1 = j * dw, min((j + 1) * dw, w)
            mask[r0:r1, c0:c1] = label
            label += 1
    return mask


def make_image(seg1, seg2, n_channels):
    """Generates Poisson noise + cell-specific intensity signals across all channels."""
    img = RNG.poisson(lam=5, size=(n_channels, IMG_H, IMG_W)).astype(np.uint16)
    
    # Add signal to seg1 (Fine cell masks)
    for label in np.unique(seg1):
        if label == 0:
            continue
        # DAPI signal in all cells
        img[0][seg1 == label] += np.uint16(RNG.integers(50, 150))
        
        # Add random differential expression for all other channels across cells
        for ch in range(1, n_channels):
            if RNG.random() > 0.4: # ~60% chance cell expresses this marker
                img[ch][seg1 == label] += np.uint16(RNG.integers(30, 180))

    # Add extra signal to seg2 (Coarse masks) for specific marker
    for label in np.unique(seg2):
        if label == 0:
            continue
        img[1][seg2 == label] += np.uint16(RNG.integers(40, 100))

    return img
###


#def make_image(seg1, seg2, n_channels):
#    """Synthetic intensities: background noise everywhere, plus a
#    per-cell signal bump so quantification isn't just noise, and a
#    stronger bump inside seg2 regions for the 2nd marker (mimics a
#    coarser "big cell" marker like a megakaryocyte stain)."""
#    img = RNG.poisson(lam=5, size=(n_channels, IMG_H, IMG_W)).astype(np.uint16)
#    for label in np.unique(seg1):
#        if label == 0:
#            continue
#        img[0][seg1 == label] += np.uint16(RNG.integers(50, 150))  # DAPI signal in every cell
#    for label in np.unique(seg2):
#        if label == 0:
#            continue
#        img[1][seg2 == label] += np.uint16(RNG.integers(80, 200))  # 127 CD61 signal in MK regions
#    # 3rd marker: random signal, unrelated to segmentation
#    img[2] += RNG.integers(0, 60, size=(IMG_H, IMG_W)).astype(np.uint16)
#    return img


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "toy_data")

    images_dir = os.path.join(out_dir, "data", "cropped_images")
    seg1_dir = os.path.join(out_dir, "data", "Cell_masks_cropped")
    seg2_dir = os.path.join(out_dir, "data", "MK_masks_cropped")
    for d in (images_dir, seg1_dir, seg2_dir):
        os.makedirs(d, exist_ok=True)

    # marker list file (one marker per line, no header - matches
    # pd.read_csv(..., header=None) in the pipeline scripts)
    marker_list_path = os.path.join(out_dir, "data", "MarkerList.txt")
    with open(marker_list_path, "w") as f:
        f.write("\n".join(MARKERS) + "\n")

    for sample in SAMPLES:
        # fine-grained single-cell mask: 4 bands
        #seg1 = make_quadrant_mask(IMG_H, IMG_W, n_labels=4)
        seg1 = make_grid_mask(IMG_H, IMG_W, n_labels=N_CELLS)
        # coarser "big cell" mask: 2 bands
        #seg2 = make_quadrant_mask(IMG_H, IMG_W, n_labels=2)
        seg2 = make_quadrant_mask(IMG_H, IMG_W, n_labels=4)

        img = make_image(seg1, seg2, n_channels=len(MARKERS))

        # filenames follow the pipeline's expected convention: after
        # replacing each suffix with "run", all three keys must match,
        # e.g. Fusion_run99_sample1 -> run99_sample1
        #      Cell_mask_99_sample1 -> run99_sample1
        #      MK_mask_99_sample1   -> run99_sample1
        tifffile.imwrite(os.path.join(images_dir, f"Fusion_run{RUN_ID}_{sample}.tif"), img)
        tifffile.imwrite(os.path.join(seg1_dir, f"Cell_mask_{RUN_ID}_{sample}.tif"), seg1)
        tifffile.imwrite(os.path.join(seg2_dir, f"MK_mask_{RUN_ID}_{sample}.tif"), seg2)

        print(f"Wrote toy data for {sample}: image {img.shape}, seg1 {seg1.shape}, seg2 {seg2.shape}")

    print(f"\nToy dataset written to: {out_dir}")
    print(f"Marker list: {marker_list_path}")
    print("Point paths.wdir in toy_config.yaml at this directory (already done by default).")


if __name__ == "__main__":
    main()
