# CODEX_pipeline

[![License: MIT](https://shields.io)](https://opensource.org)
[![Python 3.8+](https://shields.io)](https://python.org)

This repository contains a collection of scripts for the analysis and processing of high-dimensional **CODEX** (Co-detection by Indexing) multiplexed imaging data.
The pipeline enables you to run these scripts one after another on your CODEX data.
This CODEX analysis pipeline was developed while analysing data from fibrotic spleens from fibrosis mouse models published here: https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(26)00272-9 

---

## 🚀 Features

- **Single-Cell Feature Extraction**: Quantitative extraction of marker expression levels (e.g., average brightness) per cell using an image and its corresponding segmentation mask.
- **Data integration**: Integration of data from different images for comparative analyses.
- **Downstream Visualization**: Data visualization and conversion of multiplexed data into h5ad format.

## 📦 Installation

To get started, clone the repository in your project folder. The project folder where this pipeline was tested is called **CODEX_pipeline_test** 

```bash
# Clone the repository
git clone -b dev --single-branch https://github.com/CostaLab/CODEX_pipeline
# later on -> git clone https://github.com/CostaLab/CODEX_pipeline

```
## Layout
After cloning this repository, **my project folder** looks like:

```
CODEX_pipeline_test/
├── data/
|    ├── MarkerList.txt                List of marker names (sorted as in raw images)
|    ├── cell_mask/                    DAPI segmentation masks folder
|    |    |── cell_mask_run46_EV.tif
|    |    |── cell_mask_run46_Thpo.tif
|    ├── image/                        Raw images folder
|    |    |── run46_EV.tif
|    |    |── run46_Thpo.tif
|    └── MK_mask/                      Megakaryocyte masks folder
|         |── MK_mask_run46_EV.tif
|         |── MK_mask_run46_Thpo.tif
└── CODEX_pipeline/
```

And the **CODEX_pipeline** folder itself looks like:

```
CODEX_pipeline/
├── config.yaml                  <- edit this
├── run_pipeline.sh              <- submits all steps to SLURM
├── scripts/
│   ├── py_config.py             (Python config loader, imported by the 3 .py steps)
│   ├── r_config.R               (R config loader, sourced by the 3 .r steps)
│   ├── intensities.py           step 1: segmentation + quantification -> zarr
│   ├── codex_vis.py             step 2: intensity export + QC plots
│   ├── codex_vis_bgFilt.py      step 3: background-filtered spatial plots
│   ├── codex_integration_P1.r   step 4: filtered Seurat object + SCTransform
│   ├── codex_integration_P2.r   step 5: PCA/UMAP/Harmony + clustering sweep
│   └── codex_integration_P3.r   step 6: final clusters + marker plots
├── slurm/
│   ├── 01_intensities.sbatch
│   ├── 02_codex_vis.sbatch
│   ├── 03_codex_vis_bgfilt.sbatch
│   ├── 04_integration_p1.sbatch
│   ├── 05_integration_p2.sbatch
│   └── 06_integration_p3.sbatch
├── example/
│   ├── generate_toy_data.py     makes a tiny synthetic dataset (see "Toy example" below)
│   └── toy_config.yaml          config pointing at that toy dataset
└── tests/
    ├── test_py_config.py        fast unit tests for the config loader (pytest, no heavy deps)
    └── run_smoke_test.sh        runs all 6 steps locally against the toy dataset
```

## How to run the CODEX pipeline

A single `config.yaml` drives all six steps. **Make a copy of this file** and edit it to change
paths, run/condition definitions, markers, thresholds, or SLURM resources —
you should not need to touch the scripts themselves for routine reruns.

## Before first run

1. Copy and edit `config.yaml`:
   - `paths.wdir` — your project root.
   - `paths.scripts_dir` / `paths.slurm_dir` — where you copy this `pipeline/`
     folder to on the cluster (used only as documentation/reference; the
     scripts locate each other relatively, so this mostly matters if you
     want other tooling to know where things live).
   - `runs` — one entry per CODEX run, with `conditions` and `working_samples`
     matched by position (both 1-based): `conditions[i]` is the condition
     label for `working_samples[i]`'s file. **Ordering matters and is
     alphabetical**: `conditions[1]` must correspond to whichever of that
     run's sample files sorts first alphabetically, `conditions[2]` to the
     second, and so on. Bottom line: **before filling in `conditions`, list that run's
     files yourself (e.g. `ls resulting_data/quantifications/ | grep run21
     | sort`) and read off the condition for each file in that order.**
   - `integration.not_working_markers` — must be given as **R-safe names**
     (i.e. what `make.names()` turns your CSV headers into: spaces/parens
     become dots, a leading digit gets an `X` prefix, e.g.
     `"050 CD138"` → `"X050.CD138"`).
   - `slurm.py_conda_env` / `slurm.r_conda_env` — your two separate conda
     envs (e.g. one with spatialproteomics/scanpy/anndata, one with
     Seurat/harmony/sctransform/yaml). You do **not** need to unify them —
     `run_pipeline.sh` reads both from config.yaml and passes the right one
     to each job as an env var; each `.sbatch` file activates it with
     `eval "$(conda shell.bash hook)"; conda activate "$CONDA_ENV_NAME"`.
     If you submit an `.sbatch` file directly (bypassing `run_pipeline.sh`)
     with your env already active, the activation step is skipped safely.
2. Make sure `PyYAML` is installed in the Python env (`pip install pyyaml`)
   and the `yaml` package is installed in R (`install.packages("yaml")`) —
   both config loaders depend on it.

## Running

```bash
cd CODEX_pipeline
./run_pipeline.sh config.yaml   # submit all 6 steps, chained with afterok dependencies, using your config file
./run_pipeline.sh config.yaml --from 4         # re-run only steps 4-6 (e.g. after tweaking integration params)
./run_pipeline.sh config.yaml --to 3           # only run the Python steps
```

Each `.sbatch` file also runs standalone if you ever want to submit a single
step by hand: `sbatch slurm/04_integration_p1.sbatch config.yaml`.

SLURM logs land in `paths.log_dir` (default `<wdir>/pipeline_run/logs`).

## Outputs

Everything lands under `<wdir>/<pipeline_run_subdir>` (default `pipeline_run/`),
in the subfolders set by `config.yaml`. `<region>` below means a region key
like `run21_...`; `<run>`, `<condition>`, `<j>` are the run id, condition
label, and 1-based sample index for that region (e.g. `run21_d180_0`).

**Step 1 — `intensities.py`**
| File | Location |
|---|---|
| `<zarr_prefix><region>.zarr` — quantified per-cell data (image, both segmentations, merged mask, intensities, area) | `resulting_data/zarr/` |
| `segmentation_plot_pro_ind_<region>.png` — seg1 vs. seg2 masks side by side, pre-merge | `visual_output/` |
| `segmentation_plot_pro_<region>.png` — final merged segmentation mask | `visual_output/` |

**Step 2 — `codex_vis.py`**
| File | Location |
|---|---|
| `<zarr_prefix><region><intensities_file_suffix>` (CSV) — per-cell marker intensities | `resulting_data/quantifications/` |
| `allmks_hist_run<run>_<condition>_<j>.png` — histogram grid, one panel per marker | `visual_output/` |
| `run<run>_<condition>_<j><marker>_scatter.png` — one file per entry in `visualization.scatter_markers`, pairwise scatter of that marker vs. every other marker | `visual_output/` |

**Step 3 — `codex_vis_bgFilt.py`**
| File | Location |
|---|---|
| `showrun<run>_<condition>_<j>_mks_bgFilt.png` — spatial marker maps after per-marker median background subtraction (scanpy's `sc.pl.spatial(..., save=...)` prefixes the filename with `show`, at least in the scanpy version this was tested against — check your own `figures/` folder if it differs) | `figures/` next to wherever the job runs (scanpy's own convention, not `visual_output/`) |

**Step 4 — `codex_integration_P1.r`**
| File | Location |
|---|---|
| `ncells_sample.pdf` / `ncells_sample_filt.pdf` — cells per sample, before/after filtering | `visual_output/` |
| `<combination_name>_cells_filtering.pdf` — one per entry in `integration.aberrant_combinations`, showing the filtering threshold | `visual_output/` |
| `qc.pdf` — nFeature/nCount violin plots | `visual_output/` |
| `cur_data.csv` — combined, filtered intensity table across all runs/samples | `pipeline_run/` |
| `seurat_object_filt_bg_raw.rds` — Seurat object, filtered, background-subtracted, pre-SCT | `pipeline_run/` |
| `seurat_object_sct_filt_bg.rds` — same, after `SCTransform` | `pipeline_run/` |

**Step 5 — `codex_integration_P2.r`**
| File | Location |
|---|---|
| `pca_elbow_plot_filt_bg.png` | `visual_output/` |
| `umap_sample_sct_filt_bg.pdf` / `umap_run_sct_filt_bg.pdf` — pre-Harmony UMAPs colored by sample/run | `visual_output/` |
| `clustree_harmony_filt_bg.pdf` — resolution sweep tree (`clustree`) over `integration.clustering.res_init`–`res_final` | `visual_output/` |
| `seurat_object_harmony_filt_bg.rds` — after Harmony batch correction | `pipeline_run/` |
| `seurat_object_harmony_clusts_filt_bg.rds` — with all swept-resolution cluster labels attached | `pipeline_run/` |

**Step 6 — `codex_integration_P3.r`**
| File | Location |
|---|---|
| `umap_clusters_final.pdf` — final UMAP at `integration.chosen_resolution` | `visual_output/` |
| `umap_clusters_by_condition.pdf` — final UMAP colored by condition | `visual_output/` |
| `dotplot_markers_by_cluster.pdf` / `featureplots_markers.pdf` — marker expression per cluster | `visual_output/` |
| `cluster_composition_by_condition.csv` / `cluster_composition_by_sample.csv` — cell counts per cluster × condition/sample | `pipeline_run/` |
| `seurat_object_final.rds` — final annotated Seurat object | `pipeline_run/` |

Plus, for every step: SLURM `.out`/`.err` logs in `paths.log_dir`.

## Toy example & tests

`example/` has a self-contained, tiny synthetic dataset (1 run, 2 samples,
30 markers, 200x200px images) that runs through the whole pipeline in seconds.
It's useful for:

- **Sanity-checking your setup** — confirms both conda envs, the config
  structure, and the script logic all work together *before* you submit a
  real SLURM job that might sit in the queue for hours only to fail on a
  typo.
- **Regression testing** — if you modify any script, rerun the smoke test
  first. It exercises the full code path (real `spatialproteomics` calls,
  real Seurat/Harmony calls) rather than mocking anything out.

### Generate the toy dataset

```bash
conda activate <your-python-env>
python3 example/generate_toy_data.py /tmp/codex_toy
```

This writes 2 synthetic images + 2 pairs of segmentation masks + a 30-marker
list under `/tmp/codex_toy/data/...`. `example/toy_config.yaml` already
points at `/tmp/codex_toy` — edit `paths.wdir` in that file if you generate
the data somewhere else.

### Run the smoke test

```bash
./tests/run_smoke_test.sh              # runs and checks all 6 steps
./tests/run_smoke_test.sh --skip-r      # just steps 1-3, fastest way to check the Python side
```

You do **not** need to `conda activate` anything yourself first — the script
reads `slurm.py_conda_env` / `slurm.r_conda_env` out of `toy_config.yaml`
and activates each one itself for the relevant steps, the same way
`run_pipeline.sh` does for the real SLURM jobs. It just needs `conda` to be
on `PATH` (or already initialized) in the shell you run it from.

It regenerates the toy dataset fresh each run (add `--keep-data` to reuse an
existing one), runs all 6 steps locally (no SLURM — plain `python3`/`Rscript`
calls in sequence), and checks that each step produced its key output files,
printing `OK`/`MISSING` for each and exiting non-zero if anything's missing.
If the configured env name doesn't exist, it fails immediately with a clear
message rather than silently running in whatever env happened to be active.

This toy run is also how the pipeline itself was validated while building
it — steps 1-3 were run end-to-end against this exact synthetic dataset
during development, which is how a real `spatialproteomics` API version
mismatch was caught (see note below).

### Fast unit tests (no data, no heavy deps)

`tests/test_py_config.py` tests `scripts/py_config.py` in isolation — config
parsing, path resolution, run/condition-dict building — using only `pytest`
and `PyYAML`. No `spatialproteomics`/`scanpy` install needed, so these are
cheap enough to run in CI on every commit:

```bash
pip install pytest pyyaml
pytest tests/test_py_config.py -v
```

### A version pitfall this caught

While validating the toy example, `intensities.py` failed against the
latest PyPI `spatialproteomics` (0.8.1+) with:
`TypeError: merge_segmentation() got an unexpected keyword argument 'threshold'`.
Versions ≤0.8.0 take a single `threshold=` argument (what your original
scripts and this refactor both use); 0.8.1+ replaced it with
`threshold1`/`threshold2`. If you hit this on a real run, either pin
`spatialproteomics<=0.8.0` in your `py_conda_env`, or adapt
`intensities.py`'s `merge_segmentation(...)` call to the newer signature —
whichever matches what the rest of your environment already relies on.
