"""
Shared config loader for all Python steps of the CODEX pipeline.

Usage in a script:

    from py_config import load_config, resolve_paths, get_run_dicts

    cfg = load_config(config_path)
    paths = resolve_paths(cfg)
    condition_vecs, working_samples = get_run_dicts(cfg)
"""
import os
import argparse
import yaml


def load_config(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg


def resolve_paths(cfg):
    """Turn the relative paths in config.yaml into absolute paths, and
    make sure output directories exist."""
    p = cfg["paths"]
    wdir = p["wdir"]
    pipeline_run = os.path.join(wdir, p["pipeline_run_subdir"])

    seg2_subdir = p.get("seg2_masks_subdir") or None  # blank/null/missing -> single-mask mode

    paths = {
        "wdir": wdir,
        "pipeline_run": pipeline_run,
        "input_images": os.path.join(wdir, p["input_images_subdir"]),
        "seg1_masks": os.path.join(wdir, p["seg1_masks_subdir"]),
        "seg2_masks": os.path.join(wdir, seg2_subdir) if seg2_subdir else None,
        "marker_list_file": os.path.join(wdir, p["marker_list_file"]),
        "zarr_dir": os.path.join(pipeline_run, p["zarr_subdir"]),
        "quant_dir": os.path.join(pipeline_run, p["quant_subdir"]),
        "visual_output_dir": os.path.join(pipeline_run, p["visual_output_subdir"]),
        "zarr_prefix": p["zarr_prefix"],
    }

    for key in ("zarr_dir", "quant_dir", "visual_output_dir"):
        os.makedirs(paths[key], exist_ok=True)

    return paths


def get_run_dicts(cfg):
    """Build the condition_vecs / working_samples dicts, keyed like the
    original scripts expected (bare run id, e.g. '21')."""
    condition_vecs = {}
    working_samples = {}
    for run in cfg["runs"]:
        rid = str(run["id"])
        condition_vecs[rid] = run["conditions"]
        working_samples[rid] = run.get("working_samples", list(range(1, len(run["conditions"]) + 1)))
    return condition_vecs, working_samples


def get_marker_channels(marker_list_file):
    import pandas as pd
    channels = pd.read_csv(marker_list_file, header=None)
    return list(channels[0])


def parse_args(description):
    """Every script accepts a single optional positional arg: the path
    to config.yaml (defaults to config.yaml next to the scripts dir's
    parent, i.e. ../config.yaml)."""
    parser = argparse.ArgumentParser(description=description)
    default_config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config.yaml")
    parser.add_argument(
        "config", nargs="?", default=default_config,
        help="Path to config.yaml (default: %(default)s)"
    )
    return parser.parse_args()
