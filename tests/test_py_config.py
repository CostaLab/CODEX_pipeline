"""
Unit tests for scripts/py_config.py — these only need pytest + PyYAML +
pandas, not the full spatialproteomics/scanpy/Seurat stack, so they're
cheap to run anywhere (e.g. in CI) as a first line of defense before
touching real data.

Run with:
    pytest tests/test_py_config.py -v
"""
import os
import sys
import textwrap

import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from py_config import load_config, resolve_paths, get_run_dicts, get_marker_channels  # noqa: E402


SAMPLE_CONFIG = textwrap.dedent(
    """
    paths:
      wdir: "{wdir}"
      input_images_subdir: "data/cropped_images"
      seg1_masks_subdir: "data/Cell_masks_cropped"
      seg2_masks_subdir: "data/MK_masks_cropped"
      marker_list_file: "data/MarkerList.txt"
      pipeline_run_subdir: "pipeline_run"
      zarr_subdir: "resulting_data/zarr"
      quant_subdir: "resulting_data/quantifications"
      visual_output_subdir: "visual_output"
      zarr_prefix: "test_"

    runs:
      - id: "1"
        conditions: ["a", "b"]
        working_samples: [1, 2]
      - id: "2"
        conditions: ["x", "y", "z"]
        # no working_samples given on purpose - should default to [1,2,3]
    """
)


@pytest.fixture
def cfg(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(SAMPLE_CONFIG.format(wdir=str(tmp_path)))
    return load_config(str(config_path))


def test_load_config_parses_yaml(cfg):
    assert cfg["paths"]["zarr_prefix"] == "test_"
    assert len(cfg["runs"]) == 2


def test_resolve_paths_builds_absolute_paths_and_creates_dirs(cfg, tmp_path):
    paths = resolve_paths(cfg)

    assert paths["wdir"] == str(tmp_path)
    assert paths["pipeline_run"] == str(tmp_path / "pipeline_run")
    assert paths["zarr_dir"] == str(tmp_path / "pipeline_run" / "resulting_data" / "zarr")
    assert paths["quant_dir"] == str(tmp_path / "pipeline_run" / "resulting_data" / "quantifications")
    assert paths["visual_output_dir"] == str(tmp_path / "pipeline_run" / "visual_output")

    # resolve_paths should have created the output dirs
    assert os.path.isdir(paths["zarr_dir"])
    assert os.path.isdir(paths["quant_dir"])
    assert os.path.isdir(paths["visual_output_dir"])


def test_get_run_dicts_keys_are_bare_run_ids(cfg):
    condition_vecs, working_samples = get_run_dicts(cfg)

    assert set(condition_vecs.keys()) == {"1", "2"}
    assert condition_vecs["1"] == ["a", "b"]
    assert condition_vecs["2"] == ["x", "y", "z"]


def test_get_run_dicts_defaults_working_samples_when_missing(cfg):
    _, working_samples = get_run_dicts(cfg)

    assert working_samples["1"] == [1, 2]
    # run "2" didn't specify working_samples -> should default to 1..len(conditions)
    assert working_samples["2"] == [1, 2, 3]


def test_get_marker_channels_reads_headerless_csv(tmp_path):
    marker_file = tmp_path / "MarkerList.txt"
    marker_file.write_text("DAPI\n127 CD61\n010 CD14\n")

    channels = get_marker_channels(str(marker_file))

    assert channels == ["DAPI", "127 CD61", "010 CD14"]
