#!/bin/bash
# Smoke test: regenerates the toy dataset, runs all 6 pipeline steps
# locally (no SLURM — just plain python/Rscript calls, in order), and
# checks that every step produced the files it's supposed to.
#
# Good for:
#   - a fast sanity check after editing any of the scripts
#   - confirming both conda envs are set up correctly, before ever
#     touching SLURM or real data
#
# This script activates the right conda env itself for each step, reading
# slurm.py_conda_env / slurm.r_conda_env straight out of toy_config.yaml
# (the same envs run_pipeline.sh uses for the real SLURM jobs) - you don't
# need to conda activate anything yourself first.
#
# Usage:
#   ./tests/run_smoke_test.sh [--skip-r] [--keep-data]
#
#   --skip-r     only run + check steps 1-3 (Python). Useful if your R env
#                isn't set up yet, or isn't reachable from this shell.
#   --keep-data  don't regenerate the toy dataset (reuse whatever is
#                already at TOY_DIR from a previous run)

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"
EXAMPLE_DIR="${PIPELINE_DIR}/example"
CONFIG_PATH="${EXAMPLE_DIR}/toy_config.yaml"
TOY_DIR="/tmp/codex_toy_smoketest"

SKIP_R=false
KEEP_DATA=false
for arg in "$@"; do
  case "$arg" in
    --skip-r) SKIP_R=true ;;
    --keep-data) KEEP_DATA=true ;;
    *) echo "Unknown argument: $arg" >&2; exit 1 ;;
  esac
done

PASS=0
FAIL=0

check_file () {
  if [[ -e "$1" ]]; then
    echo "  OK   $1"
    PASS=$((PASS + 1))
  else
    echo "  MISSING   $1"
    FAIL=$((FAIL + 1))
  fi
}

# --- pull conda env names out of toy_config.yaml (same helper run_pipeline.sh uses) ---
yaml_get() {
  python3 - "$CONFIG_PATH" "$1" "${2:-}" <<'PYEOF'
import sys, yaml
cfg_path, dotted, default = sys.argv[1], sys.argv[2], sys.argv[3]
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)
node = cfg
for key in dotted.split('.'):
    if node is None:
        node = None
        break
    node = node.get(key)
print(node if node is not None else default)
PYEOF
}

PY_CONDA_ENV="$(yaml_get slurm.py_conda_env)"
R_CONDA_ENV="$(yaml_get slurm.r_conda_env)"

if ! command -v conda &>/dev/null; then
  echo "ERROR: 'conda' not found on PATH. Load/init conda first (e.g. 'module load anaconda3' or source your conda.sh), then re-run." >&2
  exit 1
fi
eval "$(conda shell.bash hook)"

activate_env () {
  local env_name="$1"
  if [[ -z "$env_name" || "$env_name" == "None" ]]; then
    echo "ERROR: no conda env name set for this step (check slurm.py_conda_env / slurm.r_conda_env in $CONFIG_PATH)." >&2
    exit 1
  fi
  echo "-- activating conda env: $env_name --"
  # Some envs' activate.d hooks (e.g. conda-forge's gcc_linux-64 compiler
  # activation scripts) reference variables like $SYS_SYSROOT without
  # guarding for "unset" and aren't written to be `set -u`-safe. Relax
  # nounset just for the activation call itself, then restore it.
  set +u
  conda activate "$env_name"
  set -u
}

echo "== Smoke test config: using a throwaway toy_config.yaml pointed at $TOY_DIR =="
TEST_CONFIG="$(mktemp)"
sed "s|/tmp/codex_toy|${TOY_DIR}|g" "$CONFIG_PATH" > "$TEST_CONFIG"

activate_env "$PY_CONDA_ENV"

if [[ "$KEEP_DATA" == false ]]; then
  echo "== Generating toy dataset at $TOY_DIR =="
  rm -rf "$TOY_DIR"
  python3 "${EXAMPLE_DIR}/generate_toy_data.py" "$TOY_DIR"
fi

RUN_DIR="${TOY_DIR}/pipeline_run"

echo
echo "== [1/6] intensities.py =="
python3 "${SCRIPTS_DIR}/intensities.py" "$TEST_CONFIG"
check_file "${RUN_DIR}/resulting_data/zarr/test_run99_sample1.zarr"
check_file "${RUN_DIR}/resulting_data/zarr/test_run99_sample2.zarr"
check_file "${RUN_DIR}/visual_output/segmentation_plot_pro_run99_sample1.png"

echo
echo "== [2/6] codex_vis.py =="
python3 "${SCRIPTS_DIR}/codex_vis.py" "$TEST_CONFIG"
check_file "${RUN_DIR}/resulting_data/quantifications/test_run99_sample1_asinh.csv"
check_file "${RUN_DIR}/resulting_data/quantifications/test_run99_sample2_asinh.csv"
check_file "${RUN_DIR}/visual_output/allmks_hist_run99_ctrl_0.png"

echo
echo "== [3/6] codex_vis_bgFilt.py =="
python3 "${SCRIPTS_DIR}/codex_vis_bgFilt.py" "$TEST_CONFIG"
# scanpy writes these into ./figures relative to cwd, not visual_output/
check_file "./figures/showrun99_ctrl_0_mks_bgFilt.png"

if [[ "$SKIP_R" == true ]]; then
  echo
  echo "Skipping R steps (--skip-r). Python steps: $PASS passed, $FAIL failed."
  [[ $FAIL -eq 0 ]] && exit 0 || exit 1
fi

activate_env "$R_CONDA_ENV"

echo
echo "== [4/6] codex_integration_P1.r =="
Rscript "${SCRIPTS_DIR}/codex_integration_P1.r" "$TEST_CONFIG"
check_file "${RUN_DIR}/cur_data.csv"
check_file "${RUN_DIR}/seurat_object_sct_filt_bg.rds"
check_file "${RUN_DIR}/visual_output/qc.pdf"

echo
echo "== [5/6] codex_integration_P2.r =="
Rscript "${SCRIPTS_DIR}/codex_integration_P2.r" "$TEST_CONFIG"
check_file "${RUN_DIR}/seurat_object_harmony_clusts_filt_bg.rds"
check_file "${RUN_DIR}/visual_output/clustree_harmony_filt_bg.pdf"

echo
echo "== [6/6] codex_integration_P3.r =="
Rscript "${SCRIPTS_DIR}/codex_integration_P3.r" "$TEST_CONFIG"
check_file "${RUN_DIR}/seurat_object_final.rds"
check_file "${RUN_DIR}/visual_output/umap_clusters_final.pdf"
check_file "${RUN_DIR}/cluster_composition_by_condition.csv"

echo
echo "== Summary: $PASS passed, $FAIL failed =="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
