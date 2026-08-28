#!/bin/bash
# Master pipeline submission script.
#
# Submits the 6 CODEX pipeline steps as SLURM jobs, chained with
# --dependency=afterok so each step only starts once the previous one
# finished successfully:
#
#   1. intensities.py            (segmentation + quantification -> zarr)
#   2. codex_vis.py               (intensity export + QC plots)
#   3. codex_vis_bgFilt.py         (background-filtered spatial plots)
#   4. codex_integration_P1.r      (build filtered Seurat object + SCT)
#   5. codex_integration_P2.r      (PCA/UMAP/Harmony + clustering sweep)
#   6. codex_integration_P3.r      (final clusters + marker plots)
#
# Usage:
#   ./run_pipeline.sh [config.yaml] [--from N] [--to N]
#
#   config.yaml   path to config file (default: ../config.yaml next to this script)
#   --from N      start at step N (1-6), skipping earlier steps (default 1)
#   --to N        stop after step N (1-6) (default 6)
#
# All per-step SLURM resources (time/mem/cpus/partition/account) are read
# straight out of config.yaml, so you only ever edit that one file.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="${SCRIPT_DIR}/config.yaml"
FROM_STEP=1
TO_STEP=6

# --- parse args ---
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM_STEP="$2"; shift 2 ;;
    --to)   TO_STEP="$2"; shift 2 ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done
if [[ ${#POSITIONAL[@]} -ge 1 ]]; then
  CONFIG_PATH="${POSITIONAL[0]}"
fi

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

echo "Using config: $CONFIG_PATH"
echo "Running steps ${FROM_STEP}..${TO_STEP}"

# --- pull slurm settings out of config.yaml via a tiny python helper ---
# (avoids a yq/python-yaml dependency at the shell level beyond python3+PyYAML,
#  which the pipeline already requires for the python steps)
yaml_get() {
  # yaml_get <dotted.path> [default]
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

PARTITION="$(yaml_get slurm.partition batch)"
ACCOUNT="$(yaml_get slurm.account)"
LOG_DIR="$(yaml_get paths.log_dir)"
SLURM_DIR="$(yaml_get paths.slurm_dir "${SCRIPT_DIR}/slurm")"
PY_CONDA_ENV="$(yaml_get slurm.py_conda_env)"
R_CONDA_ENV="$(yaml_get slurm.r_conda_env)"

mkdir -p "$LOG_DIR"

# submit_step <label> <sbatch_script> <resource_key> <lang: py|r> <dependency_jobid_or_empty>
submit_step () {
  local label="$1" script="$2" reskey="$3" lang="$4" dep="$5"

  local time mem cpus
  time="$(yaml_get "slurm.${reskey}.time" 24:00:00)"
  mem="$(yaml_get "slurm.${reskey}.mem" 32G)"
  cpus="$(yaml_get "slurm.${reskey}.cpus" 4)"

  local conda_env
  if [[ "$lang" == "py" ]]; then conda_env="$PY_CONDA_ENV"; else conda_env="$R_CONDA_ENV"; fi

  local sbatch_args=(
    --job-name="$label"
    --time="$time"
    --mem="$mem"
    --cpus-per-task="$cpus"
    --partition="$PARTITION"
    --output="${LOG_DIR}/${label}_%j.out"
    --error="${LOG_DIR}/${label}_%j.err"
    --export="ALL,CONDA_ENV_NAME=${conda_env}"
    --parsable
  )
  if [[ -n "$ACCOUNT" && "$ACCOUNT" != "None" ]]; then
    sbatch_args+=(--account="$ACCOUNT")
  fi
  if [[ -n "$dep" ]]; then
    sbatch_args+=(--dependency="afterok:${dep}")
  fi

  local jobid
  jobid="$(sbatch "${sbatch_args[@]}" "${SLURM_DIR}/${script}" "$CONFIG_PATH")"
  echo "  -> submitted ${label} (job ${jobid}, env: ${conda_env}, resources: ${reskey})" >&2
  echo "$jobid"
}

PREV_JOB=""

step_enabled () { [[ $1 -ge $FROM_STEP && $1 -le $TO_STEP ]]; }

if step_enabled 1; then
  echo "[1/6] intensities.py"
  PREV_JOB="$(submit_step codex_01_intensities 01_intensities.sbatch intensities py "$PREV_JOB")"
fi

if step_enabled 2; then
  echo "[2/6] codex_vis.py"
  PREV_JOB="$(submit_step codex_02_vis 02_codex_vis.sbatch codex_vis py "$PREV_JOB")"
fi

if step_enabled 3; then
  echo "[3/6] codex_vis_bgFilt.py"
  PREV_JOB="$(submit_step codex_03_vis_bgfilt 03_codex_vis_bgfilt.sbatch codex_vis_bgfilt py "$PREV_JOB")"
fi

if step_enabled 4; then
  echo "[4/6] codex_integration_P1.r"
  PREV_JOB="$(submit_step codex_04_integration_p1 04_integration_p1.sbatch integration_p1 r "$PREV_JOB")"
fi

if step_enabled 5; then
  echo "[5/6] codex_integration_P2.r"
  PREV_JOB="$(submit_step codex_05_integration_p2 05_integration_p2.sbatch integration_p2 r "$PREV_JOB")"
fi

if step_enabled 6; then
  echo "[6/6] codex_integration_P3.r"
  PREV_JOB="$(submit_step codex_06_integration_p3 06_integration_p3.sbatch integration_p3 r "$PREV_JOB")"
fi

echo "Done submitting. Track with: squeue -u \$USER"
echo "Logs will land in: $LOG_DIR"
