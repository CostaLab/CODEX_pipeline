#!/bin/bash
#SBATCH --job-name=codex_03_vis_bgfilt
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

CONFIG_PATH="${1:?Usage: sbatch 03_codex_vis_bgfilt.sbatch <config.yaml>}"

# --- environment ---
# CONDA_ENV_NAME is set by run_pipeline.sh from config.yaml
# (slurm.py_conda_env). Falls back to plain PATH if unset,
# so this also works if you call sbatch directly with your env already active.
if [[ -n "${CONDA_ENV_NAME:-}" ]]; then
  eval "$(conda shell.bash hook)"
  # Some envs' activate.d hooks (e.g. conda-forge's gcc_linux-64 compiler
  # activation scripts) reference variables like $SYS_SYSROOT without
  # guarding for "unset" and aren't written to be `set -u`-safe. Relax
  # nounset just for the activation call itself, then restore it.
  set +u
  conda activate "$CONDA_ENV_NAME"
  set -u
fi

#SCRIPT_DIR="$(dirname "$(readlink -f "$0")")/../scripts"
#python "$SCRIPT_DIR/codex_vis_bgFilt.py" "$CONFIG_PATH"

# Parse scripts_dir directly out of config.yaml
SCRIPTS_DIR="$(python3 - "$CONFIG_PATH" <<'PYEOF'
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
print(cfg.get('paths', {}).get('scripts_dir', ''))
PYEOF
)"

if [[ -z "$SCRIPTS_DIR" || ! -d "$SCRIPTS_DIR" ]]; then
  echo "Error: scripts_dir ('$SCRIPTS_DIR') not found or not set in $CONFIG_PATH" >&2
  exit 1
fi

python "$SCRIPTS_DIR/codex_vis_bgFilt.py" "$CONFIG_PATH"
