#!/usr/bin/env bash
#SBATCH --job-name=plot_summary
#SBATCH --output=logs/plot_summary_%j.out
#SBATCH --error=logs/plot_summary_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORKDIR="$SLURM_SUBMIT_DIR"
else
  WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
cd "$WORKDIR"

LOG_INPUT="${LOG:-summary_stats_log.txt}"
OUTDIR_INPUT="${OUTDIR:-plots/summary_stats}"

if [[ "$LOG_INPUT" = /* ]]; then
  LOG_PATH="$LOG_INPUT"
else
  LOG_PATH="$WORKDIR/$LOG_INPUT"
fi

if [[ "$OUTDIR_INPUT" = /* ]]; then
  OUTDIR_PATH="$OUTDIR_INPUT"
else
  OUTDIR_PATH="$WORKDIR/$OUTDIR_INPUT"
fi

if [[ ! -f "$LOG_PATH" ]]; then
  echo "Log file not found: $LOG_PATH" >&2
  echo "WORKDIR: $WORKDIR" >&2
  echo "LOG env/raw value: ${LOG:-<unset>}" >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-$WORKDIR/.venv/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || true)"
fi

if [[ -z "$PYTHON_BIN" ]]; then
  echo "No Python interpreter found (.venv/bin/python or python3)." >&2
  exit 1
fi

if ! "$PYTHON_BIN" -c "import pandas, matplotlib, seaborn" >/dev/null 2>&1; then
  echo "Missing Python packages in interpreter: $PYTHON_BIN" >&2
  echo "Required: pandas, matplotlib, seaborn" >&2
  echo "Install them in your environment before submitting." >&2
  exit 1
fi

"$PYTHON_BIN" plot_summary_stats.py --log "$LOG_PATH" --outdir "$OUTDIR_PATH"

echo "Finished plotting. Outputs in: $OUTDIR_PATH"
