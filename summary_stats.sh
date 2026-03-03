#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <assembly_fasta> [assembly_fasta ...]"
  exit 1
fi

CONTAINER="containers/qc_tools_miniconda.sif"
LOG_FILE="summary_stats_log.txt"
WORKDIR="$(pwd)"

if [[ ! -f "$CONTAINER" ]]; then
  echo "Container not found: $CONTAINER" >&2
  exit 1
fi

for assembly in "$@"; do
  if [[ ! -f "$assembly" ]]; then
    echo "Assembly file not found: $assembly" >&2
    continue
  fi

  {
    echo "============================================================"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Assembly: $assembly"
    echo "============================================================"

    singularity exec -B "$WORKDIR:$WORKDIR" "$CONTAINER" \
      stats.sh in="$assembly"

    echo
  } >> "$LOG_FILE" 2>&1

done

echo "Appended stats to $LOG_FILE"