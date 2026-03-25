#!/bin/bash
#SBATCH --job-name=preasm_stats
#SBATCH --output=logs/preassembly_stats_%j.out
#SBATCH --error=logs/preassembly_stats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONTAINER="${SCRIPT_DIR}/containers/qc_tools_miniconda.sif"
LOG_FILE="${SCRIPT_DIR}/preassembly_stats_log.txt"
WORKDIR="${SCRIPT_DIR}"

cd "${SCRIPT_DIR}"

echo "Starting pre-assembly BBTools stats job at $(date)"

if [[ ! -f "${CONTAINER}" ]]; then
  echo "Container not found: ${CONTAINER}" >&2
  exit 1
fi

declare -a read_files
if [[ $# -gt 0 ]]; then
  read_files=("$@")
else
  mapfile -t read_files < <(
    {
      find "${SCRIPT_DIR}/trimmed_reads" -type f -name '*_interleaved_trimmed_polyG_filtered.fastq.gz'
      find "${SCRIPT_DIR}/data" -type f -name '*_long_reads_filtered.fastq.gz'
    } | sort -u
  )
fi

if [[ ${#read_files[@]} -eq 0 ]]; then
  echo "No pre-assembly read files found"
  exit 0
fi

for reads in "${read_files[@]}"; do
  if [[ ! -f "${reads}" ]]; then
    echo "Read file not found: ${reads}" >&2
    continue
  fi

  {
    echo "============================================================"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Input reads: ${reads}"
    echo "============================================================"

    singularity exec -B "${WORKDIR}:${WORKDIR}" "${CONTAINER}" \
      stats.sh in="${reads}"

    echo
  } >> "${LOG_FILE}" 2>&1
done

echo "Processed ${#read_files[@]} pre-assembly read files"
echo "Appended stats to ${LOG_FILE}"
echo "Finished pre-assembly BBTools stats job at $(date)"
