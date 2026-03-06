#!/usr/bin/env bash
#SBATCH --job-name=summary_stats
#SBATCH --output=logs/summary_stats_%j.out
#SBATCH --error=logs/summary_stats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <assemblies_root_dir | assembly_fasta ...>"
  exit 1
fi

CONTAINER="containers/qc_tools_miniconda.sif"
LOG_FILE="summary_stats_log.txt"
WORKDIR="$(pwd)"

if [[ ! -f "$CONTAINER" ]]; then
  echo "Container not found: $CONTAINER" >&2
  exit 1
fi

detect_assembler() {
  local filename="$1"
  case "$filename" in
    megahit.final.contigs.fa) echo "megahit" ;;
    metaspades.contigs.fasta) echo "metaspades" ;;
    metaspades_hybrid.contigs.fasta|metaspades_hybrid.assembly.fasta) echo "metaspades_hybrid" ;;
    flye.assembly.fasta) echo "flye" ;;
    metamdbg.contigs.fasta) echo "metamdbg" ;;
    idbaud.assembly.fasta) echo "idbaud" ;;
    contigs.fasta) echo "unknown" ;;
    *)
      if [[ "$filename" == *.* ]]; then
        echo "${filename%%.*}"
      else
        echo "unknown"
      fi
      ;;
  esac
}

run_stats_for_assembly() {
  local sample="$1"
  local assembly="$2"
  local filename
  filename="$(basename "$assembly")"
  local assembler
  assembler="$(detect_assembler "$filename")"

  {
    echo "============================================================"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Sample: $sample"
    echo "Assembler: $assembler"
    echo "Assembly: $assembly"
    echo "============================================================"

    singularity exec -B "$WORKDIR:$WORKDIR" "$CONTAINER" \
      stats.sh in="$assembly"

    echo
  } >> "$LOG_FILE" 2>&1
}

if [[ $# -eq 1 && -d "$1" ]]; then
  assemblies_root="$1"
  shopt -s nullglob
  for sample_dir in "$assemblies_root"/S*; do
    [[ -d "$sample_dir" ]] || continue
    sample="$(basename "$sample_dir")"

    files=("$sample_dir"/*.fa "$sample_dir"/*.fasta)
    if [[ ${#files[@]} -eq 0 ]]; then
      echo "No FASTA files found in sample directory: $sample_dir" >&2
      continue
    fi

    for assembly in "${files[@]}"; do
      [[ -f "$assembly" ]] || continue
      run_stats_for_assembly "$sample" "$assembly"
    done
  done
  shopt -u nullglob
else
  for assembly in "$@"; do
    if [[ ! -f "$assembly" ]]; then
      echo "Assembly file not found: $assembly" >&2
      continue
    fi

    sample="$(basename "$(dirname "$assembly")")"
    run_stats_for_assembly "$sample" "$assembly"
  done
fi

echo "Appended stats to $LOG_FILE"