#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATS_SCRIPT="$SCRIPT_DIR/summary_stats.sh"
ASSEMBLY_ROOT="$SCRIPT_DIR/assemblies"

if [[ ! -x "$STATS_SCRIPT" ]]; then
  echo "Missing or non-executable stats script: $STATS_SCRIPT" >&2
  exit 1
fi

if [[ ! -d "$ASSEMBLY_ROOT" ]]; then
  echo "Assemblies directory not found: $ASSEMBLY_ROOT" >&2
  exit 1
fi

mapfile -t assembly_files < <(
  {
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.flye/assembly.fasta'
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.metamdbg/contigs.fasta'
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.idbaud/assembly.fasta'
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.metaspades/assembly.fasta'
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.megahit/final.contigs.fa'
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.metaspades_hybrid/assembly.fasta'
    find "$ASSEMBLY_ROOT" -type f -path '*/assembly.nextpolish/assembly.fasta'
  } | sort -u
)

if [[ ${#assembly_files[@]} -eq 0 ]]; then
  echo "No assembly FASTA files found under $ASSEMBLY_ROOT"
  exit 0
fi

"$STATS_SCRIPT" "${assembly_files[@]}"

echo "Processed ${#assembly_files[@]} assembly files"
