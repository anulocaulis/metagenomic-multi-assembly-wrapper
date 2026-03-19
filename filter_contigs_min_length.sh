#!/bin/bash
#SBATCH --job-name=filter_contigs
#SBATCH --output=logs/filter_contigs_%j.out
#SBATCH --error=logs/filter_contigs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-02:00:00
#SBATCH --partition=math-alderaan-gpu

# Usage:
#   sbatch filter_contigs_min_length.sh <SAMPLE> <ASSEMBLER>
# Example:
#   sbatch filter_contigs_min_length.sh S1 metaspades

set -euo pipefail

SAMPLE="${1:?ERROR: SAMPLE argument required (e.g. S1)}"
ASSEMBLER="${2:?ERROR: ASSEMBLER argument required (e.g. metaspades)}"

echo "Starting filter_contigs_min_length for ${SAMPLE} / ${ASSEMBLER} at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/qc_tools_miniconda.sif"
INPUT_CONTIGS="${WORKDIR}/assemblies/${SAMPLE}/assembly.${ASSEMBLER}/contigs.fasta"
OUTPUT_FILTERED="${WORKDIR}/assemblies/${SAMPLE}/assembly.${ASSEMBLER}/contigs.ge1000.fa"

if [[ ! -f "${INPUT_CONTIGS}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_CONTIGS}" >&2
    exit 1
fi

singularity exec \
    --bind "${WORKDIR}:${WORKDIR}" \
    "${CONTAINER_PATH}" \
    bbduk.sh \
        in="${INPUT_CONTIGS}" \
        out="${OUTPUT_FILTERED}" \
        minlength=1000

echo "Done. Output: ${OUTPUT_FILTERED}"
echo "Finished at $(date)"
