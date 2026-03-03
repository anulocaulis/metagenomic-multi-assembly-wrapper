#!/bin/bash
#SBATCH --job-name=flye_S5
#SBATCH --output=logs/flye_S5_%j.out
#SBATCH --error=logs/flye_S5_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=6-23:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

echo "Starting Flye S5 job at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/flye_assembler.sif"
INPUT_READS="${WORKDIR}/data/S5_long_reads_filtered.fastq.gz"
OUTDIR="${WORKDIR}/assemblies/S5/assembly.flye"
THREADS=16

mkdir -p "${OUTDIR}"

singularity exec \
    --bind "${WORKDIR}:${WORKDIR}" \
    "${CONTAINER_PATH}" \
    flye --nano-hq "${INPUT_READS}" --out-dir "${OUTDIR}" --threads "${THREADS}" --meta

echo "Flye S5 job finished at $(date)"
