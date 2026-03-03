#!/bin/bash
#SBATCH --job-name=metamdbg
#SBATCH --output=logs/metamdbg_%j.out
#SBATCH --error=logs/metamdbg_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=6-23:00:00
#SBATCH --partition=math-alderaan

set -euo pipefail

echo "Starting metaMDBG job at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/assembler2.sif"
THREADS=16
FORCE_RERUN=${FORCE_RERUN:-0}

cd "${WORKDIR}"

if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=(S1 S2 S5)
fi

for sample in "${SAMPLES[@]}"; do
    INPUT_READS="${WORKDIR}/data/${sample}_long_reads_filtered.fastq.gz"
    if [ "${sample}" = "S1" ]; then
        OUTDIR="${WORKDIR}/assemblies/${sample}/assembly.metamdbg_s1_retry"
    else
        OUTDIR="${WORKDIR}/assemblies/${sample}/assembly.metamdbg"
    fi

    echo ""
    echo "=== Running metaMDBG for ${sample} at $(date) ==="
    echo "Output dir: ${OUTDIR}"

    if [ ! -s "${INPUT_READS}" ]; then
        echo "Missing input reads for ${sample}: ${INPUT_READS}" >&2
        continue
    fi

    if [ "${FORCE_RERUN}" != "1" ] && [ -s "${OUTDIR}/contigs.fasta" ]; then
        echo "Existing contigs found for ${sample}; skipping. Set FORCE_RERUN=1 to rerun."
        continue
    fi

    if [ "${FORCE_RERUN}" = "1" ] && [ -d "${OUTDIR}" ]; then
        echo "FORCE_RERUN=1; removing existing output dir for ${sample}: ${OUTDIR}"
        rm -rf "${OUTDIR}"
    elif [ -d "${OUTDIR}" ] && [ ! -s "${OUTDIR}/contigs.fasta" ] && [ ! -s "${OUTDIR}/contigs.fasta.gz" ]; then
        echo "Cleaning stale incomplete output for ${sample}: ${OUTDIR}"
        rm -rf "${OUTDIR}"
    fi

    mkdir -p "${OUTDIR}"

    if ! singularity exec \
        --bind "${WORKDIR}:${WORKDIR}" \
        "${CONTAINER_PATH}" \
        metaMDBG asm --threads "${THREADS}" --in-ont "${INPUT_READS}" --out-dir "${OUTDIR}"; then
        echo "metaMDBG command failed for ${sample}" >&2
        continue
    fi

    if [ -f "${OUTDIR}/contigs.fasta.gz" ]; then
        gunzip -f "${OUTDIR}/contigs.fasta.gz"
    fi

    if [ -s "${OUTDIR}/contigs.fasta" ]; then
        echo "metaMDBG completed for ${sample}: ${OUTDIR}/contigs.fasta"
    else
        echo "metaMDBG did not produce contigs.fasta for ${sample}" >&2
        continue
    fi
done

echo "metaMDBG job finished at $(date)"
