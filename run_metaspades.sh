#!/bin/bash
#SBATCH --job-name=metaspades_miss
#SBATCH --output=logs/metaspades_%j.out
#SBATCH --error=logs/metaspades_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=900G
#SBATCH --time=2-00:00:00
#SBATCH --partition=math-alderaan-gpu

set -euo pipefail

echo "Starting metaSPAdes job at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/assembler2.sif"
THREADS=24
MEM_GB=900

cd "${WORKDIR}"

if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=(S1 S3 S4 S5 S6 S7 S8 S10 S11 S12)
fi

for sample in "${SAMPLES[@]}"; do
    INPUT_READS="${WORKDIR}/trimmed_reads/${sample}_interleaved_trimmed.fastq.gz"
    OUTDIR="${WORKDIR}/assemblies/${sample}/assembly.metaspades"

    echo ""
    echo "=== Running metaSPAdes for ${sample} at $(date) ==="
    echo "Input: ${INPUT_READS}"
    echo "Output dir: ${OUTDIR}"

    if [ ! -s "${INPUT_READS}" ]; then
        echo "Missing input reads for ${sample}: ${INPUT_READS}" >&2
        continue
    fi

    mkdir -p "${OUTDIR}"

    if [ -f "${OUTDIR}/run_spades.yaml" ] || [ -d "${OUTDIR}/pipeline_state" ]; then
        echo "Resume detected for ${sample}; trying --continue"
        if ! singularity exec -B "$PWD" "${CONTAINER_PATH}" metaspades.py --continue -o "${OUTDIR}"; then
            echo "--continue failed for ${sample}; forcing clean rerun" >&2
        fi
    fi

    if [ ! -s "${OUTDIR}/contigs.fasta" ]; then
        echo "Running fresh metaSPAdes for ${sample}"
        rm -rf "${OUTDIR}"
        singularity exec -B "$PWD" "${CONTAINER_PATH}" \
            metaspades.py --only-assembler --memory "${MEM_GB}" \
            --pe1-12 "${INPUT_READS}" -o "${OUTDIR}" -t "${THREADS}"
    fi

    if [ -s "${OUTDIR}/contigs.fasta" ]; then
        echo "metaSPAdes completed for ${sample}: ${OUTDIR}/contigs.fasta"
    else
        echo "metaSPAdes did not produce contigs.fasta for ${sample}" >&2
    fi
done

echo "metaSPAdes job finished at $(date)"
