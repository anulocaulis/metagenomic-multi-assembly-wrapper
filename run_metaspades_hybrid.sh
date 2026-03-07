#!/bin/bash
#SBATCH --job-name=metaspades_hybrid
#SBATCH --output=logs/metaspades_hybrid_%j.out
#SBATCH --error=logs/metaspades_hybrid_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=900G
#SBATCH --time=2-00:00:00
#SBATCH --partition=math-alderaan-gpu

set -euo pipefail

echo "Starting metaSPAdes hybrid job at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/assembler2.sif"
THREADS=24
MEM_GB=900

cd "${WORKDIR}"

# Default: only samples with long reads that are not yet done
# Override by passing sample names as arguments: sbatch run_metaspades_hybrid.sh S1 S2
if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=(S1 S2)
fi

for sample in "${SAMPLES[@]}"; do
    SHORT_READS="${WORKDIR}/trimmed_reads/${sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    LONG_READS="${WORKDIR}/data/${sample}_long_reads_filtered.fastq.gz"
    OUTDIR="${WORKDIR}/assemblies/${sample}/assembly.metaspades_hybrid"

    echo ""
    echo "=== Running metaSPAdes hybrid for ${sample} at $(date) ==="
    echo "Short reads: ${SHORT_READS}"
    echo "Long reads:  ${LONG_READS}"
    echo "Output dir:  ${OUTDIR}"

    if [ ! -s "${SHORT_READS}" ]; then
        echo "ERROR: Missing short reads for ${sample}: ${SHORT_READS}" >&2
        continue
    fi

    if [ ! -s "${LONG_READS}" ]; then
        echo "ERROR: Missing long reads for ${sample}: ${LONG_READS}" >&2
        continue
    fi

    mkdir -p "${OUTDIR}"

    # Resume if a partial run exists
    if [ -f "${OUTDIR}/run_spades.yaml" ] || [ -d "${OUTDIR}/pipeline_state" ]; then
        echo "Resume checkpoint detected for ${sample}; trying --continue"
        if singularity exec -B "$PWD" "${CONTAINER_PATH}" metaspades.py --continue -o "${OUTDIR}"; then
            echo "--continue succeeded for ${sample}"
        else
            echo "--continue failed for ${sample}; forcing clean rerun" >&2
            rm -rf "${OUTDIR}"
            mkdir -p "${OUTDIR}"
        fi
    fi

    # Fresh run (or after failed resume)
    if [ ! -s "${OUTDIR}/contigs.fasta" ]; then
        echo "Running fresh metaSPAdes hybrid for ${sample}"
        singularity exec -B "$PWD" "${CONTAINER_PATH}" \
            metaspades.py --only-assembler --memory "${MEM_GB}" \
                --pe1-12 "${SHORT_READS}" \
                --nanopore "${LONG_READS}" \
                -o "${OUTDIR}" -t "${THREADS}"
    fi

    # Copy to the generic assembly.fasta expected by downstream rules
    if [ -s "${OUTDIR}/contigs.fasta" ]; then
        cp "${OUTDIR}/contigs.fasta" "${OUTDIR}/assembly.fasta"
        echo "metaSPAdes hybrid completed for ${sample}: ${OUTDIR}/assembly.fasta"
    else
        echo "ERROR: metaSPAdes hybrid did not produce contigs.fasta for ${sample}" >&2
    fi
done

echo "metaSPAdes hybrid job finished at $(date)"
