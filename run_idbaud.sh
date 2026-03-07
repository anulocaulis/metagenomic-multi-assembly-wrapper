#!/bin/bash
#SBATCH --job-name=idbaud
#SBATCH --output=logs/idbaud_%j.out
#SBATCH --error=logs/idbaud_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=900G
#SBATCH --time=2-00:00:00
#SBATCH --partition=math-alderaan-gpu

set -euo pipefail

echo "Starting IDBA-UD job at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/idba-ud_151.sif"
THREADS=24

cd "${WORKDIR}"

# Pass sample names as arguments, e.g.: sbatch run_idbaud.sh S1 S2 S3
# Default: all S1-S12
if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=(S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12)
fi

for sample in "${SAMPLES[@]}"; do
    INPUT_READS="${WORKDIR}/trimmed_reads/${sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    OUTDIR="${WORKDIR}/assemblies/${sample}/assembly.idbaud"
    ASSEMBLY="${OUTDIR}/assembly.fasta"

    echo ""
    echo "=== Running IDBA-UD for ${sample} at $(date) ==="
    echo "Input: ${INPUT_READS}"
    echo "Output dir: ${OUTDIR}"

    if [ ! -s "${INPUT_READS}" ]; then
        echo "ERROR: Missing input reads for ${sample}: ${INPUT_READS}" >&2
        continue
    fi

    if [ -s "${ASSEMBLY}" ]; then
        echo "Output already exists for ${sample}, skipping: ${ASSEMBLY}"
        continue
    fi

    mkdir -p "${OUTDIR}"

    READS_FQ="${OUTDIR}/reads.interleaved.fastq"
    READS_FA="${OUTDIR}/reads.interleaved.fa"

    echo "Decompressing reads for ${sample}..."
    gunzip -c "${INPUT_READS}" > "${READS_FQ}"

    echo "Converting FASTQ -> FASTA for ${sample}..."
    singularity exec -B "$PWD" "${CONTAINER_PATH}" fq2fa --paired "${READS_FQ}" "${READS_FA}"

    echo "Running IDBA-UD for ${sample}..."
    singularity exec -B "$PWD" "${CONTAINER_PATH}" \
        idba_ud -r "${READS_FA}" -o "${OUTDIR}" --num_threads "${THREADS}"

    # Clean up large intermediates to save disk
    rm -f "${READS_FQ}" "${READS_FA}"

    if [ -s "${OUTDIR}/contig.fa" ]; then
        cp "${OUTDIR}/contig.fa" "${ASSEMBLY}"
        echo "IDBA-UD completed for ${sample}: ${ASSEMBLY}"
    else
        echo "ERROR: IDBA-UD did not produce contig.fa for ${sample}" >&2
    fi
done

echo "IDBA-UD job finished at $(date)"
