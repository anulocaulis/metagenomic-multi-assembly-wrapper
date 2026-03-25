#!/bin/bash
#SBATCH --job-name=metaconnet_miss
#SBATCH --output=logs/metaconnet_%j.out
#SBATCH --error=logs/metaconnet_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=900G
#SBATCH --time=6-18:00:00
#SBATCH --partition=math-alderaan-gpu

set -euo pipefail

echo "Starting MetaCONNET job at $(date)"

WORKDIR="/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper"
CONTAINER_PATH="${WORKDIR}/containers/metaconnet.sif"
THREADS=16
FLOWCELL="${METACONNET_FLOWCELL:-R10}"

cd "${WORKDIR}"

if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=(S1 S2 S5)
fi

for sample in "${SAMPLES[@]}"; do
    SHORT_READS="${WORKDIR}/trimmed_reads/${sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    LONG_READS="${WORKDIR}/data/${sample}_long_reads_filtered.fastq.gz"
    CONTIGS="${WORKDIR}/assemblies/${sample}/assembly.flye/assembly.fasta"
    OUTDIR="${WORKDIR}/assemblies/${sample}/assembly.metaconnet"
    OUTPUT="${OUTDIR}/${sample}_polished.fasta"

    echo ""
    echo "=== Running MetaCONNET for ${sample} at $(date) ==="
    echo "Short reads: ${SHORT_READS}"
    echo "Long reads:  ${LONG_READS}"
    echo "Contigs:     ${CONTIGS}"
    echo "Output dir:  ${OUTDIR}"

    if [ ! -s "${SHORT_READS}" ]; then
        echo "Missing short reads for ${sample}: ${SHORT_READS}" >&2
        continue
    fi
    if [ ! -s "${LONG_READS}" ]; then
        echo "Missing long reads for ${sample}: ${LONG_READS}" >&2
        continue
    fi
    if [ ! -s "${CONTIGS}" ]; then
        echo "Missing flye contigs for ${sample}: ${CONTIGS}" >&2
        continue
    fi

    mkdir -p "${OUTDIR}"

    SAMTOOLS_DIR="$PWD/${OUTDIR}/axbio_samtools"
    mkdir -p "$SAMTOOLS_DIR"
    cat > "$SAMTOOLS_DIR/samtools" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -gt 0 ] && [ "$1" = "merge" ]; then
    shift
    FILTERED_ARGS=""
    for arg in "$@"; do
        if [ "$arg" != "--write-index" ]; then
            FILTERED_ARGS="$FILTERED_ARGS $(printf '%q' "$arg")"
        fi
    done
    eval "set -- $FILTERED_ARGS"
    exec /opt/conda/envs/metaconda/bin/samtools merge "$@"
fi

exec /opt/conda/envs/metaconda/bin/samtools "$@"
EOF
    chmod +x "$SAMTOOLS_DIR/samtools"

    R1="${OUTDIR}/sr_R1.fastq.gz"
    R2="${OUTDIR}/sr_R2.fastq.gz"
    zcat "${SHORT_READS}" | \
        paste - - - - - - - - | \
        awk 'BEGIN{OFS="\\n"} {print $1,$2,$3,$4}' | gzip -c > "$R1"
    zcat "${SHORT_READS}" | \
        paste - - - - - - - - | \
        awk 'BEGIN{OFS="\\n"} {print $5,$6,$7,$8}' | gzip -c > "$R2"

    singularity exec \
        -B "$PWD" \
        -B "$SAMTOOLS_DIR:/AxBio_share/software/samtools-1.16/bin" \
        "${CONTAINER_PATH}" \
        /opt/conda/envs/metaconda/bin/metaconnet \
            --sr "$R1" "$R2" \
            --lr "$LONG_READS" \
            --c  "$CONTIGS" \
            --o  "$OUTDIR" \
            --n  "$sample" \
            --t  "$THREADS" \
            --fc "$FLOWCELL"

    if [ -s "$OUTPUT" ]; then
        echo "MetaCONNET completed for ${sample}: ${OUTPUT}"
    else
        echo "MetaCONNET did not produce expected output for ${sample}: ${OUTPUT}" >&2
    fi

    rm -f "$R1" "$R2"
done

echo "MetaCONNET job finished at $(date)"
