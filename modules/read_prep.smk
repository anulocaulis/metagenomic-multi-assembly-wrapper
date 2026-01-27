# ---
# Step 0: Filter Long Reads (BAM to FASTQ with Chopper)
# ---
# Filters long reads using Chopper to remove low-quality and short reads.
rule filter_long_reads:
    input:
        bam = config["input_reads"]["long_bam"]
    output:
        fastq = "data/{sample}_long_reads_filtered.fastq.gz"
    params:
        temp_fastq = "data/{sample}_long_reads_temp.fastq.gz", # Temporary file path
        chopper_minlength = 200,
        chopper_quality = 7,
        container_path = QC_CONTAINER
    threads: config["threads"]
    container: QC_CONTAINER
    log: "logs/filter_long_reads_{sample}.log"
    benchmark: "benchmarks/filter_long_reads_{sample}.txt"
    shell:
        """
        set -euo pipefail;

        CONTAINER="{params.container_path}"
        THREADS="{threads}"
        TEMP_FASTQ="{params.temp_fastq}"

        # 1. Convert BAM to temporary FASTQ (required by Chopper)
        echo "Converting BAM to temporary FASTQ..."
        singularity exec -B "$(pwd)" {params.container_path} reformat.sh in={input.bam} out=$TEMP_FASTQ 2>> {log}

        # 2. Filter FASTQ using Chopper
        echo "Filtering FASTQ with Chopper (minlen={params.chopper_minlength}, minq={params.chopper_quality})..."
        singularity exec -B "$(pwd)" {params.container_path} chopper --minlength {params.chopper_minlength} \
            --quality {params.chopper_quality} \
            --threads $THREADS \
            --input "$TEMP_FASTQ" | gzip > "{output.fastq}" 2>> {log}

        # 3. Clean up
        rm -f "$TEMP_FASTQ"
        """

# ---
# Step 0.5: Convert BAM to FASTQ for MetaSPAdes (Optional long-read prep)
# ---
rule bam_to_fastq:
    """
    Converts Nanopore BAM files to gzipped FASTQ format. 
    MetaSPAdes requires FASTQ for the --nanopore flag.
    """
    input:
        "test_data/{sample}_test_long_reads.bam"
    output:
        "test_data/{sample}_test_long_reads.fastq.gz"
    container:
        QC_CONTAINER
    params:
        container_path = QC_CONTAINER
    threads: 
        lambda wildcards, resources: resources.get("threads", 4)
    shell:
        """
        # samtools fastq -@ uses the thread count for parallel decompression/conversion
        singularity exec {params.container_path} samtools fastq -@ {threads} {input} | gzip > {output}
        """

# ---
# Step 1: Trim Illumina reads after fastqc
# ---
rule trim_illumina:
    """
    Trims Illumina reads for adapters and quality using bbduk.sh from BBTools.
    Takes interleaved FASTQ as input and outputs a trimmed interleaved FASTQ.
    """
    input:
        reads=config["input_reads"]["short_interleaved"]
    output:
        trimmed_reads="trimmed_reads/{sample}_interleaved_trimmed.fastq.gz"
    container:
        QC_CONTAINER
    params:
        container_path = QC_CONTAINER,
        # - ref=adapters: Uses built-in Illumina adapter list
        # - ktrim=r k=23 mink=11 hdist=1: Standard right-end adapter trimming
        # - tpe: Trim Paired Ends
        # - tbo: Trim By Overlap
        # - qtrim=rl trimq=10: Quality trimming from left and right ends
        # - minlen=50: Discards pairs where either read is < 50 bp after trimming
        bbduk_opts="ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 minlen=50"
    threads:
        config["threads"]
    log:
        "logs/trim_illumina_{sample}.log"
    shell:
        """
        singularity exec {params.container_path} bbduk.sh in={input.reads} out={output.trimmed_reads} threads={threads} {params.bbduk_opts} 2>> {log}
        """
