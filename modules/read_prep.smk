# ---
# Step 0: Adapter trimming for long reads (Porechop)
# ---
rule adapter_trim_long_reads:
    input:
        bam = config["input_reads"]["long_bam"]
    output:
        trimmed_fastq = "data/{sample}_long_reads_porechop.fastq"
    params:
        temp_fastq = "data/{sample}_long_reads_temp.fastq.gz",
        temp_fastq_uncompressed = "data/{sample}_long_reads_temp.fastq",
        temp_fastq_porechop = "data/{sample}_long_reads_porechop.fastq",
        qc_container = QC_CONTAINER,
        porechop_container = PORECHOP_CONTAINER,
        porechop_threads = config.get("porechop", {}).get("threads", 4),
        porechop_extra_args = config.get("porechop", {}).get("extra_args", "")
    threads: config["threads"]
    container: QC_CONTAINER
    log: "logs/adapter_trim_long_reads_{sample}.log"
    benchmark: "benchmarks/adapter_trim_long_reads_{sample}.txt"
    shell:
        """
        set -euo pipefail;

        TEMP_FASTQ="{params.temp_fastq}"
        TEMP_FASTQ_UNCOMPRESSED="{params.temp_fastq_uncompressed}"
        TEMP_FASTQ_PORECHOP="{params.temp_fastq_porechop}"
        WORKDIR=$(pwd)

        # 1. Convert BAM to temporary FASTQ
        echo "Converting BAM to temporary FASTQ..." >> {log}
        singularity exec -B "$WORKDIR:$WORKDIR" -W "$WORKDIR" {params.qc_container} reformat.sh in={input.bam} out=$TEMP_FASTQ 2>> {log}

        # 2. Decompress for Porechop
        echo "Decompressing FASTQ for Porechop..." >> {log}
        gunzip -c "$TEMP_FASTQ" > "$TEMP_FASTQ_UNCOMPRESSED" 2>> {log}

        # 3. Adapter trimming with Porechop
        echo "Running Porechop adapter trimming..." >> {log}
        singularity exec -B "$WORKDIR:$WORKDIR" -W "$WORKDIR" {params.porechop_container} porechop \
            -i "$TEMP_FASTQ_UNCOMPRESSED" \
            -o "$TEMP_FASTQ_PORECHOP" \
            -t {params.porechop_threads} \
            {params.porechop_extra_args} \
            2>> {log}

        # 4. Move output into place
        mv "$TEMP_FASTQ_PORECHOP" "{output.trimmed_fastq}"

        # 5. Clean up
        rm -f "$TEMP_FASTQ" "$TEMP_FASTQ_UNCOMPRESSED"
        """

# ---
# Step 0b: Filter Long Reads (FASTQ with Chopper)
# ---
# Filters long reads using Chopper to remove low-quality and short reads.
rule filter_long_reads:
    input:
        reads = "data/{sample}_long_reads_porechop.fastq"
    output:
        fastq = "data/{sample}_long_reads_filtered.fastq.gz"
    params:
        temp_fastq_uncompressed = "data/{sample}_long_reads_porechop.fastq",
        minlength = 200,
        minquality = 7,
        chopper_container = CHOPPER_CONTAINER
    threads: config["threads"]
    container: QC_CONTAINER
    log: "logs/filter_long_reads_{sample}.log"
    benchmark: "benchmarks/filter_long_reads_{sample}.txt"
    shell:
        """
        set -euo pipefail;

        THREADS={threads}
        TEMP_FASTQ_UNCOMPRESSED="{params.temp_fastq_uncompressed}"
        WORKDIR=$(pwd)

        # 1. Decompress for chopper
        echo "Decompressing FASTQ for chopper..." >> {log}
        gunzip -c "{input.reads}" > "$TEMP_FASTQ_UNCOMPRESSED" 2>> {log}

        # 2. Filter FASTQ using Chopper
        echo "Filtering FASTQ with Chopper (minlength={params.minlength}, minquality={params.minquality})..." >> {log}
        singularity exec -B "$WORKDIR:$WORKDIR" -W "$WORKDIR" {params.chopper_container} chopper \
            --minlength {params.minlength} \
            --quality {params.minquality} \
            < "$TEMP_FASTQ_UNCOMPRESSED" | gzip > "{output.fastq}" 2>> {log}

        # 3. Clean up
        rm -f "$TEMP_FASTQ_UNCOMPRESSED"
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
        WORKDIR=$(pwd)
        singularity exec -B "$WORKDIR:$WORKDIR" -W "$WORKDIR" {params.container_path} bbduk.sh in={input.reads} out={output.trimmed_reads} threads={threads} {params.bbduk_opts} 2>> {log}
        """
