# ---
# Step 0: Convert BAM to filtered FASTQ for long reads
# ---
rule bam_to_fastq_long_reads:
    input:
        bam = config["input_reads"]["long_bam"]
    output:
        fastq = "data/{sample}_long_reads_filtered.fastq.gz"
    params:
        qc_container = QC_CONTAINER,
        chopper_container = CHOPPER_CONTAINER,
        minlength = 200,
        minquality = 7
    threads: config["threads"]
    resources:
        mem_mb=360000
    container: QC_CONTAINER
    log: "logs/bam_to_fastq_long_reads_{sample}.log"
    benchmark: "benchmarks/bam_to_fastq_long_reads_{sample}.txt"
    shell:
        """
        WORKDIR=$(pwd)
        TMP_FASTQ=data/{wildcards.sample}_long_reads_uncompressed.fastq

        echo "Converting BAM to FASTQ..." >> {log}
        singularity exec -B "$WORKDIR:$WORKDIR" {params.qc_container} \
            reformat.sh in={input.bam} out=$TMP_FASTQ >> {log} 2>&1

        echo "Filtering FASTQ with Chopper (minlength={params.minlength}, minquality={params.minquality})..." >> {log}
        singularity exec -B "$WORKDIR:$WORKDIR" {params.chopper_container} chopper \
            --minlength {params.minlength} \
            --quality {params.minquality} \
            < "$TMP_FASTQ" | gzip > "{output.fastq}" 2>> {log}

        rm -f "$TMP_FASTQ"
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

# ---
# Step 2: Remove polyG tails from Illumina reads
# ---
rule trim_polyg:
    """
    Removes polyG tails from Illumina reads using BBTools polyfilter.sh
    """
    input:
        reads = "trimmed_reads/{sample}_interleaved_trimmed.fastq.gz"
    output:
        filtered_reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    container:
        QC_CONTAINER
    params:
        container_path = QC_CONTAINER
    threads:
        config["threads"]
    log:
        "logs/trim_polyg_{sample}.log"
    shell:
        """
        WORKDIR=$(pwd)
        singularity exec -B "$WORKDIR:$WORKDIR" {params.container_path} polyfilter.sh in={input.reads} out={output.filtered_reads} 2>> {log}
        """
