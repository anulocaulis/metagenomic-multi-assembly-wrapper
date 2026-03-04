# ---
# Step C1: Taxonomic classification with Kraken2
# ---
rule kraken2_classify:
    """
    Runs Kraken2 to assign taxonomic labels to trimmed short reads.
    Produces a per-read output file and a summary report.
    Requires a Kraken2 database specified in config["kraken2_db"].
    """
    input:
        reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        output_file = "classification/{sample}/kraken2/output.kraken",
        report      = "classification/{sample}/kraken2/report.txt"
    params:
        outdir         = "classification/{sample}/kraken2",
        db             = config.get("kraken2_db", "/path/to/kraken2_db"),
        container_path = CLASSIFICATION_CONTAINER
    threads: config["threads"]
    resources:
        mem_mb = 360000
    log: "logs/classification_kraken2_{sample}.log"
    shell:
        """
        set -euo pipefail
        WORKDIR=$(pwd)

        mkdir -p {params.outdir}

        echo "Running Kraken2 classification for {wildcards.sample}..." >> {log}

        singularity exec -B "$WORKDIR:$WORKDIR" {params.container_path} \
            kraken2 \
                --db {params.db} \
                --threads {threads} \
                --output {output.output_file} \
                --report {output.report} \
                --gzip-compressed \
                --paired --interleaved \
                {input.reads} 2>> {log}
        """


# ---
# Step C2: Abundance re-estimation with Bracken
# ---
rule bracken_abundance:
    """
    Runs Bracken to re-estimate species-level abundances from a Kraken2 report.
    Requires the same Kraken2 database used in kraken2_classify.
    """
    input:
        report = "classification/{sample}/kraken2/report.txt"
    output:
        bracken_output = "classification/{sample}/bracken/{sample}_bracken_species.txt",
        bracken_report = "classification/{sample}/bracken/{sample}_bracken_species.report"
    params:
        outdir         = "classification/{sample}/bracken",
        db             = config.get("kraken2_db", "/path/to/kraken2_db"),
        read_length    = config.get("bracken_read_length", 150),
        level          = config.get("bracken_level", "S"),
        threshold      = config.get("bracken_threshold", 10),
        container_path = CLASSIFICATION_CONTAINER
    threads: config["threads"]
    resources:
        mem_mb = 64000
    log: "logs/classification_bracken_{sample}.log"
    shell:
        """
        set -euo pipefail
        WORKDIR=$(pwd)

        mkdir -p {params.outdir}

        echo "Running Bracken abundance estimation for {wildcards.sample}..." >> {log}

        singularity exec -B "$WORKDIR:$WORKDIR" {params.container_path} \
            bracken \
                -d {params.db} \
                -i {input.report} \
                -o {output.bracken_output} \
                -w {output.bracken_report} \
                -r {params.read_length} \
                -l {params.level} \
                -t {params.threshold} 2>> {log}
        """
