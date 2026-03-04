# ---
# Helper function: resolve the assembled FASTA file for a given tool and sample
# ---
def get_assembly_for_binning(wildcards):
    """
    Returns the path to the assembled FASTA file for a given assembly tool and sample.
    Handles tools that use non-standard output names (megahit, metamdbg).
    """
    tool = wildcards.tool
    sample = wildcards.sample
    base = config["output_dir"]

    if tool == "megahit":
        return f"{base}/{sample}/assembly.megahit/final.contigs.fa"
    if tool == "metamdbg":
        return f"{base}/{sample}/assembly.metamdbg/contigs.fasta"
    return f"{base}/{sample}/assembly.{tool}/assembly.fasta"


# ---
# Step B1: Bin contigs with MetaWRAP (MaxBin2 + MetaBAT2 + CONCOCT)
# ---
rule metawrap_binning:
    """
    Runs MetaWRAP's binning module to produce initial bins from an assembly
    using MaxBin2, MetaBAT2, and CONCOCT in parallel.
    """
    input:
        assembly = get_assembly_for_binning,
        reads    = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        metabat2_done = "binning/{tool}_{sample}/metabat2_bins/done",
        maxbin2_done  = "binning/{tool}_{sample}/maxbin2_bins/done",
        concoct_done  = "binning/{tool}_{sample}/concoct_bins/done"
    params:
        outdir           = "binning/{tool}_{sample}",
        container_path   = METAWRAP_CONTAINER
    threads: config["threads"]
    resources:
        mem_mb = 360000
    log: "logs/binning_metawrap_{tool}_{sample}.log"
    shell:
        """
        set -euo pipefail
        WORKDIR=$(pwd)

        mkdir -p {params.outdir}

        echo "Running MetaWRAP binning for {wildcards.tool} on {wildcards.sample}..." >> {log}

        singularity exec -B "$WORKDIR:$WORKDIR" {params.container_path} \
            metawrap binning \
                -o {params.outdir} \
                -t {threads} \
                -a {input.assembly} \
                --metabat2 --maxbin2 --concoct \
                {input.reads} 2>> {log}

        touch {output.metabat2_done} {output.maxbin2_done} {output.concoct_done}
        """


# ---
# Step B2: Refine bins with MetaWRAP bin_refinement (CheckM-based)
# ---
rule metawrap_bin_refinement:
    """
    Runs MetaWRAP's bin_refinement module to produce a consolidated, high-quality
    bin set from the three initial binners using CheckM for quality assessment.
    Targets bins with >= 50% completeness and <= 10% contamination.
    """
    input:
        metabat2_done = "binning/{tool}_{sample}/metabat2_bins/done",
        maxbin2_done  = "binning/{tool}_{sample}/maxbin2_bins/done",
        concoct_done  = "binning/{tool}_{sample}/concoct_bins/done"
    output:
        stats = "binning/{tool}_{sample}/bin_refinement/metawrap_50_10_bins.stats"
    params:
        outdir         = "binning/{tool}_{sample}/bin_refinement",
        metabat2_dir   = "binning/{tool}_{sample}/metabat2_bins",
        maxbin2_dir    = "binning/{tool}_{sample}/maxbin2_bins",
        concoct_dir    = "binning/{tool}_{sample}/concoct_bins",
        container_path = METAWRAP_CONTAINER,
        completeness   = 50,
        contamination  = 10
    threads: config["threads"]
    resources:
        mem_mb = 360000
    log: "logs/binning_refinement_{tool}_{sample}.log"
    shell:
        """
        set -euo pipefail
        WORKDIR=$(pwd)

        mkdir -p {params.outdir}

        echo "Running MetaWRAP bin_refinement for {wildcards.tool} on {wildcards.sample}..." >> {log}

        singularity exec -B "$WORKDIR:$WORKDIR" {params.container_path} \
            metawrap bin_refinement \
                -o {params.outdir} \
                -t {threads} \
                -A {params.metabat2_dir} \
                -B {params.maxbin2_dir} \
                -C {params.concoct_dir} \
                -c {params.completeness} \
                -x {params.contamination} 2>> {log}
        """
