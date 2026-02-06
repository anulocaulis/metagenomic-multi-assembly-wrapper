#--
#Step 1: Long-Read Assembly with Flye (ONT-only) - STILL COMMENTED OUT
#---
rule flye_assembly:
    input:
        reads = "data/{sample}_long_reads_filtered.fastq.gz"
    output:
        assembly = f"{config['output_dir']}/{{sample}}/assembly.flye/assembly.fasta"
    params:
        outdir = f"{config['output_dir']}/{{sample}}/assembly.flye",
        container_path = FLYE_ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: FLYE_ASSEMBLY_CONTAINER
    shell:
        "singularity exec -B $PWD {params.container_path} flye --nano-raw {input.reads} --out-dir {params.outdir} --threads {threads}"


# ---
# Step 2: Long-Read Assembly with NanoMDBG (ONT-only)
# ---
rule metamdbg_assembly:
    input:
        reads = "data/{sample}_long_reads_filtered.fastq.gz"
    output:
        assembly = f"{config['output_dir']}/{{sample}}/assembly.metamdbg/contigs.fasta"
    params:
        outdir = f"{config['output_dir']}/{{sample}}/assembly.metamdbg",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: ASSEMBLY_CONTAINER
    shell:
        """
        singularity exec --bind $(pwd):$(pwd) {params.container_path} metaMDBG asm --threads {threads} --in-ont {input.reads} --out-dir {params.outdir}
        # MetaMDBG outputs contigs.fasta.gz, decompress for consistency
        gunzip -f {params.outdir}/contigs.fasta.gz
        """


# ---
# Step 2b: Run metaMDBG assembly directly
# ---
rule metamdbg_assembly_direct:
    input:
        reads = f"{config['output_dir']}/{{sample}}/assembly.metamdbg/contigs.fasta"
    output:
        assembly = f"{config['output_dir']}/{{sample}}/assembly.metamdbg/contigs.fasta"
    params:
        outdir = f"{config['output_dir']}/{{sample}}/assembly.metamdbg",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: ASSEMBLY_CONTAINER
    shell:
        """
        singularity exec --bind $(pwd):$(pwd) {params.container_path} metaMDBG asm --out-dir {params.outdir} --in-ont {input.reads} --threads {threads} && \
        gzip {params.outdir}/contigs.fasta
        """


# ---
# Step 3: Short-Read Assembly with IDBA-UD (Illumina-only)
# ---
rule idbaud_assembly:
    input:
        reads = config["input_reads"]["short_interleaved"]
    output:
        assembly = "{output_dir}/{sample}/assembly.idbaud/assembly.fasta"
    params:
        outdir = "{output_dir}/{sample}/assembly.idbaud",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: ASSEMBLY_CONTAINER
    shell:
        "singularity exec -B $PWD {params.container_path} idba_ud -r {input.reads} -o {params.outdir} --num_threads {threads}"


# ---
# Step 4: Short-Read Assembly with MetaSPAdes (Illumina-only)
# ---
rule metaspades_assembly:
    input:
        reads = config["input_reads"]["short_interleaved"]
    output:
        assembly = "{output_dir}/{sample}/assembly.metaspades/assembly.fasta"
    params:
        outdir = "{output_dir}/{sample}/assembly.metaspades",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: ASSEMBLY_CONTAINER
    shell:
        "singularity exec -B $PWD {params.container_path} metaspades.py --only-assembler --pe1-12 {input.reads} -o {params.outdir} -t {threads}"

# ---
# Step 5: Short-Read Assembly with MEGAHIT (Illumina-only)
# ---
rule megahit_assembly:
    input:
        reads = config["input_reads"]["short_interleaved"]
    output:
        assembly = "{output_dir}/{sample}/assembly.megahit/final.contigs.fa"
    params:
        outdir = "{output_dir}/{sample}/assembly.megahit",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: ASSEMBLY_CONTAINER
    shell:
        "rm -rf {params.outdir} && singularity exec -B $PWD {params.container_path} megahit --12 {input.reads} -o {params.outdir} --num-cpu-threads {threads}"


# ---
# Step 6: Hybrid Assembly with MetaSPAdes (Short-Read First)
# ---
rule metaspades_hybrid_assembly:
    input:
        short_reads = config["input_reads"]["short_interleaved"],
        long_reads = "data/{sample}_long_reads_filtered.fastq.gz"
    output:
        assembly = "{output_dir}/{sample}/assembly.metaspades_hybrid/assembly.fasta",
        spades_contigs = "{output_dir}/{sample}/assembly.metaspades_hybrid/contigs.fasta"
    params:
        outdir = "{output_dir}/{sample}/assembly.metaspades_hybrid",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    container: ASSEMBLY_CONTAINER
    shell:
        """
        # Run metaspades, interleaved reads
        singularity exec -B $PWD {params.container_path} metaspades.py --only-assembler \
            --pe1-12 {input.short_reads} \
            --nanopore {input.long_reads} \
            -o {params.outdir} -t {threads}
        # Copy SPAdes output (contigs.fasta) to the generic name expected by Snakemake (assembly.fasta)
        cp {params.outdir}/contigs.fasta {output.assembly}
        """

# ---
# Step 7: Map Short Reads to Flye Assembly with BWA MEM
# ---
rule bwa_mem_map_to_flye:
    input:
        assembly = f"{config['output_dir']}/{{sample}}/assembly.flye/assembly.fasta",
        short_reads = config["input_reads"]["short_interleaved"]
    output:
        bam = f"{config['output_dir']}/{{sample}}/assembly.flye/mapped_short_reads.bam",
        bai = f"{config['output_dir']}/{{sample}}/assembly.flye/mapped_short_reads.bam.bai"
    params:
        container_path = QC_CONTAINER,
        outdir = f"{config['output_dir']}/{{sample}}/assembly.flye"
    threads: config["threads"]
    container: QC_CONTAINER
    shell:
        """
        # 1. Index the Flye assembly with BWA
        singularity exec -B $PWD {params.container_path} bwa index {input.assembly}
        
        # 2. Map short reads to the assembly using BWA MEM
        # Output SAM, pipe to samtools to convert to sorted BAM
        singularity exec -B $PWD {params.container_path} bash -c \
            "bwa mem -t {threads} {input.assembly} {input.short_reads} | samtools sort -@ {threads} -o {output.bam} -"
        
        # 3. Index the BAM file
        singularity exec -B $PWD {params.container_path} samtools index {output.bam}
        """

# ---
# Step 8: Hybrid Assembly with Flye + NextPolish (Long-Read First)
# ---
rule nextpolish:
    input:
        long_assembly = f"{config['output_dir']}/{{sample}}/assembly.flye/assembly.fasta",
        mapped_short_reads_bam = f"{config['output_dir']}/{{sample}}/assembly.flye/mapped_short_reads.bam",
        short_reads = config["input_reads"]["short_interleaved"],
        long_reads = config["input_reads"]["long_bam"]

    output:
        assembly = f"{config['output_dir']}/{{sample}}/assembly.nextpolish/assembly.fasta",
        sgs_fofn = f"{config['output_dir']}/{{sample}}/assembly.nextpolish/sgs.fofn",
        lgs_fofn = f"{config['output_dir']}/{{sample}}/assembly.nextpolish/lgs.fofn",
        run_cfg = f"{config['output_dir']}/{{sample}}/assembly.nextpolish/run.cfg"

    params:
        outdir = f"{config['output_dir']}/{{sample}}/assembly.nextpolish",
        workdir = f"{config['output_dir']}/{{sample}}/assembly.nextpolish/01_rundir",
        container_path = "nextpolish.sif"
    threads: config["threads"]
    container: "nextpolish.sif"
    shell:
        """
        THREADS={threads}
        OUTDIR=$(readlink -f {params.outdir})
        WORKDIR=$(readlink -f {params.workdir})

        echo "Starting NextPolish preparation and assembly..."
        mkdir -p $OUTDIR $WORKDIR

        # 1. Create Short-Read FOFN (sgs.fofn)
        echo "{input.short_reads}" | tr ' ' '\n' > $OUTDIR/sgs.fofn

        # 2. Create Long-Read FOFN (lgs.fofn)
        echo "{input.long_reads}" | tr ' ' '\n' > $OUTDIR/lgs.fofn

        # 3. Create run.cfg file
        cat > $OUTDIR/run.cfg << EOL
        [General]
        job_type = local
        job_prefix = nextPolish
        task = best
        rewrite = yes
        rerun = 3
        parallel_jobs = {threads}
        multithread_jobs = {threads}
        genome = $(readlink -f {input.long_assembly})
        genome_size = auto
        workdir = $WORKDIR
        polish_options = -p {threads}

        [sgs_option]
        sgs_fofn = $OUTDIR/sgs.fofn
        sgs_options = -max_depth 100 -bwa

        [lgs_option]
        lgs_fofn = $OUTDIR/lgs.fofn
        lgs_options = -min_read_len 1k -max_depth 100
        lgs_minimap2_options = -x map-ont
        EOL

        # 4. Execute NextPolish
        singularity exec --bind $(pwd):$(pwd) {params.container_path} /opt/conda/bin/nextPolish $OUTDIR/run.cfg 2>&1 | tee -a logs/nextpolish.log
        """
