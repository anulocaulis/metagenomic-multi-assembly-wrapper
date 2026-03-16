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
    resources:
        mem_mb=512000 # Request entire node's memory for assembly rules
    container: FLYE_ASSEMBLY_CONTAINER
    shell:
        """
        if [ -f {params.outdir}/params.json ]; then
            singularity exec -B $PWD {params.container_path} flye --resume --out-dir {params.outdir} --threads {threads} --meta
        else
            singularity exec -B $PWD {params.container_path} flye --nano-raw {input.reads} --out-dir {params.outdir} --threads {threads} --meta
        fi
        """


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
    resources:
        mem_mb=512000 # Request entire node's memory for assembly rules
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
    resources:
        mem_mb=512000 # Request entire node's memory for assembly rules
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
        reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        assembly = "{output_dir}/{sample}/assembly.idbaud/assembly.fasta"
    params:
        outdir = "{output_dir}/{sample}/assembly.idbaud",
        container_path = IDBAUD_CONTAINER
    threads: config["threads"]
    resources:
        mem_mb=515000 # Request full node memory (515GB available on math-alderaan)
    container: IDBAUD_CONTAINER
    shell:
        """
        mkdir -p {params.outdir}
        READS_FQ={params.outdir}/reads.interleaved.fastq
        READS_FA={params.outdir}/reads.interleaved.fa

        gunzip -c {input.reads} > $READS_FQ
        singularity exec -B $PWD {params.container_path} fq2fa --paired $READS_FQ $READS_FA
        singularity exec -B $PWD {params.container_path} idba_ud -r $READS_FA -o {params.outdir} --num_threads {threads}

        cp {params.outdir}/contig.fa {output.assembly}
        """


# ---
# Step 4: Short-Read Assembly with MetaSPAdes (Illumina-only)
# ---
rule metaspades_assembly:
    input:
        reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        assembly = "{output_dir}/{sample}/assembly.metaspades/contigs.fasta"
    params:
        outdir = "{output_dir}/{sample}/assembly.metaspades",
        container_path = ASSEMBLY_CONTAINER
    threads: 24
    resources:
        mem_mb=900000,
        slurm_partition="math-alderaan-gpu"
    container: ASSEMBLY_CONTAINER
    shell:
        """
        if [ -f {params.outdir}/run_spades.yaml ] || [ -d {params.outdir}/pipeline_state ]; then
            singularity exec -B $PWD {params.container_path} metaspades.py --continue -o {params.outdir}
            if [ ! -s {params.outdir}/contigs.fasta ]; then
                rm -rf {params.outdir}
                singularity exec -B $PWD {params.container_path} metaspades.py --only-assembler --memory 900 --pe1-12 {input.reads} -o {params.outdir} -t {threads}
            fi
        else
            singularity exec -B $PWD {params.container_path} metaspades.py --only-assembler --memory 900 --pe1-12 {input.reads} -o {params.outdir} -t {threads}
        fi
        """

# ---
# Step 5: Short-Read Assembly with MEGAHIT (Illumina-only)
# ---
rule megahit_assembly:
    input:
        reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz"
    output:
        assembly = "{output_dir}/{sample}/assembly.megahit/final.contigs.fa"
    params:
        outdir = "{output_dir}/{sample}/assembly.megahit",
        container_path = ASSEMBLY_CONTAINER
    threads: config["threads"]
    resources:
        mem_mb=512000 # Request entire node's memory for assembly rules
    container: ASSEMBLY_CONTAINER
    shell:
        """
        if [ -f {params.outdir}/tmp/reads.lib.lib_info ]; then
            singularity exec -B $PWD {params.container_path} megahit --continue -o {params.outdir} --num-cpu-threads {threads}
        else
            rm -rf {params.outdir}
            singularity exec -B $PWD {params.container_path} megahit --12 {input.reads} -o {params.outdir} --num-cpu-threads {threads}
        fi
        """


# ---
# Step 6: Hybrid Assembly with MetaSPAdes (Short-Read First)
# ---
rule metaspades_hybrid_assembly:
    input:
        short_reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz",
        long_reads = "data/{sample}_long_reads_filtered.fastq.gz"
    output:
        assembly = "{output_dir}/{sample}/assembly.metaspades_hybrid/assembly.fasta",
        spades_contigs = "{output_dir}/{sample}/assembly.metaspades_hybrid/contigs.fasta"
    params:
        outdir = "{output_dir}/{sample}/assembly.metaspades_hybrid",
        container_path = ASSEMBLY_CONTAINER
    threads: 24
    resources:
        mem_mb=900000,
        slurm_partition="math-alderaan-gpu"
    container: ASSEMBLY_CONTAINER
    shell:
        """
        if [ -f {params.outdir}/run_spades.yaml ] || [ -d {params.outdir}/pipeline_state ]; then
            singularity exec -B $PWD {params.container_path} metaspades.py --continue -o {params.outdir}
            if [ ! -s {params.outdir}/contigs.fasta ]; then
                rm -rf {params.outdir}
                singularity exec -B $PWD {params.container_path} metaspades.py --only-assembler --memory 900 \
                    --pe1-12 {input.short_reads} \
                    --nanopore {input.long_reads} \
                    -o {params.outdir} -t {threads}
            fi
        else
            # Run metaspades, interleaved reads
            singularity exec -B $PWD {params.container_path} metaspades.py --only-assembler --memory 900 \
                --pe1-12 {input.short_reads} \
                --nanopore {input.long_reads} \
                -o {params.outdir} -t {threads}
        fi
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
    resources:
        mem_mb=512000 # Request entire node's memory for assembly rules
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
# Step 8: Long-Read Polishing / Hybrid Scaffolding with MetaCONNET
# ---
rule metaconnet:
    input:
        short_reads = "trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz",
        long_reads = "data/{sample}_long_reads_filtered.fastq.gz",
        contigs = f"{config['output_dir']}/{{sample}}/assembly.flye/assembly.fasta"
    output:
        assembly = f"{config['output_dir']}/{{sample}}/assembly.metaconnet/assembly.fasta"
    params:
        outdir = f"{config['output_dir']}/{{sample}}/assembly.metaconnet",
        container_path = METACONNET_CONTAINER,
        flowcell = config.get("metaconnet_flowcell", "R10")
    threads: config["threads"]
    resources:
        mem_mb=900000,
        slurm_partition="math-alderaan-gpu"
    container: METACONNET_CONTAINER
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        # Split interleaved paired-end reads into R1 / R2 for MetaCONNET --sr
        R1={params.outdir}/sr_R1.fastq.gz
        R2={params.outdir}/sr_R2.fastq.gz
        zcat {input.short_reads} | \
            paste - - - - - - - - | \
            awk 'BEGIN{{OFS="\\n"}} {{print $1,$2,$3,$4}}' | gzip -c > $R1
        zcat {input.short_reads} | \
            paste - - - - - - - - | \
            awk 'BEGIN{{OFS="\\n"}} {{print $5,$6,$7,$8}}' | gzip -c > $R2

        singularity exec -B $PWD {params.container_path} \
            MetaCONNET \
                --sr $R1 $R2 \
                --lr {input.long_reads} \
                --c  {input.contigs} \
                --o  {params.outdir} \
                --n  {wildcards.sample} \
                --t  {threads} \
                --fc {params.flowcell}

        # MetaCONNET names its primary output <task_name>.fasta
        cp {params.outdir}/{wildcards.sample}.fasta {output.assembly}

        # Clean up temporary split reads
        rm -f $R1 $R2
        """
