rule filter_contigs_min_length:
    input:
        contigs = "assemblies/{sample}/assembly.{assembler}/contigs.fasta"
    output:
        filtered = "assemblies/{sample}/assembly.{assembler}/contigs.ge1000.fa"
    params:
        container_path = QC_CONTAINER
    threads: config["threads"]
    resources:
        mem_mb = 32000,
        slurm_partition = "math-alderaan-gpu"
    container: QC_CONTAINER
    shell:
        """
        singularity exec -B $PWD {params.container_path} reformat.sh in={input.contigs} out={output.filtered} minlength=1000
        """