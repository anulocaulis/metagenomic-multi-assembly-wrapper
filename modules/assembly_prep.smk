import os


def _resolve_contigs_input(wildcards):
    base_dir = f"assemblies/{wildcards.sample}/assembly.{wildcards.assembler}"
    if wildcards.assembler == "metamdbg":
        native = f"{base_dir}/metamdbg.contigs.fasta"
        legacy = f"{base_dir}/contigs.fasta"
        return native if os.path.exists(native) else legacy
    return f"{base_dir}/contigs.fasta"


rule filter_contigs_min_length:
    input:
        contigs = _resolve_contigs_input
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