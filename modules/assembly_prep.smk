rule filter_contigs_min_length:
    input:
        contigs="assemblies/{sample}/assembly.{assembler}/contigs.fasta"
    output:
        filtered="assemblies/{sample}/assembly.{assembler}/contigs.ge1000.fa"
    container:
        "qc_tools_miniconda.sif"
    shell:
        """
        bbduk.sh in={input.contigs} out={output.filtered} minlength=1000
        """