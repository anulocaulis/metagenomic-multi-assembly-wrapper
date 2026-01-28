# Metagenomic Multi-Assembly Wrapper - Main Snakefile
# Author: Mike Beitner
# Date: January 28, 2026

import os

# ===========================
# Configuration
# ===========================

configfile: "config.yaml"

# Container Paths
ASSEMBLY_CONTAINER = "containers/assembler2.sif"
FLYE_ASSEMBLY_CONTAINER = "containers/flye_assembler.sif"
QC_CONTAINER = "containers/qc_tools_miniconda.sif"

# ===========================
# Include Module Rules
# ===========================

include: "modules/read_prep.smk"
include: "modules/assembly.smk"

# ===========================
# Target Rule
# ===========================

rule all:
    input:
        # Read prep outputs
        expand("data/{sample}_long_reads_filtered.fastq.gz", sample=config["samples"]),
        expand("trimmed_reads/{sample}_interleaved_trimmed.fastq.gz", sample=config["samples"]),
        
        # Assembly outputs - choose which assemblies you want to run
        # Uncomment the ones you want to produce:
        
        # Long-read assemblies
        expand("{output_dir}/{sample}/assembly.flye/assembly.fasta", 
               output_dir=config["output_dir"], sample=config["samples"]),
        expand("{output_dir}/{sample}/assembly.metamdbg/contigs.fasta", 
               output_dir=config["output_dir"], sample=config["samples"]),
        
        # Short-read assemblies
        expand("{output_dir}/{sample}/assembly.megahit/final.contigs.fa", 
               output_dir=config["output_dir"], sample=config["samples"]),
        expand("{output_dir}/{sample}/assembly.metaspades/assembly.fasta", 
               output_dir=config["output_dir"], sample=config["samples"]),
        # expand("{output_dir}/{sample}/assembly.idbaud/assembly.fasta", 
        #        output_dir=config["output_dir"], sample=config["samples"]),
        
        # Hybrid assemblies
        expand("{output_dir}/{sample}/assembly.metaspades_hybrid/assembly.fasta", 
               output_dir=config["output_dir"], sample=config["samples"]),
        # expand("{output_dir}/{sample}/assembly.nextpolish/assembly.fasta", 
        #        output_dir=config["output_dir"], sample=config["samples"]),
