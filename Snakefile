# Metagenomic Multi-Assembly Wrapper - Main Snakefile
# Author: Mike Beitner
# Date: January 28, 2026

import os

# ===========================
# Configuration
# ===========================

configfile: "config.yaml"

# Determine which samples to run
if config.get("test_mode", False):
    LONG_READ_SAMPLES = ["S5"]  # Only S5 for testing (has both long and short reads)
    ALL_SAMPLES = ["S5"]
    print("=" * 60)
    print("TEST MODE ENABLED - Running on S5 only")
    print("=" * 60)
else:
    LONG_READ_SAMPLES = config.get("long_read_samples", [])
    ALL_SAMPLES = config.get("all_samples", [])
    print("=" * 60)
    print(f"PRODUCTION MODE - Running all {len(ALL_SAMPLES)} samples")
    print(f"Long-read samples: {LONG_READ_SAMPLES}")
    print(f"Short-read only: {[s for s in ALL_SAMPLES if s not in LONG_READ_SAMPLES]}")
    print("=" * 60)

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
        # Read prep outputs (all samples get trimmed short reads)
        expand("trimmed_reads/{sample}_interleaved_trimmed.fastq.gz", sample=ALL_SAMPLES),
        
        # Long-read filtering (only long-read samples)
        expand("data/{sample}_long_reads_filtered.fastq.gz", sample=LONG_READ_SAMPLES),
        
        # Long-read assemblies (only on samples with long reads)
        expand("{output_dir}/{sample}/assembly.flye/assembly.fasta", 
               output_dir=config["output_dir"], sample=LONG_READ_SAMPLES),
        expand("{output_dir}/{sample}/assembly.metamdbg/contigs.fasta", 
               output_dir=config["output_dir"], sample=LONG_READ_SAMPLES),
        
        # Short-read assemblies (all samples)
        expand("{output_dir}/{sample}/assembly.megahit/final.contigs.fa", 
               output_dir=config["output_dir"], sample=ALL_SAMPLES),
        expand("{output_dir}/{sample}/assembly.metaspades/assembly.fasta", 
               output_dir=config["output_dir"], sample=ALL_SAMPLES),
        # expand("{output_dir}/{sample}/assembly.idbaud/assembly.fasta", 
        #        output_dir=config["output_dir"], sample=ALL_SAMPLES),
        
        # Hybrid assemblies (only for long-read samples)
        expand("{output_dir}/{sample}/assembly.metaspades_hybrid/assembly.fasta", 
               output_dir=config["output_dir"], sample=LONG_READ_SAMPLES),
        # expand("{output_dir}/{sample}/assembly.nextpolish/assembly.fasta", 
        #        output_dir=config["output_dir"], sample=LONG_READ_SAMPLES),
