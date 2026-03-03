# metagenomic-multi-assembly-wrapper

This repository is a Snakemake workflow for generating metagenomic assemblies from:

- Illumina short reads (all samples)
- ONT long reads (subset of samples)

It is designed as a standalone assembly module in a larger assembly comparison pipeline.

## What this repo does

### 1) Read preparation

Implemented in `modules/read_prep.smk`.

- Convert ONT BAM to FASTQ (`bam_to_fastq_long_reads`)
- Filter ONT reads with Chopper (`filter_long_reads`)
- Trim Illumina reads with BBDuk (`trim_illumina`)
- Remove polyG tails from Illumina reads (`trim_polyg`)
- Includes an optional/test BAM→FASTQ conversion rule (`bam_to_fastq`)

Outputs include:

- `data/{sample}_long_reads_filtered.fastq.gz`
- `trimmed_reads/{sample}_interleaved_trimmed.fastq.gz`
- `trimmed_reads/{sample}_interleaved_trimmed_polyG_filtered.fastq.gz`

### 2) Assembly workflows

Implemented in `modules/assembly.smk`.

- **Flye** long-read assembly (`flye_assembly`)
- **metaMDBG / nanoMDBG** long-read assembly (`metamdbg_assembly`)
- **IDBA-UD** short-read assembly (`idbaud_assembly`)
- **metaSPAdes** short-read assembly (`metaspades_assembly`)
- **MEGAHIT** short-read assembly (`megahit_assembly`)
- **Hybrid metaSPAdes** short+long assembly (`metaspades_hybrid_assembly`)
- **Flye + NextPolish** polishing workflow (`bwa_mem_map_to_flye`, `nextpolish`)

Primary outputs are written under `assemblies/{sample}/...`.

### 3) Containerized execution

Container variables are defined in `Snakefile`.

- `ASSEMBLY_CONTAINER = containers/assembler2.sif`
- `IDBAUD_CONTAINER = containers/idba-ud_151.sif`
- `FLYE_ASSEMBLY_CONTAINER = containers/flye_assembler.sif`
- `QC_CONTAINER = containers/qc_tools_miniconda.sif`
- `CHOPPER_CONTAINER = containers/chopper_0.7.0--hdcf5f25_0.sif`
- `PORECHOP_CONTAINER = containers/porechop_0.2.4--py39h2de1943_9.sif`

`idbaud_assembly` is configured to run with `containers/idba-ud_151.sif`.

## Repository layout

- `Snakefile`: workflow entrypoint, sample selection, global container paths, target rule
- `config.yaml`: sample lists, input path templates, thread defaults, output dir
- `modules/read_prep.smk`: preprocessing rules
- `modules/assembly.smk`: assembly and polishing rules
- `containers/`: Singularity images and definition files
- `logs/`, `benchmarks/`: rule logs and runtime benchmarks
- `assemblies/`: per-sample assembly outputs

## Configuration

Edit `config.yaml` for:

- `test_mode` (`true` runs S5-only)
- `long_read_samples`
- `all_samples`
- `input_reads.long_bam`
- `input_reads.short_interleaved`
- `output_dir`
- `threads`

## Running the workflow

Make sure your conda env is activated first (so `snakemake` is on PATH), then run from repo root.

### Dry run

```bash
snakemake -n
```

Dry-run one target example:

```bash
snakemake -n assemblies/S1/assembly.idbaud/assembly.fasta
```

### Execute

```bash
snakemake --cores 16
```

Use your cluster profile when needed:

```bash
snakemake --profile slurm_profile
```

## Notes

- Rules currently request high-memory resources for assembly jobs.
- Some targets in `rule all` may be toggled/commented during development.
- Existing files in `assemblies/` can satisfy targets and cause jobs to be skipped.
