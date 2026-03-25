# Subsampling module for S1 long reads
# Creates 12 subsamples at 0%, 8.33%, 16.67%, ..., 91.67% of total nucleotides
# Uses BBTools reformat.sh with nucleotide-based sampling
# Operates on already-filtered S1 reads from read_prep.smk

# S1 long reads configurations
S1_LONG_READS_FILTERED = "data/S1_long_reads_filtered.fastq.gz"
SUBSAMPLE_RATIOS = list(range(12))  # [0, 1, 2, ..., 11] for 0/12, 1/12, ..., 11/12
SUBSAMPLE_SEED = 42  # Fixed seed for reproducibility
CONTAINER = "containers/qc_tools_miniconda.sif"

rule decompress_s1_fastq:
    """Decompress S1 filtered FASTQ for subsampling (reformat.sh needs uncompressed input)."""
    input:
        fastq_gz=S1_LONG_READS_FILTERED,
    output:
        fastq=temp("subsampling/S1_long_reads_filtered.fastq"),
    log:
        "logs/subsampling_decompress.log",
    shell:
        """
        zcat {input.fastq_gz} > {output.fastq} 2> {log}
        """

rule count_total_nucleotides:
    """Count total nucleotides in S1 filtered FASTQ."""
    input:
        fastq="subsampling/S1_long_reads_filtered.fastq",
    output:
        count_file="subsampling/total_nucleotides.txt",
    log:
        "logs/subsampling_count_nucleotides.log",
    shell:
        """
        awk 'NR%4==2 {{sum+=length($0)}} END {{print sum}}' {input.fastq} > {output.count_file}
        """

rule subsample_s1_single:
    """Create a single subsample of S1 filtered reads at a specified nucleotide percentage."""
    input:
        fastq="subsampling/S1_long_reads_filtered.fastq",
        count_file="subsampling/total_nucleotides.txt",
    output:
        fastq="subsampling/S1_long_reads_subsample_{ratio_idx}.fastq.gz",
    params:
        ratio_idx="{ratio_idx}",
        total_ratio=12,
        subsample_seed=SUBSAMPLE_SEED,
        container=CONTAINER,
    log:
        "logs/subsampling_subsample_{ratio_idx}.log",
    run:
        import subprocess
        
        # Read total nucleotides
        with open(input.count_file) as f:
            total_nucleotides = int(f.read().strip())
        
        # Calculate target for this subsample
        ratio_num = int(params.ratio_idx)
        target_bases = int((total_nucleotides * ratio_num) / params.total_ratio)
        
        output_file = str(output.fastq).replace('.gz', '')
        
        # Handle edge case: 0% subsample
        if target_bases == 0:
            # Create empty FASTQ file
            open(output_file, 'w').close()
        else:
            # Use reformat.sh to subsample
            cmd = (
                f"singularity exec -B $(pwd):$(pwd) {params.container} "
                f"reformat.sh "
                f"in={input.fastq} "
                f"out={output_file} "
                f"samplereadstarget={target_bases} "
                f"sampleseed={params.subsample_seed}"
            )
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        
        # Compress the output
        subprocess.run(f"gzip -f {output_file}", shell=True, check=True)
        
        # Log the subsample details
        with open(log.out, 'w') as log_f:
            log_f.write(
                f"Subsample {ratio_num}/12 (ratio={ratio_num/params.total_ratio:.4f}, "
                f"total_bases={total_nucleotides}, target_bases={target_bases})\n"
            )

# Create individual subsample rules for each ratio
for ratio_idx in SUBSAMPLE_RATIOS:
    pass  # Rules are dynamically expanded via expand() in main Snakefile

rule subsample_summary:
    """Generate summary of subsampling results."""
    input:
        expand("subsampling/S1_long_reads_subsample_{ratio_idx}.fastq.gz", ratio_idx=SUBSAMPLE_RATIOS),
        count_file="subsampling/total_nucleotides.txt",
    output:
        "subsampling/subsample_summary.txt",
    log:
        "logs/subsampling_summary.log",
    shell:
        """
        {{
            echo "S1 Filtered Reads Subsampling Summary"
            echo "======================================"
            total_bases=$$(cat {input.count_file})
            echo "Total nucleotides (filtered): ${{total_bases}}"
            echo ""
            echo "Subsample\tRatio\t\tTarget_Bases\tFile_Size"
            echo "================================================================"
            for i in {{0..11}}; do
                ratio=$$(echo "scale=4; $i / 12" | bc)
                target=$$(echo "scale=0; ${{total_bases}} * $i / 12" | bc)
                file="subsampling/S1_long_reads_subsample_${{i}}.fastq.gz"
                if [[ -f "${{file}}" ]]; then
                    size=$$(stat -c%s "${{file}}" 2>/dev/null || stat -f%z "${{file}}" 2>/dev/null)
                    printf "S1_subsample_%d\t%s\t\t%d\t%d\n" "$$i" "$$ratio" "$$target" "$$size"
                fi
            done
        }} > {output}
        cat {output}
        """
