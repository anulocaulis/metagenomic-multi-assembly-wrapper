#####Container Documentation: Porechop
**Date:** February 6, 2026  
**Author:** Mike Beitner

Porechop is used to remove ONT adapter sequences from long-read FASTQ files prior to quality/length filtering and assembly. We pull the container directly from the Biocontainers project so we can pin a specific build tag and avoid rebuilding a local definition. This is reproducible because the exact tag is encoded in the image name.

The instruction to create this container is:

```bash
singularity pull docker://quay.io/biocontainers/porechop:0.2.4--py39h2de1943_9
```

If you need to build the container on your own machine (without direct pull access on the HPC), run the same pull command locally and then copy the resulting .sif to the repo’s containers directory:

```bash
# On your local machine
singularity pull docker://quay.io/biocontainers/porechop:0.2.4--py39h2de1943_9

# Copy to the HPC repo (example)
scp porechop_0.2.4--py39h2de1943_9.sif <user>@<host>:/storage/biology/projects/miller-lowry/beitner/metagenomic-multi-assembly-wrapper/containers/
```
