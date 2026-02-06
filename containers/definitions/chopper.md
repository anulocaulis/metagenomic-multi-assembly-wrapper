#####Container Documentation: Chopper
**Date:** February 1, 2026  
**Author:** Mike Beitner

To make a container for Chopper, we pull directly from the Biocontainers project. The "recipe" for this container effectively exists upstream in the Bioconda Github Repository. This is as reproducible as rebuilding the sifs used in this workflow from the provided .def files, as the version 0.7.0 is specified in the image tag below below.


The instruction to create this container is:


```bash
singularity pull docker://quay.io/biocontainers/chopper:0.7.0--hdcf5f25_0
```