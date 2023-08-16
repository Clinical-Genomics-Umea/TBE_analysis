# TBE_analysis

## Analysis of TBE virus from nanopore sequencing

## Quick start:
### Build singularity image:
```console
sudo singularity build singularity/tbe_image.sif singularity/tbe_container.def
```
### Run pipeline
```console
singularity exec tbe_image.sif snakemake -c <CORES> -s Snakefile.py
```


## Log

### 230816
- Started the repo

- metaFlye does only produce contigs about 3000 in length.
- tried Raven, but that has a minimum contig of 10k
- canu doesnt work either...


- map the reads to the NC_001672 genome to see how the pcr tiling went 




