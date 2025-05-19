<p align="center">
  <img src="WEPP_logo.svg" width="300">
</p>

<h1 align="center">
  Wastewater-based Epidemiology using Phylogenetic Placements
</h1>


## Installation (Follow steps 2-3 if you want to use docker)
1. `git clone --recurse-submodules https://github.com/TurakhiaLab/SARS2-WBE.git`. Switch the branch if needed and check if the `src/Freyja` is not empty. If it is empty then go inside `src/Freyja` and use
```
git pull --recurse 
```
2. Create a docker image by going to the `docker` and running,
```
docker build -t {image_name} .
```
3. Return to the main folder and run 
```
docker run -it -v "$PWD":/workspace -w /workspace {image_name} /bin/bash
```

## Running
We assume that all the different wastewater samples are stored in the `data` folder under different `DIR` names. Each wastewater `DIR` should have the following files:
1. Sequencing Reads: Ending with `*R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Reference Genome (REF)
3. Mutation Annotated Tree (TREE)
4. [OPTIONAL] Genome Masking File: `mask.bed` whose third column specifies sites to be excluded from analysis.

WEPP's snakemake pipeline requires the following parameters passed either through the config file (`config/config.yaml`) or command line using the `--config` argument. The command line arguments take precedence over the config file.
1. `DIR` - Folder with the samples 
2. `FILE_PREFIX` - Prefix for all intermediate files 
3. `REF` - Reference Genome
4. `TREE` - MAT
5. `SEQUENCING_TYPE` - Sequencing read type (Illumina single-ended, Illumina double-ended, ONT long reads)
6. `PRIMER_BED` - BED file for primers in the `database` folder
7. `MIN_AF` - Alleles in reads below this Allele Frequency threshold get masked 
8. `MIN_Q` - Alleles below this Phred score get masked
9. `MAX_READS` - Maximum number of reads considered for analysis. Helpful for reducing run time
10. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT

Example:
1. Using parameters from the config file
```
snakemake --cores 32 --use-conda
```

2. Overriding MIN_Q and PRIMER_BED through command line
```
snakemake --config PRIMER_BED=none.bed MIN_Q=10 --cores 32 --use-conda
```

**WEPP** is a novel phylogenetic method for detecting pathogen variants from the wastewater. Since, WEPP is based on a Phylogentic method, it can be used to detect the variants at the resolution of haplotypes. 

