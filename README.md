<div align="center">
    
# Wastewater-Based Epidemiology using Phylogenetic Placements

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/WEPP/blob/main/LICENSE

[![License][license-badge]][license-link]
[<img src="https://img.shields.io/badge/Build with-CMake-green.svg?logo=snakemake">](https://cmake.org)
[<img src="https://img.shields.io/badge/Made with-Snakemake-aquamarine.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)

<div align="center">
  <img src="images/WEPP_logo.svg" width="300"/>
</div>

</div>

## Table of Contents
- [Introduction](#intro) ([Wiki](https://turakhia.ucsd.edu/WEPP))
- [Installation](#install)
  - [Summary](#summary) 
  - [Using Install Script](#script)
  - [Using Dockerfile](#docker)
- [Run WEPP](#run)
  - [Default mode](#default)
  - [Iterative mode](#iterative)
- [Contributions](#contribution)
- [Citing WEPP](#cite)

<br>


## <a name="intro"></a> Introduction

WEPP (**T**all and **Wi**de A**lig**nments at **H**igh **T**hroughput) 

By default, TWILIGHT requires an unaligned sequence file in FASTA format and an input guide tree in Newick format to generate the output alignment in FASTA format (Fig. 1a, <a name="default"></a>**default mode**). 

<div align="center">
    <div><b>Figure 1: Overview of WEPP alogorithm</b></div>
    <img src="docs/WEPP_overview.svg" width="800"/>
</div>


## <a name="install"></a> Installation
### <a name="summary"></a> Summary (choose your installation method)

WEPP offers multiple installation methods:
- Install script for directly WEPP running on your system
- Docker (built from the provided Dockerfile) is recommended to prevents any conflict with existing packages

### <a name="script"></a> Using installation script (requires sudo access if certain common libraries are not already installed)  

Users without sudo access are advised to install WEPP via [Docker](#docker).

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
cd WEPP
```
**Step 2:** Install dependencies (might require sudo access)
WEPP depends on the following common system libraries, which are typically pre-installed on most development environments:
```bash
- wget
- curl
- pip
- build-essential 
- python3-pandas
- pkg-config
- zip
- cmake 
- libtbb-dev
- libprotobuf-dev
- protobuf-compiler
- snakemake
```

For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
```bash
sudo apt-get install -y wget pip curl python3-pip build-essential python3-pandas pkg-config zip cmake libtbb-dev libprotobuf-dev protobuf-compiler snakemake
```

### <a name="docker"></a> Using Dockerfile
The Dockerfile installs all the dependencies and tools for WEPP. 

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
cd WEPP
```
**Step 2:** Build a docker image (ensure Docker is installed first)
```
cd docker
docker build -t wepp .
cd ..
```
**Step 3:** Start and run docker container
```
docker run -it -v "$PWD":/workspace -w /workspace wepp /bin/bash
```

## <a name="run"></a> Running WEPP
### <a name="default"></a> Organizing Data
We assume that all the different wastewater samples are stored in the `data` folder under different `DIR` names. Each wastewater `DIR` should have the following files:
1. Sequencing Reads: Ending with `*R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Reference Genome (REF)
3. Mutation Annotated Tree (TREE)
4. [OPTIONAL] Genome Masking File: `mask.bed` whose third column specifies sites to be excluded from analysis.

### <a name="default"></a> Snakemake Arguments
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

### <a name="default"></a> Run Command
Example:
1. Using parameters from the config file
```
snakemake --cores 32 --use-conda
```

2. Overriding MIN_Q and PRIMER_BED through command line
```
snakemake --config PRIMER_BED=none.bed MIN_Q=10 --cores 32 --use-conda
```