<div align="center">
    
# Wastewater-Based Epidemiology using Phylogenetic Placements

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/WEPP/blob/main/LICENSE

[![License][license-badge]][license-link]
[<img src="https://img.shields.io/badge/Build with-CMake-green.svg?logo=snakemake">](https://cmake.org)
[<img src="https://img.shields.io/badge/Made with-Snakemake-aquamarine.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)

<div align="center">
  <img src="docs/images/WEPP_logo.svg" width="300"/>
</div>

</div>

## Table of Contents
- [Introduction](#intro) ([Wiki](https://turakhia.ucsd.edu/WEPP))
- [Installation](#install)
  - [Option-1: Install via DockerHub](#dockerhub)
  - [Option-2: Install via Dockerfile](#dockerfile)
  - [Option-3: Install via Shell Commands](#script)
- [Quick Start](#example)
     - [Example-1: RSV-A Dataset (Runs Quickly)](#rsv_a_example)
     - [Example-2: SARS-CoV-2 Dataset (Longer Runtime)](#sars_cov_2_example)
- [User Guide](#guide)
  - [Organizing Data](#data)
  - [WEPP Arguments](#arguments)
  - [Run Command](#snakemake)
  - [Getting Mutation-Annotated Trees](#mat)
- [Contributions](#contribution)
- [Citing WEPP](#cite)

<br>


## <a name="intro"></a> Introduction

WEPP (**W**astewater-Based **E**pidemiology using **P**hylogenetic **P**lacements) is a pathogen-agnostic pipeline that enhances wastewater surveillance by leveraging the pathogen’s full phylogeny. It reports haplotype and lineage abundances, maps reads parsimoniously to selected haplotypes, and flags *Unaccounted Alleles* — those observed in the sample but unexplained by selected haplotypes, potentially indicating novel variants (Figure 1A). An interactive dashboard enables visualization of haplotypes in the global phylogenetic tree (Figure 1B(i)) and read-level analysis by selecting a haplotype (Figure 1B(ii)). Additional information about individual reads or haplotypes can be accessed by clicking on their corresponding objects (Figures 1B(iii-iv)), respectively. 

WEPP performs parsimonious read placement on the mutation-annotated tree (MAT) to select a subset of haplotypes and adds their neighbors to form an initial candidate pool, which is passed to a deconvolution algorithm to estimate their relative abundances. WEPP retains haplotypes above an abundance threshold, and iteratively adds their neighbors and recomputes abundances until convergence or a maximum iteration count. An outlier detection algorithm flags Unaccounted Alleles from the deconvolution residue (Figure 1C).


<div style="width: 600px; margin: 0 auto;">
   <img src="docs/images/WEPP_Overview_Display.svg" style="display: block; margin-left: auto; margin-right: auto; width: 100%;">
  <div style="text-align: justify; margin-top: 10px;">
    <b>Figure 1: Overview of the WEPP pipeline.</b> (A) WEPP input and output. (B) Features of the interactive Dashboard: (i) Phylogenetic view of WEPP-inferred haplotypes with their proportions, associated lineages, and uncertain haplotypes. Unaccounted alleles and their possible haplotype sources are shown in a separate panel; (ii) Read analysis panel highlighting accounted and unaccounted alleles contained in reads mapped to a selected haplotype; (iii) Read information panel displaying all possible haplotypes and unaccounted alleles for a selected read; (iv) Haplotype information panel listing the possible unaccounted alleles associated with the selected haplotype. (C) Key stages of WEPP’s phylogenetic algorithm for haplotype detection and abundance estimation.
    </div>
</div>


## <a name="install"></a> Installation
WEPP offers multiple installation methods. Using a [Docker](https://docs.docker.com/engine/install/) is recommended to prevent any conflict with existing packages.
1. Docker image from DockerHub
2. Dockerfile 
3. Shell Commands

⚠️ The Docker image is currently built for the `linux/amd64` platform. While it can run on `arm64` systems (e.g., Apple Silicon or Linux aarch64) via emulation, this may lead to reduced performance.

### <a name="dockerhub"></a> Option-1: Install via DockerHub
The Docker image includes all dependencies required to run WEPP.

**Step 1:** Get the image from DockerHub 
```bash
docker pull pranavgangwar/wepp:latest
```
**Step 2:** Start and run Docker container. The command below will take you inside the docker container with WEPP already installed.
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 80 or 8080)
# Use this command if your datasets can be downloaded from the Web
docker run -it -p 80:80 pranavgangwar/wepp:latest

# Use this command if your datasets are present in your current directory
docker run -it -p 80:80 -v "$PWD":/WEPP -w /WEPP pranavgangwar/wepp:latest
```
**Step 3:** Confirm proper working by running the following command. This should print WEPP's help menu.
```bash
run-wepp help --cores 1 --use-conda
```

All set to try the [examples](#example).


### <a name="dockerfile"></a> Option-2: Install via Dockerfile
The Dockerfile contains all dependencies required to run WEPP.

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git 
cd WEPP
chmod +x run-wepp
```
**Step 2:** Build a Docker Image
```bash
cd docker
docker build -t wepp . 
cd ..
```
**Step 3:** Start and run Docker container. The command below will take you inside the docker container with the view of the current directory.
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 80 or 8080)
docker run -it -p 80:80 -v "$PWD":/workspace -w /workspace wepp
```

All set to try the [examples](#example).


### <a name="script"></a> Option-3: Install via Shell Commands (requires sudo access)  
Users without sudo access are advised to install WEPP via [Docker Image](#dockerhub).

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
cd WEPP
chmod +x run-wepp
```
**Step 2:** Install dependencies (might require sudo access)
WEPP depends on the following common system libraries, which are typically pre-installed on most development environments:
```text
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
- conda
- nodejs(v18+)
- nginx
```

For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
```bash
sudo apt-get update
sudo apt-get install -y wget pip curl python3-pip build-essential python3-pandas pkg-config zip cmake libtbb-dev libprotobuf-dev protobuf-compiler snakemake nginx
```

Note: WEPP expects the `python` command to be available. If your system only provides python3, you can optionally set up a symlink:
```bash
update-alternatives --install /usr/bin/python python /usr/bin/python3 1
```

If you do not have Node.js v18 or higher installed, follow these steps to install Node.js v22:
```bash
# Update and install prerequisites
apt-get install -y curl gnupg ca-certificates

# Add NodeSource Node.js 22 repo
curl -fsSL https://deb.nodesource.com/setup_22.x | bash -

# Install Node.js 22
apt-get install -y nodejs
```

```bash
# Install Yarn package manager globally
npm install -g yarn
```

If your system doesn't have Conda, you can install it with:
```bash
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
```

All set to try the [examples](#example).


##  <a name="example"></a> Quick Start
The following steps will download real wastewater datasets and analyze them using WEPP.

### <a name="rsv_a_example"></a> Example - 1: RSV-A Dataset (Runs Quickly: Under 10 minutes on 32 cores)
**Step 1:** Download the RSV-A test dataset
```bash
mkdir -p data/RSVA_real
cd data/RSVA_real
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/011/ERR14763711/ERR14763711_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/011/ERR14763711/ERR14763711_2.fastq.gz https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz
gunzip GCF_002815475.1_ASM281547v1_genomic.fna.gz 
mv ERR14763711_1.fastq.gz ERR14763711_R1.fastq.gz
mv ERR14763711_2.fastq.gz ERR14763711_R2.fastq.gz
cd ../../
```
This will save the datasets on a separate data/RSVA_real folder within the repository.

**Step 2:**  Run the pipeline
```bash
run-wepp --config DIR=RSVA_real FILE_PREFIX=test_run TREE=rsvA.2025-04-25.pb.gz REF=GCF_002815475.1_ASM281547v1_genomic.fna CLADE_LIST=annotation_1 CLADE_IDX=0 DASHBOARD_ENABLED=True --cores 32 --use-conda
```

**Step 3:**  Analyze Results

All results generated by WEPP are available in the `results/RSVA_real` directory. These include haplotype and lineage abundances, associated uncertain haplotypes, and the potential haplotypes corresponding to each detected unaccounted allele.

⚠️ Make sure port forwarding is enabled when accessing services on external servers.


### <a name="sars_cov_2_example"></a> Example - 2: SARS-CoV-2 Dataset (Longer Runtime: ~20 minutes on 32 cores)
**Step 1:** Download the SARS-CoV-2 test dataset
```bash
mkdir -p data/SARS_COV_2_real
cd data/SARS_COV_2_real
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR185/041/SRR18541041/SRR18541041_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR185/041/SRR18541041/SRR18541041_2.fastq.gz https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
mv SRR18541041_1.fastq.gz SRR18541041_R1.fastq.gz
mv SRR18541041_2.fastq.gz SRR18541041_R2.fastq.gz
cp ../../NC_045512v2.fa .
cd ../../
```
This will save the datasets on a separate data/SARS_COV_2_real folder within the repository.

**Step 2:**  Run the pipeline
```bash
run-wepp --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=public-2021-12-05.all.masked.pb.gz REF=NC_045512v2.fa DASHBOARD_ENABLED=True --cores 32 --use-conda
```

**Step 3:**  Analyze Results

All results generated by WEPP are available in the `results/SARS_COV_2_real` directory. These include haplotype and lineage abundances, associated uncertain haplotypes, and the potential haplotypes corresponding to each detected unaccounted allele.

⚠️ Make sure port forwarding is enabled when accessing services on external servers.

## <a name="guide"></a> User Guide
### <a name="data"></a> Organizing Data
We assume that all wastewater samples are organized in the `data` directory, each within its own subdirectory given by `DIR` argument (see Run Command). For each sample, WEPP generates intermediate and output files in corresponding subdirectories under `intermediate` and `result`, respectively. 
Each created `DIR` inside `data` is expected to contain the following files:
1. Sequencing Reads: Ending with `*_R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Reference Genome in FASTA format
3. Mutation-Annotated Tree (MAT)
4. [OPTIONAL] Genome Masking File: `mask.bed`, whose third column specifies sites to be excluded from analysis.
5. [OPTIONAL] Taxonium `.jsonl` file to be used for visualizing results in the WEPP dashboard. 

Visualization of WEPP's workflow directories
```text
📁 WEPP
└───📁data                                   # [User Created] Contains data to analyze 
    ├───📁SARS_COV_2_real                    # SARS-CoV-2 run wastewater samples - 1
         ├───sars_cov_2_reads_R1.fastq.gz    # Paired-ended reads
         ├───sars_cov_2_reads_R2.fastq.gz
         ├───sars_cov_2_reference.fa 
         ├───mask.bed                        # OPTIONAL 
         ├───sars_cov_2_taxonium.jsonl.gz    # OPTIONAL 
         └───sars_cov_2_mat.pb.gz

└───📁intermediate                           # [WEPP Generated] Contains intermediate stage files 
    ├───📁SARS_COV_2_real                
         ├───file_1
         └───file_2

└───📁results                                # [WEPP Generated] Contains final WEPP results
    ├───📁SARS_COV_2_real                
         ├───file_1
         └───file_2
```

### <a name="arguments"></a> WEPP Arguments
The WEPP Snakemake pipeline requires the following arguments, which can be provided either via the configuration file (`config/config.yaml`) or passed directly on the command line using the `--config` argument. The command line arguments take precedence over the config file.
1. `DIR` - Folder name containing the wastewater reads.
2. `FILE_PREFIX` - File Prefix for all intermediate files. 
3. `REF` - Reference Genome in fasta.
4. `TREE` - Mutation-Annotated Tree.
5. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads).
6. `PRIMER_BED` - BED file argument for primers, with few primers provided in the `primers` folder. Requires path to the file.
7. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked (Illumina: 0.5%, Ion Torrent: 1.5%, ONT: 2%).
8. `MIN_DEPTH` - Sites with read depth below this threshold will be masked. 
9. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked.
10. `MIN_PROP` - Minimum Proportion of haplotypes (Wastewater Samples: 0.5%, Clinical Samples: 5%).
11. `MIN_LEN` - Minimum read length to be considered after ivar trim (Deafult: 80).
12. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime.
13. `CLADE_LIST` - List the clade annotation schemes used in the MAT. SARS-CoV-2 MAT uses both nextstrain and pango lineage naming systems, so use "nextstrain,pango" for it. 
14. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT. Use '1' for Pango naming and '0' for Nextstrain naming for SARS-CoV-2. Other pathogens usually follow a single lineage annotation system, so work with '0'. In case of NO lineage annotations, use '-1'. Lineage Annotations could be checked by running: "matUtils summary -i {TREE} -C {FILENAME}" -> Use '0' for annotation_1 and '1' for annotation_2. 
15. `DASHBOARD_ENABLED` - Set to `True` to enable the interactive dashboard for viewing WEPP results, or `False` to disable it.
16. `TAXONIUM_FILE` [Optional] - Name of the user-provided Taxonium `.jsonl` file for visualization. If specified, this file will be used instead of generating a new one from the given MAT. Ensure that the provided Taxonium file corresponds to the same MAT used for WEPP.

### <a name="snakemake"></a> Run Command
WEPP's snakemake workflow requires `DIR` and `FILE_PREFIX` as config arguments through the command line, while the remaining ones can be taken from the config file. It also requires `--cores` from the command line, which specifies the number of threads used by the workflow.

Examples:
1. Using all the parameters from the config file.
```bash
run-wepp --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=sars_cov_2_mat.pb.gz REF=sars_cov_2_reference.fa --cores 32 --use-conda
```

2. Overriding MIN_Q and CLADE_IDX through command line.
```bash
run-wepp --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=sars_cov_2_mat.pb.gz REF=sars_cov_2_reference.fa MIN_Q=25 CLADE_IDX=1 --cores 32 --use-conda
```

3. To visualize results from a previous WEPP analysis that was run without the dashboard, set `DASHBOARD_ENABLED` to `True` and re-run only the dashboard components, without reanalyzing the dataset.
```bash
run-wepp --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=sars_cov_2_mat.pb.gz REF=sars_cov_2_reference.fa MIN_Q=25 CLADE_IDX=1 DASHBOARD_ENABLED=True --cores 32 --use-conda --forcerun dashboard_serve
```
⚠️ Use the same configuration parameters (DIR, FILE_PREFIX, etc.) as were used for the specific project. This ensures the dashboard serves the correct results for your chosen dataset.

⚠️ Make sure port forwarding is enabled when running on external servers to view results on your personal machine.

###  <a name="mat"></a> Getting Mutation-Annotated Trees
Mutation-annotated trees (MAT) for different pathogens are maintained by the UShER team, which can be found [here](https://dev.usher.bio). You can also create your own MAT for any pathogen from the consensus genome assemblies using [viral_usher](https://github.com/AngieHinrichs/viral_usher).

##  <a name="contribution"></a> Contributions
We welcome contributions from the community to enhance the capabilities of **WEPP**. If you encounter any issues or have suggestions for improvement, please open an issue on [WEPP GitHub page](https://github.com/TurakhiaLab/WEPP/issues). For general inquiries and support, reach out to our team.

##  <a name="cite"></a> Citing WEPP
If you use WEPP in your research or publications, please cite the following paper:<br>
* Pranav Gangwar, Pratik Katte, Manu Bhatt, Yatish Turakhia, "<i>WEPP: Phylogenetic Placement Achieves Near-Haplotype Resolution in Wastewater-Based Epidemiology</i>", medRxiv 2025.06.09.25329287; doi: [10.1101/2025.06.09.25329287](https://doi.org/10.1101/2025.06.09.25329287)
