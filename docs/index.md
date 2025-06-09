# <b>Welcome to WEPP Wiki</b>
<div align="center">
    <img src="images/WEPP_logo.svg" width="300">
</div>

## <b>Introduction</b> 
### <b>Overview</b><a name="overview"></a>
WEPP (**W**astewater-Based **E**pidemiology using **P**hylogenetic **P**lacements) is a phylogeny-based pipeline that estimates haplotype proportions from wastewater sequencing reads using a mutation-annotated tree (MAT) (Figure 1A). By improving the resolution of pathogen variant detection, WEPP enables critical epidemiological applications previously feasible only through clinical sequencing. It also flags potential novel variants via *Unaccounted Mutations*, which can be examined at the read level using the interactive dashboard (Figure 1B).

WEPP begins by placing reads on the mutation-annotated tree (MAT) and identifying an initial set of candidate haplotypes. It expands this set by including neighbors around each selected haplotype to form a candidate pool, which is passed to a deconvolution algorithm to estimate haplotype abundances. Haplotypes above a frequency threshold are retained, and their neighbors are again added to form a new candidate pool. This process is repeated iteratively until the haplotype set stabilizes or the maximum number of iterations is reached (Figure 1C).


<div align="center">
    <img src="images/WEPP_Overview_Display.svg" width="600">
    <div><b>Figure 1: Overview of WEPP</b></div>
</div>

### <b>Key Features</b>

#### <b>Haplotype Proportion Estimation</b>  
WEPP's *Phylogenetic Placement* of reads enables accurate estimation of haplotype proportions from wastewater samples. These estimates can be interactively explored using the integrated dashboard (Figure 1B(i)), which displays each haplotype’s abundance, associated lineage, and phylogenetic uncertainty via *Uncertain Haplotypes* - neighboring haplotypes that cannot be confidently disambiguated.

#### <b>Lineage Proportion Estimation</b>  
WEPP infers lineage proportions by combining abundances of haplotypes belonging to each lineage. This approach accounts for intra-lineage diversity, resulting in more accurate and robust estimates.

#### <b>Unccounted Alleles </b>
WEPP reports a list of *Unaccounted Alleles* - alleles observed in wastewater that are not explained by the selected haplotypes, along with the inferred haplotype(s) they are most likely associated with (Figure 1B(i)). These *Unaccounted Alleles* can serve as early indicators of novel variants and often resemble the 'cryptic' mutations described in previous studies.

#### <b>Read-Level Analysis </b>
WEPP supports detailed analysis of sequencing reads in the context of selected haplotypes (Figure 1B(ii)). It also facilitates interpretation of *Unaccounted Alleles* by examining their presence in reads relative to the haplotypes they are mapped to. Additional information about individual reads or haplotypes can be accessed by selecting them within the interactive panel (Figure 1B(iii) and Figure 1B(iv)).
 

## <b>Installation</b> <a name="install"></a>

WEPP offers multiple installation methods. Using a [Docker](https://docs.docker.com/engine/install/) is recommended to prevent any conflict with existing packages.

1. Docker image from DockerHub
2. Dockerfile 
3. Shell Commands 

!!!Note
     ⚠️The Docker image is currently built for the `linux/amd64` platform. While it can run on `arm64` systems (e.g., Apple Silicon or Linux aarch64) via emulation, this may lead to reduced performance.


### **Option-1: Install via DockerHub** <a name=dockerhub></a> 
The Docker image includes all dependencies required to run WEPP.

**Step 1:** Get the image from DockerHub 
```bash
docker pull pranavgangwar/wepp:latest
```
**Step 2:** Start and run Docker container
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 80 or 8080)
# Use this command if your datasets can be downloaded from the Web
docker run -it -p 80:80 pranavgangwar/wepp:latest

# Use this command if your datasets are present in your current directory
docker run -it -p 80:80 -v "$PWD":/WEPP -w /WEPP pranavgangwar/wepp:latest
```
**Step 3:** Confirm proper working by running 
```bash
snakemake test --cores 1 --use-conda
```

### **Option-2: Install via Dockerfile** <a name=dockerfile></a> 
The Dockerfile contains all dependencies required to run WEPP.

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git 
cd WEPP
```
**Step 2:** Build a Docker Image
```bash
cd docker
docker build -t wepp . 
cd ..
```
**Step 3:** Start and run Docker container
```bash
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 80 or 8080)
# Use this command if your datasets can be downloaded from the Web
docker run -it -p 80:80 wepp

# Run this command if your datasets are in the current directory
docker run -it -p 80:80 -v "$PWD":/workspace -w /workspace wepp
```

### **Option-3: Install via Shell Commands (requires sudo access)** <a name=script></a>
Users without sudo access are advised to install WEPP via [Docker Image](#dockerhub).

**Step 1:** Clone the repository
```bash
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
cd WEPP
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
- nodejs
- npm 
- nginx
```

For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
```bash
sudo apt-get install -y wget pip curl python3-pip build-essential python3-pandas pkg-config zip cmake libtbb-dev libprotobuf-dev protobuf-compiler snakemake nodejs npm nginx
```

```bash
# Install Yarn package manager globally
npm install -g yarn
# Install TaxoniumTools Python package
pip install taxoniumtools
```

If your system doesn't have Conda, you can install it with:
```bash
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
```

## <b>Quick Start</b> <a name="example"></a>
The following steps will download a real wastewater RSVA dataset and analyze it with WEPP.

**Step 1:** Download the test dataset
```bash
mkdir -p data/RSVA_real
cd data/RSVA_real
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/011/ERR14763711/ERR14763711_*.fastq.gz https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz
gunzip GCF_002815475.1_ASM281547v1_genomic.fna.gz 
mv ERR14763711_1.fastq.gz ERR14763711_R1.fastq.gz
mv ERR14763711_2.fastq.gz ERR14763711_R2.fastq.gz
cd ../../
```
This will save the datasets on a separate data/RSVA_real folder within the repository.

**Step 2:**  Run the pipeline
```bash
snakemake --config DIR=RSVA_real FILE_PREFIX=test_run PRIMER_BED=RSVA_all_primers_best_hits.bed TREE=rsvA.2025-04-25.pb.gz REF=GCF_002815475.1_ASM281547v1_genomic.fna CLADE_IDX=0 DASHBOARD_ENABLED=True --cores 32 --use-conda
```

**Step 3:**  Analyze Results

All results generated by WEPP can be found in the `results/RSVA_real` directory.


## <b>User Guide</b> <a name="guide"></a>
### <b>Organizing Data</b> <a name="data"></a>
We assume that all wastewater samples are organized in the `data` directory, each within its own subdirectory given by `DIR` argument (see Run Command). For each sample, WEPP generates intermediate and output files in corresponding subdirectories under `intermediate` and `result`, respectively. 

Each created `DIR` inside `data` is expected to contain the following files:

1. Sequencing Reads: Ending with `*R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Reference Genome fasta
3. Mutation-Annotated Tree (MAT)
4. [OPTIONAL] Genome Masking File: `mask.bed`, whose third column specifies sites to be excluded from analysis.
5. [OPTIONAL] Taxonium `.jsonl` file to be used for visualizing results in the WEPP dashboard. 

Visualization of WEPP's workflow directories
```text
📁 WEPP
└───📁data                                # [User Created] Contains data to analyze 
    ├───📁SARS-CoV-2_test_1               # SARS-CoV-2 run wastewater samples
         ├───sars_cov_2_reads.fastq.gz    # Single-ended reads 
         ├───sars_cov_2_reference.fa
         ├───mask.bed                     # OPTIONAL 
         └───sars_cov_2_mat.pb.gz
    ├────📁RSVA_test_1                    # RSVA run wastewater samples 
         ├───rsva_reads_R1.fastq.gz       # Paired-ended reads
         ├───rsva_reads_R2.fastq.gz       # Paired-ended reads
         ├───rsva_reference.fa 
         └───rsva_mat.pb.gz

└───📁intermediate                        # [WEPP Generated] Contains intermediate stage files 
    ├───📁SARS-CoV-2_test_1                
         ├───file_1
         └───file_2
    ├────📁RSVA_test_1                      
         ├───file_1
         └───file_2

└───📁results                             # [WEPP Generated] Contains final WEPP results
    ├───📁SARS-CoV-2_test_1                
         ├───file_1
         └───file_2
    ├────📁RSVA_test_1                      
         ├───file_1
         └───file_2
```

### <b>WEPP Arguments</b><a name="arguments"></a>
The WEPP Snakemake pipeline requires the following arguments, which can be provided either via the configuration file (`config/config.yaml`) or passed directly on the command line using the `--config` argument. The command line arguments take precedence over the config file.

1. `DIR` - Folder name containing the wastewater reads
2. `FILE_PREFIX` - File Prefix for all intermediate files 
3. `REF` - Reference Genome in fasta
4. `TREE` - Mutation-Annotated Tree
5. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads).
6. `PRIMER_BED` - BED file for primers from the `primers` folder.
7. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked. 
8. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked.
9. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime.
10. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT. Generally '1' for SARS-CoV-2 MATs and '0' for others. Could be checked by running: "matUtils summary -i {TREE} -C {FILENAME}" -> Use '0' for annotation_1 and '1' for annotation_2
11. `DASHBOARD_ENABLED` - Set to `True` to enable the interactive dashboard for viewing WEPP results, or `False` to disable it.
12. `TAXONIUM_FILE` [Optional] - Name of the user-provided Taxonium `.jsonl` file for visualization. If specified, this file will be used instead of generating a new one from the given MAT. Ensure that the provided Taxonium file corresponds to the same MAT used for WEPP. 

### <b>Run Command</b> <a name="snakemake"></a>
WEPP's snakemake workflow requires `DIR` and `FILE_PREFIX` as config arguments through the command line, while the remaining ones can be taken from the config file. It also requires `--cores` from the command line, which specifies the number of threads used by the workflow.

Examples:

1. Using all the parameters from the config file
```bash
snakemake --config DIR=SARS-CoV-2_test_1 FILE_PREFIX=test_run --cores 32 --use-conda
```

2. Overriding MIN_Q, PRIMER_BED, and DASHBOARD_ENABLED through command line.
```bash
snakemake --config DIR=RSVA_test_1 FILE_PREFIX=test_run MIN_Q=25 PRIMER_BED=none.bed DASHBOARD_ENABLED=True --cores 32 --use-conda
```
3. To visualize results from a previous WEPP analysis that was run without the dashboard, set `DASHBOARD_ENABLED` to `True` and re-run only the dashboard components, without reanalyzing the dataset.
```bash
snakemake --config DIR=SARS-CoV-2_test_1 FILE_PREFIX=test_run DASHBOARD_ENABLED=True --cores 32 --use-conda --forcerun dashboard_serve
```

!!!Note
     ⚠️ Use the same configuration parameters (DIR, FILE_PREFIX, etc.) as were used for the specific project. This ensures the dashboard serves the correct results for your chosen dataset.

## <b>Contributions</b>
We welcome contributions from the community to enhance the capabilities of WEPP. If you encounter any issues or have suggestions for improvement, please open an issue on [WEPP GitHub page](https://github.com/TurakhiaLab/WEPP). For general inquiries and support, reach out to our team.

## <b>Citing WEPP</b>
TBA.