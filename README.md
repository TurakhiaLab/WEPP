<p align="center">
  <img src="WEPP_logo.svg" width="300" height="300">
</p>

<h1 align="center">
  Wastewater-based Epidemiology using Phylogenetic Placement
</h1>

## Installation and running
0. `git clone --recurse-submodules https://github.com/TurakhiaLab/SARS2-WBE.git`.
1. Place the reads under the `./data/{dataset}/` folder ending with the name `*R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Create a config/config.yaml file that contains the name of the MAT under the key "TREE". The MAT should also be placed in `./data/{dataset}` folder. Pick the primer bed file from the database folder.
3. Run snakemake. Do not forget the --use-conda flag. Optionally, you can specific configuration options here instead.
```
snakemake ./results/{dataset}/{file_prefix}_run.txt --cores 16 --use-conda
```

**WEPP** is a novel phylogenetic method for detecting the SARS CoV-2 variants from the wastewater. Since, WEPP is based on a Phylogentic method, it can be used to detect the variants at the resolution of haplotypes. We have two version of WEPP - one is based on MAT called WEPP-MAT and the other uses PANMAT called WEPP-PANMAT.  

This repository contains the complete workflow, which requires the following inputs to generate results:
1. Sequencing reads in `fastq.gz` format
2. A phylogenetic tree in `MAT` format
3. A primer `BED` file
