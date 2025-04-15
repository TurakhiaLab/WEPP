<p align="center">
  <img src="WEPP_logo.svg" width="200" height="200">
</p>

<h1 align="center">
  Wastewater-based Epidemiology using Phylogenetic Placement
</h1>

## Installation and running
0. `git clone --recurse-submodules https://github.com/TurakhiaLab/SARS2-WBE.git`.
1. Place the reads under the `./data/{dataset}/` folder ending with the name `*R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Create a config/config.yaml file that contains the name of the MAT under the key "TREE" and the reference sequence file for the MAT under the key "REF". The MAT and reference sequence should also be placed in `./data/{dataset}` folder. Pick the primer bed file from the database folder.
3. Run snakemake. Do not forget the --use-conda flag. Optionally, you can specific configuration options here instead.
```
snakemake ./results/{dataset}/{file_prefix}_run.txt --cores 16 --use-conda
```

**WEPP** is a novel phylogenetic method for detecting the SARS CoV-2 variants from the wastewater. Since, WEPP is based on a Phylogentic method, it can be used to detect the variants at the resolution of haplotypes. We have two version of WEPP - one is based on MAT called WEPP and the other uses PANMAT called WEPP-PANMAT.  

This repository consists of the entire workflow, which only requires the following things to generate the output. It uses C-WAP pipeline internally to filter and convert fastq to bam format, which is then used by WEPP.
1. Reads in fastq format
2. Phylogenetic Tree as MAT
