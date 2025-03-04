# Wastewater based Epidemiology using Phylogenetic Placement (WEPP)

## Snakemake running
Do not forget the --use-conda flag
```
snakemake ./results/manuscript_swampy_dec_2022/my_vcf_run.txt --cores 16 --use-conda
```

<div align="center">
<img src="images/WBE_gif.gif" style="margin: 0px 0px -20px 0px;"/>
</div>

**WEPP** is a novel phylogenetic method for detecting the SARS CoV-2 variants from the wastewater. Since, WEPP is based on a Phylogentic method, it can be used to detect the variants at the resolution of haplotypes. We have two version of WEPP - one is based on MAT called WEPP and the other uses PANMAT called WEPP-PANMAT.  

This repository consists of the entire workflow, which only requires the following things to generate the output. It uses C-WAP pipeline internally to filter and convert fastq to bam format, which is then used by WEPP.
1. Reads in fastq format
2. Phylogenetic Tree as MAT


## Installation
1. git clone --recurse-submodules https://github.com/TurakhiaLab/SARS2-WBE.git
2. cd SARS2-WBE/install
3. conda env create -f environment.yml
4. conda activate usher
5. cd ..
6. mkdir build
7. cd build
8. wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
9. tar -xvzf 2019_U9.tar.gz
10. cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake ..
11. make -j
12. conda deactivate
13. cd ../src/Freyja
14. conda create -n freyja-env
15. conda config --add channels defaults
16. conda config --add channels bioconda
17. conda config --add channels conda-forge
18. conda install freyja
19. make dev
20. conda activate freyja-env
21. pip install clarabel
22. conda deactivate
23. C-WAP installation steps can be found: https://github.com/pgangwar-ucsd/C-WAP
