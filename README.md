# Wastewater based Epidemiology using Phylogenetic Placement (WEPP)

[license-badge]:
[license-link]:


## Installation
1. git clone --recurse-submodules https://github.com/TurakhiaLab/SARS2-WBE.git
2. cd SARS2-WBE/install
3. conda env create -f environment.yml
4. conda activate SARS2-WBE
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