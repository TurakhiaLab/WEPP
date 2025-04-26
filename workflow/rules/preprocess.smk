from os import listdir 
from os.path import join

build_inps = []
for root, _, files in os.walk("src/WBE"):
    for f in files:
        build_inps.append(join(root, f))

rule install:
    input:
        "CMakeLists.txt"
    output:
        "build/Makefile"
    conda:
        "../envs/wbe.yml"
    shell:
        "./workflow/scripts/install.sh"

rule build_wbe:
    input:
        "build/Makefile",
        build_inps
    output:
        "build/wbe"
    conda:
        "../envs/wbe.yml"
    threads:
        workflow.cores
    params:
        af_thresh=lambda wildcards: config["AF"],
        conda_path=lambda wildcards: config["CONDA_PATH"],
        minq=lambda wildcards: config["Min_Q"],
        clade_idx=lambda wildcards: config["CLADE_IDX"]
    shell:
        """
        # Update FREQ_READ_THRESHOLD in config.hpp
        sed -i 's#FREQ_READ_THRESHOLD.*$#FREQ_READ_THRESHOLD = {params.af_thresh};#' src/WBE/config.hpp
        
        # Update CONDA_PATH in config.hpp
        sed -i 's#CONDA_PATH.*$#CONDA_PATH = "{params.conda_path}/etc/profile.d/conda.sh";#' src/WBE/config.hpp

        # Update PHRED_SCORE_THRESHOLD in config.hpp
        sed -i 's#PHRED_SCORE_THRESHOLD.*$#PHRED_SCORE_THRESHOLD = {params.minq};#' src/WBE/config.hpp

        # Update CLADE_IDX in config.hpp
        sed -i 's#CLADE_IDX.*$#CLADE_IDX = {params.clade_idx};#' src/WBE/config.hpp
        
        echo "Starting build..."
        cd build && make -j || exit 1
        """
