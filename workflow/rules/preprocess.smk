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
        conda_path=lambda wildcards: config["CONDA_PATH"]
    shell:
        """
        # Update FREQ_READ_THRESHOLD in config.hpp
        sed -i 's#FREQ_READ_THRESHOLD.*$#FREQ_READ_THRESHOLD = {params.af_thresh};#' src/WBE/config.hpp
        
        # Update CONDA_PATH in config.hpp
        sed -i 's#CONDA_PATH.*$#CONDA_PATH = "{params.conda_path}/etc/profile.d/conda.sh";#' src/WBE/config.hpp
        
        echo "Starting build..."
        cd build && make -j || exit 1
        """
