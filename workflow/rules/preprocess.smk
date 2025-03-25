from os import listdir 
from os.path import join

build_inps = []
excluded_dirs = {"src/C-WAP", "src/Freyja"}
for root, _, files in os.walk("src"):
    if any(root.startswith(excluded) for excluded in excluded_dirs):
        continue

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
    shell:
        "cd build && make -j"
