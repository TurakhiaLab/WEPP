from os.path import join

build_inps = []
for root, _, files in os.walk("src/WEPP"):
    for f in files:
        build_inps.append(join(root, f))

rule install:
    input:
        "CMakeLists.txt"
    output:
        "build/Makefile"
    conda:
        "../envs/wepp.yml"
    shell:
        "./workflow/scripts/install.sh"

rule build_wepp:
    input:
        "build/Makefile",
        build_inps
    output:
        "build/wepp"
    conda:
        "../envs/wepp.yml"
    threads:
        workflow.cores
    shell:
        """
        echo "Starting build..."
        cd build && make -j || exit 1
        """
