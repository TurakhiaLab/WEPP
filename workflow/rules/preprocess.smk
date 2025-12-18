from os.path import join

src_dir = str(BASE_DIR / "src/WEPP")
build_inps = []
if os.path.exists(src_dir):
    for root, _, files in os.walk(src_dir):
        for f in files:
            build_inps.append(os.path.join(root, f))

rule install:
    input:
        str(BASE_DIR / "CMakeLists.txt")
    output:
        str(BASE_DIR / "build/Makefile")
    conda:
        wepp_env_file
    params:
        source_dir = str(BASE_DIR),
        script = str(BASE_DIR / "workflow/scripts/install.sh")
    shell:
        """
        workspace=$(pwd)

        cd {params.source_dir}
        bash {params.script}

        cd $workspace
        """

rule build_wepp:
    input:
        str(BASE_DIR / "build/Makefile"),
        src = build_inps
    output:
        str(BASE_DIR / "build/wepp")
    conda:
        wepp_env_file
    threads:
        workflow.cores
    params: 
        build_dir = str(BASE_DIR / "build")
    shell:
        """
        echo "Starting build..."
        cd {params.build_dir}
        make -j || exit 1
        """
