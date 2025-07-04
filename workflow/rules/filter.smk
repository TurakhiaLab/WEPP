rule filter:
    input:
        "build/wepp",
        "data/{DIR}/" + config["TREE"],
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb",
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv"
    output:
        # other results too, but probably should add all at some point
        "intermediate/{DIR}/{FILE_PREFIX}_run_tmp.txt"
    conda:
        "../envs/wepp.yml"
    threads:
        workflow.cores
    params:
        tmp_file=lambda wildcards: f"intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_run_tmp.txt"
    shell:
        "mkdir -p results/{wildcards.DIR} && "
        "./build/wepp detectPeaks -T {threads} -i " + config["TREE"] + " -p '{wildcards.FILE_PREFIX}' -f " + config["REF"] + " -d '{wildcards.DIR}'" + " -a " + str(config["MIN_AF"]) + " -c " + str(config["CLADE_IDX"]) +
        " | tee {params.tmp_file}"