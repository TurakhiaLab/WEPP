rule filter:
    input:
        "build/wepp",
        "data/{DIR}/" + config["TREE"],
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb",
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv"
    output:
        # other results too, but probably should add all at some point
        "results/{DIR}/{FILE_PREFIX}_run.txt"
    conda:
        "../envs/wepp.yml"
    threads:
        workflow.cores
    params:
        tmp_file=lambda wildcards: f"intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_run_tmp.txt"
    shell:
        "mkdir -p results/{wildcards.DIR} && "
        "./build/wepp detectPeaks -T {threads} -i " + config["TREE"] + " -p '{wildcards.FILE_PREFIX}' -d '{wildcards.DIR}'" + " -a " + str(config["MIN_AF"]) +
        " | tee {params.tmp_file} && mv {params.tmp_file} {output}"