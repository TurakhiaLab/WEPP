configfile: config.get("config_path","config/config.yaml")

rule filter:
    input:
        "build/wbe",
        "data/{dataset}/" + config["TREE"],
        "intermediate/{dataset}/{file_prefix}_reads.pb",
        "intermediate/{dataset}/{file_prefix}_corrected_variants.tsv",
        "intermediate/{dataset}/{file_prefix}_depth.tsv"
    output:
        # other results too, but probably should add all at some point
        "results/{dataset}/{file_prefix}_run.txt"
    conda:
        "../envs/wbe.yml"
    threads:
        workflow.cores
    params:
        tmp_file=lambda wildcards: f"intermediate/{wildcards.dataset}/{wildcards.file_prefix}_run_tmp.txt"
    shell:
        "mkdir -p results/{wildcards.dataset} && "
        "./build/wbe detectPeaks -T {threads} -i " + config["TREE"] + " -v '{wildcards.file_prefix}' -d '{wildcards.dataset}'"
        " | tee {params.tmp_file} && mv {params.tmp_file} {output}"