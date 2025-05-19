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
        "./build/wbe detectPeaks -T {threads} -i " + config["TREE"] + " -p '{wildcards.file_prefix}' -f " + config["REF"] + " -d '{wildcards.dataset}'" + " -a " + config["Min_AF"] + " -c " + config["CLADE_IDX"] +
        " | tee {params.tmp_file} && mv {params.tmp_file} {output}"