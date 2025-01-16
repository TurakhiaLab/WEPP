configfile: config.get("config_path","config/config.yaml")

rule filter:
    input:
        "build/wbe",
        "data/{dataset}/" + config["TREE"],
        "intermediate/{dataset}/{file_prefix}_reads.pb",
        "intermediate/{dataset}/{file_prefix}_corrected_variants.tsv",
        "intermediate/{dataset}/{file_prefix}_depth.tsv"
    output:
        "results/{dataset}/{file_prefix}_run.txt"
    conda:
        "../envs/wbe.yml"
    threads:
        workflow.cores
    shell:
        "mkdir -p results/{wildcards.dataset} && "
        "./build/wbe detectPeaks -T {threads} -i " + config["TREE"] + " -p '{wildcards.file_prefix}' -f " + config["REF"] + " -d '{wildcards.dataset}'"
        " | tee {output}"