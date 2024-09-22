configfile: config.get("config_path","config/config.yaml")

rule post_filter:
    input:
        "build/wbe",
        "data/{dataset}/{file_prefix}_first_checkpoint.txt",
        "data/{dataset}/" + config["TREE"],
        "data/{dataset}/{file_prefix}_variants.tsv",
        "data/{dataset}/{file_prefix}_depth.tsv"
    output:
        "data/{dataset}/{file_prefix}_run.txt"
    shell:
        "./build/wbe post_filter -T {workflow.cores} -i " + config["TREE"] + " -v {wildcards.file_prefix} -f " + config["REF"] + " -o data/{wildcards.dataset} |"
        "tee {output}"