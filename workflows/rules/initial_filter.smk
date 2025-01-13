configfile: config.get("config_path","config/config.yaml")

rule initial_filter:
    input:
        "build/wbe",
        "data/{dataset}/{file_prefix}_reads.pb",
        "data/{dataset}/" + config["TREE"]
    output:
        "data/{dataset}/{file_prefix}_first_checkpoint.txt"
    shell:
        "./build/wbe initial_filter -T {workflow.cores} -i " + config["TREE"] + " -v {wildcards.file_prefix} -f " + config["REF"] + " -o data/{wildcards.dataset}"