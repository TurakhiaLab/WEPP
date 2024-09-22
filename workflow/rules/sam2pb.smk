configfile: config.get("config_path","config/config.yaml")

rule sam2pb:
    input:
        "build/wbe",
        "data/{dataset}/{file_prefix}_alignment.sam"
    output:
        "data/{dataset}/{file_prefix}_reads.pb"
    shell:
        "./build/wbe sam2PB -T {workflow.cores} -i " + config["TREE"] + " -v '{wildcards.file_prefix}' -f " + config["REF"] + " -s {wildcards.file_prefix}_alignment.sam -o 'data/{wildcards.dataset}'"