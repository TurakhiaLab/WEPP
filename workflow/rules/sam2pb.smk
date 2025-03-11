configfile: config.get("config_path","config/config.yaml")

rule sam2pb:
    input:
        "build/wbe",
        "intermediate/{dataset}/{file_prefix}_alignment.sam"
    output:
        "intermediate/{dataset}/{file_prefix}_reads.pb"
    conda:
        "../envs/wbe.yml"
    threads:
        # sam2pb is early so we may want to allow other 
        # jobs to run as well
        max(workflow.cores - 1, 1)
    shell:
        "mkdir -p intermediate/{wildcards.dataset} && "
        "./build/wbe sam2PB -T {threads} -i " + config["TREE"] + " -p '{wildcards.file_prefix}' -f " + config["REF"] + " -d '{wildcards.dataset}'" + " -m  " + config.get("MAX_READS", 1e9)