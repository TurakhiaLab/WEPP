configfile: config.get("config_path","config/config.yaml")

rule sorted_sam:
    input:
        "intermediate/{dataset}/{file_prefix}_resorted.bam"
    output:
        "intermediate/{dataset}/{file_prefix}_alignment.sam"
    conda:
        "../envs/wbe.yml"
    shell:
        "samtools view -h -o intermediate/{wildcards.dataset}/{wildcards.file_prefix}_alignment.sam intermediate/{wildcards.dataset}/{wildcards.file_prefix}_resorted.bam"

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
        "./build/wbe sam2PB -T {threads} -i " + config["TREE"] + " -v '{wildcards.file_prefix}' -d '{wildcards.dataset}'"