rule sorted_sam:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_alignment.sam"
    conda:
        "../envs/wepp.yml"
    shell:
        "samtools view -h -o intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_alignment.sam intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_resorted.bam"

rule sam2pb:
    input:
        "build/wepp",
        "intermediate/{DIR}/{FILE_PREFIX}_alignment.sam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb"
    conda:
        "../envs/wepp.yml"
    threads:
        workflow.cores
    shell:
        "mkdir -p intermediate/{wildcards.DIR} && "
        "./build/wepp sam2PB -T {threads} -i " + config["TREE"] + " -p '{wildcards.FILE_PREFIX}' -d '{wildcards.DIR}'" + " -m " + str(config.get("MAX_READS", str(int(1e9)))) + " -a " + str(config["MIN_AF"]) + " -q " + str(config["MIN_Q"]) + " -c " + str(config["MIN_DEPTH"])