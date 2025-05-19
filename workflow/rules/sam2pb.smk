rule sorted_sam:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_alignment.sam"
    conda:
        "../envs/wbe.yml"
    shell:
        "samtools view -h -o intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_alignment.sam intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_resorted.bam"

rule sam2pb:
    input:
        "build/wbe",
        "intermediate/{DIR}/{FILE_PREFIX}_alignment.sam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb"
    conda:
        "../envs/wbe.yml"
    threads:
        workflow.cores
    shell:
        "mkdir -p intermediate/{wildcards.DIR} && "
        "./build/wbe sam2PB -T {threads} -i " + config["TREE"] + " -p '{wildcards.FILE_PREFIX}' -f " + config["REF"] + " -d '{wildcards.DIR}'" + " -m " + str(config.get("MAX_READS", str(int(1e9)))) + " -a " + str(config["MIN_AF"]) + " -q " + str(config["MIN_Q"]) + " -c " + str(config["CLADE_IDX"])