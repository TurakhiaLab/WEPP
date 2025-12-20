rule sorted_sam:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_alignment.sam"
    conda:
        wepp_env_file
    shell:
        "samtools view -h -o {output} {input}"

rule sam2pb:
    input:
        binary = str(BASE_DIR / "build/wepp"),
        sam = "intermediate/{DIR}/{FILE_PREFIX}_alignment.sam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_reads.pb"
    conda:
        wepp_env_file
    threads:
        workflow.cores
    params:
        tree = config["TREE"],
        ref = config["REF"],
        max_reads = int(config.get("MAX_READS", 5e6)),
        min_af = config["MIN_AF"],
        min_depth = config["MIN_DEPTH"],
        min_q = config["MIN_Q"]
    shell:
        """
        mkdir -p intermediate/{wildcards.DIR} && \
        {input.binary} sam2PB \
            -T {threads} \
            -i {params.tree} \
            -p '{wildcards.FILE_PREFIX}' \
            -f {params.ref} \
            -d '{wildcards.DIR}' \
            -m {params.max_reads} \
            -a {params.min_af} \
            -c {params.min_depth} \
            -q {params.min_q}
        """