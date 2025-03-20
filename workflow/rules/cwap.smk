rule cwap
    input:
        "data/{dataset}/{prefix}_{ext}.fastq.gz"
    output:
        "intermediate/{dataset}/{prefix}_resorted.bam"
    conda:
        "../envs/wbe.yml"
    shell:
        echo "TODO"