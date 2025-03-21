rule cwap:
    input:
        "data/{dataset}/{prefix}_reads.fastq.gz"
    output:
        "intermediate/{dataset}/{prefix}_resorted.bam"
    conda:
        "../envs/wbe.yml"
    shell:
        """
        in_path=$(realpath data/{wildcards.dataset}/)
        out_path=$(realpath intermediate/{wildcards.dataset})
        cd ./src/C-WAP
        source run.sh $in_path $out_path """ + config.get("SEQUENCING_TYPE", "i") + " " + config.get("PRIMER_BED", "none.bed") + " {wildcards.prefix}"