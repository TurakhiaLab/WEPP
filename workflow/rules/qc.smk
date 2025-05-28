import glob

def get_fastqs(wildcards):
    files = glob.glob(f"data/{wildcards.DIR}/*.fastq.gz")
    if len(files) == 1:
        return files[0]  # Single-end
    elif len(files) == 2:
        return files     # Paired-end
    else:
        raise ValueError(f"Unexpected number of FASTQ files found for {wildcards.DIR}: {files}")

rule qc:
    input:
        get_fastqs
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam"
    conda:
        "../envs/qc.yml"
    params:
        seq_type=lambda wildcards: config["SEQUENCING_TYPE"],
        primer_bed=lambda wildcards: config.get("PRIMER_BED", "none.bed"),
        ref=lambda wildcards: config["REF"]
    threads:
        workflow.cores
    shell:
        """
        in_path=$(realpath data/{wildcards.DIR}/)
        out_path=$(realpath intermediate/{wildcards.DIR})

        python src/qc_preprocess.py --platform {params.seq_type} --primers {params.primer_bed} --in $in_path --out $out_path --threads {threads} --reference {params.ref} --prefix {wildcards.FILE_PREFIX}
        """