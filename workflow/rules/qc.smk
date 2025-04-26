import glob

def get_fastqs(wildcards):
    files = glob.glob(f"data/{wildcards.dataset}/*.fastq.gz")
    if len(files) == 1:
        return files[0]  # Single-end
    elif len(files) == 2:
        return files     # Paired-end
    else:
        raise ValueError(f"Unexpected number of FASTQ files found for {wildcards.dataset}: {files}")

rule qc:
    input:
        get_fastqs
    output:
        "intermediate/{dataset}/{file_prefix}_resorted.bam"
    conda:
        "../envs/qc.yml"
    params:
        seq_type=lambda wildcards: config.get("SEQUENCING_TYPE", "s"),
        primer_bed=lambda wildcards: config.get("PRIMER_BED", "none.bed"),
        ref=lambda wildcards: config["REF"]
    threads:
        workflow.cores
    shell:
        """
        in_path=$(realpath data/{wildcards.dataset}/)
        out_path=$(realpath intermediate/{wildcards.dataset})

        python src/WBE/qc_preprocess.py --platform {params.seq_type} --primers {params.primer_bed} --in $in_path --out $out_path --threads {threads} --reference {params.ref} --prefix {wildcards.file_prefix}
        """