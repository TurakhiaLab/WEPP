rule freyja:
    input:
        "intermediate/{dataset}/{file_prefix}_resorted.bam"
    output:
        "intermediate/{dataset}/{file_prefix}_corrected_variants.tsv",
        "intermediate/{dataset}/{file_prefix}_depth.tsv"
    conda:
        "../envs/wbe.yml"
    shell:
        """
        cd ./src/Freyja
        make dev
        cd ../..
        freyja variants ./intermediate/{wildcards.dataset}/{wildcards.file_prefix}_resorted.bam --variants ./intermediate/{wildcards.dataset}/{wildcards.file_prefix}_variants.tsv --depths ./intermediate/{wildcards.dataset}/{wildcards.file_prefix}_depth.tsv
        # Correct ivar deletion errors
        python3 src/WBE/ivar_correction.py '{wildcards.dataset}' '{wildcards.file_prefix}'
        """
