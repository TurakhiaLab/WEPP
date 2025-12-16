rule freyja:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv"
    conda:
        "../envs/freyja.yml"
    params:
        ref=lambda wildcards: config["REF"],
        minq=lambda wildcards: config["MIN_Q"]
    shell:
        """
        cd ./src/Freyja
        make dev
        cd ../..
        freyja variants ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_resorted.bam --variants ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_variants.tsv --depths ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_depth.tsv --ref ./data/{wildcards.DIR}/{params.ref} --minq {params.minq}
        # Correct ivar deletion errors
        python3 src/WEPP/ivar_correction.py '{wildcards.DIR}' '{wildcards.FILE_PREFIX}'
        """
