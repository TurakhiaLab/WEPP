rule freyja:
    input:
        "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam"
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv"
    conda:
        str(BASE_DIR / "workflow/envs/freyja.yml")
    params:
        ref=lambda wildcards: config["REF"],
        minq=lambda wildcards: config["MIN_Q"],
        freyja_dir = str(BASE_DIR / "src/Freyja"),
        ivar_script = str(BASE_DIR / "src/WEPP/ivar_correction.py")
    shell:
        """
        (cd {params.freyja_dir} && make dev)
        freyja variants {input} \
            --variants ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_variants.tsv \
            --depths ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_depth.tsv \
            --ref ./data/{wildcards.DIR}/{params.ref} \
            --minq {params.minq}
        python3 {params.ivar_script} '{wildcards.DIR}' '{wildcards.FILE_PREFIX}'
        """
