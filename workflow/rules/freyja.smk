rule install_modified_freyja:
    output:
        touch(".snakemake/.modified_freyja_installed")
    conda:
        str(BASE_DIR / "workflow/envs/freyja.yml")
    params:
        freyja_dir = str(BASE_DIR / "src/Freyja"),
        installer  = str(BASE_DIR / "workflow/scripts/install_modified_freyja.py")
    shell:
        "python {params.installer} {params.freyja_dir}"


rule freyja:
    input:
        bam = "intermediate/{DIR}/{FILE_PREFIX}_resorted.bam",
        freyja_installed = ancient(".snakemake/.modified_freyja_installed")
    output:
        "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv"
    conda:
        str(BASE_DIR / "workflow/envs/freyja.yml")
    params:
        ref=lambda wildcards: config["REF"],
        minq=lambda wildcards: config["MIN_Q"],
        ivar_script = str(BASE_DIR / "src/WEPP/ivar_correction.py")
    shell:
        """
        freyja variants {input.bam} \
            --variants ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_variants.tsv \
            --depths ./intermediate/{wildcards.DIR}/{wildcards.FILE_PREFIX}_depth.tsv \
            --ref ./data/{wildcards.DIR}/{params.ref} \
            --minq {params.minq}
        python3 {params.ivar_script} '{wildcards.DIR}' '{wildcards.FILE_PREFIX}'
        """
