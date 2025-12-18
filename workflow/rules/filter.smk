rule filter:
    input:
        binary = str(BASE_DIR / "build/wepp"),
        tree = lambda w: f"data/{w.DIR}/{config['TREE']}",
        pb = "intermediate/{DIR}/{FILE_PREFIX}_reads.pb",
        variants = "intermediate/{DIR}/{FILE_PREFIX}_corrected_variants.tsv",
        depth = "intermediate/{DIR}/{FILE_PREFIX}_depth.tsv"
    output:
        tmp = "intermediate/{DIR}/{FILE_PREFIX}_run_tmp.txt",
        unaccounted = "results/{DIR}/{FILE_PREFIX}_unaccounted_alleles.txt"
    conda:
        str(BASE_DIR / "workflow/envs/freyja.yml")
    threads:
        workflow.cores
    params:
        tree = config["TREE"],
        ref = config["REF"],
        min_af = config["MIN_AF"],
        min_prop = config["MIN_PROP"],
        clade_idx = config["CLADE_IDX"]
    shell:
        """
        mkdir -p results/{wildcards.DIR} && \
        {input.binary} detectPeaks \
            -T {threads} \
            -i {params.tree} \
            -p '{wildcards.FILE_PREFIX}' \
            -f {params.ref} \
            -d '{wildcards.DIR}' \
            -a {params.min_af} \
            -r {params.min_prop} \
            -n {params.clade_idx} \
        | tee {output.tmp} && \
        mv intermediate/{wildcards.DIR}/residual_mutations.txt {output.unaccounted}
        """