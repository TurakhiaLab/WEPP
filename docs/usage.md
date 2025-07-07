# <b>User Guide</b> <a name="guide"></a>
## <b>Organizing Data</b> <a name="data"></a>
We assume that all wastewater samples are organized in the `data` directory, each within its own subdirectory given by `DIR` argument (see Run Command). For each sample, WEPP generates intermediate and output files in corresponding subdirectories under `intermediate` and `result`, respectively. 

Each created `DIR` inside `data` is expected to contain the following files:

1. Sequencing Reads: Ending with `*_R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
2. Reference Genome in FASTA format
3. Mutation-Annotated Tree (MAT)
4. [OPTIONAL] Genome Masking File: `mask.bed`, whose third column specifies sites to be excluded from analysis.
5. [OPTIONAL] Taxonium `.jsonl` file to be used for visualizing results in the WEPP dashboard. 

Visualization of WEPP's workflow directories
```text
📁 WEPP
└───📁data                                   # [User Created] Contains data to analyze 
    ├───📁SARS_COV_2_real                    # SARS-CoV-2 run wastewater samples - 1
         ├───sars_cov_2_reads_R1.fastq.gz    # Paired-ended reads
         ├───sars_cov_2_reads_R2.fastq.gz
         ├───sars_cov_2_reference.fa 
         ├───mask.bed                        # OPTIONAL 
         ├───sars_cov_2_taxonium.jsonl.gz    # OPTIONAL 
         └───sars_cov_2_mat.pb.gz

└───📁intermediate                           # [WEPP Generated] Contains intermediate stage files 
    ├───📁SARS_COV_2_real                
         ├───file_1
         └───file_2

└───📁results                                # [WEPP Generated] Contains final WEPP results
    ├───📁SARS_COV_2_real                
         ├───file_1
         └───file_2
```

## <b>WEPP Arguments</b><a name="arguments"></a>
The WEPP Snakemake pipeline requires the following arguments, which can be provided either via the configuration file (`config/config.yaml`) or passed directly on the command line using the `--config` argument. The command line arguments take precedence over the config file.

1. `DIR` - Folder name containing the wastewater reads
2. `FILE_PREFIX` - File Prefix for all intermediate files 
3. `REF` - Reference Genome in fasta
4. `TREE` - Mutation-Annotated Tree
5. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)
6. `PRIMER_BED` - BED file for primers from the `primers` folder
7. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked. 
8. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked.
9. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime
10. `CLADE_LIST` - List the clade annotation schemes used in the MAT. SARS-CoV-2 MAT uses both nextstrain and pango lineage naming systems, so use "nextstrain,pango" for it. 
11. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT. Use '1' for Pango naming and '0' for Nextstrain naming for SARS-CoV-2. Other pathogens usually follow a single lineage annotation system, so work with '0'. In case of NO lineage annotations, use '-1'. Lineage Annotations could be checked by running: "matUtils summary -i {TREE} -C {FILENAME}" -> Use '0' for annotation_1 and '1' for annotation_2. 
12. `DASHBOARD_ENABLED` - Set to `True` to enable the interactive dashboard for viewing WEPP results, or `False` to disable it.
13. `TAXONIUM_FILE` [Optional] - Name of the user-provided Taxonium `.jsonl` file for visualization. If specified, this file will be used instead of generating a new one from the given MAT. Ensure that the provided Taxonium file corresponds to the same MAT used for WEPP.

## <b>Run Command</b> <a name="snakemake"></a>
WEPP's snakemake workflow requires `DIR` and `FILE_PREFIX` as config arguments through the command line, while the remaining ones can be taken from the config file. It also requires `--cores` from the command line, which specifies the number of threads used by the workflow.

Examples:

1. Using all the parameters from the config file.
```bash
snakemake --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=sars_cov_2_mat.pb.gz REF=sars_cov_2_reference.fa --cores 32 --use-conda
```

2. Overriding MIN_Q and CLADE_IDX through command line.
```bash
snakemake --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=sars_cov_2_mat.pb.gz REF=sars_cov_2_reference.fa MIN_Q=25 CLADE_IDX=1 --cores 32 --use-conda
```

3. To visualize results from a previous WEPP analysis that was run without the dashboard, set `DASHBOARD_ENABLED` to `True` and re-run only the dashboard components, without reanalyzing the dataset.
```bash
snakemake --config DIR=SARS_COV_2_real FILE_PREFIX=test_run TREE=sars_cov_2_mat.pb.gz REF=sars_cov_2_reference.fa MIN_Q=25 CLADE_IDX=1 DASHBOARD_ENABLED=True --cores 32 --use-conda --forcerun dashboard_serve
```

!!!Note
     ⚠️ Use the same configuration parameters (DIR, FILE_PREFIX, etc.) as were used for the specific project. This ensures the dashboard serves the correct results for your chosen dataset.

## <b>MAT Download</b> <a name="mat"></a>
Mutation-annotated trees (MAT) for different pathogens are maintained by the UShER team, which can be found [here](https://dev.usher.bio). You can also create your own MAT for any pathogen from the consensus genome assemblies using [viral_usher](https://github.com/AngieHinrichs/viral_usher).