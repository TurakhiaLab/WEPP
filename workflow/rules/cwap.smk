rule cwap:
    input:
        "data/{dataset}/{prefix}_R1.fastq.gz"
    output:
        "intermediate/{dataset}/{prefix}_resorted.bam"
    conda:
        "../envs/cwap.yml"
    params:
        seq_type=lambda wildcards: config.get("SEQUENCING_TYPE", "i"),
        primer_bed=lambda wildcards: config.get("PRIMER_BED", "none.bed")
    shell:
        """
        in_path=$(realpath data/{wildcards.dataset}/)
        out_path=$(realpath intermediate/{wildcards.dataset})
        
        # Navigate to C-WAP directory
        cd ./src/C-WAP

        # Prepare environment
        ./prepare_envs.sh

        # Update primerBedFile in nextflow.config
        substitution='s#primerBedFile.*$#primerBedFile = "$projectDir/covidRefSequences/{params.primer_bed}"#'
        sed -i "$substitution" nextflow.config

        # Clean up previous runs
        rm -rf work
        rm .nextflow.log
        
        # Run nextflow workflow
        nextflow run startWorkflow.nf --platform {params.seq_type} --in $in_path --out point_loma --verbose -wait
        
        # Locate and copy the resulting BAM file
        bam_file=$(find work -name resorted.bam | head -n 1)
        cp $bam_file $out_path/{wildcards.prefix}_resorted.bam
        """