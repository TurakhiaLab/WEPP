#!/bin/bash

# Check that the correct number of arguments were passed in
if [ $# -eq 7 ]; then
    genomes_fasta=$1
    genome_abundances=$2
    reference_fasta=$3
    output_filename_prefix=$4
    output_path=$5
    n_reads=$6
    read_length=$7
    
    # Check that the files corresponding to the required input arguments exist
    if [ ! -f "$1" ]; then
        echo "Genome fasta file does not exist"
        exit 1
    fi
    if [ ! -f "$2" ]; then
        echo "Genome abundance file does not exist"
        exit 1
    fi
    if [ ! -f "$3" ]; then
        echo "Reference fasta file does not exist"
        exit 1
    fi

    mkdir -p ${output_path}/temp/
    output_temp=${output_path}/temp/
    output_swampy=${output_path}

    #SWAMPy    
    echo -e "SWAMPy Amplicon generating"
    python ../SWAMPy/src/simulate_metagenome.py \
        --genomes_file ${genomes_fasta} \
        --temp_folder ${output_temp} \
        --genome_abundances ${genome_abundances} \
        --primer_set n2\
        --output_folder ${output_swampy} \
        --output_filename_prefix ${output_filename_prefix} \
        --n_reads  ${n_reads} \
        --read_length ${read_length} \
        --amplicon_pseudocounts 1000 \
        --autoremove

    # concatenate r1 and r2 to get final reads
    cat ${output_swampy}/${output_filename_prefix}_R1.fastq ${output_swampy}/${output_filename_prefix}_R2.fastq > ${output_swampy}/${output_filename_prefix}_R1+R2.fastq

    #ALIGNMENT
    echo -e "Running alignment"
    bowtie2-build ${reference_fasta} ${output_temp}/ref_index
    bowtie2 -x ${output_temp}/ref_index -U ${output_swampy}/${output_filename_prefix}_R1+R2.fastq -S ${output_swampy}/${output_filename_prefix}_alignment.sam


elif [ $# -eq 4 ]; then
    reads_fastq=$1
    reference_fasta=$2
    output_filename_prefix=$3
    output_path=$4
    
    # Check that the files corresponding to the required input arguments exist
    if [ ! -f "$1" ]; then
        echo "Reads fastq file does not exist"
        exit 1
    fi
    if [ ! -f "$2" ]; then
        echo "Reference fasta file does not exist"
        exit 1
    fi
    
    mkdir -p ${output_path}/temp/
    output_temp=${output_path}/temp/
    
    echo -e "Running alignment"
    bowtie2-build ${reference_fasta} ${output_temp}/ref_index
    bowtie2 -x ${output_temp}/ref_index -U ${reads_fastq} -S ${output_path}/${output_filename_prefix}_alignment.sam
fi


###python WBE/align_mod.py ${output_swampy}/${output_filename_prefix}_alignment.sam
