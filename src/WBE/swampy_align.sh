#!/bin/bash

# Check that the correct number of arguments were passed in
if [ $# -lt 7 ]; then
    echo "Usage: $0 <genomes_fasta> <genome_abundances> <reference_fasta> <output_filename_prefix> <output_path> <n_reads> <read_length>"
    exit 1
fi

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

# Read in the alignments
genomes_fasta=$1
genome_abundances=$2
reference_fasta=$3
output_filename_prefix=$4
output_path=$5
n_reads=$6
read_length=$7

mkdir -p ${output_path}/temp/
output_temp=${output_path}/temp/
output_swampy=${output_path}

#SWAMPy    
echo -e "SWAMPy Amplicon generating"
python ../SWAMPy/src/simulate_metagenome.py \
    --genomes_file ${genomes_fasta} \
    --temp_folder ${output_temp} \
    --genome_abundances ${genome_abundances} \
    --primer_set PointLoma\
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

###python WBE/align_mod.py ${output_swampy}/${output_filename_prefix}_alignment.sam
