#!/bin/bash

# This script implements the following steps:
# Inputs: path to VCF file and path to corresponding reference genome (fasta)
# 1. Convert the VCF file to a fasta file, using the reference genome
# 2. Align the resulting fasta file to the reference genome using bowtie2, and output the alignment to a SAM file
# 3. Convert the SAM file to a BAM file
# 4. Convert the BAM file back to a FASTA file
# 5. Convert the FASTA file to a VCF file --> output

# Usage: bash get_aligned_vcf.sh <path to VCF file> <path to reference genome (fasta)> <output file name>
# Example: bash get_aligned_vcf.sh ../data/chr22_maf0.05_miss0.05.vcf ../data/chr22.fa

# Check the user uses the syntax correctly
if [ $# -ne 3 ]
then
    echo "Usage: bash get_aligned_vcf.sh <path to VCF file> <path to reference genome (fasta)> <output file name>"
    exit 1
fi

VCF=$1
REF=$2
OUT=$3

mkdir -p ./intermediate_files
mkdir -p ./intermediate_files/outputs
mkdir -p ./intermediate_files/ref

SAM=./intermediate_files/outputs/output.sam
BAM=./intermediate_files/outputs/output.bam
FASTA1=./intermediate_files/outputs/output.fasta
FASTA2=./intermediate_files/outputs/output2.fasta
REF_INDEX=./intermediate_files/ref/ref_index

# Convert the VCF file to a fasta file, using the reference genome
python vcf_to_fasta.py $VCF $REF $FASTA1

# The result of the above program is stored in output.fasta

# Align the resulting fasta file to the reference genome using bowtie2, and output the alignment to a SAM file
bowtie2-build -f $REF $REF_INDEX
bowtie2 -x $REF_INDEX -f $FASTA1 -S $SAM

# Convert the SAM file to a sorted BAM file
samtools view -bS $SAM | samtools sort -o $BAM

# Convert the BAM file back to a FASTA file
samtools fasta $BAM > $FASTA2

# Convert the FASTA file to a VCF file
python fasta_to_vcf.py $FASTA2 $REF $OUT