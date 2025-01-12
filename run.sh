#!/bin/bash

#/home/swalia@AD.UCSD.EDU/work/panman/build/panman/viralmsa_8M.panman

set -e

#This script assumes that you have WBE, SWAMPy, and Freyja setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_path="panmat_exp_aug_2023"
file_prefix="my_vcf"
MAT="sars_8M_annotated_nextclade.panman"
full_file_path=$(realpath "$file_path")
fastq_file=$(find ${file_path} -name *fastq* | head -n 1)

##minimap2 -a --sam-hit-only -2 -x sr ./NC_045512v2.fa ${file_path}/${file_prefix}_R1+R2.fastq -t 32 -o ${file_path}/${file_prefix}_alignment.sam 

# Run C-WAP
cd ./src/C-WAP
source run.sh $full_file_path

# Run Freyja
cd ../Freyja
source run.sh $full_file_path
cd ../../

#Run WEPP
python src/ivar_correction.py ${file_path}
samtools view -h -o ${file_path}/${file_prefix}_alignment.sam ${file_path}/resorted.bam 
wbe sam2PB -i ${MAT} -v ${file_prefix} -s ${file_prefix}_alignment.sam -o ${file_path}
wbe initial_filter -T 56 -i ${MAT} -v ${file_prefix} -o ${file_path}
#wbe post_filter -T 56 -i ${MAT} -v ${file_prefix} -o ${file_path}
# gdb --args wbe detectPeaks -T 32 -i ${MAT} -v ${file_prefix} -o ${file_path}

#usher_to_taxonium -i debug/gisaidAndPublic.2022-02-08.masked.pb.gz -o debug/gisaidAndPublic.2022-02-08.masked.jsonl.gz --name_internal_nodes --clade_types nextstrain,pango