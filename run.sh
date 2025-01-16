#!/bin/bash
set -e

#This script assumes that you have WBE, SWAMPy, and C-WAP setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_path="data/manuscript_swampy_dec_2022"
#file_path_cmp=
file_prefix="my_vcf"
#file_prefix_cmp=
#MAT="updated_gisaidAndPublic.2023-12-15.masked.pb.gz"
#MAT="pruned_public-2023-12-25.all.masked.pb.gz"
MAT="pruned_public-2023-12-25.all.masked.pb.gz"
#SRA="SRR29616810"
full_file_path=$(realpath "$file_path")

## Create Directory and Download reads
#mkdir -p ${file_path}
#cd ${file_path}
#cp ../${MAT} .
#cp ../test/NC_045512v2.fa .
#cd ../../sratoolkit.3.1.1-ubuntu64/
#prefetch ${SRA}
#fasterq-dump ${SRA}
##gzip ${SRA}.fastq
#gzip ${SRA}_1.fastq
#gzip ${SRA}_2.fastq
#rm -rf ${SRA}
##mv ${SRA}.fastq.gz ../SARS2-WBE/${file_path}/
#mv ${SRA}_1.fastq.gz ../SARS2-WBE/${file_path}/
#mv ${SRA}_2.fastq.gz ../SARS2-WBE/${file_path}/
#cd ../SARS2-WBE
#
#
# Conver Freyja files for Quick run
##sed -i 's/NC_045512v2/NC_045512.2/g' ${file_path}/my_vcf_alignment.sam 
##samtools sort -O BAM -o ${file_path}/resorted.bam ${file_path}/my_vcf_alignment.sam
#cd ./src/Freyja
#source run.sh $full_file_path
#cd -

#### Simulating reads from SWAMPy
##conda activate SWAMPy
##source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_samples.fa ${file_path}/${file_prefix}_samples.tsv ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path} 800000 149
##conda deactivate
#
# Run C-WAP
# cd ./src/C-WAP
# source run.sh $full_file_path

# Run Freyja
cd ./src/Freyja
source run.sh $full_file_path

## Run WEPP Pipeline
cd ../../
python src/WBE/ivar_correction.py ${file_path}
samtools view -h -o ${file_path}/${file_prefix}_alignment.sam ${file_path}/resorted.bam 

## Selecting subset of reads if needed
#mv ${file_path}/${file_prefix}_alignment.sam ${file_path}/${file_prefix}_alignment.sam_orig
#python src/WBE/select_subset_reads.py ${file_path}/${file_prefix}_alignment.sam_orig ${file_path}/${file_prefix}_alignment.sam

wbe sam2PB -v ${file_prefix} -f NC_045512v2.fa -s ${file_prefix}_alignment.sam -o ${file_path}
wbe detectPeaks -T 56 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}


##Analysis Scripts
#wbe analyzePeaks -T 1 -i ${MAT} -r ${MAT} -v ${file_prefix} -w ${file_prefix_cmp} -o ${file_path} -p ${file_path_cmp} -f NC_045512v2.fa 
#usher_to_taxonium -i debug/gisaidAndPublic.2022-02-08.masked.pb.gz -o debug/gisaidAndPublic.2022-02-08.masked.jsonl.gz --name_internal_nodes --clade_types nextstrain,pango
