#!/bin/bash
set -e

#This script assumes that you have WBE, SWAMPy, and C-WAP setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_path="NWSS_request"
file_prefix="SRR30906372"
MAT="public-2024-11-05.all.masked.pb.gz"

## Copy files to Freyja for Quick run
###cp ${file_path}/${file_prefix}_alignment.sam ../Freyja/${file_prefix}_alignment.sam
#cp ${file_path}/resorted.bam src/Freyja/
#cd src/Freyja
#source run.sh
#cd -

#### Simulating reads from SWAMPy
##conda activate SWAMPy
##source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_samples.fa ${file_path}/${file_prefix}_samples.tsv ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path} 800000 149
##conda deactivate

# Run C-WAP
fastq_path=$PWD/$file_path
cd ../C-WAP
source run.sh $fastq_path

# Run Freyja
cd ../SARS2-WBE/src/Freyja
source run.sh

# Run WEPP Pipeline
cd ../../
samtools fastq ${file_path}/resorted.bam > ${file_path}/${file_prefix}_reads.fastq
source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_reads.fastq ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path}

### Selecting subset of reads if needed
##mv ${file_path}/${file_prefix}_alignment.sam ${file_path}/${file_prefix}_alignment.sam_orig
##python src/WBE/select_subset_reads.py ${file_path}/${file_prefix}_alignment.sam_orig ${file_path}/${file_prefix}_alignment.sam
##
##cp ${file_path}/${file_prefix}_alignment.sam_orig ./src/Freyja/${file_prefix}_alignment.sam
##cd ./src/Freyja
##source run.sh
##cd -

wbe sam2PB -v ${file_prefix} -f NC_045512v2.fa -s ${file_prefix}_alignment.sam -o ${file_path}
wbe detectPeaks -T 56 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
##wbe detectPeaks -T 56 -i ${MAT} -r ${REF_MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}

##Analysis Scripts
#wbe analyzePeaks -T 1 -i ${MAT} -r ${MAT} -v ${file_prefix} -w ${file_prefix_cmp} -o ${file_path} -p ${file_path_cmp} -f NC_045512v2.fa 
#usher_to_taxonium -i debug/gisaidAndPublic.2022-02-08.masked.pb.gz -o debug/gisaidAndPublic.2022-02-08.masked.jsonl.gz --name_internal_nodes --clade_types nextstrain,pango