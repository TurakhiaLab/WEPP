#!/bin/bash

set -e

#This script assumes that you have WBE, SWAMPy, and Freyja setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_path="houston_as_nov_29_2021"
file_prefix="H_as_nov_29"
#MAT="updated_gisaidAndPublic.2022-05-15.masked.pb"
MAT="gisaidAndPublic.2021-11-15.masked.pb.gz"
#REF_MAT="gisaidAndPublic.2022-05-01.masked.pb.gz"

##cp ${file_path}/${file_prefix}_alignment.sam ../Freyja/${file_prefix}_alignment.sam
cp ${file_path}/resorted.bam src/Freyja/
cd src/Freyja
source run.sh
cd -


#wbe analyzePeaks -T 1 -i ${MAT} -r ${MAT} -v ${file_prefix} -w ${file_prefix_cmp} -o ${file_path} -p ${file_path_cmp} -f NC_045512v2.fa 

#### Simulating reads from SWAMPy
##conda activate SWAMPy
##source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_samples.fa ${file_path}/${file_prefix}_samples.tsv ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path} 800000 149
##conda deactivate
##
##mv ${file_path}/${file_prefix}_R1.fastq ${file_path}/${file_prefix}_R1.orig_fastq
##mv ${file_path}/${file_prefix}_R2.fastq ${file_path}/${file_prefix}_R2.orig_fastq
##gzip ${file_path}/${file_prefix}_R1+R2.fastq

## Running C-WAP
#fastq_path=$PWD/$file_path
#cd ../C-WAP
#source run.sh $fastq_path
#
## Prepare files for Freyja
#cd ../Freyja
#source run.sh
#
## Run WEPP Pipeline
#cd ../SARS2-WBE
#samtools fastq ${file_path}/resorted.bam > ${file_path}/${file_prefix}_reads.fastq
#
#source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_reads.fastq ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path}
#
### Prepare Files
##mv ${file_path}/${file_prefix}_alignment.sam ${file_path}/${file_prefix}_alignment.sam_orig
##python select_subset_reads.py ${file_path}/${file_prefix}_alignment.sam_orig ${file_path}/${file_prefix}_alignment.sam
##
##cp ${file_path}/${file_prefix}_alignment.sam_orig ../Freyja/${file_prefix}_alignment.sam
#####cp ${file_path}/resorted.bam ../Freyja/
##cd ../Freyja
##source run.sh
##cd -
#
#wbe sam2PB -v ${file_prefix} -f NC_045512v2.fa -s ${file_prefix}_alignment.sam -o ${file_path}

#DETECTING PEAKS
wbe detectPeaks -T 56 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
##wbe detectPeaks -T 56 -i ${MAT} -r ${REF_MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
#### gdb --args wbe detectPeaks -T 32 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}

#CALCULATING MUTATION DISTANCE
# wbe refinePeaks -T 32 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}

#FREYJA
# rm ../Freyja/my_output_latest.txt ../Freyja/${file_prefix}_reads_freyja.*
# cp ${file_path}/${file_prefix}_reads_freyja.* ../Freyja/
# cd ../Freyja
# conda activate freyja-env
# make dev
# mkdir -p data
# freyja update --buildlocal --outdir data
# freyja demix ${file_prefix}_reads_freyja.vcf ${file_prefix}_reads_freyja.depth --barcodes data/usher_barcodes.csv --meta data/curated_lineages.json --output my_output_latest.txt
# conda deactivate
# cd -
# python src/WBE/freyja_correct_format.py my_output_latest.txt ${file_prefix} ${file_path}


#usher_to_taxonium -i debug/gisaidAndPublic.2022-02-08.masked.pb.gz -o debug/gisaidAndPublic.2022-02-08.masked.jsonl.gz --name_internal_nodes --clade_types nextstrain,pango