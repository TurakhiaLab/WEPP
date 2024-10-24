#!/bin/bash

#/home/swalia@AD.UCSD.EDU/work/panman/build/panman/viralmsa_8M.panman

set -e

#This script assumes that you have WBE, SWAMPy, and Freyja setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_prefix="my_vcf"
# file_prefix="PL_2023_03_01"

# ont (15)
# MAT="updated_1_public-2023-04-10.all.masked.pb.gz"
# file_path="golden_mixture1_v41_control"

# ont (8)/6
# MAT="updated_public-2023-04-10.all.masked.pb.gz"
# file_path="golden_mixture6_v41_control"

# ont (8)/5
# MAT="updated_1_public-2023-04-10.all.masked.pb.gz"
# file_path="mixture-05"

# illumina
MAT="sars_8M.panman"
#MAT="sars_20_viralmsa.panman"
file_path="single_hap_simulated"

# point loma
# MAT="updated_public-2023-04-10.all.masked.pb.gz"
# file_path="point_loma"

#REF="test/NC_045512v2.fa"

##Setting up directory
#rm -r ${file_path}
#mkdir -p ${file_path}
#cp ${MAT} ${file_path}
#cp ${REF} ${file_path}

##SWAMPy + ALIGNMENT
#conda activate SWAMPy
#source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_samples.fa ${file_path}/${file_prefix}_samples.tsv ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path} 800000 149
# source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_reads.fastq ${REF} ${file_prefix} ${file_path}
#conda deactivate

#wbe sam2PB -i ${MAT} -v ${file_prefix} -s ${file_prefix}_alignment.sam -o ${file_path}
#gdb --args wbe sam2PB -v ${file_prefix} -s ${file_prefix}_alignment.sam -o ${file_path}

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

#DETECTING PEAKS
wbe initial_filter -T 56 -i ${MAT} -v ${file_prefix} -o ${file_path}
#wbe post_filter -T 56 -i ${MAT} -v ${file_prefix} -o ${file_path}
# gdb --args wbe detectPeaks -T 32 -i ${MAT} -v ${file_prefix} -o ${file_path}

#usher_to_taxonium -i debug/gisaidAndPublic.2022-02-08.masked.pb.gz -o debug/gisaidAndPublic.2022-02-08.masked.jsonl.gz --name_internal_nodes --clade_types nextstrain,pango