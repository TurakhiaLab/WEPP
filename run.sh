#!/bin/bash

#This script assumes that you have WBE, SWAMPy, and Freyja setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_prefix="my_vcf"
MAT="public-2023-08-17.all.masked.nextclade.pangolin.pb"
#MAT="public-2023-12-31.all.masked.pb.gz"
REF="test/NC_045512v2.fa"
file_path="output_files"

##Setting up directory
#rm -r ${file_path}
#mkdir -p ${file_path}
#cp ${MAT} ${file_path}
#cp ${REF} ${file_path}
#
##HAPLOTYPE SELECTION
#wbe selectHaplotypes -i ${MAT} -f NC_045512v2.fa -l BA.2.10,EG.5.1,HV.1,XBB.1,XBB.1.5.1,XBB.1.9.1,XBB.1.16.1,XBB.2.3.2 -d 0.03,0.2,0.12,0.05,0.02,0.3,0.25,0.03 -v ${file_prefix} -w 100 -o ${file_path}
#python src/WBE/haplotype_abundance.py ${file_prefix} ${file_path}
#
##SWAMPy + ALIGNMENT
#conda activate SWAMPy
#source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_samples.fa ${file_path}/${file_prefix}_samples.tsv ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path} 800000 116
#conda deactivate
#wbe sam2VCF -v ${file_prefix} -f NC_045512v2.fa -s ${file_prefix}_alignment.sam -o ${file_path}
#
##FREYJA
#rm ../Freyja/my_output_latest.txt ../Freyja/${file_prefix}_reads_freyja.*
#cp ${file_path}/${file_prefix}_reads_freyja.* ../Freyja/
#cd ../Freyja
#conda activate freyja-env
#make dev
#mkdir -p data
#freyja update --buildlocal --outdir data
#freyja demix ${file_prefix}_reads_freyja.vcf ${file_prefix}_reads_freyja.depth --barcodes data/usher_barcodes.csv --meta data/curated_lineages.json --output my_output_latest.txt
#conda deactivate
#cd -
#python src/WBE/freyja_correct_format.py my_output_latest.txt ${file_prefix} ${file_path}

#FILTERING LINEAGES
wbe filterLineages -T 48 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}

#DETECTING PEAKS
wbe detectPeaks -T 48 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path} -p BA.2.10,EG.5.1,HV.1,XBB.1,XBB.1.5.1,XBB.1.9.1,XBB.1.16.1,XBB.2.3.2

#CALCULATING MUTATION DISTANCE
wbe refinePeaks -T 48 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
