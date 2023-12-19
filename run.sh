#!/bin/bash

#This script assumes that you have WBE, SWAMPy, and Freyja setup in parallel directories 
export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#Variables
file_prefix="my_vcf"
MAT="public-2021-05-31.all.masked.nextclade.pangolin.pb"
REF="test/NC_045512v2.fa"
file_path="output_files"

#Setting up directory
rm -r ${file_path}
mkdir -p ${file_path}
cp ${MAT} ${file_path}
cp ${REF} ${file_path}

#HAPLOTYPE SELECTION
##wbe selectHaplotypes -i public-2023-08-17.all.masked.nextclade.pangolin.pb -f NC_045512v2.fa -l AY.100,B.1.1.529,B.1.2,B.1.526,B.1.582,BA.1.18,XBB.1.5,XBB.1.9.1,P.1 -d 0.1,0.1,0.15,0.15,0.2,0.05,0.1,0.05,0.1 -v ${file_prefix} -w 20
wbe selectHaplotypes -i ${MAT} -f NC_045512v2.fa -l B.1.160,B.1.177.7,B.1.429,P.1,B.42,R.1,B.33 -d 0.2,0.15,0.15,0.15,0.1,0.05,0.2 -v ${file_prefix} -w 20 -o ${file_path}
python src/WBE/haplotype_abundance.py ${file_prefix} ${file_path}

#SWAMPy + ALIGNMENT
conda activate SWAMPy
source src/WBE/swampy_align.sh ${file_path}/${file_prefix}_samples.fa ${file_path}/${file_prefix}_samples.tsv ${file_path}/NC_045512v2.fa ${file_prefix} ${file_path} 200000 150
conda deactivate
wbe sam2VCF -v ${file_prefix} -f NC_045512v2.fa -s ${file_prefix}_alignment.sam -o ${file_path}

#FREYJA
rm ../Freyja/my_output_latest.txt ../Freyja/${file_prefix}_reads_freyja.*
cp ${file_path}/${file_prefix}_reads_freyja.* ../Freyja/
cd ../Freyja
conda activate freyja-env
#make dev
#mkdir -p data
#freyja update --buildlocal --outdir data
freyja demix ${file_prefix}_reads_freyja.vcf ${file_prefix}_reads_freyja.depth --barcodes data/usher_barcodes.csv --meta data/curated_lineages.json --output my_output_latest.txt
conda deactivate
cd -
python src/WBE/freyja_correct_format.py my_output_latest.txt ${file_prefix} ${file_path}

#FILTERING LINEAGES
wbe filterLineages -T 48 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
echo -e "\nFREYJA - LINEAGE FILTER"
python src/WBE/peaks_filtering.py ${file_prefix} ${file_path}

#DETECTING PEAKS
wbe detectPeaks -T 48 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
echo -e "\nFREYJA - PEAKS FILTER"
python src/WBE/peaks_filtering.py ${file_prefix} ${file_path}

#CALCULATING MUTATION DISTANCE
wbe refinePeaks -T 48 -i ${MAT} -v ${file_prefix} -f NC_045512v2.fa -o ${file_path}
