export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

###wbe selectHaplotypes -i public-2023-08-17.all.masked.nextclade.pangolin.pb -l AY.100,B.1.1.529,B.1.2,B.1.526,B.1.582,BA.1.18,XBB.1.5,XBB.1.9.1,P.1 -d 0.1,0.1,0.15,0.15,0.2,0.05,0.1,0.05,0.1 -v my_vcf -w 20

wbe filterLineages -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf -f test/NC_045512v2.fa
echo -e "\nFREYJA - LINEAGE FILTER"
python src/WBE/peaks_filtering.py my_vcf
mv my_vcf_haplotype_abundance.csv my_vcf_haplotypeLINEAGE_abundance.csv
mv my_vcf_haplotypes.vcf my_vcf_haplotypesLINEAGE.vcf
mv my_vcf_peaks.vcf my_vcf_peaksLINEAGE.vcf
mv my_vcf_barcode.csv my_vcf_barcodeLINEAGE.csv

wbe detectPeaks -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf -f test/NC_045512v2.fa
echo -e "\nFREYJA - PEAKS FILTER"
python src/WBE/peaks_filtering.py my_vcf
wbe refinePeaks -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf -f test/NC_045512v2.fa