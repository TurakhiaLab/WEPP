export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

wbe filterLineages -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf -f test/NC_045512v2.fa
echo -e "\nFREYJA - LINEAGE FILTER"
python src/WBE/peaks_filtering.py my_vcf
wbe detectPeaks -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf -f test/NC_045512v2.fa
echo -e "\nFREYJA - PEAKS FILTER"
python src/WBE/peaks_filtering.py my_vcf
wbe refinePeaks -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf -f test/NC_045512v2.fa

###matUtils haplotype_pruning -i public-2021-05-31.all.masked.nextclade.pangolin.pb -w public-2023-08-17.all.masked.nextclade.pangolin_pruned.pb -s remove_sample.txt -d 5
#
###Estimating using regression based approach
#####Getting the clades from haplotypes using USHeR
#####usher -i public-2023-08-17.all.masked.nextclade.pangolin_pruned.pb -v my_vcf_haplotypes.vcf --no-add
####usher -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf_haplotypes.vcf --no-add
####echo "Haplotype,Clade" > my_vcf_hap_clade.csv
####awk -F'\t' '{print $1 "," $NF}' clades.txt >> my_vcf_hap_clade.csv
##