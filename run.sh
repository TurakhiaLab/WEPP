export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

#matUtils haplotype_pruning -i public-2023-08-17.all.masked.nextclade.pangolin.pb -w public-2023-08-17.all.masked.nextclade.pangolin_pruned.pb -s remove_sample.txt -d 5
matUtils place_read -T 48 -i public-2023-08-17.all.masked.nextclade.pangolin.pb -j public-2023-08-17.all.masked.nextclade.pangolin.pb -l AY.100,B.1.1.529,B.1.2,B.1.526,B.1.582,BA.1.18,XBB.1.5,XBB.1.9.1,P.1 -d 0.1,0.1,0.15,0.15,0.2,0.05,0.1,0.05,0.1 -v my_vcf -r 150 -w 20 -e 0 -s 100 -f test/NC_045512v2.fa

#Estimating using regression based approach
echo -e "\nPYTHON-FREYJA"
python regression_abundance_estimate.py

#Getting the clades from haplotypes using USHeR
#usher -i public-2023-08-17.all.masked.nextclade.pangolin_pruned.pb -v my_vcf_haplotypes.vcf --no-add
usher -i public-2023-08-17.all.masked.nextclade.pangolin.pb -v my_vcf_haplotypes.vcf --no-add
echo "Haplotype,Clade" > my_vcf_hap_clade.csv
awk -F'\t' '{print $1 "," $NF}' clades.txt >> my_vcf_hap_clade.csv

echo -e "\nPOST PROCESSING"
matUtils post_processing -T 48 -i public-2023-08-17.all.masked.nextclade.pangolin.pb -v my_vcf

#######################################
#Estimating using EM approach
#echo -e "\nPYTHON-EM"
#Rscript abundance.r 

#cp my_vcf_reads_freyja.* ../Freyja/
#cd ../Freyja/
#conda activate freyja-env
#source run.sh 
#conda deactivate
#cd -

#usher -i public-2023-08-17.all.masked.nextclade.pangolin.pb -v my_vcf_haplotypes.vcf --no-add --write-parsimony-scores-per-node
#awk -F'\t' '$5 == "y"' parsimony-scores.tsv 
