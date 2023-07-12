export PATH=$PATH:$PWD/build
cd build
make -j
cd ../
matUtils place_read -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l B.1.160,B.1.177.7,B.1.429,P.1,B.42,R.1,B.33 -d 0.2,0.15,0.15,0.15,0.1,0.05,0.2 -v my_vcf -r 150 -w 20 -e 0 -s 100 -f test/NC_045512v2.fa

#Estimating using regression based approach
echo "PYTHON-FREYJA"
python regression_abundance_estimate.py

#Getting the clades from haplotypes using USHeR
echo "USHeR RESULTS"
usher -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf_haplotypes.vcf --no-add
echo "Haplotype,Clade" > my_vcf_hap_clade.csv
awk '{ print $1 "," $3}' clades.txt >> my_vcf_hap_clade.csv

echo "POST PROCESSING"
matUtils post_processing -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf

#Estimating using EM approach
echo "PYTHON-EM"
Rscript abundance.r 

#cp my_vcf_reads_freyja.* ../Freyja/
#cd ../Freyja/
#conda activate freyja-env
#source run.sh 
#conda deactivate
#cd -

#usher -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf_haplotypes.vcf --no-add --write-parsimony-scores-per-node
#awk -F'\t' '$5 == "y"' parsimony-scores.tsv 