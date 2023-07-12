export PATH=$PATH:$PWD/build
cd build
make -j
cd ../
matUtils place_read -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l B.1.160,B.1.177.7,B.1.429,P.1,B.42,R.1,B.33 -d 0.2,0.15,0.15,0.15,0.1,0.05,0.2 -v my_vcf -r 150 -w 20 -e 0 -s 100 -f test/NC_045512v2.fa

#Estimating using regression based approach
python regression_abundance_estimate.py
matUtils post_processing -T 48 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf

#Estimating using EM approach
Rscript abundance.r 

#cp my_vcf_reads_freyja.* ../Freyja/
#cd ../Freyja/
#conda activate freyja-env
#source run.sh 
#conda deactivate
#cd -

#usher -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v my_vcf_haplotypes.vcf --no-add --write-parsimony-scores-per-node
#awk -F'\t' '$5 == "y"' parsimony-scores.tsv 