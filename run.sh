export PATH=$PATH:$PWD/build
cd build
make -j
cd ../
matUtils place_read -T 8 -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l B.1.1.198 -d 1 -v my_vcf -r 150 -w 20 -e 0 -s 2 -f test/NC_045512v2.fa
#grep "Clade:" read_info | cut -d ' ' -f 4 | sort -h | uniq


#matUtils place_read -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l B.1.1.117 -d 1 -v my_vcf -r 150 -w 20 -e 0 -s 2 -f test/NC_045512v2.fa
#vi -d  my_vcf_samples.vcf my_vcf_reads.vcf


#matUtils place_read -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l B.1.1.207 -d 1 -v my_vcf -r 150 -w 1500 -s 2 -e 0 -f test/NC_045512v2.fa
#matUtils extract -i public-2021-05-31.all.masked.nextclade.pangolin.pb -c 'B.1.1.207' -v ref_vcf.vcf
###vi -d my_vcf_samples.vcf my_vcf_reads.vcf
#vi -d ref_vcf.vcf my_vcf_samples.vcf
###matUtils place_read -i public-2021-05-31.all.masked.nextclade.pangolin.pb -l AP.1 -d 1  -v my_vcf.vcf 

#matUtils extract -i public-2021-05-31.all.masked.nextclade.pangolin.pb -v samples.vcf
