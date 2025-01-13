export PATH=$PATH:$PWD/build
cd build
make -j
cd ../
#matUtils mutation_paths public-2022-08-15.all.masked.nextclade.pangolin.pb.gz llm_mutation_leaves_aug2022_full_path.tsv
#matUtils mutation_paths public-2022-08-15.all.masked.nextclade.pangolin.pb.gz llm_mutations_leaves_aa_full_path_2022_08_15.tsv
matUtils mutation_paths gisaidAndPublic.2023-12-15.masked.pb.gz llm_mutations_leaves_aa_full_path_2023_12_15.tsv