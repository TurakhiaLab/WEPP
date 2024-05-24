#!/bin/bash

MAT="updated_public-2023-04-10.all.masked.pb.gz"
# List of input FASTA files
fasta_files=("Control_15" "Control_17" "Control_19" "Control_50" "Control_23" "Control_48" "Control_51" "Control_62")

# Loop through each FASTA file
for fasta_file in "${fasta_files[@]}"
do
    # Perform alignment using MAFFT
    mafft --thread 10 --auto --keeplength --addfragments "./$fasta_file.fasta" "../test/NC_045512v2.fa" > "./aligned_$fasta_file"

    # Convert alignment to SAM format
    ../faToVcf "./aligned_$fasta_file" "./$fasta_file.vcf"

    # Add the new vcf to the tree
    usher -i ${MAT} -v "./$fasta_file.vcf" -o ${MAT}
done
