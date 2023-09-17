# This uses bowtie2 for local alignment

# Take in 3 required arguments + 1 optional argument
# 1. Path to the reference FASTA
# 2. Path to the input fastq
# 3. Path to output VCF
# 4. (Optional) --compile: whether to compile the C++ code or not

# Check that the correct number of arguments were passed in
if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference_fasta> <input_fastq> <output_vcf> [--compile]"
    exit 1
fi

# Exit if any command fails
set -e

# Check if the optional argument was passed in
if [ $# -eq 4 ]; then
    if [ "$4" == "--compile" ]; then
        echo "Compiling C++ code"
        g++ -o fasta_to_vcf fasta_to_vcf.cpp
        g++ -o sort_vcf sort_vcf.cpp
        g++ -o group_vcf group_vcf.cpp
        g++ -o generate_freyja_files generate_freyja_files.cpp
    else
        echo "Usage: $0 <reference_fasta> <input fastq> <output_vcf> [--compile]"
        exit 1
    fi
fi

# Check that the files corresponding to the required input arguments exist
if [ ! -f "$1" ]; then
    echo "Reference FASTA file does not exist"
    exit 1
fi
if [ ! -f "$2" ]; then
    echo "Input Fastq file does not exist"
    exit 1
fi

# Remove the output VCF if it already exists
rm -f $3

# Remove intermediate files if they already exist
rm -f reads.fasta reads_aligned.fasta
mkdir -p intermediate_files

# Read in the arguments
reference_fasta=$1
input_fastq=$2
output_vcf=$3

output_vcf_no_ext=${output_vcf%.*}
output_vcf_freyja=${output_vcf_no_ext}_freyja.vcf
output_vcf_depth=${output_vcf_no_ext}_freyja.depth

bowtie2-build $reference_fasta intermediate_files/ref_index
bowtie2 -x intermediate_files/ref_index -U $input_fastq -S intermediate_files/alignment.sam
python process_rc.py $input_fastq
# Append _processed to the input fastq name
input_fastq_processed=${input_fastq%.*}_processed.fastq
bowtie2 -x intermediate_files/ref_index -U $input_fastq_processed -S intermediate_files/alignment.sam

# bwa index $reference_fasta
# bwa mem $reference_fasta $input_fastq > intermediate_files/alignment.sam
python new_alignment_positions.py
python new_vcf_generator.py intermediate_files/alignment_modified.sam my_test_output.vcf
./sort_vcf
./group_vcf my_test_output_2.vcf $output_vcf
./generate_freyja_files my_test_output_2.vcf $reference_fasta $output_vcf_freyja $output_vcf_depth
python plotting_data.py $output_vcf > mutation_counts_grouped.txt