# This uses bowtie2 for local alignment

# Take in 3 required arguments + 1 optional argument
# 1. Path to the reference FASTA
# 2. Path to the input fastq
# 3. Path to output VCF
# 4. (Optional) --compile: whether to compile the C++ code or not

# Check that the correct number of arguments were passed in
if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference_fasta> <input_fastq> <output_vcf> [--compile] [--no-indels]"
    exit 1
fi

# Exit if any command fails
set -e

# Create a boolean variable to indicate whether indels should be included
include_indels=true

# Check if the optional argument was passed in
if [ $# -gt 3 ]; then

    if [ "$4" == "--compile" ] || [ "$5" == "--compile" ]; then
        echo "Compiling C++ code"
        g++ -o sort_vcf sort_vcf.cpp
        g++ -o group_vcf group_vcf.cpp
        g++ -o generate_freyja_files generate_freyja_files.cpp
    fi

    if [ "$4" == "--no-indels" ] || [ "$5" == "--no-indels" ]; then
        echo "Running without indels"
        include_indels=false
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
python process_rc.py intermediate_files/alignment.sam $input_fastq
input_fastq_processed=${input_fastq%.*}_processed.fastq
bowtie2 -x intermediate_files/ref_index -U $input_fastq_processed -S intermediate_files/alignment_2.sam
python new_alignment_positions.py intermediate_files/alignment_2.sam intermediate_files/alignment_modified.sam

# Check if include_indels is true
if [ "$include_indels" = true ]; then
    ./sam_to_vcf intermediate_files/alignment_modified.sam $reference_fasta intermediate_files/vcf_unsorted.vcf 4
else
    python sam_to_vcf_no_indels.py intermediate_files/alignment_modified.sam $reference_fasta intermediate_files/vcf_unsorted.vcf
fi

./sort_vcf intermediate_files/vcf_unsorted.vcf intermediate_files/vcf_sorted.vcf
./group_vcf intermediate_files/vcf_sorted.vcf $output_vcf
./generate_freyja_files intermediate_files/vcf_sorted.vcf $reference_fasta $output_vcf_freyja $output_vcf_depth $output_vcf

python plotting_data.py $output_vcf > mutation_counts_grouped.txt