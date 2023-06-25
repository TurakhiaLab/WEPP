# Take in 3 required arguments + 1 optional argument
# 1. Path to the reference FASTA
# 2. Path to the reads, VCF
# 3. Path to output VCF
# 4. (Optional) --compile: whether to compile the C++ code or not

# Check that the correct number of arguments were passed in
if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference_fasta> <reads_vcf> <output_vcf> [--compile]"
    exit 1
fi

# Exit if any command fails
set -e

# Check if the optional argument was passed in
if [ $# -eq 4 ]; then
    if [ "$4" == "--compile" ]; then
        echo "Compiling C++ code"
        g++ -o vcf_to_fasta vcf_to_fasta.cpp
        g++ -o local_alignment local_alignment.cpp
        g++ -o fasta_to_vcf fasta_to_vcf.cpp
    else
        echo "Usage: $0 <reference_fasta> <reads_vcf> <output_vcf> [--compile]"
        exit 1
    fi
fi

# Check that the files corresponding to the required input arguments exist
if [ ! -f "$1" ]; then
    echo "Reference FASTA file does not exist"
    exit 1
fi
if [ ! -f "$2" ]; then
    echo "Reads VCF file does not exist"
    exit 1
fi

# Remove the output VCF if it already exists
rm -f $3

# Read in the arguments
reference_fasta=$1
reads_vcf=$2
output_vcf=$3

# Run the pipeline
echo "Running the pipeline"

echo "Converting VCF to FASTA"
time ./vcf_to_fasta $reads_vcf $reference_fasta reads.fasta

echo "Aligning reads to reference"
time ./local_alignment $reference_fasta reads.fasta reads_aligned.fasta

echo "Converting FASTA to VCF"
time ./fasta_to_vcf reads_aligned.fasta $reference_fasta $output_vcf

# Remove the intermediate files
rm -f reads.fasta reads_aligned.fasta

# Print location where output VCF was written
echo "Output VCF written to $output_vcf"

