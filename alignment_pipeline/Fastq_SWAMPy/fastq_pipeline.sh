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

# Run the pipeline
echo "Running the pipeline"

echo "Converting FASTQ to FASTA"
time python fastq_to_fasta.py $input_fastq reads.fasta
echo ""

echo "Aligning reads to reference"

bowtie2-build $reference_fasta intermediate_files/ref_index
bowtie2 -x intermediate_files/ref_index -f -U reads.fasta --local -S intermediate_files/alignment.sam

cat intermediate_files/alignment.sam | python process_sam.py > intermediate_files/alignment_forward.sam
sed '1s/^\xEF\xBB\xBF//' intermediate_files/alignment_forward.sam > intermediate_files/alignment.sam
python new_alignment_positions.py
samtools view -S -b intermediate_files/alignment_modified.sam > intermediate_files/alignment.bam
samtools sort intermediate_files/alignment.bam -o intermediate_files/sorted_alignment.bam
samtools index intermediate_files/sorted_alignment.bam
samtools fasta intermediate_files/sorted_alignment.bam > reads_aligned.fasta
# samtools view -F 4 intermediate_files/sorted_alignment.bam | samtools fasta - > reads_aligned.fasta

# bwa index $reference_fasta
# bwa mem $reference_fasta reads.fasta > intermediate_files/alignment.sam
# python new_alignment_positions.py
# samtools view -S -b intermediate_files/alignment_modified.sam > intermediate_files/alignment.bam
# samtools sort intermediate_files/alignment.bam -o intermediate_files/sorted_alignment.bam
# samtools index intermediate_files/sorted_alignment.bam
# samtools fasta intermediate_files/sorted_alignment.bam > reads_aligned.fasta

echo ""

echo "Converting FASTA to VCF"
python process_aligned_fasta.py
python indel_fasta_to_vcf.py
time ./fasta_to_vcf reads_aligned.fasta $reference_fasta $output_vcf
echo ""

# Print location where output VCF was written
echo "Output VCF written to $output_vcf"
echo ""

# Remove the intermediate files
# rm -f reads.fasta reads_aligned.fasta

