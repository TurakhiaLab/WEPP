import sys

# Check that correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python fasta_to_vcf.py <input_fasta_file> <reference_fasta> <output_vcf>")
    sys.exit(1)

reference = sys.argv[2]
input_file = sys.argv[1]
output_file = sys.argv[3]

# Open the reference and read files
with open(reference, "r") as ref_file:
    ref_lines = ref_file.readlines()

num_samples = 0
aligned_reads = []

sample_names = []
with open("./intermediate_files/sample_names.txt","r") as f:
    sample_names = f.readline().split(",")
    f.close()

with open(input_file, "r") as read_file:
    read_lines = read_file.readlines()
    # Count the number of samples by counting the number of lines that start with ">"
    for line in read_lines:
        if line.startswith(">"):
            num_samples += 1
        else:
            aligned_reads.append(line.strip())

# Parse the reference sequence and chromosome name
ref_seq = ""
chromosome_name = ""
for line in ref_lines:
    if line.startswith(">"):
        chromosome_name = line.strip().lstrip(">")
    else:
        ref_seq += line.strip()

# Open the output VCF file for writing
with open(output_file, "w") as vcf_file:
    # Format it with boolean entires to indicate whether or not a given mutation is present in a sample
    vcf_file.write("##fileformat=VCFv4.2\n")
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")

    for i in range(num_samples):
        vcf_file.write(sample_names[i] + "\t")

    vcf_file.write("\n")

    # Iterate through the reference sequence and aligned fasta in parallel
    # Under the SampleX header, write a 1 if the reference and read sequences differ at that position, and a 0 otherwise
    for i in range(len(ref_seq)):
        if not all([ref_seq[i] == aligned_reads[k][i] for k in range(num_samples)]):
            mismatch_indices = [k for k in range(num_samples) if ref_seq[i] != aligned_reads[k][i]]
            mutations = set([aligned_reads[k][i] for k in mismatch_indices])

            for mutation in mutations:
                vcf_file.write(chromosome_name + "\t")
                vcf_file.write(str(i+1) + "\t")
                ref = ref_seq[i]
                alt = mutation
                id = ref + str(i+1) + alt
                vcf_file.write(id + "\t")
                vcf_file.write(ref + "\t")
                vcf_file.write(alt + "\t")
                vcf_file.write(".\t")
                vcf_file.write(".\t")
                vcf_file.write(".\t")
                vcf_file.write(".\t")


                for j in range(num_samples):
                    if aligned_reads[j][i] == mutation:
                        vcf_file.write("1\t")
                    else:
                        vcf_file.write("0\t")

                vcf_file.write("\n")

    # Remove the last line in the file, which is a newline
    vcf_file.seek(0, 2)
    vcf_file.seek(vcf_file.tell() - 2, 0)
    vcf_file.truncate()
    vcf_file.close()



