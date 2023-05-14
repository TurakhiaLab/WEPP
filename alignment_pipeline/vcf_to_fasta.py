import sys

# Check that correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python vcf_to_fasta.py <vcf_file> <reference_fasta> <output_fasta>")
    sys.exit(1)

reference = sys.argv[2]
vcf_file = sys.argv[1]
output_file = sys.argv[3]

# Read reference fasta file into memory
ref_fasta = open(reference, "r")
ref_lines = ref_fasta.readlines()
ref_fasta.close()
ref_seq = "".join(ref_lines[1:]).replace("\n", "")

# Read VCF file into memory
vcf_file = open(vcf_file, "r")
vcf_lines = vcf_file.readlines()
vcf_file.close()

# Get list of sample names from VCF file header
header_line = vcf_lines[3].strip()
headers = header_line.split("\t")
sample_names = headers[9:]

# Create dictionary to store mutated sequences
mutated_seqs = {}

# Loop through each sample in VCF file
for sample_index in range(len(sample_names)):
    sample_name = sample_names[sample_index]
    mutated_seq = list(ref_seq)

    # Loop through each variant row in VCF file for the current sample
    for vcf_line in vcf_lines[4:]:
        vcf_fields = vcf_line.strip().split("\t")
        chrom = vcf_fields[0]
        pos = int(vcf_fields[1])
        ref_base = vcf_fields[3]
        alt_base = vcf_fields[4]
        gt = vcf_fields[9+sample_index]

        # Check if sample is mutated at current position
        if gt == "1":
            # Check that reference base at current position matches expected base
            if mutated_seq[pos-1] == ref_base:
                mutated_seq[pos-1] = alt_base
            else:
                print(f"Reference base ({ref_base}) does not match at position {pos} for sample {sample_name}")
        
    # Create mutated sequence and store in dictionary
    mutated_seq_str = "".join(mutated_seq)
    mutated_seqs[sample_name] = mutated_seq_str

# Write mutated sequences to Fasta file
mutated_fasta = open(output_file, "w")
for sample_name, mutated_seq in mutated_seqs.items():
    mutated_fasta.write(f">{sample_name}\n")
    mutated_fasta.write(f"{mutated_seq}\n")
mutated_fasta.close()
