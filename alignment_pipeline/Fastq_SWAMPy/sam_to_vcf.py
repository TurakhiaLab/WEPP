import sys
import re

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python sam_to_vcf.py <sam_file> <reference_file> <vcf_file>")
    sys.exit(1)

reference_file = sys.argv[2]

ref_genome = ""
with open(reference_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            continue
        ref_genome += line.strip()

def parse_mdz(mdz, cigar, seq, pos):
    mismatches = []
    
    # split the MD string into tokens
    tokens = re.findall(r'(\d+|\^[A-Z]+|[A-Z])', mdz)

    read_pos = 0
    ref_pos = pos

    cigar_tokens = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    for token in tokens:
        if token.startswith("^"): # deletion
            del_len = len(token) - 1
            ref_pos += del_len
            mismatches.append((ref_pos, token[1:], "-"))
        elif token.isdigit(): # matching bases
            read_pos += int(token)
            ref_pos += int(token)
        else: # mismatch
            ref_base = ref_genome[ref_pos-1]
            alt_base = seq[read_pos]
            mismatches.append((ref_pos, ref_base, alt_base))
            read_pos += 1
            ref_pos += 1

    for count, operation in cigar_tokens:
        count = int(count)
        if operation == 'I': # Insertion
            alt_base = seq[read_pos:read_pos+count]
            mismatches.append((ref_pos, "-", alt_base))
            read_pos += count
        elif operation in 'DN': # Deletion
            ref_pos += count

    return mismatches

def parse_sam_line(line):
    parts = line.strip().split("\t")
    read_name, flag, ref_name, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = parts[:11]
    tags = dict((tag[:2], tag[5:]) for tag in parts[11:] if tag.startswith("MD:Z"))
    return read_name, ref_name, int(pos), seq, tags.get("MD", ""), cigar

variants_by_pos = {}
reads = set()

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('@'): # Skip header lines
            continue
        
        read_name, ref_name, pos, seq, mdz, cigar = parse_sam_line(line)
        mismatches = parse_mdz(mdz, cigar, seq, pos)
        
        reads.add(read_name)
        for position, ref_base, alt_base in mismatches:
            key = (ref_name, position, ref_base, alt_base)
            if key not in variants_by_pos:
                variants_by_pos[key] = set()
            variants_by_pos[key].add(read_name)

out_lines = []
out_lines.append("##fileformat=VCFv4.2")
out_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sorted(reads)))

for (ref_name, position, ref_base, alt_base), variant_reads in variants_by_pos.items():
    genotype_data = ["0" if read not in variant_reads else "1" for read in sorted(reads)]
    out_lines.append("\t".join([ref_name, str(position), ref_base + str(position) + alt_base , ref_base, alt_base, ".", ".", ".", ".", *genotype_data]))


with open(sys.argv[3], 'w') as f:
    f.write("\n".join(out_lines))
