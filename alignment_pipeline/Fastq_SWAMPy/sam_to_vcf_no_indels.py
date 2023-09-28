import sys
import re

if len(sys.argv) != 4:
    print("Usage: python sam_to_vcf_no_indels.py <sam_file> <reference_file> <vcf_file>")
    sys.exit(1)

sam_file = sys.argv[1]
reference_file = sys.argv[2]
vcf_file = sys.argv[3]

ref_genome = ""
with open(reference_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            continue
        ref_genome += line.strip()

def parse_cigar_and_mdz(cigar, mdz, seq, pos):
    mismatches = []
    read_pos = 0
    ref_pos = pos
    
    mdz_tokens = iter(re.findall(r'(\d+|\^[A-Z]+|[A-Z])', mdz))
    cigar_tokens = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    for count, operation in cigar_tokens:
        count = int(count)
        
        if operation == 'M':
            for _ in range(count):
                md_token = next(mdz_tokens, '')
                if md_token.isdigit():
                    read_pos += int(md_token)
                    ref_pos += int(md_token)
                elif md_token.startswith("^"):
                    continue  # Ignore deletions
                else:
                    ref_base = ref_genome[ref_pos - 1]
                    alt_base = seq[read_pos]
                    mismatches.append((ref_pos, ref_base, alt_base))
                    read_pos += 1
                    ref_pos += 1
        elif operation == 'I':
            read_pos += count  # Skip insertions
        elif operation == 'D':
            ref_pos += count  # Skip deletions
        else:
            # For other operations like N, S, H, etc., do not alter read_pos or ref_pos
            continue

    return mismatches

def parse_sam_line(line):
    parts = line.strip().split("\t")
    read_name, flag, ref_name, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = parts[:11]
    tags = dict((tag[:2], tag[5:]) for tag in parts[11:] if tag.startswith("MD:Z"))
    return read_name, ref_name, int(pos), seq, tags.get("MD", ""), cigar

variants_by_pos = {}
reads = set()

with open(sam_file, 'r') as f:
    for line in f:
        if line.startswith('@'):
            continue
        read_name, ref_name, pos, seq, mdz, cigar = parse_sam_line(line)
        mismatches = parse_cigar_and_mdz(cigar, mdz, seq, pos)
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
    out_lines.append("\t".join([ref_name, str(position), ".", ref_base, alt_base, ".", ".", ".", ".", *genotype_data]))

with open(vcf_file, 'w') as f:
    f.write("\n".join(out_lines))
