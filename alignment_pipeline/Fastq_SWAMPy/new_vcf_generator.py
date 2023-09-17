import sys


ref_genome = ""
with open("NC_045512v2.fa", 'r') as f:
    for line in f:
        if line.startswith('>'):
            continue
        ref_genome += line.strip()

def parse_mdz(mdz, seq, pos):
    mismatches = []
    
    # split the MD string into tokens, where each token is either a number of matching bases or a mismatch
    tokens = []
    token = ""
    for char in mdz:
        if char.isdigit():
            token += char
        else:
            if token:
                tokens.append(int(token))
                token = ""
            tokens.append(char)


    read_pos = 1
    ref_pos = pos
   
    for token in tokens:
        if token == "^": # skip over the next character, which is the reference base
            ref_pos += 1
            read_pos += 1
            del tokens[tokens.index(token) + 1] # Delete the next token, which is the reference base
        
        elif isinstance(token, int): # matching bases
            read_pos += token
            ref_pos += token
        else: # mismatch
            ref_base = ref_genome[ref_pos-1]
            if read_pos > len(seq):
                break
            alt_base = seq[read_pos - 1]
            if ref_base != alt_base:
                mismatches.append((ref_pos, ref_base, alt_base))
            read_pos += 1
            ref_pos += 1

    return mismatches

def parse_sam_line(line):
    parts = line.strip().split("\t")
    read_name, flag, ref_name, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = parts[:11]
    tags = dict((tag[:2], tag[5:]) for tag in parts[11:] if tag.startswith("MD:Z"))
    return read_name, ref_name, int(pos), seq, tags.get("MD", "")

variants_by_pos = {}
reads = set()

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('@'): # Skip header lines
            continue
        
        read_name, ref_name, pos, seq, mdz = parse_sam_line(line)
        mismatches = parse_mdz(mdz, seq, pos)
        
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


with open(sys.argv[2], 'w') as f:
    f.write("\n".join(out_lines))
