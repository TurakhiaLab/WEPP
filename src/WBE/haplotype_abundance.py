import sys

def read_fasta(path: str) -> dict:
    """Reads fasta file and returns a dictionary of sequences with their count of appearances."""
    sequences = {}
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                sequences[name] = sequences.get(name, 0) + 1
    return sequences

def calculate_abundance(seqs: dict) -> dict:
    """Calculates the abundance of each sequence based on its appearances."""
    total_count = sum(seqs.values())
    abundance = {name: (count / total_count) * 100 for name, count in seqs.items()}
    return abundance

#Check arguments
if len(sys.argv) != 3:
    printf("USAGE: python abundance_generation.py <file_prefix> <file_path>")
    sys.exit(1)

# Paths
file_prefix = sys.argv[1]
file_path = sys.argv[2]

# Reading and calculating abundance
seqs = read_fasta(file_path + "/" + file_prefix + "_samples.fa")
seqs_with_abundance = calculate_abundance(seqs)

# Writing to file
with open(file_path + "/" + file_prefix + "_samples.tsv", 'w') as f:
    for name, abundance in seqs_with_abundance.items():
        f.write(f'{name}\t{abundance:.1f}\n')
