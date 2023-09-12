import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

if len(sys.argv) != 2:
    print("Usage: python plotting_data.py <input_vcf>")
    sys.exit(1)

input_vcf = sys.argv[1]
f = open(input_vcf, 'r')
lines = f.readlines()
f.close()

# Remove header lines
lines = [line.strip() for line in lines]
lines = [line for line in lines if not line.startswith('#')]

line_length = len(lines[0].split('\t'))

# for each line, split on tabs
lines = [line.split('\t')[9:] for line in lines]


data = np.array(lines, dtype=int)

# Find number of non-zero elements in each column
data = np.count_nonzero(data, axis=0)

# Find the number of reads with each number of mutations
num_mutations = np.bincount(data)
num_reads = np.arange(len(num_mutations))

print(num_mutations)
print(num_reads)

# Plot a bar plot of number of reads vs number of mutations
plt.bar(num_reads, num_mutations)
plt.xlabel('Number of mutations')
plt.ylabel('Number of reads')
plt.show()
plt.savefig('mutation_counts.png')