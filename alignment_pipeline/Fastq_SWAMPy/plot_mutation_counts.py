f = open("mutation_counts.txt", 'r')
lines = f.readlines()
f.close()

lines = [line.strip() for line in lines]
num_reads = [int(line.split(' ')[0]) for line in lines]
num_mutations = [int(line.split(' ')[1]) for line in lines]

import matplotlib.pyplot as plt

# Plot a bar plot of number of reads vs number of mutations
plt.bar(num_mutations, num_reads)
plt.xlabel('Number of mutations')
plt.ylabel('Number of reads')
plt.show()