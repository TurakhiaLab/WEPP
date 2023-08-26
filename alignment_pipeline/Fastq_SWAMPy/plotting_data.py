import numpy as np
import matplotlib.pyplot as plt

input_vcf = "my_indel_output.vcf"
f = open(input_vcf, 'r')
lines = f.readlines()
f.close()

# Remove header lines
lines = [line for line in lines if not line.startswith('#')]
# for each line, split on tabs
lines = [line.split('\t')[9:] for line in lines]

# Convert to numpy array, integers
lines = np.array(lines, dtype=int)

# Find the number of non-zero entries in each column (sample)
num_non_zero = np.count_nonzero(lines, axis=0)

# Plot histogram
plt.hist(num_non_zero, bins=100)
plt.show()
