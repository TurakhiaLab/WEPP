import sys
import random

#Check arguments
if len(sys.argv) != 3:
    print("USAGE: python select_subset_reads.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Open input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Keep the first three lines
    for _ in range(3):
        line = infile.readline()
        outfile.write(line)
    
    # Read the remaining lines into a list
    remaining_lines = infile.readlines()

# Shuffle the remaining lines and select 1M
random.shuffle(remaining_lines)
selected_lines = remaining_lines[:1000000]

# Write the selected lines to the output file
with open(output_file, 'a') as outfile:
    outfile.writelines(selected_lines)

