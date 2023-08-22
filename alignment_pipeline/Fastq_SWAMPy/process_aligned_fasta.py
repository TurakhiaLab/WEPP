f = open("reads_aligned.fasta")
lines = f.readlines()

lines_sequences = lines[1::2]
lines_sequences = [line.strip() for line in lines_sequences]

#  Find the indices for which lines_sequences[i] contains the a character other than A,C,G,T
indices = [i for i in range(len(lines_sequences)) if not all(c in "ACGT" for c in lines_sequences[i])]
print("Indices of reads with no ACTG:")
print(indices)
f.close()

# Find the first line which ends with _READ_0_149
for j in range(len(lines)):
    if lines[j].endswith("_READ_0_149\n"):
        break

# Remove this line all the way to the end of the file
lines = lines[:j]

# Write the lines to the same file
f = open("reads_aligned.fasta", "w")
f.writelines(lines)
f.close()