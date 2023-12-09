# f = open("swampy.vcf")
# lines = f.readlines()
# f.close()

# for j, line in enumerate(lines):
#     if line[0] != "#":
#         line = line.split("\t")
#         line = line[:-1] # Remove the last element, which is a "\n"
#         line = line[9:]
#         # Find the number of '0' s
#         num_zeros = line.count("0")
#         # Find the number of '1' s
#         num_ones = line.count("1")
#         # Find the number of '2' s
#         num_twos = line.count("2")
#         # Find the number of '3' s
#         num_threes = line.count("3")
#         #print(num_zeros, num_ones, num_twos, num_threes)
#         # Check if there's anything else other than 0, 1, 2, 3, and if so print what it is
#         for i in line:
#             if i not in ["0", "1", "2", "3"]:
#                 print("Something else: ", i, "in line number: ", j)
#                 print(line)

# print("Done")

f = open("redone_swampy.vcf")
lines = f.readlines()
f.close()

# first_read_mutations = []

# # Print out the mutations that are present in the first read
# for j, line in enumerate(lines):
#     if line[0] != "#":
#         line = line.split("\t")
#         line = line[:-1] # Remove the last element, which is a "\n"
#         if line[9] != "0":
#             first_read_mutations.append(j)
#             print("Mutation in first read: ", line[2] + ":" + line[9])


f2 = open("redone_swampy_parallel.vcf")
lines2 = f2.readlines()
f2.close()

# Check if every line is the same in both files
for i in range(len(lines)):
    if lines[i] != lines2[i]:
        print("Line number: ", i, " is different")
        print(lines[i][:100])
        print(lines2[i][:100])
        break
    else:
        print("Line number: ", i, " is the same")