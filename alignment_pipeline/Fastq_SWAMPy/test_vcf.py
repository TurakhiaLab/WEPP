f = open("swampy.vcf")
lines = f.readlines()
f.close()

for j, line in enumerate(lines):
    if line[0] != "#":
        line = line.split("\t")
        line = line[:-1] # Remove the last element, which is a "\n"
        line = line[9:]
        # Find the number of '0' s
        num_zeros = line.count("0")
        # Find the number of '1' s
        num_ones = line.count("1")
        # Find the number of '2' s
        num_twos = line.count("2")
        # Find the number of '3' s
        num_threes = line.count("3")
        #print(num_zeros, num_ones, num_twos, num_threes)
        # Check if there's anything else other than 0, 1, 2, 3, and if so print what it is
        for i in line:
            if i not in ["0", "1", "2", "3"]:
                print("Something else: ", i, "in line number: ", j)
                print(line)

print("Done")