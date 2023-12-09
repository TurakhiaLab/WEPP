f = open("redone_swampy.vcf")
lines = f.readlines()
f.close()

read_dict = {} # key: read name, value: list of mutations
reads_list = []

for line in lines:
    if line.startswith("#CHROM"):
        line = line.strip().split("\t")
        reads_list = line[9:]
        break

for line in lines:
    line = line.strip()
    if line[0] == "#":
        continue
    line = line.split("\t")
    mutations = line[2].split(",")
    for mut in mutations:
        for read in reads_list:
            if read not in read_dict:
                read_dict[read] = []
            read_dict[read].append(mut)

for read in read_dict:
    numbers = read.split("_")
    start_position = int(numbers[-3])
    end_position = int(numbers[-2])
    mdz = "MD:Z:"
    previous_position = start_position
    for mut in read_dict[read]:
        ref_base = mut[0]
        alt_base = mut[-1]
        position = int(mut[1:-1])
        mdz += str(position - previous_position) + ref_base
        previous_position = position + 1
    
    mdz += str(end_position - previous_position + 1)
    print("Read: " + read, "MD:Z: " + mdz)
