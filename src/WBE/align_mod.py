import sys

if len(sys.argv) != 2:
    print("Usage: python3 align_mod.py <sam_file>")
    exit()

unknown_count = 0
sam_file = sys.argv[1]

#Read the sam file
with open(sam_file, 'r') as f:
    sam_lines = f.readlines()

#Modify the sam file
with open(sam_file, 'w') as output_file:
    for i in range(len(sam_lines)):
        if i < 3:
            output_file.write(sam_lines[i])
        else:
            lines = sam_lines[i].strip().split('\t')
            #Check if all reads have a name
            read_name = lines[0]
            try:
                flag = int(lines[1])
            except ValueError:
                unknown_count += 1
                read_name = 'Unknown-' + str(unknown_count)
                flag = int(lines[2])
            #Modify the read name
            lines[0] = read_name

            #Join the elements into a tab-separated string
            modified_line = '\t'.join(lines)
            # Write the modified line to the output file
            output_file.write(modified_line + '\n')