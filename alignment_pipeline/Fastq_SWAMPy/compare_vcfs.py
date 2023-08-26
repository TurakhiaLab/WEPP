def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Remove lines starting with a hashtag
        lines = [line.strip() for line in lines if not line.startswith('#')]

        # Group every 2 lines together
        lines = [lines[i:i+2] for i in range(0, len(lines), 2)]

        # Order the groups by the first line in alphabetical order
        lines = sorted(lines, key=lambda x: x[0])

        # Flatten the list of lines to make a list of strings
        lines = [line for group in lines for line in group]

        return lines

def print_differences(file1, file2):
    lines1 = read_file(file1)
    lines2 = read_file(file2)
    count = 0

    if lines1 == lines2:
        print("The files are identical.")
    else:
        print("The files are not identical.")
        print("Differences:")

        for line1, line2 in zip(lines1, lines2):
            if line1 != line2:
                # print(f"Line {count}")
                count += 1
                # line1_list = line1.split('\t')[9:]
                # line2_list = line2.split('\t')[9:]

                # # Find which elements are different
                # for i, (elem1, elem2) in enumerate(zip(line1_list, line2_list)):
                #     if elem1 != elem2:
                #         print(f"Element {i}: {elem1} != {elem2}")
            
                print(f"Local Alignment: {line1}")
                print(f"Global Alignment: {line2}")
                print()
                print(count)

        # Check if there are additional lines in either file
        print("Additional Lines")
        if len(lines1) > len(lines2):
            for line in lines1[len(lines2):]:
                print(f"File 1: {line}")

        if len(lines2) > len(lines1):
            for line in lines2[len(lines1):]:
                print(f"File 2: {line}")

    print(f"Total differences: {count}")

# Usage
file1 = "./reads_aligned.fasta"
file2 = "./reads_aligned_normal.fasta"
print_differences(file1, file2)
