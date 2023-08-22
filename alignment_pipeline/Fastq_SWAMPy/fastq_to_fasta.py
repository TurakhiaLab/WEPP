import sys

def filter_lines(input_file, output_file):
    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()

            f2 = open("invalid_lines.txt", "w")
            
            extracted_lines = []
            for i in range(0, len(lines), 4):  # iterate over lines in steps of 4
                extracted_lines.extend(lines[i:i+2])  # add the first 2 lines of each set of 4 to the result

            # for lines starting with @, remove the @ and replace with >
            extracted_lines = [line.replace('@', '>') if line.startswith('@') else line for line in extracted_lines]

            for i in range(1, len(extracted_lines), 2):  # iterate over lines in steps of 2)
                # check if the line contains anything except A, C, G, T
                if any([char not in ['A', 'C', 'G', 'T'] for char in extracted_lines[i].strip()]):
                    f2.write(extracted_lines[i-1])
                    f2.write(extracted_lines[i])

            f2.close()
            
            with open(output_file, 'w') as outfile:
                outfile.writelines(extracted_lines)

        print(f"Successfully filtered lines from '{input_file}' and saved to '{output_file}'.")
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python fastq_to_fasta.py <input_file> <output_file>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    filter_lines(input_file_path, output_file_path)
