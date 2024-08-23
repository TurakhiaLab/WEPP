from Bio import SeqIO

# Input files
fasta_file = "/data/pgangwar/PANMAT_data/sequences.fa"  # Your FASTA file
sequence_names_file = "junk1"  # File containing sequence names

# Output file
output_fasta_file = "output_sequences.fasta"

with open(sequence_names_file, 'r') as names_file:
    sequence_names = [line.strip() for line in names_file]

# Create a hash table (dictionary) to store the number of occurrences of each sequence name
sequence_name_counts = {}

# Populate the hash table with counts of each sequence name
for name in sequence_names:
    if name in sequence_name_counts:
        sequence_name_counts[name] += 1
    else:
        sequence_name_counts[name] = 1

print("Read Sequence Names")

# List to store the line numbers where sequences start
sequence_start_positions = {}

# Read the FASTA file line by line to find the sequence headers and track line numbers
with open(fasta_file, 'r') as fasta:
    current_line_number = 0  # Initialize line number
    for line in fasta:
        current_line_number += 1  # Increment line number for each line read
        
        if line.startswith(">"):  # Check if it's a header line
            # Extract the sequence name from the header line
            sequence_name = line[1:].split()[0]  # Assuming the name is the first part after '>'
            
            # Check if the sequence name is in the sequence_name_counts dictionary
            if sequence_name in sequence_name_counts:
                # Store the current line number and how many times to print the sequence
                sequence_start_positions[current_line_number] = sequence_name_counts[sequence_name]

print("Created Fasta List")

# Now iterate through the FASTA file again to print sequences based on sequence_start_positions
with open(fasta_file, 'r') as fasta, open(output_fasta_file, 'w') as output_file:
    current_line_number = 0
    sequence_buffer = []  # Buffer to store the current sequence lines
    write_flag = False
    print_count = 0
    
    for line in fasta:
        current_line_number += 1
        
        # Check if the current line number is in sequence_start_positions
        if current_line_number in sequence_start_positions:
            # If we were already writing a sequence, flush the buffer to the output file
            if sequence_buffer:
                for _ in range(print_count):
                    output_file.write("".join(sequence_buffer))
                sequence_buffer = []  # Clear the buffer

            write_flag = True
            print_count = sequence_start_positions[current_line_number]  # Set the number of times to print the sequence

        # If write_flag is True, add the line to the buffer
        if write_flag:
            sequence_buffer.append(line)
        
        # Check if the current line is a header (">"), indicating the start of a new sequence
        if line.startswith(">") and len(sequence_buffer) > 1:
            for _ in range(print_count):
                output_file.write("".join(sequence_buffer[:-1]))
            sequence_buffer = [] 
            write_flag = False  # Turn off writing until next valid sequence start

    # After finishing the loop, flush any remaining sequences in the buffer
    if sequence_buffer:
        for _ in range(print_count):
            output_file.write("".join(sequence_buffer))

print("Completed Fasta file Dump")