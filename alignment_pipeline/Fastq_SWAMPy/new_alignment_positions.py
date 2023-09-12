# Open the SAM file
with open("intermediate_files/alignment.sam", "r") as samfile:
    # Open the output file
    with open("intermediate_files/alignment_modified.sam", "w") as outfile:
        # Process each line in the SAM file
        for line in samfile:
            # Skip header lines
            if line.startswith('@'):
                outfile.write(line)
                continue

            # Split the line into fields
            fields = line.rstrip().split('\t')

            # Calculate the end position
            # Note: This is a simplification and may not be accurate if the alignment includes insertions, deletions, or soft clipping
            start = int(fields[3])  # Start position is the 4th field in SAM format
            seq = fields[9]  # Sequence is the 10th field
            end = start + len(seq) - 1  # End position = start position + length of the sequence - 1

            # If fields[0] contains '_amplicon_', then remove the '_amplicon_' and everything after it
            # if '_amplicon_' in fields[0]:
            #     fields[0] = fields[0].split('_amplicon_')[0]

            if start == 0: # Skip unmapped reads
                continue

            # Append positions to read name (the 1st field)
            fields[0] += "_READ_" + str(start) + "_" + str(end)

            # Write the modified line to the output file
            outfile.write('\t'.join(fields) + '\n')
