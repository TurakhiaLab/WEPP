def correct_fastq(file_path, output_path):
    with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
        lines = infile.readlines()
        
        i = 0
        while i < len(lines):
            if lines[i].startswith('@'):
                # Write header line
                outfile.write(lines[i])
                
                # Process sequence line
                sequence = ''
                i += 1
                while i < len(lines) and not lines[i].startswith('+'):
                    sequence += lines[i].strip()
                    i += 1
                outfile.write(sequence + '\n')
                
                # Write plus line
                if i < len(lines):
                    outfile.write(lines[i])
                
                # Process quality score line
                quality = ''
                i += 1
                while i < len(lines) and not lines[i].startswith('@'):
                    quality += lines[i].strip()
                    i += 1
                outfile.write(quality + '\n')
            else:
                i += 1

# Usage
input_file = 'point_loma/PL_2023_03_01_reads.fastq_orig'  # Replace with the path to your input FASTQ file
output_file = 'point_loma/PL_2023_03_01_reads.fastq'  # Replace with the desired path for the corrected output FASTQ file

correct_fastq(input_file, output_file)
