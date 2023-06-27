# Add two required arguments for the file name of the reads_aligned.fasta file and differences.txt file
# The arguments should be required

import sys

# Check the correct number of arguments have been provided
if len(sys.argv) != 3:
    print("Usage: python calculate_alignment_rate.py reads_aligned.fasta differences.txt")
    sys.exit(1)

aligned_reads = sys.argv[1]
differences = sys.argv[2]

file = open(aligned_reads, "r")
output_file = open(differences,"w")
output_file.write("Start_Position_Before End_Position_Before Start_Position_After End_Position_After max(abs(Start_position_before - Start_position_after),abs(end_position_before - end_position_after))\n")
lines = file.readlines()
file.close()

lines = lines[0::2]
lines = [line.strip() for line in lines]
same_alignments = 0
for line in lines:
    if line.split('_')[-1] == line.split('_')[-3] and line.split('_')[-2] == line.split('_')[-4]:
        same_alignments += 1

    else:
        output_file.write(line.split('_')[-4] + ' ' + line.split('_')[-3]  + ' ' + line.split('_')[-2]  + ' ' + line.split('_')[-1] + ' ' + str(max(abs(int(line.split('_')[-4]) - int(line.split('_')[-2])),abs(int(line.split('_')[-3]) - int(line.split('_')[-1])))) + '\n')

total_number_of_alignments = len(lines)
alignment_rate = same_alignments / total_number_of_alignments

print("Alignment rate: " + str(alignment_rate))
file.close()
output_file.close()