import argparse
from Bio import SeqIO
import pysam
import bisect
import time

# Function to read a FASTA file and store sequences in a dictionary
def read_fasta(file):
    record = next(SeqIO.parse(file, "fasta"))
    sequence_id = record.id
    sequence = str(record.seq)
    return sequence

# Function to read a SAM file and store alignments in a list
def read_sam(file):
    alignments = []
    with pysam.AlignmentFile(file, "r") as samfile:
        for read in samfile.fetch():
            # Filtering out unaligned reads
            if not read.is_unmapped:
                alignments.append(read)
    return alignments

# Function to find the highest block start point less than or equal to the given val
def index_update(blocks, val):
    start_points = [block[0] for block in blocks]
    pos = bisect.bisect_right(start_points, val) - 1
    if pos >= 0:
        return val + blocks[pos][1]
    else:
        return val

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Takes MSA Reference FASTA, Aligned Read SAM file, and output SAM file")
    parser.add_argument("fasta_file", help="Path to the input FASTA file")
    parser.add_argument("sam_file", help="Path to the input SAM file")
    parser.add_argument("output_file", help="Path to the output file")

    # Parse arguments
    args = parser.parse_args()

    # Start timer
    start_time = time.time()

    # Read FASTA file
    ref_sequence = read_fasta(args.fasta_file)

    # Get non-gapped blocks from ref_sequence
    non_gap_blocks = []
    gap_count = 0
    in_block = False
    for i, char in enumerate(ref_sequence):
        if char != '-':
            if not in_block:
                non_gap_blocks.append((i + 1 - gap_count, gap_count))
                in_block = True
        else:
            gap_count += 1
            in_block = False

    # Read SAM file
    aligned_reads = read_sam(args.sam_file)

    # Write SAM file
    with open(args.output_file, 'w') as outfile:
        outfile.write("Read_Name\tReference_Start\tRead_Sequence\tRead_Qualities\n")
        for alignment in aligned_reads:
            new_index = index_update(non_gap_blocks, alignment.reference_start + 1)
            read_quality = ''.join(chr(q + 33) for q in alignment.query_qualities)
            outfile.write(
                    f"{alignment.query_name}\t"
                    f"{new_index}\t"
                    f"{alignment.query_sequence}\t"
                    f"{read_quality}\n"
            )

    print(f"{args.output_file} written in {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()
