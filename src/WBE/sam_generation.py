import sys
import csv
import os
import time
import subprocess
from collections import defaultdict


def read_csv_file(file):
    read_haps = defaultdict(set)
    haplotypes = {}
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for idx, row in enumerate(reader):
            value, *keys = row
            haplotypes[value] = idx+1
            for key in keys:
                read_haps[key].add(value)   
    return read_haps, haplotypes

def write_sam_files(input_sam_file):
    write_file_sam = os.path.join(directory, f"{file_prefix}_haplotype_reads.sam")
    write_file_bam = os.path.join(directory, f"{file_prefix}_haplotype_reads.bam")
    
    with open(write_file_sam, mode='w', newline='') as w_file:
        add_RG = True
        with open(input_sam_file, 'r') as r_file:
            for line in r_file:
                line = line.strip()
                if line.startswith('@'):
                    w_file.write(line+"\n")
                else:
                    if add_RG:
                        add_RG = False
                        for hap, idx in haplotypes.items():
                            w_file.write(f"@RG\tID:group{idx}\tDS:Node:{hap}\n")
                    else:
                        tokens = line.split()
                        w_file.write(line+"\tXG:Z")
                        for idx, hap in enumerate(read_hap_map[tokens[0]]):
                            if idx == 0:
                                w_file.write(f":group{haplotypes[hap]}")
                            else:
                                w_file.write(f",group{haplotypes[hap]}")
                        w_file.write("\n")

        # Convert sam to bam
        command_1 = f"samtools view -Sb {write_file_sam} | samtools sort -o {write_file_bam}"
        command_2 = f"samtools index {write_file_bam}"
        commands = [command_1, command_2]
        # Execute sam to bam conversion
        for command in commands:
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                print(f"{command} failed")
                print("Error message:", result.stderr.decode())
                sys.exit(1)      

#Check arguments
if len(sys.argv) != 3:
    print("USAGE: python src/WBE/sam_generation.py <directory> <file_prefix>")
    sys.exit(1)

# Start time
start_time = time.time()
file_prefix = sys.argv[2]
directory = sys.argv[1]

# Reading File
read_hap_map, haplotypes = read_csv_file(directory + "/" + file_prefix + "_haplotype_reads.csv")

# Writing File
write_sam_files(directory + "/" + file_prefix + "_alignment.sam")

print(f"\nElapsed time: {time.time() - start_time} seconds")