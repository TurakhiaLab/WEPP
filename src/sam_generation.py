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

def read_csv_abundance_file(file):
    hap_abun = {}
    hap_lineage = {}
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            hap, lineage, abun = row[0], row[1], row[2] 
            hap_abun[hap] = float(abun)
            hap_lineage[hap] = lineage
    return hap_abun, hap_lineage

def read_csv_uncertainty_file(file):
    hap_nodes = defaultdict(set)
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            hap, *nodes = row
            hap_nodes[hap].update(nodes)
    return hap_nodes

def read_tsv_file(file):
    # Increase CSV field size limit
    max_int = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_int)
            break
        except OverflowError:
            max_int = int(max_int / 10)

    hap_sams = defaultdict(set)
    with open(file, newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for idx, row in enumerate(reader):
            value, *keys = row
            hap_sams[value] = keys
    return hap_sams

def write_sam_files(input_sam_file):
    write_file_sam_reads = os.path.join(results_directory, f"{file_prefix}_haplotype_reads.sam")
    write_file_bam_reads = os.path.join(results_directory, f"{file_prefix}_haplotype_reads.bam")
    write_file_sam_haps = os.path.join(results_directory, f"{file_prefix}_haplotypes.sam")
    write_file_bam_haps = os.path.join(results_directory, f"{file_prefix}_haplotypes.bam")
    
    with open(write_file_sam_reads, mode='w', newline='') as w_file_reads:
        add_RG = True
        ref_name = ""
        with open(input_sam_file, 'r') as r_file_reads:
            for line in r_file_reads:
                line = line.strip()
                if line.startswith('@'):
                    w_file_reads.write(line+"\n")
                    # Find reference name
                    if line.startswith('@SQ'):
                        fields = line.split()
                        for field in fields:
                            if field.startswith('SN:'):
                                ref_name = field[3:]
                else:
                    if add_RG:
                        add_RG = False
                        # Write unaccounted groups of mutations
                        for mut, idx in mutations.items():
                            w_file_reads.write(f"@CO\tUM:unaccounted{idx}\tUS:{mut}\n")
                        # Write read groups of haplotypes
                        for hap, idx in haplotypes.items():
                            w_file_reads.write(f"@RG\tID:group{idx}\tDS:Node:{hap}")
                            # Write unaccounted groups of haplotypes
                            if len(hap_muts[hap]):
                                w_file_reads.write(f"\tUM:Z")
                                for m_idx, mut in enumerate(hap_muts[hap]):
                                    if mut in mutations:
                                        if m_idx == 0:
                                            w_file_reads.write(f":unaccounted{mutations[mut]}")
                                        else:
                                            w_file_reads.write(f",unaccounted{mutations[mut]}")
                            w_file_reads.write(f"\tHS:Z:{hap_abun[hap]}\tHL:Z:{hap_lineage[hap]}")
                            if (len(hap_nodes[hap])):
                                w_file_reads.write(f"\tUH:Z:{','.join(hap_nodes[hap])}")
                            w_file_reads.write("\n")
                    else:
                        tokens = line.split()
                        # Write haplotypes mapping to this read
                        w_file_reads.write(line+"\tRG:Z")
                        for idx, hap in enumerate(read_haps[tokens[0]]):
                            if idx == 0:
                                w_file_reads.write(f":group{haplotypes[hap]}")
                            else:
                                w_file_reads.write(f",group{haplotypes[hap]}")
                        w_file_reads.write(f"\tEP:i:{len(read_haps[tokens[0]])}")

                        # Write unaccounted mutations present in this read
                        if (len(read_muts[tokens[0]])):
                            w_file_reads.write("\tUM:Z")
                            for idx, mut in enumerate(read_muts[tokens[0]]):
                                if idx == 0:
                                    w_file_reads.write(f":unaccounted{mutations[mut]}")
                                else:
                                    w_file_reads.write(f",unaccounted{mutations[mut]}")

                        w_file_reads.write("\n")

        # Convert sam to bam
        command_1 = f"samtools view -Sb {write_file_sam_reads} | samtools sort -o {write_file_bam_reads}"
        command_2 = f"samtools index {write_file_bam_reads}"
        commands = [command_1, command_2]
        # Execute sam to bam conversion
        for command in commands:
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                print(f"{command} failed")
                print("Error message:", result.stderr.decode())
                sys.exit(1)      

    with open(write_file_sam_haps, mode='w', newline='') as w_file_haps:
        first_hap, first_row = next(iter(hap_sams.items()))
        w_file_haps.write(f"@SQ\tSN:{ref_name}\tLN:{len(first_row[-3])}\n@CO\tHP_SEQ\n")
        # Write unaccounted groups of mutations
        for mut, idx in mutations.items():
            w_file_haps.write(f"@CO\tUM:unaccounted{idx}\tUS:{mut}\n")
        # Write read groups of haplotypes
        for hap, idx in haplotypes.items():
            w_file_haps.write(f"@RG\tID:group{idx}\tDS:Node:{hap}")
            # Write unaccounted groups of haplotypes
            if len(hap_muts[hap]):
                w_file_haps.write(f"\tUM:Z")
                for m_idx, mut in enumerate(hap_muts[hap]):
                    if mut in mutations:
                        if m_idx == 0:
                            w_file_haps.write(f":unaccounted{mutations[mut]}")
                        else:
                            w_file_haps.write(f",unaccounted{mutations[mut]}")
            w_file_haps.write(f"\tHS:Z:{hap_abun[hap]}\tHL:Z:{hap_lineage[hap]}")
            if (len(hap_nodes[hap])):
                w_file_haps.write(f"\tUH:Z:{','.join(hap_nodes[hap])}")
            w_file_haps.write("\n")
        
        #Write haplotypes now
        for hap, idx in haplotypes.items():
            w_file_haps.write(hap + "\t" + "\t".join(hap_sams[hap]) + str(idx))
            # Write unaccounted groups of haplotypes
            if len(hap_muts[hap]):
                w_file_haps.write(f"\tUM:Z")
                for m_idx, mut in enumerate(hap_muts[hap]):
                    if mut in mutations:
                        if m_idx == 0:
                            w_file_haps.write(f":unaccounted{mutations[mut]}")
                        else:
                            w_file_haps.write(f",unaccounted{mutations[mut]}")
            w_file_haps.write("\n")

        # Convert sam to bam
        command_1 = f"samtools view -Sb {write_file_sam_haps} | samtools sort -o {write_file_bam_haps}"
        command_2 = f"samtools index {write_file_bam_haps}"
        commands = [command_1, command_2]
        # Execute sam to bam conversion
        for command in commands:
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                print(f"{command} failed")
                print("Error message:", result.stderr.decode())
                sys.exit(1)      

#Check arguments
if len(sys.argv) != 4:
    print("USAGE: python src/WBE/sam_generation.py <results_directory> <intermediate_directory> <file_prefix>")
    sys.exit(1)

# Start time
start_time = time.time()
results_directory = sys.argv[1]
intermediate_directory = sys.argv[2]
file_prefix = sys.argv[3]

# Reading haplotype_reads File
read_haps, haplotypes = read_csv_file(results_directory + "/" + file_prefix + "_haplotype_reads.csv")

# Reading mutations_reads File
read_muts, mutations = read_csv_file(results_directory + "/" + file_prefix + "_mutation_reads.csv")

# Reading mutations_haplotypes File
hap_muts, _ = read_csv_file(results_directory + "/" + file_prefix + "_mutation_haplotypes.csv")

# Reading haplotype_abundance File
hap_abun, hap_lineage = read_csv_abundance_file(results_directory + "/" + file_prefix + "_haplotype_abundance.csv")

# Reading haplotype_uncertainty File
hap_nodes = read_csv_uncertainty_file(results_directory + "/" + file_prefix + "_haplotype_uncertainty.csv")

# Reading haplotypes File
hap_sams = read_tsv_file(results_directory + "/" + file_prefix + "_haplotypes.tsv")

# Writing File
write_sam_files(intermediate_directory + "/" + file_prefix + "_alignment.sam")

print(f"\nElapsed time: {time.time() - start_time} seconds")