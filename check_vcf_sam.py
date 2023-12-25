import re 

def read_fasta(file_path):
    sequence = ""
    with open(file_path, "r") as file:
        # Skip the header line (starts with ">")
        header = file.readline().strip()
        
        # Read the sequence lines and concatenate them
        for line in file:
            sequence += line.strip()

    return header, sequence

def mismatch_cigar_file(file_path):
    with open(file_path, "r") as file:
        # Read each line and strip newline characters
        lines = [line.strip() for line in file]

    return lines

def compare_files(file1_path, file2_path, cigar_list, ref_seq):
    mismatch_count = 0

    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        line2 = next(file2, None)
        for line1_count, line1 in enumerate(file1, start=1):
            parts1 = line1.strip().split('\t')
            parts2 = line2.strip().split('\t')
            split_parts = parts2[0].rsplit('_', 3)
            parts2[0] = split_parts[0]
            start = int(split_parts[2])
            end = int(split_parts[3])
            try:
                if parts1[0] != parts2[0] or parts1[1] != parts2[1]:
                    #print(f"{line1_count} -> File 1: {line1.strip()}")
                    #print(f"File 2: {line2.strip()}")
                    if parts1[0] == parts2[0]:
                        ##print(f"ERROR -> {line1.strip()}")
                        cigar = re.findall(r'(\d+)([A-Za-z])', cigar_list[mismatch_count])
                        ins_offset = 0
                        del_offset = 0
                        idx_offset = 0
                        seq_mismatch = 0
                        
                        #Comparing with original SAM for "M" and "N", skipping insertions and comparing deletions with ref
                        for c in cigar:
                            if c[1] == "M" or c[1] == "N":
                                for i in range(int(c[0])):
                                    if (parts2[1][i + idx_offset + del_offset] != parts1[1][i + idx_offset + ins_offset]):
                                        seq_mismatch += 1
                                idx_offset += int(c[0])
                            elif c[1] == "I":
                                ins_offset += int(c[0]) 
                            elif c[1] == "D":
                                for i in range(int(c[0])):
                                    if (parts2[1][i + idx_offset + del_offset] != ref_seq[start + i + idx_offset + del_offset - 1]):
                                        seq_mismatch += 1
                                del_offset += int(c[0]) 
                        if (seq_mismatch):
                            print(seq_mismatch, cigar_list[mismatch_count])

                        line2 = next(file2)
                    
                    else:
                        print("Read mismatch")
                    
                    mismatch_count += 1
                else:
                    line2 = next(file2)
            except StopIteration:
                mismatch_count += 1
                #print(f"{line1_count} -> File 1: {line1.strip()}")

    return mismatch_count

# Replace 'file1.txt' and 'file2.txt' with your actual file paths
file1_path = 'align_sam_cmp'
file2_path = 'align_vcf_cmp'
fasta_file = "test/NC_045512v2.fa"
cigar_file = "mismatch_cigars"

_, ref_seq = read_fasta(fasta_file)
cigar_list = mismatch_cigar_file(cigar_file)

mismatch_count = compare_files(file1_path, file2_path, cigar_list, ref_seq)
print(f"Total lines with mismatch: {mismatch_count}")