import pysam
import os
import sys
import subprocess
from collections import defaultdict, Counter
from scipy.stats import entropy

def parse_sam_and_compute_af(sam_file, output_file):
    # Dictionary to store counts of alleles at different reference positions
    allele_dict = defaultdict(Counter)

    # Open the SAM file
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam.fetch():
            if read.is_unmapped:  # Skip unmapped reads
                continue

            # Get the query sequence and the CIGAR string
            query_seq = read.query_sequence
            cigar = read.cigartuples
            ref_start = read.reference_start  # 0-based reference start

            # Process the CIGAR string to map query positions to reference positions
            ref_pos = ref_start  # Start at the 0-based reference position
            query_pos = 0  # Query sequence is 0-based

            for op, length in cigar:
                if op == 0:  # Match or mismatch (M)
                    for i in range(length):
                        allele_dict[ref_pos + 1][query_seq[query_pos]] += 1  # Convert to 1-based
                        ref_pos += 1
                        query_pos += 1
                elif op == 1:  # Insertion (I)
                    query_pos += length  # Skip these positions in the query
                elif op == 2:  # Deletion (D)
                    for i in range(length):
                        allele_dict[ref_pos + 1]["-"] += 1  # Store deletion as "-" at the reference position
                        ref_pos += 1  # Skip these positions in the reference
                elif op == 4:  # Soft clipping (S)
                    query_pos += length  # Skip soft-clipped bases
                elif op == 5:  # Hard clipping (H)
                    continue  # Hard clipping does not consume query sequence
                else:
                    # Other CIGAR operations (e.g., padding) can be ignored here
                    continue

    # Compute allele frequencies and depth
    af_dict = {}
    depth_dict = {}

    for ref_pos, alleles in allele_dict.items():
        total_counts = sum(alleles.values())
        af_dict[ref_pos] = {allele: count / total_counts for allele, count in alleles.items()}
        depth_dict[ref_pos] = total_counts 

    # Write allele frequencies and depth to a file
    with open(output_file, "w") as f:
        f.write("Position\tAllele\tFrequency\tDepth\n")
        for ref_pos, afs in sorted(af_dict.items()):
            depth = depth_dict[ref_pos]
            for allele, freq in afs.items():
                f.write(f"{ref_pos}\t{allele}\t{freq:.10f}\t{depth}\n")

    print(f"'{output_file}' generated.")

def compute_kl_divergence(ref_file, subsampled_file):
    # Read allele frequencies from files
    def read_allele_frequencies(file):
        af_dict = {}
        depth_dict = {} 
        with open(file, "r") as f:
            next(f)  # Skip header
            for line in f:
                pos, allele, freq, depth = line.strip().split()
                pos = int(pos)
                freq = float(freq)
                dep = int(depth)
                if pos not in af_dict:
                    af_dict[pos] = {}
                af_dict[pos][allele] = freq
                depth_dict[pos] = dep
        return af_dict, depth_dict

    af1, dep1 = read_allele_frequencies(ref_file)
    af2, dep2 = read_allele_frequencies(subsampled_file)

    total_depth = sum(dep1.values())

    # Aggregate allele frequencies across all positions
    p_aggregated = []
    q_aggregated = []

    for pos, alleles in af1.items():  # Iterate through positions and alleles in af1
        for allele in alleles.keys():  # Iterate through all alleles in af1
            freq1 = af1[pos][allele]  # Frequency in ref_file
            freq2 = af2.get(pos, {}).get(allele, 1e-10)  # Frequency in subsampled_file with pseudocounts

            # Weight the frequencies by the depth at this position
            depth = dep1.get(pos, 0)
            #p_aggregated.append(freq1 * depth / total_depth)
            #q_aggregated.append(freq2 * depth / total_depth)
            p_aggregated.append(freq1)
            q_aggregated.append(freq2)

    # Compute the overall KL divergence
    total_weighted_kl = entropy(p_aggregated, q_aggregated)

    return total_weighted_kl

def check_and_process_files(folder, prefix):
    ref_sam_file = os.path.join(folder, f"{prefix}_alignment.sam_orig")
    ref_af_file = os.path.join(folder, "ref_allele_frequencies.txt")
    sampled_sam_file = os.path.join(folder, f"{prefix}_alignment.sam")
    sampled_af_file = os.path.join(folder, "sampled_allele_frequencies.txt")

    # Check if the reference output file exists
    if not os.path.exists(ref_af_file):
        parse_sam_and_compute_af(ref_sam_file, ref_af_file)
    else:
        print(f"'{ref_af_file}' already exists. Skipping processing.")

    # Loop to generate subsampled files and compute KL divergence
    min_kl_divergence = float("inf")
    best_sampled_file = None
    
    kl_divergence = compute_kl_divergence(ref_af_file, sampled_af_file)
    for i in range(50):
        # Run the subsampling script
        command = [
            "python",
            "src/WBE/select_subset_reads.py",
            ref_sam_file,
            sampled_sam_file,
        ]
        subprocess.run(command, check=True)
        
        # Compute allele frequencies for the newly generated subsampled file
        parse_sam_and_compute_af(sampled_sam_file, sampled_af_file)

        # Compute KL divergence
        kl_divergence = compute_kl_divergence(ref_af_file, sampled_af_file)
        print(f"KL divergence for iteration {i + 1}: {kl_divergence}")

        # Update the best KL divergence and corresponding file
        if kl_divergence < min_kl_divergence:
            min_kl_divergence = kl_divergence
            best_sampled_file = f"{sampled_sam_file}_best"
            # Save the best subsampled SAM file
            os.rename(sampled_sam_file, best_sampled_file)

    print(f"Best KL divergence: {min_kl_divergence}")
    print(f"Best subsampled SAM file saved as: {best_sampled_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <folder> <prefix>")
        sys.exit(1)

    folder = sys.argv[1]
    prefix = sys.argv[2]

    check_and_process_files(folder, prefix)
