#!/usr/bin/env python3

import os
import re
import sys
import glob
import subprocess
import argparse

all_sample_names = []

def print_usage():
    print("Example usage: ")
    print("\tpython3 qc_preprocess.py --platform i --primers path/to/bed --in path/to/fastq/ --out path/to/output --threads n --reference ref --prefix pre")
    print("\nPlatform options:")
    print("\t-h\t\t: Show this help message and exit.")
    print("\t-i\t\t: Illumina platform, paired-end mode")
    print("\t-s\t\t: Illumina platform, single read mode")
    print("\t-n\t\t: ONT platform")
    sys.exit()

def parse_args():
    parser = argparse.ArgumentParser(description='QC workflow')
    parser.add_argument('--platform', required=True)
    parser.add_argument('--primers', required=True)
    parser.add_argument('--in', dest='input_dir', required=True)
    parser.add_argument('--out', dest='output_dir', required=True)
    parser.add_argument('--threads', type=int, required=False, default=8)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--prefix', required=True)
    return parser.parse_args()

def resolve_path(path, launch_dir):
    return path if path.startswith('/') or path.startswith('~') else os.path.join(launch_dir, path)

def get_sample_name(filename):
    name = os.path.basename(filename).split('.')[0].split('_')[0]
    if name in all_sample_names:
        raise ValueError(f"Duplicate sample name detected: {name}")
    all_sample_names.append(name)
    return name

def find_fastq_files(input_dir, paired):
    fastq_files = sorted(glob.glob(os.path.join(input_dir, "*.fastq*")))
    if paired:
        if len(fastq_files) != 2:
            raise ValueError(f"Expected exactly 2 FASTQ files for paired-end data, found {len(fastq_files)}")
        fq1 = fq2 = None
        for fq in fastq_files:
            fname = os.path.basename(fq)
            if re.search(r"_R1", fname):
                fq1 = fq
            elif re.search(r"_R2", fname):
                fq2 = fq
        if not fq1 or not fq2:
            raise ValueError("Both R1 and R2 files must be present and correctly named (e.g., *_R1.fastq, *_R2.fastq)")
        sample = get_sample_name(fq1)
        return {sample: (fq1, fq2)}
    
    else:
        if len(fastq_files) != 1:
            raise ValueError(f"Expected exactly 1 FASTQ file for single-end data, found {len(fastq_files)}")
        fq = fastq_files[0]
        fname = os.path.basename(fq)
        if re.search(r"_R[12]", fname):
            raise ValueError("Single-end FASTQ should not be named with _R1/_R2")
        sample = get_sample_name(fq)
        return {sample: (fq, None)}

def reference_alignment(output_dir, r1, r2, platform, reference, threads, prefix):
    ref_base = reference.replace(".fa", "")
    out_sam = os.path.join(output_dir, prefix + "_aligned_reads.sam")

    cmd = ["minimap2", "-a", "--sam-hit-only", "--MD", "-2"]
    if platform == "Illumina":
        cmd += ["-x", "sr"]
        if r2:
            cmd += [f"{ref_base}.mmi", r1, r2]
        else:
            cmd += [f"{ref_base}.mmi", r1]
    elif platform == "ONT":
        cmd += ["-x", "map-ont", f"{ref_base}.mmi", r1]
    cmd += ["-t", str(threads), "-o", out_sam]

    subprocess.run(cmd, check=True)
    return out_sam

def trimming(output_dir, sam_file, primer_bed, threads, prefix):
    # Sort the SAM file to get sorted.bam
    sorted_bam = os.path.join(output_dir, prefix + "_sorted.bam")
    subprocess.run(["samtools", "sort", sam_file, "-o", sorted_bam, "-@", str(threads)], check=True)

    # Trim with ivar and output directly to trimmed.bam
    trimmed_path = os.path.join(output_dir, prefix + "_trimmed")
    trim_cmd = [
        "ivar", "trim", "-e", "-b", primer_bed, "-p", trimmed_path, "-i", sorted_bam, "-q", "1", "-m", "80", "-x", "3"
    ]
    subprocess.run(trim_cmd, check=True)

    # Sort the trimmed.bam to get resorted.bam
    resorted_bam = os.path.join(output_dir, prefix + "_resorted.bam")
    subprocess.run(["samtools", "sort", trimmed_path + ".bam", "-o", resorted_bam, "-@", str(threads)], check=True)

    # Clean up intermediate files
    os.remove(sorted_bam)
    os.remove(trimmed_path + ".bam")

def main():
    args = parse_args()

    platform_map = {
        "d": ("Illumina", True),
        "s": ("Illumina", False),
        "n": ("ONT", False),
    }

    if args.platform == "h":
        print_usage()

    platform_info = platform_map.get(args.platform)
    if not platform_info:
        print("ERROR: unknown sequencing modality")
        print_usage()
        return 1

    platform, is_paired = platform_info
    print(f"{platform} mode")

    primer_bed_file = resolve_path(args.primers, "./database")
    reference_file = resolve_path(args.reference, "./database")
    fastq_files = find_fastq_files(args.input_dir, is_paired)

    sample_name, (r1, r2) = next(iter(fastq_files.items()))
    sam_file = reference_alignment(args.output_dir, r1, r2, platform, reference_file, args.threads, args.prefix)
    trimming(args.output_dir, sam_file, primer_bed_file, args.threads, args.prefix)
    print(f"QC done on {sample_name}")


if __name__ == "__main__":
    main()
