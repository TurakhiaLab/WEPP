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
    print("\tpython3 src/WEPP/qc_preprocess.py --platform i --primers bed_name --in path_to_fastq --out path_to_output --min_len l --threads n --reference ref --prefix prefix")
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
    parser.add_argument('--min_len', type=int, required=True)
    parser.add_argument('--threads', type=int, required=False, default=16)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--prefix', required=True)
    return parser.parse_args()

def resolve_path(path, launch_dir):
    return path if path.startswith('/') or path.startswith('~') else os.path.join(launch_dir, path)

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
        return fq1, fq2
    
    else:
        if len(fastq_files) != 1:
            raise ValueError(f"Expected exactly 1 FASTQ file for single-end data, found {len(fastq_files)}")
        fq = fastq_files[0]
        fname = os.path.basename(fq)
        if re.search(r"_R[12]", fname):
            raise ValueError("Single-end FASTQ should not be named with _R1/_R2")
        return fq, None

def reference_alignment(output_dir, r1, r2, platform, reference, threads, prefix):
    ref_base = os.path.splitext(reference)[0]
    ref_mmi = f"{ref_base}.mmi"
    ref_mmi_path = os.path.join(output_dir, os.path.basename(ref_mmi))
    out_sam = os.path.join(output_dir, prefix + "_aligned_reads.sam")

    print(f"Building Index file in {output_dir} ...")
    subprocess.run(["minimap2", "-d", ref_mmi_path, reference], check=True)

    cmd = ["minimap2", "-a", "--sam-hit-only", "--MD", "-2"]
    if platform == "Illumina":
        cmd += ["-x", "sr"]
        if r2:
            cmd += [ref_mmi_path, r1, r2]
        else:
            cmd += [ref_mmi_path, r1]
    elif platform == "ONT":
        cmd += ["-x", "map-ont", ref_mmi_path, r1]
    cmd += ["-t", str(threads), "-o", out_sam]

    subprocess.run(cmd, check=True)
    return out_sam

def trimming(output_dir, sam_file, primer_bed, min_len, threads, prefix):
    # Sort the SAM file to get sorted.bam
    sorted_bam = os.path.join(output_dir, prefix + "_sorted.bam")
    subprocess.run(["samtools", "sort", sam_file, "-o", sorted_bam, "-@", str(threads)], check=True)

    # Trim with ivar and output directly to trimmed.bam
    trimmed_path = os.path.join(output_dir, prefix + "_trimmed")
    trim_cmd = [
        "ivar", "trim", "-e", "-b", primer_bed, "-p", trimmed_path, "-i", sorted_bam, "-q", "1", "-m", str(min_len), "-x", "3"
    ]
    subprocess.run(trim_cmd, check=True)

    # Sort the trimmed.bam to get resorted.bam
    resorted_bam = os.path.join(output_dir, prefix + "_resorted.bam")
    subprocess.run(["samtools", "sort", trimmed_path + ".bam", "-o", resorted_bam, "-@", str(threads)], check=True)

    # Clean up intermediate files
    os.remove(sam_file)
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

    primer_bed_file = resolve_path(args.primers, "./primers")
    reference_file = resolve_path(args.reference, args.input_dir)
    r1, r2 = find_fastq_files(args.input_dir, is_paired)

    sam_file = reference_alignment(args.output_dir, r1, r2, platform, reference_file, args.threads, args.prefix)
    trimming(args.output_dir, sam_file, primer_bed_file, args.min_len, args.threads, args.prefix)
    print(f"QC Completed")


if __name__ == "__main__":
    main()
