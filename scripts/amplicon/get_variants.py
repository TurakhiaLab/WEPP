#!/bin/usr/python3
# Script to generate VCF files for given FASTQ amplicon reads
import argparse
import errno
import os
import os.path
import sys
import subprocess

def options():
    print("usage: get_variants.py [-h] -f READS [-a ALIGN_WITH] [-r]")
    print()
    print("optional arguments:")
    print("-h, --help            show this help message and exit")
    print("-a, --align_with")
    print("               Specify aligner to use. default = minimap2")
    print("-r             Specify -r flag to collapse duplicate amplicon reads with same starting position.")
    print()
    print(" Required:")
    print("-f READS, --reads READS          Amplicon reads fastq file.")
    print()

# Downloaded reference files
mask_sites = "problematic_sites_sarsCov2.vcf"
covid_ref = "wuhCor1.fa"

# Output files
output_sam = "data/amplicons.sam"
aligned_amplicons = "data/amplicons.fa"
output_vcf = "data/amplicons.vcf"

# Option parsing
parser = argparse.ArgumentParser()
req_grp = parser.add_argument_group(title='Required')
req_grp.add_argument("-f", "--reads", required=True, help="Amplicon reads fastq file.")
parser.add_argument("-a", "--align_with",  default="minimap2", help="Specify aligner to use. default = minimap2")
parser.add_argument("-r", "--remove", action='store_true', help="Specify -r flag to collapse duplicate amplicon reads with same starting position.")
args = parser.parse_args()

amplicon_reads = args.reads
aligner = args.align_with
collapse_duplicates = args.remove

print("FASTQ reads file given: ", amplicon_reads)
print("Using {} aligner".format(aligner))
print("Collapse reads = ", collapse_duplicates)

if (aligner != "minimap2" and aligner != "bwamem2"):
    # Print options function
    print()
    print("ERROR: Invalid aligner given. See options below.")
    options()
    print()
    raise argparse.ArgumentTypeError("Invalid aligner specified, please see options for specifying aligner with -a flag above.")

# Check correct amplicon reads fastq file given
if (os.path.exists(amplicon_reads) == False):
    print()
    print("ERROR: Amplicon read FASTQ file not found.  Make sure the amplicon read file is in current directory.")
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), amplicon_reads)

# Check correct masking sites file downloaded
if (os.path.exists(mask_sites) == False):
    print()
    print("ERROR: 'problematic_sites_sarsCov2.vcf' file not found. Make sure to download this file before running this pipeline.")
    print("This file can be downloaded here:")
    print()
    print("wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf")
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mask_sites)

# Check correct covid reference genome downloaded
if (os.path.exists(covid_ref) == False):
    print()
    print("ERROR: 'wuhCor1.fa' SARS-CoV2 reference file not found. Make sure to download this file before running this pipeline.")
    print("This file can be downloaded and uncompressed here:")
    print()
    print("wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz")
    print("gunzip wuhCor1.fa.gz")
    print()
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), covid_ref)

# Create data dir to dump intermediate files
subprocess.run(["mkdir", "-p", "data"])

# Make sure alignment script has permissions 
subprocess.run(["chmod", "a+x", "align.sh"])

#TODO: Add read filtering steps
#if (collapse_duplicates == True):
    

# Call alignment script: ./aligner <aligner> <amplicon_reads>
run_aligner = [ "./align.sh", aligner, amplicon_reads ]
print(run_aligner)
subprocess.run(run_aligner)
print("Alignment with {} finished and alignment file generated here: {}".format(aligner, output_sam))


# Check correct amplicon reads fastq file given
if (os.path.exists(output_sam) == False):
    print("ERROR: Alignment file, {} not found.".format(output_sam))
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_sam)

# Generate fasta alignment file 
generate_msa = ["python3", "samtomsa.py"]
subprocess.run(generate_msa)

# Check correct aligned amplicons file was generated
if (os.path.exists(aligned_amplicons) == False):
    print("ERROR: MSA file, {} not found.".format(aligned_amplicons))
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), aligned_amplicons)

# Generate variant call file using amplicon multifasta
generate_msa = ["faToVcf", "-maskSites=problematic_sites_sarsCov2.vcf", aligned_amplicons, output_vcf]
subprocess.run(generate_msa)

if (os.path.exists(output_vcf) == False):
    print("ERROR: Variant call file, {} not found.".format(output_vcf))
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_vcf)

print("SUCCESS: Variant call file generated for input amplicon reads: {}".format(output_vcf))

