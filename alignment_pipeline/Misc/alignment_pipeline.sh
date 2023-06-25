# Extract Reads and Generate a FASTA file
python3 vcf_to_fasta.py my_vcf_reads.vcf NC_045512v2.fa reads.fasta

# Align each read to the reference genome
bowtie2-build NC_045512v2.fa ref
bowtie2 -x ref -f -U reads.fasta -S align.sam --local

# Get the alignment positions
awk '{print $1"\t"$4"\t"($4+length($10))}' align.sam > align_positions.txt

# Convert SAM to BAM
samtools view -b align.sam > align.bam

# Sort BAM file
samtools sort align.bam -o align_sorted.bam

# Index the sorted BAM file
samtools index align_sorted.bam

python3 -c "
import pandas as pd
from Bio import SeqIO

# load the alignment positions into a dictionary
df = pd.read_csv('align_positions.txt', sep='\t', header=None, names=['id', 'start', 'end'])
align_positions = df.set_index('id').T.to_dict('list')

# load the sequences and update the identifiers
reads = SeqIO.parse('reads.fasta', 'fasta')
updated_reads = []
for read in reads:
    read_id = read.id.split('_READ_')[0]
    if read_id in align_positions:
        read.id += f'_READ_{align_positions[read_id][0]}_{align_positions[read_id][1]}'
    updated_reads.append(read)

# write the updated sequences to a new file
SeqIO.write(updated_reads, 'updated_reads.fasta', 'fasta')
"

# Align each read to the reference genome
# Build Bowtie2 index
bowtie2-build NC_045512v2.fa reference

# Perform the alignment
bowtie2 -x reference -f -U reads.fasta -S align.sam

# Convert SAM to BAM
samtools view -b align.sam > align.bam

# Sort BAM file
samtools sort align.bam -o align_sorted.bam

# Index the sorted BAM file
samtools index align_sorted.bam

# Generate a VCF file
bcftools mpileup -Ou -f NC_045512v2.fa align_sorted.bam | bcftools call -vmO z -o new.vcf.gz

# Index the VCF file
tabix new.vcf.gz

# Generate the final VCF file
python3 -c "
import pandas as pd
import pysam

# Load the alignment positions
positions = pd.read_csv('align_positions.txt', sep='\t', header=None, names=['id', 'start', 'end'])

# Load the original VCF file
orig_vcf = pd.read_csv('my_vcf_reads.vcf', sep='\t', comment='#', header=None)

# Load the new VCF file
new_vcf = pysam.VariantFile('new.vcf.gz')

# Initialize a dictionary to store the new sample names
samples_dict = dict(zip(orig_vcf.iloc[:, 9:].columns, positions['id']))

# Create a new DataFrame with the same structure as the original VCF file
new_vcf_df = orig_vcf.copy()

# Update the sample names in the new DataFrame
new_vcf_df.rename(columns=samples_dict, inplace=True)

# Iterate over the new VCF records
for rec in new_vcf:
    # Find the corresponding row in the new DataFrame
    row_idx = new_vcf_df[new_vcf_df[1] == rec.pos].index[0]

    # Update the mutation information in the DataFrame
    new_vcf_df.at[row_idx, 3] = rec.ref
    new_vcf_df.at[row_idx, 4] = ','.join(rec.alts)

    # Update the sample information in the DataFrame
    for sample in rec.samples:
        if sample in new_vcf_df.columns:
            new_vcf_df.at[row_idx, sample] = rec.samples[sample]['GT'][0]

# Write the new DataFrame to a VCF file
new_vcf_df.to_csv('updated.vcf', sep='\t', index=False, header=False)
"
