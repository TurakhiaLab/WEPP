#!/bin/bash

# Script to split BAM file by read groups
# Usage: ./split_bam_by_readgroups.sh input.bam

# Check if input file is provided
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 input.bam [output_directory]"
    echo "This script splits a BAM file by read groups into separate files"
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_DIR="${2:-.}"  # Use current directory if not specified

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if input file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input file '$INPUT_BAM' not found"
    exit 1
fi

# Get base filename without extension
BASE_NAME=$(basename "$INPUT_BAM" .bam)

echo "Processing BAM file: $INPUT_BAM"
echo "Base filename: $BASE_NAME"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Extract read group IDs from BAM header
echo "Extracting read group IDs from BAM header..."
READ_GROUPS=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | sed 's/.*ID:\([^[:space:]]*\).*/\1/')

# Check if any read groups were found
if [ -z "$READ_GROUPS" ]; then
    echo "No read groups found in BAM file header"
    exit 1
fi

echo "Found read groups:"
echo "$READ_GROUPS"
echo ""

# Counter for progress
TOTAL_GROUPS=$(echo "$READ_GROUPS" | wc -l)
CURRENT=0

# Process each read group
for GROUP_ID in $READ_GROUPS; do
    CURRENT=$((CURRENT + 1))
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    
    echo "[$CURRENT/$TOTAL_GROUPS] Processing read group: $GROUP_ID"
    echo "  Output file: $OUTPUT_FILE"
    
    # Filter BAM by read group using grep method
    samtools view -h "$INPUT_BAM" | grep -E "(^@|$GROUP_ID)" | samtools view -b > "$OUTPUT_FILE"
    
    # Check if output file was created successfully
    if [ -f "$OUTPUT_FILE" ]; then
        # Get read count for verification
        READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
        echo "  Created successfully with $READ_COUNT reads"
        
        # Index the output BAM file
        echo "  Indexing $OUTPUT_FILE..."
        samtools index "$OUTPUT_FILE"
        echo "  Indexed successfully"
    else
        echo "  Error: Failed to create $OUTPUT_FILE"
    fi
    echo ""
done

echo "Processing complete!"
echo "Created BAM files in $OUTPUT_DIR:"
ls -la "${OUTPUT_DIR}/${BASE_NAME}"-*.bam

# Summary
echo ""
echo "Summary:"
for GROUP_ID in $READ_GROUPS; do
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    if [ -f "$OUTPUT_FILE" ]; then
        READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
        FILE_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
        echo "  $OUTPUT_FILE: $READ_COUNT reads, $FILE_SIZE"
    fi
done