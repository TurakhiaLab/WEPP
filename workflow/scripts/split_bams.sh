# #!/bin/bash

# # Script to split BAM file by read groups
# # Usage: ./split_bam_by_readgroups.sh input.bam

# # Check if input file is provided
# if [ $# -lt 1 ] || [ $# -gt 2 ]; then
#     echo "Usage: $0 input.bam [output_directory]"
#     echo "This script splits a BAM file by read groups into separate files"
#     exit 1
# fi

# INPUT_BAM="$1"
# OUTPUT_DIR="${2:-.}"  # Use current directory if not specified

# # Create output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# # Check if input file exists
# if [ ! -f "$INPUT_BAM" ]; then
#     echo "Error: Input file '$INPUT_BAM' not found"
#     exit 1
# fi

# # Get base filename without extension
# BASE_NAME=$(basename "$INPUT_BAM" .bam)

# echo "Processing BAM file: $INPUT_BAM"
# echo "Base filename: $BASE_NAME"
# echo "Output directory: $OUTPUT_DIR"
# echo ""

# # Extract read group IDs from BAM header
# echo "Extracting read group IDs from BAM header..."
# READ_GROUPS=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | sed 's/.*ID:\([^[:space:]]*\).*/\1/')

# # Check if any read groups were found
# if [ -z "$READ_GROUPS" ]; then
#     echo "No read groups found in BAM file header"
#     exit 1
# fi

# echo "Found read groups:"
# echo "$READ_GROUPS"
# echo ""

# # Counter for progress
# TOTAL_GROUPS=$(echo "$READ_GROUPS" | wc -l)
# CURRENT=0

# # Process each read group
# for GROUP_ID in $READ_GROUPS; do
#     CURRENT=$((CURRENT + 1))
#     OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    
#     echo "[$CURRENT/$TOTAL_GROUPS] Processing read group: $GROUP_ID"
#     echo "  Output file: $OUTPUT_FILE"
    
#     # Filter BAM by read group using grep method
#     samtools view -h "$INPUT_BAM" | grep -E "(^@|$GROUP_ID)" | samtools view -b > "$OUTPUT_FILE"
    
#     # Check if output file was created successfully
#     if [ -f "$OUTPUT_FILE" ]; then
#         # Get read count for verification
#         READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
#         echo "  Created successfully with $READ_COUNT reads"
        
#         # Index the output BAM file
#         echo "  Indexing $OUTPUT_FILE..."
#         samtools index "$OUTPUT_FILE"
#         echo "  Indexed successfully"
#     else
#         echo "  Error: Failed to create $OUTPUT_FILE"
#     fi
#     echo ""
# done

# echo "Processing complete!"
# echo "Created BAM files in $OUTPUT_DIR:"
# ls -la "${OUTPUT_DIR}/${BASE_NAME}"-*.bam

# # Summary
# echo ""
# echo "Summary:"
# for GROUP_ID in $READ_GROUPS; do
#     OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
#     if [ -f "$OUTPUT_FILE" ]; then
#         READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
#         FILE_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
#         echo "  $OUTPUT_FILE: $READ_COUNT reads, $FILE_SIZE"
#     fi
# done


# #!/bin/bash

# # Script to split BAM file by read groups
# # Usage: ./split_bam_by_readgroups.sh input.bam

# # Check if input file is provided
# if [ $# -lt 1 ] || [ $# -gt 2 ]; then
#     echo "Usage: $0 input.bam [output_directory]"
#     echo "This script splits a BAM file by read groups into separate files"
#     exit 1
# fi

# INPUT_BAM="$1"
# OUTPUT_DIR="${2:-.}"  # Use current directory if not specified

# # Create output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# # Check if input file exists
# if [ ! -f "$INPUT_BAM" ]; then
#     echo "Error: Input file '$INPUT_BAM' not found"
#     exit 1
# fi

# # Get base filename without extension
# BASE_NAME=$(basename "$INPUT_BAM" .bam)

# echo "Processing BAM file: $INPUT_BAM"
# echo "Base filename: $BASE_NAME"
# echo "Output directory: $OUTPUT_DIR"
# echo ""

# # Extract read group IDs from BAM header
# echo "Extracting read group IDs from BAM header..."
# READ_GROUPS=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | sed 's/.*ID:\([^[:space:]]*\).*/\1/')

# # Check if any read groups were found
# if [ -z "$READ_GROUPS" ]; then
#     echo "No read groups found in BAM file header"
#     exit 1
# fi

# echo "Found read groups:"
# echo "$READ_GROUPS"
# echo ""

# # Extract non-@RG header lines (these will be common to all output files)
# echo "Extracting common header information..."
# COMMON_HEADER=$(samtools view -H "$INPUT_BAM" | grep -v "^@RG")

# # Counter for progress
# TOTAL_GROUPS=$(echo "$READ_GROUPS" | wc -l)
# CURRENT=0

# # Process each read group
# for GROUP_ID in $READ_GROUPS; do
#     CURRENT=$((CURRENT + 1))
#     OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    
#     echo "[$CURRENT/$TOTAL_GROUPS] Processing read group: $GROUP_ID"
#     echo "  Output file: $OUTPUT_FILE"
    
#     # Extract the specific @RG line for this group
#     SPECIFIC_RG_HEADER=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | grep "ID:$GROUP_ID")
    
#     # Create temporary header file with common headers + specific RG header
#     TEMP_HEADER=$(mktemp)
#     {
#         echo "$COMMON_HEADER"
#         echo "$SPECIFIC_RG_HEADER"
#     } > "$TEMP_HEADER"
    
#     # Filter reads by read group and combine with custom header
#     echo "  Creating BAM with group-specific header..."
#     {
#         # First output the custom header
#         cat "$TEMP_HEADER"
#         # Then output only the reads for this specific read group
#         samtools view "$INPUT_BAM" | grep -F "$GROUP_ID"
#     } | samtools view -b > "$OUTPUT_FILE"
    
#     # Clean up temporary file
#     rm "$TEMP_HEADER"
    
#     # Check if output file was created successfully
#     if [ -f "$OUTPUT_FILE" ]; then
#         # Get read count for verification
#         READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
#         echo "  Created successfully with $READ_COUNT reads"
        
#         # Verify the header contains only one @RG line
#         RG_COUNT=$(samtools view -H "$OUTPUT_FILE" | grep -c "^@RG")
#         echo "  Header contains $RG_COUNT @RG line(s) (should be 1)"
        
#         # Index the output BAM file
#         echo "  Indexing $OUTPUT_FILE..."
#         samtools index "$OUTPUT_FILE"
#         echo "  Indexed successfully"
#     else
#         echo "  Error: Failed to create $OUTPUT_FILE"
#     fi
#     echo ""
# done

# echo "Processing complete!"
# echo "Created BAM files in $OUTPUT_DIR:"
# ls -la "${OUTPUT_DIR}/${BASE_NAME}"-*.bam

# # Summary
# echo ""
# echo "Summary:"
# for GROUP_ID in $READ_GROUPS; do
#     OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
#     if [ -f "$OUTPUT_FILE" ]; then
#         READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
#         FILE_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
#         RG_COUNT=$(samtools view -H "$OUTPUT_FILE" | grep -c "^@RG")
#         echo "  $OUTPUT_FILE: $READ_COUNT reads, $FILE_SIZE, $RG_COUNT @RG header(s)"
#     fi
# done

# echo ""
# echo "Header verification:"
# echo "Each output BAM should contain:"
# echo "  - All common headers (@HD, @SQ, @PG, etc.)"
# echo "  - Only one @RG line (specific to that read group)"


#!/bin/bash

# Parallel script to split BAM file by read groups
# Usage: ./parallel_split_bam_by_readgroups.sh input.bam [output_directory] [max_parallel_jobs]

# Check if input file is provided
if [ $# -lt 1 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 input.bam [output_directory] [max_parallel_jobs]"
    echo "This script splits a BAM file by read groups into separate files using parallel processing"
    echo "  max_parallel_jobs: default is number of CPU cores"
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_DIR="${2:-.}"  # Use current directory if not specified
MAX_JOBS="${3:-$(nproc)}"  # Use number of CPU cores if not specified

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
echo "Max parallel jobs: $MAX_JOBS"
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

# Extract non-@RG header lines (these will be common to all output files)
echo "Extracting common header information..."
COMMON_HEADER=$(samtools view -H "$INPUT_BAM" | grep -v "^@RG")

# Save common header to temporary file for reuse
COMMON_HEADER_FILE=$(mktemp)
echo "$COMMON_HEADER" > "$COMMON_HEADER_FILE"

# Counter for progress
TOTAL_GROUPS=$(echo "$READ_GROUPS" | wc -l)
echo "Total read groups to process: $TOTAL_GROUPS"
echo "Starting parallel processing..."
echo ""

# Function to process a single read group
process_read_group() {
    local GROUP_ID="$1"
    local INPUT_BAM="$2"
    local OUTPUT_DIR="$3"
    local BASE_NAME="$4"
    local COMMON_HEADER_FILE="$5"
    
    local OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    
    echo "[$$] Processing read group: $GROUP_ID"
    echo "[$$]   Output file: $OUTPUT_FILE"
    
    # Extract the specific @RG line for this group
    local SPECIFIC_RG_HEADER=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | grep "ID:$GROUP_ID")
    
    # Create temporary header file with common headers + specific RG header
    local TEMP_HEADER=$(mktemp)
    {
        cat "$COMMON_HEADER_FILE"
        echo "$SPECIFIC_RG_HEADER"
    } > "$TEMP_HEADER"
    
    # Filter reads by read group and combine with custom header
    echo "[$$]   Creating BAM with group-specific header..."
    {
        # First output the custom header
        cat "$TEMP_HEADER"
        # Then output only the reads for this specific read group
        samtools view "$INPUT_BAM" | grep -F "$GROUP_ID"
    } | samtools view -b > "$OUTPUT_FILE"
    
    # Clean up temporary file
    rm "$TEMP_HEADER"
    
    # Check if output file was created successfully
    if [ -f "$OUTPUT_FILE" ]; then
        # Get read count for verification
        local READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
        echo "[$$]   Created successfully with $READ_COUNT reads"
        
        # Verify the header contains only one @RG line
        local RG_COUNT=$(samtools view -H "$OUTPUT_FILE" | grep -c "^@RG")
        echo "[$$]   Header contains $RG_COUNT @RG line(s) (should be 1)"
        
        # Index the output BAM file
        echo "[$$]   Indexing $OUTPUT_FILE..."
        samtools index "$OUTPUT_FILE"
        echo "[$$]   Indexed successfully"
        
        # Write success status to temp file for summary
        echo "$GROUP_ID:$READ_COUNT:$(stat -c%s "$OUTPUT_FILE"):$RG_COUNT" >> "${OUTPUT_DIR}/.processing_results"
    else
        echo "[$$]   Error: Failed to create $OUTPUT_FILE"
        echo "$GROUP_ID:ERROR" >> "${OUTPUT_DIR}/.processing_results"
    fi
    
    echo "[$$] Completed read group: $GROUP_ID"
}

# Export the function so it can be used by background processes
export -f process_read_group

# Initialize results file
> "${OUTPUT_DIR}/.processing_results"

# Method 1: Using background processes with job control
echo "Starting parallel processing with job control..."

# Process read groups in parallel, limiting concurrent jobs
job_count=0
for GROUP_ID in $READ_GROUPS; do
    # Wait if we've reached the maximum number of jobs
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 0.1
    done
    
    # Start background job
    process_read_group "$GROUP_ID" "$INPUT_BAM" "$OUTPUT_DIR" "$BASE_NAME" "$COMMON_HEADER_FILE" &
    job_count=$((job_count + 1))
    
    echo "Started job $job_count/$TOTAL_GROUPS for read group: $GROUP_ID (PID: $!)"
done

# Wait for all background jobs to complete
echo ""
echo "Waiting for all jobs to complete..."
wait

# Clean up common header file
rm "$COMMON_HEADER_FILE"

echo ""
echo "Processing complete!"
echo "Created BAM files in $OUTPUT_DIR:"
ls -la "${OUTPUT_DIR}/${BASE_NAME}"-*.bam

# Generate summary from results
echo ""
echo "Summary:"
if [ -f "${OUTPUT_DIR}/.processing_results" ]; then
    while IFS=':' read -r group_id read_count file_size rg_count; do
        if [ "$read_count" = "ERROR" ]; then
            echo "  ${OUTPUT_DIR}/${BASE_NAME}-${group_id}.bam: FAILED"
        else
            file_size_human=$(numfmt --to=iec-i --suffix=B $file_size 2>/dev/null || echo "${file_size} bytes")
            echo "  ${OUTPUT_DIR}/${BASE_NAME}-${group_id}.bam: $read_count reads, $file_size_human, $rg_count @RG header(s)"
        fi
    done < "${OUTPUT_DIR}/.processing_results"
    
    # Clean up results file
    rm "${OUTPUT_DIR}/.processing_results"
else
    # Fallback to original method if results file doesn't exist
    for GROUP_ID in $READ_GROUPS; do
        OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
        if [ -f "$OUTPUT_FILE" ]; then
            READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
            FILE_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
            RG_COUNT=$(samtools view -H "$OUTPUT_FILE" | grep -c "^@RG")
            echo "  $OUTPUT_FILE: $READ_COUNT reads, $FILE_SIZE, $RG_COUNT @RG header(s)"
        fi
    done
fi

echo ""
echo "Header verification:"
echo "Each output BAM should contain:"
echo "  - All common headers (@HD, @SQ, @PG, etc.)"
echo "  - Only one @RG line (specific to that read group)"