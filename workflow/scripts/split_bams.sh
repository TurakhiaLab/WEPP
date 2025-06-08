#!/bin/bash

# Smart parallel script to split BAM file by read groups
# Usage: ./smart_split_bam_by_readgroups.sh input.bam [output_directory] [max_parallel_jobs] [--force]

# Check if input file is provided
if [ $# -lt 1 ] || [ $# -gt 4 ]; then
    echo "Usage: $0 input.bam [output_directory] [max_parallel_jobs] [--force]"
    echo "This script splits a BAM file by read groups into separate files using parallel processing"
    echo "  max_parallel_jobs: default is number of CPU cores"
    echo "  --force: overwrite existing files"
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_DIR="${2:-.}"  # Use current directory if not specified
MAX_JOBS="${3:-$(nproc)}"  # Use number of CPU cores if not specified
FORCE_OVERWRITE=false

# Check for --force flag in any position
for arg in "$@"; do
    if [ "$arg" = "--force" ]; then
        FORCE_OVERWRITE=true
        echo "Force overwrite mode enabled"
        break
    fi
done

# If third argument is --force, adjust MAX_JOBS
if [ "$3" = "--force" ]; then
    MAX_JOBS="$(nproc)"
fi

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
echo "Force overwrite: $FORCE_OVERWRITE"
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

# Check existing files
echo "Checking for existing output files..."
EXISTING_FILES=()
MISSING_FILES=()
INCOMPLETE_FILES=()

for GROUP_ID in $READ_GROUPS; do
    OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    INDEX_FILE="${OUTPUT_FILE}.bai"
    
    if [ -f "$OUTPUT_FILE" ]; then
        # Check if file is complete (not empty and has valid BAM header)
        if [ -s "$OUTPUT_FILE" ] && samtools view -H "$OUTPUT_FILE" >/dev/null 2>&1; then
            # Check if index exists
            if [ -f "$INDEX_FILE" ]; then
                EXISTING_FILES+=("$GROUP_ID")
                echo "  ✓ $GROUP_ID: BAM and index exist"
            else
                INCOMPLETE_FILES+=("$GROUP_ID")
                echo "  ⚠ $GROUP_ID: BAM exists but index missing"
            fi
        else
            INCOMPLETE_FILES+=("$GROUP_ID")
            echo "  ⚠ $GROUP_ID: BAM exists but appears corrupt/incomplete"
        fi
    else
        MISSING_FILES+=("$GROUP_ID")
        echo "  ✗ $GROUP_ID: Missing"
    fi
done

echo ""
echo "File status summary:"
echo "  Complete files: ${#EXISTING_FILES[@]}"
echo "  Incomplete/corrupt files: ${#INCOMPLETE_FILES[@]}"
echo "  Missing files: ${#MISSING_FILES[@]}"
echo ""

# Determine what needs to be processed
GROUPS_TO_PROCESS=()

if [ "$FORCE_OVERWRITE" = true ]; then
    echo "Force mode: Processing all read groups"
    GROUPS_TO_PROCESS=($READ_GROUPS)
else
    # Only process missing and incomplete files
    GROUPS_TO_PROCESS=(${MISSING_FILES[@]} ${INCOMPLETE_FILES[@]})
    
    if [ ${#EXISTING_FILES[@]} -gt 0 ]; then
        echo "Skipping ${#EXISTING_FILES[@]} existing complete files:"
        printf "  %s\n" "${EXISTING_FILES[@]}"
        echo ""
    fi
    
    if [ ${#GROUPS_TO_PROCESS[@]} -eq 0 ]; then
        echo "All files already exist and are complete!"
        echo "Use --force to regenerate all files"
        exit 0
    fi
fi

echo "Processing ${#GROUPS_TO_PROCESS[@]} read groups..."

# Interactive confirmation for large number of existing files
if [ ${#EXISTING_FILES[@]} -gt 10 ] && [ "$FORCE_OVERWRITE" = false ]; then
    echo ""
    echo "Found ${#EXISTING_FILES[@]} existing files. Options:"
    echo "1. Process only missing/incomplete files (default)"
    echo "2. Process all files (overwrite existing)"
    echo "3. Exit"
    read -p "Choose option [1]: " choice
    
    case $choice in
        2)
            echo "Processing all files..."
            GROUPS_TO_PROCESS=($READ_GROUPS)
            FORCE_OVERWRITE=true
            ;;
        3)
            echo "Exiting..."
            exit 0
            ;;
        *)
            echo "Processing only missing/incomplete files..."
            ;;
    esac
fi

# Extract non-@RG header lines (these will be common to all output files)
echo "Extracting common header information..."
COMMON_HEADER=$(samtools view -H "$INPUT_BAM" | grep -v "^@RG")

# Save common header to temporary file for reuse
COMMON_HEADER_FILE=$(mktemp)
echo "$COMMON_HEADER" > "$COMMON_HEADER_FILE"

# Counter for progress
TOTAL_GROUPS=${#GROUPS_TO_PROCESS[@]}
echo "Total read groups to process: $TOTAL_GROUPS"

if [ $TOTAL_GROUPS -eq 0 ]; then
    echo "Nothing to process!"
    rm "$COMMON_HEADER_FILE"
    exit 0
fi

echo "Starting parallel processing..."
echo ""

# Function to process a single read group
process_read_group() {
    local GROUP_ID="$1"
    local INPUT_BAM="$2"
    local OUTPUT_DIR="$3"
    local BASE_NAME="$4"
    local COMMON_HEADER_FILE="$5"
    local FORCE_OVERWRITE="$6"
    
    local OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
    local INDEX_FILE="${OUTPUT_FILE}.bai"
    local PID=$$
    
    echo "[${PID}] Processing read group: $GROUP_ID"
    echo "[${PID}]   Output file: $OUTPUT_FILE"
    
    # Final check before processing (in case of race conditions)
    if [ "$FORCE_OVERWRITE" = false ] && [ -f "$OUTPUT_FILE" ] && [ -f "$INDEX_FILE" ]; then
        if samtools view -H "$OUTPUT_FILE" >/dev/null 2>&1; then
            local READ_COUNT=$(samtools view -c "$OUTPUT_FILE" 2>/dev/null || echo "0")
            echo "[${PID}]   File already exists and is valid ($READ_COUNT reads), skipping"
            echo "$GROUP_ID:$READ_COUNT:$(stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo "0"):1:SKIPPED" >> "${OUTPUT_DIR}/.processing_results"
            return 0
        fi
    fi
    
    # Remove existing files if we're overwriting
    if [ "$FORCE_OVERWRITE" = true ]; then
        rm -f "$OUTPUT_FILE" "$INDEX_FILE"
    fi
    
    # Extract the specific @RG line for this group
    local SPECIFIC_RG_HEADER=$(samtools view -H "$INPUT_BAM" | grep "^@RG" | grep "ID:$GROUP_ID")
    
    if [ -z "$SPECIFIC_RG_HEADER" ]; then
        echo "[${PID}]   Error: No @RG header found for group $GROUP_ID"
        echo "$GROUP_ID:ERROR:No RG header" >> "${OUTPUT_DIR}/.processing_results"
        return 1
    fi
    
    # Create temporary header file with common headers + specific RG header
    local TEMP_HEADER=$(mktemp)
    {
        cat "$COMMON_HEADER_FILE"
        echo "$SPECIFIC_RG_HEADER"
    } > "$TEMP_HEADER"
    
    # Filter reads by read group and combine with custom header
    echo "[${PID}]   Creating BAM with group-specific header..."
    
    # Use temporary output file to ensure atomic writes
    local TEMP_OUTPUT="${OUTPUT_FILE}.tmp.${PID}"
    
    {
        # First output the custom header
        cat "$TEMP_HEADER"
        # Then output only the reads for this specific read group
        samtools view "$INPUT_BAM" | grep -F "$GROUP_ID"
    } | samtools view -b > "$TEMP_OUTPUT"
    
    # Clean up temporary header file
    rm "$TEMP_HEADER"
    
    # Check if temporary file was created successfully
    if [ -f "$TEMP_OUTPUT" ] && [ -s "$TEMP_OUTPUT" ]; then
        # Verify the BAM file is valid
        if samtools view -H "$TEMP_OUTPUT" >/dev/null 2>&1; then
            # Move temporary file to final location
            mv "$TEMP_OUTPUT" "$OUTPUT_FILE"
            
            # Get read count for verification
            local READ_COUNT=$(samtools view -c "$OUTPUT_FILE")
            echo "[${PID}]   Created successfully with $READ_COUNT reads"
            
            # Verify the header contains only one @RG line
            local RG_COUNT=$(samtools view -H "$OUTPUT_FILE" | grep -c "^@RG")
            echo "[${PID}]   Header contains $RG_COUNT @RG line(s) (should be 1)"
            
            # Index the output BAM file
            echo "[${PID}]   Indexing $OUTPUT_FILE..."
            if samtools index "$OUTPUT_FILE"; then
                echo "[${PID}]   Indexed successfully"
                
                # Write success status to temp file for summary
                echo "$GROUP_ID:$READ_COUNT:$(stat -c%s "$OUTPUT_FILE"):$RG_COUNT:SUCCESS" >> "${OUTPUT_DIR}/.processing_results"
            else
                echo "[${PID}]   Warning: Indexing failed"
                echo "$GROUP_ID:$READ_COUNT:$(stat -c%s "$OUTPUT_FILE"):$RG_COUNT:INDEX_FAILED" >> "${OUTPUT_DIR}/.processing_results"
            fi
        else
            echo "[${PID}]   Error: Created BAM file is invalid"
            rm -f "$TEMP_OUTPUT" "$OUTPUT_FILE"
            echo "$GROUP_ID:ERROR:Invalid BAM" >> "${OUTPUT_DIR}/.processing_results"
            return 1
        fi
    else
        echo "[${PID}]   Error: Failed to create $OUTPUT_FILE"
        rm -f "$TEMP_OUTPUT"
        echo "$GROUP_ID:ERROR:Creation failed" >> "${OUTPUT_DIR}/.processing_results"
        return 1
    fi
    
    echo "[${PID}] Completed read group: $GROUP_ID"
}

# Export the function so it can be used by background processes
export -f process_read_group

# Initialize results file
> "${OUTPUT_DIR}/.processing_results"

# Process read groups in parallel, limiting concurrent jobs
echo "Starting parallel processing with job control..."

job_count=0
for GROUP_ID in "${GROUPS_TO_PROCESS[@]}"; do
    # Wait if we've reached the maximum number of jobs
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 0.1
    done
    
    # Start background job
    process_read_group "$GROUP_ID" "$INPUT_BAM" "$OUTPUT_DIR" "$BASE_NAME" "$COMMON_HEADER_FILE" "$FORCE_OVERWRITE" &
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
ls -la "${OUTPUT_DIR}/${BASE_NAME}"-*.bam 2>/dev/null || echo "No BAM files found"

# Generate summary from results
echo ""
echo "Summary:"
if [ -f "${OUTPUT_DIR}/.processing_results" ]; then
    SUCCESS_COUNT=0
    ERROR_COUNT=0
    SKIP_COUNT=0
    
    while IFS=':' read -r group_id read_count file_size rg_count status; do
        case "$status" in
            "SUCCESS")
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                file_size_human=$(numfmt --to=iec-i --suffix=B $file_size 2>/dev/null || echo "${file_size} bytes")
                echo "  ✓ ${BASE_NAME}-${group_id}.bam: $read_count reads, $file_size_human, $rg_count @RG header(s)"
                ;;
            "SKIPPED")
                SKIP_COUNT=$((SKIP_COUNT + 1))
                file_size_human=$(numfmt --to=iec-i --suffix=B $file_size 2>/dev/null || echo "${file_size} bytes")
                echo "  → ${BASE_NAME}-${group_id}.bam: $read_count reads, $file_size_human (skipped - already exists)"
                ;;
            "INDEX_FAILED")
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                file_size_human=$(numfmt --to=iec-i --suffix=B $file_size 2>/dev/null || echo "${file_size} bytes")
                echo "  ⚠ ${BASE_NAME}-${group_id}.bam: $read_count reads, $file_size_human (index failed)"
                ;;
            *)
                ERROR_COUNT=$((ERROR_COUNT + 1))
                echo "  ✗ ${BASE_NAME}-${group_id}.bam: FAILED ($read_count)"
                ;;
        esac
    done < "${OUTPUT_DIR}/.processing_results"
    
    echo ""
    echo "Final summary:"
    echo "  Successfully processed: $SUCCESS_COUNT"
    echo "  Skipped (already existed): $SKIP_COUNT"
    echo "  Failed: $ERROR_COUNT"
    echo "  Total files in output directory: $(ls -1 "${OUTPUT_DIR}/${BASE_NAME}"-*.bam 2>/dev/null | wc -l)"
    
    # Clean up results file
    rm "${OUTPUT_DIR}/.processing_results"
else
    echo "No results file found - using fallback summary"
    # Fallback to original method if results file doesn't exist
    for GROUP_ID in $READ_GROUPS; do
        OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}-${GROUP_ID}.bam"
        if [ -f "$OUTPUT_FILE" ]; then
            READ_COUNT=$(samtools view -c "$OUTPUT_FILE" 2>/dev/null || echo "0")
            FILE_SIZE=$(ls -lh "$OUTPUT_FILE" 2>/dev/null | awk '{print $5}' || echo "unknown")
            RG_COUNT=$(samtools view -H "$OUTPUT_FILE" 2>/dev/null | grep -c "^@RG" || echo "0")
            echo "  $OUTPUT_FILE: $READ_COUNT reads, $FILE_SIZE, $RG_COUNT @RG header(s)"
        fi
    done
fi

echo ""
echo "Header verification:"
echo "Each output BAM should contain:"
echo "  - All common headers (@HD, @SQ, @PG, etc.)"
echo "  - Only one @RG line (specific to that read group)"