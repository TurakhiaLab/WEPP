#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 8 ]; then
    echo "Usage: treePrune.sh -i <path to input .mat file> -o <path to output .mat file> -s <selected sample> -d <distance>"
    exit 1
fi

while getopts ":i:o:s:d:" opt; do
  case $opt in
    i) input_file="$OPTARG"
    ;;
    o) output_file="$OPTARG"
    ;;
    s) selected_sample="$OPTARG"
    ;;
    d) distance="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Step 1: Extract the samples that are within 'd' distance of the selected sample
matUtils extract -i $input_file -dist $selected_sample:$distance -u samples_to_prune.txt

# Step 2: Prune the identified samples from the original tree
matUtils extract -i $input_file --prune samples_to_prune.txt -o $output_file

# Cleanup the temporary file
# rm samples_to_prune.txt

echo "Pruning complete. The pruned tree can be found at: $output_file"
