#!/bin/bash

# Script to run Fastp on single-end FASTQ files
# Author: Pranav
# Usage: ./fastp.sh

# Directories
DATA_DIR="/home/bit/project/yeast_stress/data"
OUTPUT_DIR="/home/bit/project/yeast_stress/results/trimmed"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all single-end FASTQ files
for file in ${DATA_DIR}/*.fastq; do
    # Extract the base filename
    base=$(basename "$file" .fastq)

    echo "Running Fastp for sample: $base"

    # Run Fastp for single-end reads
    fastp \
      -i "${file}" \
      -o "${OUTPUT_DIR}/${base}_trimmed.fastq" \
      --html "${OUTPUT_DIR}/${base}_fastp_report.html" \
      --json "${OUTPUT_DIR}/${base}_fastp_report.json" \
      --thread 4

    # Check if Fastp succeeded
    if [ $? -eq 0 ]; then
        echo "Fastp completed successfully for ${base}"
    else
        echo "Error running Fastp for ${base}" >&2
    fi
done

echo "All samples processed. Check the output in ${OUTPUT_DIR}"
