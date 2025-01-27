#!/bin/bash

# Script to align trimmed single-end reads using HISAT2
# Usage: ./hisat2.sh

# Directories
DATA_DIR="/home/bit/project/yeast_stress/results/trimmed"
OUTPUT_DIR="/home/bit/project/yeast_stress/results/alignment"
GENOME_INDEX="/home/bit/project/yeast_stress/data/genome/sacCer3"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Align each trimmed FASTQ file
for file in ${DATA_DIR}/*_trimmed.fastq; do
    # Extract the base name
    base=$(basename "$file" _trimmed.fastq)

    echo "Aligning sample: $base"

    # Run HISAT2
    hisat2 -x "$GENOME_INDEX" \
           -U "$file" \
           -S "${OUTPUT_DIR}/${base}_aligned.sam" \
           --summary-file "${OUTPUT_DIR}/${base}_alignment_summary.txt"

    # Check success
    if [ $? -eq 0 ]; then
        echo "Alignment completed for ${base}"
    else
        echo "Error during alignment for ${base}" >&2
    fi
done

echo "All alignments completed. SAM files are saved in ${OUTPUT_DIR}"
