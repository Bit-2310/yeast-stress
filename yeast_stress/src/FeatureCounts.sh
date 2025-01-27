#!/bin/bash

# Script to run FeatureCounts on aligned BAM files
# Usage: ./FeatureCounts.sh

# Directories
ALIGN_DIR="/home/bit/project/yeast_stress/results/alignment"
OUTPUT_DIR="/home/bit/project/yeast_stress/results/counts"
GTF_FILE="/home/bit/project/yeast_stress/data/Saccharomyces_cerevisiae.R64-1-1.59.gtf"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FeatureCounts
featureCounts -T 4 -a "$GTF_FILE" -o "${OUTPUT_DIR}/gene_counts.txt" \
    ${ALIGN_DIR}/*_sorted.bam

echo "FeatureCounts completed. Gene counts saved to ${OUTPUT_DIR}/gene_counts.txt"
