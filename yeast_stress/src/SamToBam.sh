#!/bin/bash

# Script to convert SAM to BAM and sort the files using SAMtools
# Usage: ./SamToBam.sh

# Directories
ALIGN_DIR="/home/bit/project/yeast_stress/results/alignment"
OUTPUT_DIR="/home/bit/project/yeast_stress/results/alignment"

# Loop through SAM files and convert to sorted BAM
for file in ${ALIGN_DIR}/*_aligned.sam; do
    base=$(basename "$file" _aligned.sam)

    echo "Converting $base to BAM and sorting..."

    # Convert to BAM and sort
    samtools view -bS "$file" | samtools sort -o "${OUTPUT_DIR}/${base}_sorted.bam"

    # Index the BAM file
    samtools index "${OUTPUT_DIR}/${base}_sorted.bam"

    echo "Completed: ${base}_sorted.bam"
done

echo "All SAM files converted to sorted BAM."
