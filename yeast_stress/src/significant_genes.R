# Load required libraries
library("dplyr")

# File paths
results_file <- "/home/bit/project/yeast_stress/results/edger/edger_results.csv"
output_file <- "/home/bit/project/yeast_stress/results/edger/significant_genes.txt"

# Load results
res <- read.csv(results_file, row.names = 1)

# Filter significant genes
significant_genes <- res %>%
    filter(PValue < 0.05, abs(logFC) > 1) %>%
    rownames()

# Save gene list (one gene per line)
writeLines(significant_genes, output_file)

cat("Significant gene list saved to:", output_file, "\n")
