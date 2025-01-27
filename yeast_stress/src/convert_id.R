# Load required libraries
library("biomaRt")
library("dplyr")

# File paths
input_file <- "/home/bit/project/yeast_stress/results/edger/significant_genes.txt"
output_file <- "/home/bit/project/yeast_stress/results/edger/converted_genes.txt"

# Load gene list
gene_list <- readLines(input_file)

# Connect to Ensembl BioMart for Yeast
mart <- useMart("fungi_mart", dataset = "scerevisiae_eg_gene")

# Convert gene symbols to systematic IDs
converted <- getBM(
    attributes = c("external_gene_name", "systematic_name"),
    filters = "external_gene_name",
    values = gene_list,
    mart = mart
)

# Save converted gene list
write.table(
    converted$systematic_name,
    file = output_file,
    row.names = FALSE, col.names = FALSE, quote = FALSE
)

cat("Converted gene list saved to:", output_file, "\n")
