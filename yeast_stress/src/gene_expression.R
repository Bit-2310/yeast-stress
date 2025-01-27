# Load required libraries
library("edgeR")
library("ggplot2")
library("dplyr")
library("ggrepel")

# Define file paths
counts_file <- "/home/bit/project/yeast_stress/results/counts/gene_counts.txt"
output_dir <- "/home/bit/project/yeast_stress/results/edger"
dir.create(output_dir, showWarnings = FALSE)

# Load gene counts
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, -(1:5)]  # Remove metadata columns if present

# Check data structure
if (ncol(counts) < 2) stop("Error: Counts file has insufficient columns.")

# Define sample information
group <- factor(c("control", "control", "treated", "treated", "treated", "treated"))

# Create DGEList object
dge <- DGEList(counts = counts, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize data using TMM (Trimmed Mean of M-values)
dge <- calcNormFactors(dge)

# Estimate dispersion
dge <- estimateDisp(dge)

# Fit the model and perform exact test
et <- exactTest(dge, pair = c("control", "treated"))

# Get results table
res <- topTags(et, n = Inf)$table
res <- res %>% arrange(PValue)

# Add significance categories
res$significant <- "Not Significant"
res$significant[res$PValue < 0.05 & res$logFC > 1] <- "Upregulated"
res$significant[res$PValue < 0.05 & res$logFC < -1] <- "Downregulated"

# Save results to CSV
write.csv(res, file = file.path(output_dir, "edger_results.csv"))

# Volcano Plot with Labels for Top 10 DEGs
top_genes <- res %>%
    filter(PValue < 0.05) %>%
    head(10)

volcano <- ggplot(res, aes(x = logFC, y = -log10(PValue), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                    size = 4, max.overlaps = 10) +
    scale_color_manual(values = c("Not Significant" = "gray50",
                                  "Upregulated" = "#E41A1C",
                                  "Downregulated" = "#377EB8")) +
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dotted") +
    labs(title = "Volcano Plot of Differential Expression (edgeR)",
         x = "Log2 Fold Change",
         y = "-Log10 P-Value",
         color = "Significance") +
    theme_minimal(base_size = 14)

# Save Volcano Plot
ggsave(file.path(output_dir, "volcano_plot.png"), plot = volcano, width = 8, height = 6)

# DEG Summary
deg_summary <- res %>%
    group_by(significant) %>%
    summarise(count = n())

# Save DEG summary
write.csv(deg_summary, file = file.path(output_dir, "deg_summary.csv"))

# Save session info
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

# Completion message
cat("edgeR analysis completed successfully!\n")
cat("Results saved to:", file.path(output_dir, "edger_results.csv"), "\n")
cat("Volcano plot saved to:", file.path(output_dir, "volcano_plot.png"), "\n")
cat("DEG summary saved to:", file.path(output_dir, "deg_summary.csv"), "\n")
cat("Session info saved to:", file.path(output_dir, "session_info.txt"), "\n")
