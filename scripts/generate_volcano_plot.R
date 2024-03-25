# if (!requireNamespace("ggplot2", quietly = TRUE)) {
#     install.packages("ggplot2")
# }

library(ggplot2)

# Load the DESeq2 results
results_path <- 'outputs/deseq2/deseq2_results.csv'
res <- read.csv(results_path)

# Calculate the -log10 of the p-value
res$negLogPvalue <- -log10(res$pvalue)

# Define significance thresholds
log2FC_threshold <- 1
pvalue_threshold <- 0.05

# Create a basic volcano plot
ggplot(res, aes(x=log2FoldChange, y=negLogPvalue)) +
  geom_point(aes(color=(abs(log2FoldChange) >= log2FC_threshold & pvalue < pvalue_threshold)), alpha=0.5) +
  scale_color_manual(values=c("FALSE"="black", "TRUE"="red")) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 p-value") +
  theme_minimal() +
  theme(legend.position="none") +
  geom_vline(xintercept=c(-log2FC_threshold, log2FC_threshold), linetype="dashed", color="blue") +
  geom_hline(yintercept=-log10(pvalue_threshold), linetype="dashed", color="blue")

# Save the plot
ggsave('outputs/deseq2/volcano_plot.png', width=10, height=8)
