# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DESeq2")
# if (!requireNamespace("ggplot2", quietly = TRUE)) {
#     install.packages("ggplot2")
# }

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  all install.packages("BiocManager")
#
# BiocManager::install("apeglm")

library(DESeq2)

count_data <- read.csv('outputs/deseq2/count_matrix.csv', row.names=1)
col_data <- read.csv('outputs/deseq2/column_data.csv')

# Run DESeq
col_data$condition <- factor(col_data$condition, levels = c("HBR", "UHRR"))
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
resLFC <- lfcShrink(dds, coef="condition_UHRR_vs_HBR", type="apeglm")

# Order shrunken LFC results by p-value
resOrderedLFC <- resLFC[order(resLFC$pvalue),]

# Write results to file
write.csv(as.data.frame(resOrdered), file='outputs/deseq2/deseq2_results.csv')

# # MA plot
# png(file=snakemake@output[["deseq2_ma_plot"]])
# plotMA(res, main="MA Plot", ylim=c(-2,2))
# dev.off()

# # Volcano plot
# png(file=snakemake@output[["deseq2_volcano_plot"]])
# with(resOrdered, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2,2)))
# abline(h = -log10(0.05), col="red")
# dev.off()