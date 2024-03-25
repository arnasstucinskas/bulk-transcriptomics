library(DESeq2)
library(ggplot2)

count_data <- read.csv("outputs/deseq2/count_matrix.csv", row.names = 1)
metadata <- read.csv("outputs/deseq2/column_data.csv")

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = DataFrame(metadata), design = ~ condition)

# Variance stabilizing transformation on the full dataset
vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)

# Run PCA on VST-transformed data
pca <- prcomp(t(assay(vst_data)))

# Prepare a data frame for plotting
pca_data <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])
metadata$Sample <- gsub("-", ".", metadata$Sample)
pca_data <- merge(pca_data, metadata, by = "Sample")

# Generate PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point() +
  theme_minimal() +
  ggtitle("PCA of DEGs by sample preparation method") +
  xlab(paste("PC1 - variance explained:", round(summary(pca)$importance[2,1]*100, 2), "%")) +
  ylab(paste("PC2 - variance explained:", round(summary(pca)$importance[2,2]*100, 2), "%")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))

# Save the plot
ggsave("outputs/deseq2/pca_plot_de_genes.png", width = 10, height = 8)