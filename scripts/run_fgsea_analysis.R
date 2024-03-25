library(fgsea)
library(data.table)
library(ggplot2)
library(org.Hs.eg.db)

deseq2Results <- read.csv("outputs/deseq2/deseq2_results.csv", header = TRUE)
ranks <- with(deseq2Results, setNames(log2FoldChange, X))
ranks <- ranks[is.finite(ranks)]

# Load the Reactome pathways
ensembl_ids <- names(ranks)
ensembl_ids_no_version <- gsub("\\..*$", "", ensembl_ids)
entrez_ids <- select(org.Hs.eg.db, keys=ensembl_ids_no_version, keytype="ENSEMBL", column="ENTREZID")
valid_entrez_ids <- entrez_ids[!is.na(entrez_ids$ENTREZID), ]
pathways <- reactomePathways(valid_entrez_ids$ENTREZID)
set.seed(42)

# Map Ensembl IDs to Entrez IDs
ranks_ensembl <- names(ranks)
version <- gsub(".*\\.", "", ranks_ensembl)
ensembl_ids <- gsub("\\..*", "", ranks_ensembl)
matched_ids <- match(ensembl_ids, valid_entrez_ids$ENSEMBL)
matched_ids <- matched_ids[!is.na(matched_ids)]
mapped_entrez_ids <- valid_entrez_ids$ENTREZID[matched_ids]
mapped_data <- data.frame(ENSEMBL = ensembl_ids[matched_ids], ENTREZID = mapped_entrez_ids, Version = version[matched_ids], Rank = ranks[ranks_ensembl[matched_ids]])
rank_vector <- setNames(mapped_data$Rank, as.character(mapped_data$ENTREZID))

# Run fgsea
fgseaRes <- fgsea(pathways = pathways, stats = rank_vector, minSize  = 5, maxSize  = 1000)

# Save the enrichment results to a CSV file
enrichment_results <- fgseaRes[, .(pathway, pval, padj, log2err, ES, NES, size)]
write.csv(enrichment_results, "outputs/deseq2/enrichment_results.csv", row.names = FALSE)

# Order the results by p-value and preview
head(fgseaRes[order(pval), ])

# Plot the enrichment results
plot <- plotEnrichment(pathways[["Neutrophil degranulation"]], rank_vector) + labs(title="Neutrophil degranulation")

# Save the plot image
ggsave("outputs/deseq2/fgsea_analysis_plot.png", plot, width = 8, height = 6)