# ==========================================================
# Libraries
# ==========================================================
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)
library(genefilter)

# ==========================================================
# Load raw count data
# ==========================================================
countdata <- read.csv(
  "FeatureCounts - Counts.csv",
  header = TRUE,
  row.names = "Geneid",
  check.names = FALSE
)
countdata <- as.matrix(countdata)

# ==========================================================
# Experimental design
# ==========================================================
condition <- factor(
  c(rep("Control", 3), rep("Ab42", 3)),
  levels = c("Control", "Ab42")
)

coldata <- data.frame(
  row.names = colnames(countdata),
  condition
)

# ==========================================================
# Create DESeq2 dataset
# ==========================================================
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ condition
)

# ==========================================================
# Filter low-count genes 
# ==========================================================
keep <- rowSums(counts(dds) >= 5) >= 2
dds <- dds[keep, ]

# ==========================================================
# Run DESeq2
# ==========================================================
dds <- DESeq(dds)
res <- results(dds)

# ==========================================================
# SAVE UNFILTERED DESEQ2 OUTPUT
# ==========================================================
write.csv(
  as.data.frame(res),
  file = "ControlvsAb42_Deseq_output_unfiltered.csv"
)
# ==========================================================
# DEG filtering
# padj < 0.05 AND |log2FC| > 1
# ==========================================================
deg_idx <- which(
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1
)

res_DEG <- res[deg_idx, ]

# Save DEG table
write.csv(
  as.data.frame(res_DEG),
  file = "ControlvsAb42_DEGs.csv"
)

# ==========================================================
# Extract RAW COUNTS for DEGs only
# ==========================================================
deg_genes <- rownames(res_DEG)

raw_counts_DEG <- countdata[deg_genes, ]

write.csv(
  raw_counts_DEG,
  file = "ControlvsAb42_DEG_raw_counts.csv"
)

# ==========================================================
# Load DEG raw counts
# ==========================================================
countdata <- read.csv(
  "ControlvsAb42_DEG_raw_counts.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

countdata <- as.matrix(countdata)

# ==========================================================
# Experimental design
# ==========================================================
condition <- factor(
  c(rep("Control", 3), rep("Ab42", 3)),
  levels = c("Control", "Ab42")
)

coldata <- data.frame(
  row.names = colnames(countdata),
  condition
)

# ==========================================================
# Create DESeq2 object
# ==========================================================
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ condition
)

# ==========================================================
# Transformations (visualization only)
# ==========================================================
rld <- rlog(dds, blind = FALSE)

# ==========================================================
# PCA
# ==========================================================
plotPCA(rld, intgroup = "condition") +
  ggtitle("PCA (rlog, DEGs only)")

# ==========================================================
# Sample distance heatmaps
# ==========================================================
colors <- colorRampPalette(
  rev(brewer.pal(9, "Blues"))
)(255)

# rlog
sampleDistsRlog <- dist(t(assay(rld)))
pheatmap(
  as.matrix(sampleDistsRlog),
  clustering_distance_rows = sampleDistsRlog,
  clustering_distance_cols = sampleDistsRlog,
  col = colors,
  main = "Sample distances (rlog, DEGs)"
)

########################################################
# Export rlog matrix for downstream analysis (Python)
########################################################
rlog_counts <- assay(rld)

write.csv(
  rlog_counts,
  file = "ControlvsAb42_DEG_rlog_counts.csv"
)
