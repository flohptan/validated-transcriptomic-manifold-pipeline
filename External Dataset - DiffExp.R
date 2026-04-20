#--------------------------------------------
# EXTERNAL DATASET VALIDATION
#--------------------------------------------

library(dplyr)
library(readr)
library(purrr)
library(tidyr)


# 1. Identify the specific GSM files (picking up 238 to 243, including -a and -b)
target_pattern <- "GSM297824[0-3]|GSM297823[8-9]"
files_to_merge <- list.files(path = path_to_files, pattern = target_pattern, full.names = TRUE)

# 2. Read, Label, and Sum
external_data <- files_to_merge %>%
  map_df(function(f) {
    # Read the file
    df <- read_tsv(f, col_names = c("GeneID", "Count"), show_col_types = FALSE)
    
    # Get just the filename 
    filename <- basename(f)
    # Extract just the first 10 characters (the GSM ID)
    gsm_id <- substr(filename, 1, 10) 
    
    df %>% mutate(SampleID = gsm_id)
  }) %>%
  group_by(GeneID, SampleID) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = TotalCount)

# 3. Save the final result
write_csv(external_data, file.path(path_to_files, "GSE110135_validation_counts.csv"))


#--------------------------------------------
# PIPELINE VALIDATION
#--------------------------------------------

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)

# ==========================================================
# Load merged count data
# ==========================================================
countdata <- read.csv(
  "GSE110135_validation_counts.csv",
  row.names = 1,
  check.names = FALSE
)

countdata <- as.matrix(countdata)

# ==========================================================
# Experimental design 
# ==========================================================
condition <- factor(
  c(rep("Young", 3), rep("Old", 3)),
  levels = c("Young", "Old")
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
# Save unfiltered results
# ==========================================================
write.csv(
  as.data.frame(res),
  file = "YoungvsOld_Deseq_output_unfiltered.csv"
)

# ==========================================================
# DEG filtering
# ==========================================================
deg_idx <- which(
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > 1
)

res_DEG <- res[deg_idx, ]

write.csv(
  as.data.frame(res_DEG),
  file = "YoungvsOld_DEGs.csv"
)

deg_genes <- rownames(res_DEG)

raw_counts_DEG <- countdata[deg_genes, ]

write.csv(
  raw_counts_DEG,
  file = "YoungvsOld_DEG_raw_counts.csv"
)

# Load DEG counts
countdata <- read.csv(
  "YoungvsOld_DEG_raw_counts.csv",
  row.names = 1,
  check.names = FALSE
)

countdata <- as.matrix(countdata)

# Recreate coldata
condition <- factor(
  c(rep("Young", 3), rep("Old", 3)),
  levels = c("Young", "Old")
)

coldata <- data.frame(
  row.names = colnames(countdata),
  condition
)

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ condition
)

# rlog
rld <- rlog(dds, blind = FALSE)

plotPCA(rld, intgroup = "condition") +
  ggtitle("PCA (rlog, external dataset)")

colors <- colorRampPalette(
  rev(brewer.pal(9, "Blues"))
)(255)

sampleDistsRlog <- dist(t(assay(rld)))

pheatmap(
  as.matrix(sampleDistsRlog),
  clustering_distance_rows = sampleDistsRlog,
  clustering_distance_cols = sampleDistsRlog,
  col = colors,
  main = "Sample distances (external dataset)"
)

rlog_counts <- assay(rld)

write.csv(
  rlog_counts,
  file = "YoungvsOld_DEG_rlog_counts.csv"
)