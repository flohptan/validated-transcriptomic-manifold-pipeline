# Project Title
**Integrating Theorem-Guided Unsupervised Learning and Transcriptomics to Identify Neurodegenerative Modules in an Aβ42 Drosophila Model**

Brief description: This project analyzes differential gene expression in Drosophila melanogaster, performs clustering, and quantifies transcriptomic complexity using the LMC framework.

---

## Directory Structure
/Project_Folder
│-- 01_Differential_Expression.R
│-- 02_DimRed_Clustering_BioRelev.ipynb
│-- 03_LMC_Complexity_Analysis.R
│-- FeatureCounts - Counts.csv
│-- dmel_human_orthologs_disease_fb_2025_02.tsv
│-- Outputs/ (folder containing CSVs, figures, etc.)
│-- README.md


---

## Input Files
- `FeatureCounts - Counts.csv`: Raw read counts for gene expression.  
- `dmel_human_orthologs_disease_fb_2025_02.tsv`: Ortholog/disease mapping for annotation.  

**Note:** These files should be in the same directory as the scripts for the code to run correctly.

---

## Code Files
1. `01_Differential_Expression.R`  
   - Performs differential expression analysis using rlog-transformed counts.  

2. `02_DimRed_Clustering_BioRelev.ipynb`  
   - Performs dimensionality reduction (UMAP) and K-means clustering.  
   - Visualizes clusters and annotates biological relevance.  

3. `03_LMC_Complexity_Analysis.R`  
   - Computes LMC statistical complexity (C_JSD and C_euc) for each cluster.  
   - Performs size-preserving permutation tests for significance.  
   - Generates Figures 6a–c (entropy vs disequilibrium, size vs C_JSD, rank-ordered C_JSD).

---

## Outputs
- Output files are stored in the `Outputs/` folder.  
- Includes CSVs of cluster assignments, complexity metrics with p-values, and figures.

---

## How to Run
1. Place the input files in the same directory as the scripts.  
2. Run the scripts in order:  
   1. `01_Differential_Expression.R`  
   2. `02_DimRed_Clustering_BioRelev.ipynb`  
   3. `03_LMC_Complexity_Analysis.R`  
3. All outputs (CSVs and figures) will be generated in the same directory as the scripts.

---

## Notes
- The R scripts use packages: `DESeq2`, `pheatmap`, `RColorBrewer`, `EnhancedVolcano`, `genefilter`, `dplyr`, `ggplot2`, `tibble`, `purrr`, `readr`.  
- The Python notebook uses `pandas`, `numpy`, `matplotlib`, `seaborn`, `umap-learn`, `scikit-learn`, `gprofiler-official`, `networkx`, `requests`, `itertools`.  
- C_JSD is emphasized due to robustness for compositional transcriptomic data.
