# ==========================================================
# Libraries
# ==========================================================
# Check what folders actually exist
list.dirs("C:/Users/Admin/Desktop/Flo - sem 3/CDS590/Paper/Codes")

# Then try again
setwd("C:/Users/Admin/Desktop/Flo - sem 3/CDS590/Paper/Codes/Project_Folder")
getwd()

library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)

# ==========================================================
# 1. LOAD DATA
# ==========================================================
rlog <- read_csv("ControlvsAb42_DEG_rlog_counts.csv")
clusters <- read_csv("kmeans_final_gene_clusters.csv")

# Ensure column names match
colnames(rlog)[1] <- "GeneID"

# Join expression data with cluster assignments
df <- rlog %>% left_join(clusters, by = "GeneID")
df <- df %>%
  rename(KMeans_Cluster = Cluster)

expr_cols <- colnames(rlog)[-1]  # all expression columns
# ==========================================================
# 2. UTILITY FUNCTIONS
# ==========================================================
# Convert vector to probability (positive + small offset)
to_prob <- function(vec) {
  v <- vec - min(vec) + 1e-9
  p <- v / sum(v)
  return(p)
}

# Shannon entropy
compute_entropy <- function(vec) {
  p <- to_prob(vec)
  -sum(p * log(p))
}

# Jensen-Shannon divergence (symmetric)
compute_jsd <- function(P, Q) {
  P <- P / sum(P)
  Q <- Q / sum(Q)
  M <- 0.5 * (P + Q)
  
  kl_div <- function(A, B) {
    valid <- A > 0
    sum(A[valid] * log(A[valid] / B[valid]))
  }
  
  0.5 * kl_div(P, M) + 0.5 * kl_div(Q, M)
}

# Disequilibrium
compute_diseq <- function(vec) {
  p <- to_prob(vec)
  n <- length(p)
  u <- rep(1/n, n)
  list(
    euc = sqrt(sum((p - u)^2)),
    jsd = sqrt(compute_jsd(p, u))
  )
}

# Cluster statistics
cluster_stats <- function(cluster_df) {
  mean_expr <- rowMeans(cluster_df[, expr_cols, drop = FALSE])
  H <- compute_entropy(mean_expr)
  D <- compute_diseq(mean_expr)
  
  tibble(
    size = nrow(cluster_df),
    H = H,
    D_euc = D$euc,
    D_jsd = D$jsd,
    C_euc = H * D$euc,
    C_jsd = H * D$jsd
  )
}

# ==========================================================
# BENCHMARKING: START TIMER
# ==========================================================
total_start <- Sys.time()

# ==========================================================
# 3. RAW CLUSTER COMPLEXITY
# ==========================================================
# Track time for raw calculation
raw_start <- Sys.time()
cluster_complexity <- df %>%
  group_by(KMeans_Cluster) %>%
  group_modify(~ cluster_stats(.x)) %>%
  ungroup() %>%
  mutate(
    KMeans_Cluster = as.integer(KMeans_Cluster),
    rank_C_jsd = rank(-C_jsd)
  )
raw_end <- Sys.time()
raw_time <- as.numeric(difftime(raw_end, raw_start, units = "secs"))

# ==========================================================
# 4. COMPUTE PERMUTATIONS (Most Intensive Step)
# ==========================================================
perm_start <- Sys.time()
cluster_sizes <- df %>% count(KMeans_Cluster) %>% pull(n)

permute_once <- function() {
  shuffled <- df %>%
    mutate(perm_cluster = sample(rep(seq_along(cluster_sizes), times = cluster_sizes)))
  
  shuffled %>%
    group_by(perm_cluster) %>%
    group_modify(~ cluster_stats(.x)) %>%
    ungroup() %>%
    pull(C_jsd)
}

set.seed(42)
Nperm <- 500
# replicate is the core computational loop
perm_matrix <- replicate(Nperm, permute_once())

real_C <- cluster_complexity$C_jsd
p_values <- map_dbl(seq_along(real_C), function(i) mean(perm_matrix[i, ] >= real_C[i]))

cluster_complexity <- cluster_complexity %>%
  mutate(p_jsd = p_values)
perm_end <- Sys.time()
perm_time <- as.numeric(difftime(perm_end, perm_start, units = "secs"))

# ==========================================================
# 4B. PERMUTATION STABILITY TEST (N=500 vs N=1000)
# ==========================================================

cat("\n--- PERMUTATION STABILITY TEST ---\n")

# Save ranking from N = 500
rank_500 <- cluster_complexity %>%
  arrange(desc(C_jsd)) %>%
  select(KMeans_Cluster, C_jsd) %>%
  mutate(rank_500 = row_number())

# Run N = 1000 permutations
set.seed(42)
Nperm_1000 <- 1000

perm_matrix_1000 <- replicate(Nperm_1000, permute_once())

real_C <- cluster_complexity$C_jsd

p_values_1000 <- sapply(seq_along(real_C), function(i) {
  mean(perm_matrix_1000[i, ] >= real_C[i])
})

# Compute ranking for N = 1000
rank_1000 <- cluster_complexity %>%
  mutate(p_jsd_1000 = p_values_1000) %>%
  arrange(desc(C_jsd)) %>%
  select(KMeans_Cluster, C_jsd) %>%
  mutate(rank_1000 = row_number())

# Merge rankings
rank_compare <- rank_500 %>%
  inner_join(rank_1000, by = c("KMeans_Cluster", "C_jsd"))

# Spearman correlation
cor_test_perm <- cor.test(rank_compare$rank_500, rank_compare$rank_1000, method = "spearman")

rho_perm <- round(cor_test_perm$estimate, 4)

cat("Spearman correlation (500 vs 1000):", rho_perm, "\n")
col_main <- "#1f77b4"
# Optional visualization (VERY strong for supplement)
p_stability <- ggplot(rank_compare, aes(rank_500, rank_1000)) +
  geom_point(color = col_main, size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = col_trend) +
  labs(
    title = "Permutation Stability (N=500 vs N=1000)",
    subtitle = paste("Spearman rho =", rho_perm),
    x = "Rank (N = 500)",
    y = "Rank (N = 1000)"
  ) +
  base_theme

ggsave("permutation_stability.png", p_stability, width = 7, height = 5)

cat("----------------------------------------\n")
# ==========================================================
# BENCHMARKING: END & PRINT PROOF
# ==========================================================
total_end <- Sys.time()
total_time_min <- as.numeric(difftime(total_end, total_start, units = "mins"))

cat("\n--- COMPUTATIONAL BENCHMARKS (PROOF) ---\n")
cat("System Info:", Sys.info()["sysname"], "| CPU:", Sys.info()["machine"], "\n")
cat("1. Raw Complexity Calc:", round(raw_time, 2), "seconds\n")
cat("2. 500 Permutations:", round(perm_time, 2), "seconds (", round(perm_time/60, 2), "mins )\n")
cat("TOTAL RUNTIME:", round(total_time_min, 2), "minutes\n")
cat("----------------------------------------\n")

write_csv(cluster_complexity, "cluster_complexity_with_pvalues.csv")

# ==========================================================
# GLOBAL FIGURE STYLE (CONSISTENT ACROSS ALL FIGURES)
# ==========================================================
col_main <- "#1f77b4"      # Blue (all clusters)
col_highlight <- "#d62728" # Red (Cluster 32)
col_trend <- "#555555"     # Grey (trend lines)

base_theme <- theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

# ==========================================================
# 5. FIGURE: ETROPTY VS DISEQUILLIBRIUM
# ==========================================================
p1 <- ggplot(cluster_complexity, aes(H, D_jsd)) +
  geom_point(color = col_main, size = 3, alpha = 0.7) +
  geom_point(
    data = cluster_complexity %>% filter(KMeans_Cluster == 32),
    aes(H, D_jsd),
    color = col_highlight, size = 4
  ) +
  geom_text(
    data = cluster_complexity %>% filter(KMeans_Cluster == 32),
    aes(label = "Cluster 32"),
    hjust = -0.2, vjust = -0.5
  ) +
  labs(
    title = "Entropy vs Disequilibrium of Clusters",
    x = "Shannon Entropy (H)",
    y = "Disequilibrium (√Jensen–Shannon Divergence)"
  ) +
  base_theme

ggsave("entropy_vs_disequilibrium.png", p1, width = 7, height = 5)
# ==========================================================
# 6. FIGURE: CLUSTER SIZE VS C_JSD
# ==========================================================
p2 <- ggplot(cluster_complexity, aes(size, C_jsd)) +
  geom_point(color = col_main, size = 3, alpha = 0.7) +
  geom_point(
    data = cluster_complexity %>% filter(KMeans_Cluster == 32),
    aes(size, C_jsd),
    color = col_highlight, size = 4
  ) +
  labs(
    title = "Cluster Size vs LMC Complexity",
    x = "Cluster Size (Number of Genes)",
    y = "LMC Complexity (H × √JSD)"
  ) +
  base_theme

ggsave("size_vs_Cjsd.png", p2, width = 7, height = 5)
# ==========================================================
# 7. FIGURE: RANK ORDER C_JSD
# ==========================================================
cluster_complexity <- cluster_complexity %>%
  arrange(desc(C_jsd)) %>%
  mutate(rank = row_number())

p3 <- ggplot(cluster_complexity, aes(rank, C_jsd)) +
  geom_line(color = col_main, size = 1) +
  geom_point(
    data = cluster_complexity %>% filter(KMeans_Cluster == 32),
    aes(rank, C_jsd),
    color = col_highlight, size = 4
  ) +
  labs(
    title = "Rank-ordered LMC Complexity Across Clusters",
    x = "Cluster Rank",
    y = "LMC Complexity (H × √JSD)"
  ) +
  base_theme

ggsave("rank_ordered_Cjsd.png", p3, width = 7, height = 5)
# ==========================================================
# 8. CORRELATION WITH BIOLOGICAL RELEVANCE (UPDATED)
# ==========================================================
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# 1. Export gprofiler results
enrichment_dir <- "gprofiler_results_kmeans_clusters"

# 2. List the files
files <- list.files(path = enrichment_dir, 
                    pattern = "kmeans_cluster_.*_enrichment_significant.csv", 
                    full.names = TRUE)

if(length(files) == 0) {
  stop("No files found! Check folder path.")
}

# 3. Read and extract the strongest signal (-log10 p-value)
biorelev_list <- lapply(files, function(f) {
  c_id <- as.integer(str_extract(basename(f), "\\d+"))
  df_temp <- read_csv(f, show_col_types = FALSE)
  
  if (nrow(df_temp) > 0) {
    min_p <- min(df_temp$p_value, na.rm = TRUE)
    return(tibble(KMeans_Cluster = c_id, neg_log10_p = -log10(max(min_p, 1e-300))))
  } else {
    return(NULL)
  }
}) %>% bind_rows()

# 4. Merge with LMC data
cor_data <- cluster_complexity %>%
  inner_join(biorelev_list, by = "KMeans_Cluster")

# 5. Spearman Correlation
cor_test <- cor.test(cor_data$C_jsd, cor_data$neg_log10_p, method = "spearman")
rho_val <- round(cor_test$estimate, 3)
c32_data <- filter(cor_data, KMeans_Cluster == 32)

# 6. Final Plot with Enhanced Labeling
p4 <- ggplot(cor_data, aes(x = C_jsd, y = neg_log10_p)) +
  geom_smooth(
    method = "lm",
    linetype = "dashed",
    color = col_trend,
    se = TRUE
  ) +
  geom_point(color = col_main, size = 3, alpha = 0.6) +
  
  geom_vline(
    xintercept = c32_data$C_jsd,
    linetype = "dotted",
    color = col_highlight,
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = c32_data$neg_log10_p,
    linetype = "dotted",
    color = col_highlight,
    alpha = 0.5
  ) +
  
  geom_point(
    data = c32_data,
    color = col_highlight,
    size = 5
  ) +
  
  geom_label(
    data = c32_data,
    aes(label = paste("Cluster 32\nRank:", rank_C_jsd)),
    nudge_y = 0.8,
    fontface = "bold",
    fill = "white"
  ) +
  
  labs(
    title = "LMC Complexity vs Biological Signal",
    subtitle = paste("Spearman rho =", rho_val),
    x = "LMC Complexity (H × √JSD)",
    y = "Biological Significance (-log10 adjusted p-value)"
  ) +
  base_theme

ggsave("complexity_vs_enrichment_final.png", p4, width = 7, height = 5)

# 7. Print Cluster 32 details 
print(paste("Cluster 32 Complexity:", round(c32_data$C_jsd, 4)))
print(paste("Cluster 32 -log10(p):", round(c32_data$neg_log10_p, 4)))
############################################################
# END
############################################################
