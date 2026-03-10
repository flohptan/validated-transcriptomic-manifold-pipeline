# ==========================================================
# Libraries
# ==========================================================

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
# 3. RAW CLUSTER COMPLEXITY
# ==========================================================
cluster_complexity <- df %>%
  group_by(KMeans_Cluster) %>%
  group_modify(~ cluster_stats(.x)) %>%
  ungroup() %>%
  mutate(
    KMeans_Cluster = as.integer(KMeans_Cluster),
    rank_C_jsd = rank(-C_jsd)
  )

write_csv(cluster_complexity, "cluster_complexity_full_final.csv")

# ==========================================================
# 4. COMPUTE PERMUTATIONS
# ==========================================================
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
perm_matrix <- replicate(Nperm, permute_once())

real_C <- cluster_complexity$C_jsd
p_values <- map_dbl(seq_along(real_C), function(i) mean(perm_matrix[i, ] >= real_C[i]))

cluster_complexity <- cluster_complexity %>%
  mutate(p_jsd = p_values)

write_csv(cluster_complexity, "cluster_complexity_with_pvalues.csv")

# ==========================================================
# 5. FIGURE: ETROPTY VS DISEQUILLIBRIUM
# ==========================================================
p1 <- ggplot(cluster_complexity, aes(H, D_jsd)) +
  geom_point(color = "#1f78b4", size = 3, alpha = 0.8) +
  geom_point(data = cluster_complexity %>% filter(KMeans_Cluster == 32),
             aes(H, D_jsd), color = "#e31a1c", size = 4) +
  geom_text(data = cluster_complexity %>% filter(KMeans_Cluster == 32),
            aes(label = "Cluster 32"), hjust = -0.2, vjust = -0.5) +
  labs(
    title = "Entropy vs Disequilibrium of Clusters",
    x = "Shannon Entropy (H)",
    y = "Jensen-Shannon Disequilibrium (D_jsd)"
  ) +
  theme_bw(base_size = 14)

ggsave("entropy_vs_disequilibrium.png", p1, width = 7, height = 5)

# ==========================================================
# 6. FIGURE: CLUSTER SIZE VS C_JSD
# ==========================================================
p2 <- ggplot(cluster_complexity, aes(size, C_jsd)) +
  geom_point(color = "#33a02c", size = 3, alpha = 0.8) +
  geom_point(data = cluster_complexity %>% filter(KMeans_Cluster == 32),
             aes(size, C_jsd), color = "#e31a1c", size = 4) +
  labs(
    title = "Cluster Size vs LMC Complexity (C_jsd)",
    x = "Cluster Size",
    y = "C_jsd"
  ) +
  theme_bw(base_size = 14)

ggsave("size_vs_Cjsd.png", p2, width = 7, height = 5)

# ==========================================================
# 7. FIGURE: RANK ORDER C_JSD
# ==========================================================
cluster_complexity <- cluster_complexity %>% arrange(desc(C_jsd)) %>% mutate(rank = row_number())

p3 <- ggplot(cluster_complexity, aes(rank, C_jsd)) +
  geom_line(color = "#1f78b4", size = 1) +
  geom_point(data = cluster_complexity %>% filter(KMeans_Cluster == 32),
             aes(rank, C_jsd), color = "#e31a1c", size = 4) +
  labs(
    title = "Rank-ordered LMC Complexity Across Clusters",
    x = "Rank",
    y = "C_jsd"
  ) +
  theme_bw(base_size = 14)

ggsave("rank_ordered_Cjsd.png", p3, width = 7, height = 5)

############################################################
# END
############################################################
