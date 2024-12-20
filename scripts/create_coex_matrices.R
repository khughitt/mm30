#
# compute_coex_matrices
#
library(arrow)
library(tidyverse)

snek <- snakemake

feat_name <- ifelse(snek@wildcards$feat_level == "gene", "symbol", "gene_set")

# generate expr subsets as before
ranked_feats <- read_feather(snek@input[[1]]) %>%
  pull(feat_name)

top100 <- head(ranked_feats, 100)
top500 <- head(ranked_feats, 500)
top1000 <- head(ranked_feats, 1000)
top5000 <- head(ranked_feats, 5000)

# compute co-expression matrix for the largest subset, and then use that as the basis for the
# smaller subsets

# 1) unscaled
expr <- read_feather(snek@input[[2]]) %>%
  column_to_rownames(feat_name)

cor_mat <- cor(t(expr), use='pairwise.complete') %>%
  as.data.frame()

# determine indices to use when subsetting coex matrix
ind100 <- rownames(cor_mat) %in% top100
ind500 <- rownames(cor_mat) %in% top500
ind1000 <- rownames(cor_mat) %in% top100
ind5000 <- rownames(cor_mat) %in% top5000

cor_mat[ind100, ind100] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[1]])

cor_mat[ind500, ind500] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[2]])

cor_mat[ind1000, ind1000] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[3]])

cor_mat[ind5000, ind5000] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[4]])

# 2) scaled
expr_scaled <- read_feather(snek@input[[3]]) %>%
  column_to_rownames(feat_name)

cor_mat <- cor(t(expr_scaled), use='pairwise.complete') %>%
  as.data.frame()

cor_mat[ind100, ind100] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[5]])

cor_mat[ind500, ind500] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[6]])

cor_mat[ind1000, ind1000] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[7]])

cor_mat[ind5000, ind5000] %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[8]])
