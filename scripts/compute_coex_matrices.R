#
# compute_coex_matrices
#
library(arrow)
library(tidyverse)

snek <- snakemake

feat_name <- ifelse(snek@wildcards$feat_level == "gene", "symbol", "gene_set")

expr <- read_feather(snek@input[[1]]) %>%
  column_to_rownames(feat_name)

cor_mat <- cor(t(expr), use='pairwise.complete')

cor_mat %>%
  as.data.frame() %>%
  rownames_to_column(feat_name) %>%
  write_feather(snek@output[[1]])
