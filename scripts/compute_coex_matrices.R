#
# compute_coex_matrices
#
library(arrow)
library(tidyverse)

snek <- snakemake

expr <- read_feather(snek@input[[1]]) %>%
  column_to_rownames("symbol")

cor_mat <- cor(t(expr), use='pairwise.complete')

cor_mat %>%
  as.data.frame() %>%
  rownames_to_column('symbol') %>%
  write_feather(snek@output[[1]])
