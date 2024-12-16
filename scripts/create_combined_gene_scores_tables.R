#
# create_combined_gene_scores_tables.R
#
library(arrow)
library(tidyverse)

snek <- snakemake

category <- snek@wildcards[["category"]]

metap <- read_feather(snek@input[[1]])

metafor <- read_feather(snek@input[[2]]) %>%
  select(-num_present, -num_missing)

metap %>%
  left_join(metafor) %>%
  select(symbol, ends_with("pval"), starts_with("num_sig"), everything()) %>%
  write_feather(snek@output[[1]])
