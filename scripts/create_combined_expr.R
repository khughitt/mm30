#!/bin/env Rscript
#
# Combines expression data from individual experiments making up MM30
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

# exclude genes below a specified threshold of significance;
# current approach: keep all genes with >= 1 P-values < 0.001.
gene_scores <- read_feather(snek@input[[1]])

to_keep <- gene_scores %>%
  filter(num_sig_p001 >= 1) %>%
  pull(symbol)

# create a list of inidividual dataframes
infiles <- c(
  Sys.glob(snek@config$expr_data$geo),
  snek@config$expr_data$mmrf
)

# load data and merge into a single dataframe
dat <- lapply(infiles, read_feather) %>%
  purrr::reduce(full_join, by = "symbol")

# create a size factor "scaled" version prior to filtering
dat_scaled <- dat
dat_scaled[, -1] <- sweep(dat_scaled[, -1], 2, colSums(dat_scaled[, -1], na.rm=TRUE), "/") * 1E6

# exclude low significance genes
dat <- dat %>%
  filter(symbol %in% to_keep)

dat_scaled <- dat_scaled %>%
  filter(symbol %in% to_keep)

# store combined datasets
write_feather(dat, snek@output[[1]])
write_feather(dat_scaled, snek@output[[2]])
