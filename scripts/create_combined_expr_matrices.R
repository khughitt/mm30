#!/bin/env Rscript
#
# Combines expression data from individual experiments making up MM30
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

feat_level <- snek@wildcards$feat_level
feat_name <- ifelse(feat_level == "gene", "symbol", "gene_set")

ranked_feats <- read_feather(snek@input[[1]]) %>%
  pull(feat_name)

top100 <- head(ranked_feats, 100)
top500 <- head(ranked_feats, 500)
top1000 <- head(ranked_feats, 1000)
top5000 <- head(ranked_feats, 5000)

# load individual expression matrices
infiles <- c(
  Sys.glob(snek@config$expr_data[[feat_level]][["geo"]]),
  snek@config$expr_data[[feat_level]][["mmrf"]]
)

# merge into a single df
dat <- lapply(infiles, read_feather) %>%
  purrr::reduce(full_join, by=feat_name)

# create a size factor scaled version
dat_scaled <- dat
dat_scaled[, -1] <- sweep(dat_scaled[, -1], 2, colSums(dat_scaled[, -1], na.rm=TRUE), "/") * 1E6

# 1) store unscaled expression matrices
dat %>%
  write_feather(snek@output[[1]])

dat %>%
  filter(!!rlang::sym(feat_name) %in% top100) %>%
  write_feather(snek@output[[2]])

dat %>%
  filter(!!rlang::sym(feat_name) %in% top500) %>%
  write_feather(snek@output[[3]])

dat %>%
  filter(!!rlang::sym(feat_name) %in% top1000) %>%
  write_feather(snek@output[[4]])

dat %>%
  filter(!!rlang::sym(feat_name) %in% top5000) %>%
  write_feather(snek@output[[5]])

# 2) store scaled expression matrices
dat_scaled %>%
  write_feather(snek@output[[6]])

dat_scaled %>%
  filter(!!rlang::sym(feat_name) %in% top100) %>%
  write_feather(snek@output[[7]])

dat_scaled %>%
  filter(!!rlang::sym(feat_name) %in% top500) %>%
  write_feather(snek@output[[8]])

dat_scaled %>%
  filter(!!rlang::sym(feat_name) %in% top1000) %>%
  write_feather(snek@output[[9]])

dat_scaled %>%
  filter(!!rlang::sym(feat_name) %in% top5000) %>%
  write_feather(snek@output[[10]])
