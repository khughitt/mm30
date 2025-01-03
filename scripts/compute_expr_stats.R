#!/bin/env Rscript
#
# Computes basic gene- or gene set-level statistics within each dataset:
#  1. mean expr
#  2. median expr
#  3. expr var
#  4. ratio non-zero expr
#  5. coefficient of variation
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

feat_level <- snek@wildcards$feat_level
feat_name <- ifelse(feat_level == "gene", "symbol", "gene_set")

# lists to store stats
expr_mean <- list()
expr_median <- list()
expr_var <- list()
expr_cv <- list()
expr_ratio_nonzero <- list()
expr_ratio_missing <- list()

geo_infiles <- Sys.glob(snek@config$expr_data[[feat_level]][["geo"]])

ratio_nonzero <- function(x) {
  sum(x != 0, na.rm=TRUE) / length(x)
}

ratio_missing <- function(x) {
  sum(is.na(x)) / length(x)
}

# GEO
for (infile in geo_infiles) {
  acc <- stringr::str_extract(infile, "GSE[0-9]+")

  df <- read_feather(infile)

  df <- df %>%
    column_to_rownames(feat_name)

  dataset_means <- apply(df, 1, mean, na.rm=TRUE)
  dataset_medians <- apply(df, 1, median, na.rm=TRUE)
  dataset_vars <- apply(df, 1, var, na.rm=TRUE)
  dataset_cvs <- sqrt(dataset_vars) / dataset_means
  dataset_ratio_nonzero <- apply(df, 1, ratio_nonzero)
  dataset_ratio_missing <- apply(df, 1, ratio_missing)

  expr_mean[[acc]]   <- enframe(dataset_means, name=feat_name, value=acc)
  expr_median[[acc]] <- enframe(dataset_medians, name=feat_name, value=acc)
  expr_var[[acc]]    <- enframe(dataset_vars, name=feat_name, value=acc)
  expr_cv[[acc]]     <- enframe(dataset_cvs, name=feat_name, value=acc)
  expr_ratio_nonzero[[acc]] <- enframe(dataset_ratio_nonzero, name=feat_name, value=acc)
  expr_ratio_missing[[acc]] <- enframe(dataset_ratio_missing, name=feat_name, value=acc)
}

# MMRF
infile <- snek@config$expr_data[[feat_level]][["mmrf"]]

df <- read_feather(infile) %>%
  column_to_rownames(feat_name)

mmrf_means <- apply(df, 1, mean, na.rm=TRUE)
mmrf_medians <- apply(df, 1, median, na.rm=TRUE)
mmrf_vars <- apply(df, 1, var, na.rm=TRUE)
mmrf_ratio_nonzero <- apply(df, 1, ratio_nonzero)
mmrf_ratio_missing <- apply(df, 1, ratio_missing)
mmrf_cvs <- sqrt(mmrf_vars) / mmrf_means

expr_mean[["mmrf"]]     <- enframe(mmrf_means, name=feat_name, value="mmrf")
expr_median[["mmrf"]]   <- enframe(mmrf_medians, name=feat_name, value="mmrf")
expr_var[["mmrf"]]      <- enframe(mmrf_vars, name=feat_name, value="mmrf")
expr_cv[["mmrf"]]       <- enframe(mmrf_cvs, name=feat_name, value="mmrf")
expr_ratio_nonzero[["mmrf"]] <- enframe(mmrf_ratio_nonzero, name=feat_name, value="mmrf")
expr_ratio_missing[["mmrf"]] <- enframe(mmrf_ratio_missing, name=feat_name, value="mmrf")

# save results
expr_mean %>%
  purrr::reduce(full_join, by=feat_name) %>%
  write_feather(snek@output[[1]])

expr_median %>%
  purrr::reduce(full_join, by=feat_name) %>%
  write_feather(snek@output[[2]])

expr_var %>%
  purrr::reduce(full_join, by=feat_name) %>%
  write_feather(snek@output[[3]])

expr_cv %>%
  purrr::reduce(full_join, by=feat_name) %>%
  write_feather(snek@output[[4]])

expr_ratio_nonzero %>%
  purrr::reduce(full_join, by=feat_name) %>%
  write_feather(snek@output[[5]])

expr_ratio_missing %>%
  purrr::reduce(full_join, by=feat_name) %>%
  write_feather(snek@output[[6]])
