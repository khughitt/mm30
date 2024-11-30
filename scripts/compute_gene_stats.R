#!/bin/env Rscript
#
# Computes basic gene-level statistics within each dataset:
#  1. mean expr
#  2. median expr
#  3. expr var
#  4. ratio non-zero expr
#  5. coefficient of variation
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

# lists to store stats
expr_mean <- list()
expr_median <- list()
expr_var <- list()
expr_cv <- list()
expr_ratio_nonzero <- list()

infiles <- Sys.glob(snek@config$expr_data$geo)

ratio_nonzero <- function(x) {
  sum(x != 0, na.rm=TRUE) / length(x)
}

# GEO
for (infile in infiles) {
  acc <- basename(dirname(infile))

  df <- read_feather(infile) %>%
    column_to_rownames('symbol')

  gene_means <- apply(df, 1, mean, na.rm=TRUE)
  gene_medians <- apply(df, 1, median, na.rm=TRUE)
  gene_vars <- apply(df, 1, var, na.rm=TRUE)
  gene_ratio_nonzeros <- apply(df, 1, ratio_nonzero)
  gene_cvs <- sqrt(gene_vars) / gene_means

  expr_mean[[acc]] <- enframe(gene_means, name="symbol", value=acc)
  expr_median[[acc]] <- enframe(gene_medians, name="symbol", value=acc)
  expr_var[[acc]] <- enframe(gene_vars, name="symbol", value=acc)
  expr_ratio_nonzero[[acc]] <- enframe(gene_ratio_nonzeros, name="symbol", value=acc)
  expr_cv[[acc]] <- enframe(gene_cvs, name="symbol", value=acc)
}

# MMRF
infile <- snek@config$expr_data$mmrf

df <- read_feather(infile) %>%
  column_to_rownames('symbol')

gene_means <- apply(df, 1, mean, na.rm=TRUE)
gene_medians <- apply(df, 1, median, na.rm=TRUE)
gene_vars <- apply(df, 1, var, na.rm=TRUE)
gene_ratio_nonzeros <- apply(df, 1, ratio_nonzero)
gene_cvs <- sqrt(gene_vars) / gene_means

expr_mean[["mmrf"]] <- enframe(gene_means, name="symbol", value="mmrf")
expr_median[["mmrf"]] <- enframe(gene_medians, name="symbol", value="mmrf")
expr_var[["mmrf"]] <- enframe(gene_vars, name="symbol", value="mmrf")
expr_ratio_nonzero[["mmrf"]] <- enframe(gene_ratio_nonzeros, name="symbol", value="mmrf")
expr_cv[["mmrf"]] <- enframe(gene_cvs, name="symbol", value="mmrf")

# save results
expr_mean %>%
  purrr::reduce(full_join, by="symbol") %>%
  write_feather(snakemake@output[[1]])

expr_median %>%
  purrr::reduce(full_join, by="symbol") %>%
  write_feather(snakemake@output[[2]])

expr_var %>%
  purrr::reduce(full_join, by="symbol") %>%
  write_feather(snakemake@output[[3]])

expr_cv %>%
  purrr::reduce(full_join, by="symbol") %>%
  write_feather(snakemake@output[[4]])

expr_ratio_nonzero %>%
  purrr::reduce(full_join, by="symbol") %>%
  write_feather(snakemake@output[[5]])
