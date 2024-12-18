#
# create_combined_score_tables.R
#
library(arrow)
library(tidyverse)

snek <- snakemake

feat_level <- snek@wildcards$feat_level
feat_name <- ifelse(feat_level == "gene", "symbol", "gene_set")

metap <- read_feather(snek@input[[1]])

metafor <- read_feather(snek@input[[2]]) %>%
  select(-num_present, -num_missing)

means <- read_feather(snek@input[[3]])
medians <- read_feather(snek@input[[4]])
vars <- read_feather(snek@input[[5]])
cvs <- read_feather(snek@input[[6]])
ratio_nonzeros <- read_feather(snek@input[[7]])

ratio_na <- apply(means[, -1], 1, function(x) {
  sum(is.na(x)) / length(x)
})

df <- data.frame(
  mean = apply(means[, -1], 1, mean, na.rm=TRUE),
  median = apply(medians[, -1], 1, mean, na.rm=TRUE),
  var = apply(vars[, -1], 1, mean, na.rm=TRUE),
  cv = apply(cvs[, -1], 1, mean, na.rm=TRUE),
  ratio_nonzero = apply(ratio_nonzeros[, -1], 1, mean, na.rm=TRUE),
  ratio_missing = ratio_na
)

df[[feat_name]] <- means %>%
  pull(feat_name)

metap %>%
  left_join(metafor) %>%
  left_join(df) %>%
  select(1, ends_with("pval"), starts_with("num_sig"), everything()) %>%
  write_feather(snek@output[[1]])
