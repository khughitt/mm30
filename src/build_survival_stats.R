#!/bin/env Rscript
#
# Combines feature associations from multiple sources
#
# Note: for some useful notes and diagrams on choosing a method to combine p-values,
# refer to figures 4-5 in the [metap compare
# vignette](https://cran.r-project.org/web/packages/metap/vignettes/compare.pdf).
# ----------
#
suppressMessages(library(arrow))
suppressMessages(library(metap))
suppressMessages(library(tidyverse))

# load dataset gene- or pathway-level test statistics and model coefficients
stats <- read_feather(snakemake@input[['stats']])
coefs <- read_feather(snakemake@input[['coefs']])

# "genes" or "gene sets"
id_field <- colnames(stats)[1]

# load feature-phenotype association metadata
mdata <- read_feather(snakemake@input[['mdata']])

# get a list of associations of the desired category
surv_subset <- mdata %>%
  filter(method == 'survival')

# remove covariates that are not in the specified category
cols_to_keep <- sprintf("%s_%s", surv_subset$dataset, surv_subset$phenotype)
cols_to_keep <- c(id_field, cols_to_keep)

mask <- colnames(stats) %in% cols_to_keep

stats <- stats[, mask]
coefs <- coefs[, mask]

# drop any features that no longer have any non-missing values after filtering
num_non_na <- apply(stats, 1, function(x) {
  sum(!is.na(x))
})

stats <- stats[num_non_na > 1, ]
coefs <- coefs[num_non_na > 1, ]

########################

# normalize contributions from each dataset, if enabled;
# note: for some metap methods, weights can also be specified for each p-value..
if (snakemake@config$normalize_dataset_contributions) {
  max_stat_list <- list(pull(stats, id_field))
  max_coef_list <- list(pull(coefs, id_field))

  dataset_ids <- unique(str_split(colnames(stats)[-1], '_', simplify = TRUE)[, 1])

  for (id_ in dataset_ids) {
    # collapse columns for dataset into a single column
    mask <- startsWith(colnames(stats), id_)

    # get maximum statistics for each dataset
    max_stats <- suppressWarnings(apply(stats[, mask, drop = FALSE], 1, function (x) {
      x[which.max(abs(x))]
    }))

    # convert nulls to na's so that they are preserved during call to "unlist" and
    # flatten
    is.na(max_stats) <- lapply(max_stats, length) == 0
    max_stats <- unlist(max_stats)

    max_stats[is.infinite(max_stats)] <- NA
    max_stat_list <- c(max_stat_list, list(c(max_stats)))

    # get maximum coefficients for each dataset
    max_coefs <- suppressWarnings(apply(coefs[, mask, drop = FALSE], 1, function (x) {
      x[which.max(abs(x))]
    }))

    # convert nulls to na's so that they are preserved during call to "unlist" and
    # flatten
    is.na(max_coefs) <- lapply(max_coefs, length) == 0
    max_coefs <- unlist(max_coefs)

    max_coefs[is.infinite(max_coefs)] <- NA
    max_coef_list <- c(max_coef_list, list(c(max_coefs)))
  } 

  # convert back to a tibble
  names(max_stat_list) <- c(id_field, dataset_ids)
  names(max_coef_list) <- c(id_field, dataset_ids)

  stats <- as_tibble(max_stat_list)
  coefs <- as_tibble(max_coef_list)
}

# create matrix versions of the p-values without the id column
stat_mat <- stats %>%
  select(-all_of(id_field)) %>%
  as.matrix()

coef_mat <- coefs %>%
  select(-all_of(id_field)) %>%
  as.matrix()

# count number of missing values for each gene or gene set
num_missing <- apply(stat_mat, 1, function(x) {
  sum(is.na(x))
})
num_present <- ncol(stat_mat) - num_missing

# summary dataframe containing aggregated test stats / coefs
res <- data.frame(
  pull(stats, id_field),
  mean_stat    = apply(stat_mat, 1, mean, na.rm = TRUE),
  median_stat  = apply(stat_mat, 1, median, na.rm = TRUE),
  mean_coef    = apply(coef_mat, 1, mean, na.rm = TRUE),
  median_coef  = apply(coef_mat, 1, median, na.rm = TRUE),
  num_present,
  num_missing
)
colnames(res)[1] <- id_field

# store results
write_feather(res, snakemake@output[[1]])
