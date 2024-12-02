#!/bin/env Rscript
#
# Aggregates and summarizes survival effect sizes and standard errors from multiple sources
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

# load dataset gene- or pathway-level test effect sizes and standard errors
effects <- read_feather(snek@input[["effects"]])
errors <- read_feather(snek@input[["errors"]])

# "genes" or "gene sets"
id_field <- colnames(effects)[1]

# load feature-phenotype association metadata
mdata <- read_feather(snek@input[["mdata"]])

# get a list of associations of the desired category
surv_subset <- mdata %>%
  filter(category == "survival")

# remove covariates that are not in the specified category
cols_to_keep <- sprintf("%s_%s", surv_subset$dataset, surv_subset$phenotype)
cols_to_keep <- c(id_field, cols_to_keep)

mask <- colnames(effects) %in% cols_to_keep

effects <- effects[, mask]
errors <- errors[, mask]

# drop any features that no longer have any non-missing values after filtering
num_non_na <- apply(effects, 1, function(x) {
  sum(!is.na(x))
})

effects <- effects[num_non_na > 1, ]
errors <- errors[num_non_na > 1, ]

########################

# normalize contributions from each dataset, if enabled;
# note: for some metap methods, weights can also be specified for each p-value..
if (snek@config$normalize_dataset_contributions) {
  max_effect_list <- list(pull(effects, id_field))
  max_error_list <- list(pull(errors, id_field))

  dataset_ids <- unique(str_split(colnames(effects)[-1], "_", simplify = TRUE)[, 1])

  for (id_ in dataset_ids) {
    # collapse columns for dataset into a single column
    mask <- startsWith(colnames(effects), id_)

    # get maximum effectistics for each dataset
    max_effects <- suppressWarnings(apply(effects[, mask, drop = FALSE], 1, function (x) {
      x[which.max(abs(x))]
    }))

    # convert nulls to na's so that they are preserved during call to "unlist" and
    # flatten
    is.na(max_effects) <- lapply(max_effects, length) == 0
    max_effects <- unlist(max_effects)

    max_effects[is.infinite(max_effects)] <- NA
    max_effect_list <- c(max_effect_list, list(c(max_effects)))

    # get maximum std errors for each dataset
    max_errors <- suppressWarnings(apply(errors[, mask, drop = FALSE], 1, function (x) {
      x[which.max(abs(x))]
    }))

    # convert nulls to na's so that they are preserved during call to "unlist" and
    # flatten
    is.na(max_errors) <- lapply(max_errors, length) == 0
    max_errors <- unlist(max_errors)

    max_errors[is.infinite(max_errors)] <- NA
    max_error_list <- c(max_error_list, list(c(max_errors)))
  } 

  # convert back to a tibble
  names(max_effect_list) <- c(id_field, dataset_ids)
  names(max_error_list) <- c(id_field, dataset_ids)

  effects <- as_tibble(max_effect_list)
  errors <- as_tibble(max_error_list)
}

# create matrix versions of the effect/error matrices without the id column
effect_mat <- effects %>%
  select(-all_of(id_field)) %>%
  as.matrix()

error_mat <- errors %>%
  select(-all_of(id_field)) %>%
  as.matrix()

# count number of missing values for each gene or gene set
num_missing <- apply(effect_mat, 1, function(x) {
  sum(is.na(x))
})
num_present <- ncol(effect_mat) - num_missing

# summary dataframe containing aggregated test effects / errors
res <- data.frame(
  pull(effects, id_field),
  mean_effect    = apply(effect_mat, 1, mean, na.rm = TRUE),
  median_effect  = apply(effect_mat, 1, median, na.rm = TRUE),
  mean_error     = apply(error_mat, 1, mean, na.rm = TRUE),
  median_error   = apply(error_mat, 1, median, na.rm = TRUE),
  num_present,
  num_missing
)
colnames(res)[1] <- id_field

# store results
write_feather(res, snek@output[[1]])
