#!/bin/env Rscript
#
# compute_metafor_scores.R
#
suppressMessages(library(arrow))
suppressMessages(library(metafor))
suppressMessages(library(tidyverse))

snek <- snakemake

# load dataset gene- or pathway-level effect sizes, and standard errors calculated by fassoc, along
# with the experiment metadata
effects <- read_feather(snek@input[["effects"]])
errors <- read_feather(snek@input[["errors"]])
mdata <- read_feather(snek@input[["mdata"]])

id_field <- colnames(effects)[1]

# normalize contributions from each dataset?
# if enabled, the max effect size associated with each dataset will be used
if (snek@config$normalize_dataset_contributions) {
  feature_ids <- pull(effects, id_field)

  max_effect_list <- list(feature_ids)
  max_effect_error_list <- list(feature_ids)

  dataset_ids <- unique(str_split(colnames(effects)[-1], "_", simplify=TRUE)[, 1])

  for (id_ in dataset_ids) {
    # get columns corresponding to dataset
    mask <- startsWith(colnames(effects), id_)

    if (sum(mask) == 1) {
      max_effect_list <- c(max_effect_list, list(pull(effects[, mask], 1)))
      max_effect_error_list <- c(max_effect_error_list, list(pull(errors[, mask], 1)))
    } else {
      # get max effect size + corresponding std error
      effects_subset <- as.matrix(effects[, mask, drop=FALSE])
      errors_subset <- as.matrix(errors[, mask, drop=FALSE])

      max_ind <- max.col(abs(effects_subset), "first")

      max_effects <- c()
      max_effect_errors <- c()

      for (i in seq_len(nrow(effects))) {
        if (is.na(max_ind[i])) {
          max_effects <- c(max_effects, NA)
          max_effect_errors <- c(max_effect_errors, NA)
        } else {
          max_effects <- c(max_effects, as.numeric(effects_subset[i, max_ind[i]]))
          max_effect_errors <- c(max_effect_errors, as.numeric(errors_subset[i, max_ind[i]]))
        }
      }
      max_effect_list <- c(max_effect_list, list(c(max_effects)))
      max_effect_error_list <- c(max_effect_error_list, list(c(max_effect_errors)))
    }
  }

  # convert back to tibbles
  names(max_effect_list) <- c(id_field, dataset_ids)
  names(max_effect_error_list) <- c(id_field, dataset_ids)

  effects <- as_tibble(max_effect_list)
  errors <- as_tibble(max_effect_error_list)
}

# create matrix versions of each without the id columns
effect_mat <- effects %>%
  select(-all_of(id_field)) %>%
  as.matrix()

error_mat <- errors %>%
  select(-all_of(id_field)) %>%
  as.matrix()

# replace infinite values with "NA"
effect_mat[is.infinite(effect_mat)] <- NA
error_mat[is.infinite(error_mat)] <- NA

# replace "0" effect sizes with "NA"
effect_mat[effect_mat == 0] <- NA

# sanity check
if (!all(colnames(effect_mat) %in% mdata$dataset)) {
  stop("Missing metadata for one or more datasets!")
}

# count number of missing effect sizes for each gene or gene set
num_missing <- apply(effect_mat, 1, function(x) {
  sum(is.na(x))
})
num_present <- ncol(effect_mat) - num_missing

# effect size aggregation
metafor_pvals <- c()

for (i in seq_len(nrow(effect_mat))) {
  gene_effects <- as.numeric(effect_mat[i, ])
  gene_std_errors <- as.numeric(error_mat[i, ])

  # log-transform hazard ratios
  if (snakemake@wildcards$category %in% c("survival_os", "survival_pfs")) {
    gene_effects <- log(gene_effects)
  }

  pval <- tryCatch({
    # coef(summary(fit))
    #         estimate       se     zval      pval     ci.lb    ci.ub
    # intrcpt 17.10576 16.90056 1.012141 0.3114705 -16.01874 50.23026
    fit <- rma(gene_effects, sei=gene_std_errors)
    coefs <- coef(summary(fit))
    coefs[, "pval"]
  }, error=function(e) {
    NA
  })

  metafor_pvals <- c(metafor_pvals, pval)
}

# construct summary dataframe
res <- data.frame(
  pull(effects, id_field),
  mean_effect   = apply(effect_mat, 1, mean, na.rm=TRUE),
  median_effect = apply(effect_mat, 1, median, na.rm=TRUE),
  mean_error    = apply(error_mat, 1, mean, na.rm=TRUE),
  median_error  = apply(error_mat, 1, median, na.rm=TRUE),
  metafor_pval  = metafor_pvals,
  num_present,
  num_missing
)
colnames(res)[1] <- id_field

# drop any genes with missing values for the aggregated scores
res <- res[complete.cases(res), ]

# reorder and store results
res <- res %>%
  arrange(metafor_pval)

# store results
write_feather(res, snek@output[[1]])
