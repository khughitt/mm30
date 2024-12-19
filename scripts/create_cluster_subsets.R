#!/bin/env Rscript
#
# Creates subsetted versions of the full MM30 p-value / test effect size results which
# include only covariates associated with a particular cluster
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# load dataset gene- or gene set-level p-values calculcated by fassoc
pvals <- read_feather(snakemake@input[['pvals']])
effects <- read_feather(snakemake@input[['effects']])
errors <- read_feather(snakemake@input[['errors']])

# "genes" or "gene sets"
id_field <- colnames(pvals)[1]

# load covariate clustering
clusters <- read_feather(snakemake@input[['clusters']])

# get covariates in cluster
covariates <- clusters %>%
  filter(cluster == snakemake@wildcards$cluster_num) %>%
  pull(covariate)

mask <- colnames(pvals) %in% c(id_field, covariates)

pvals <- pvals[, mask]
effects <- effects[, mask]
errors <- errors[, mask]

# drop any features that no longer have any non-missing values after filtering
num_non_na <- apply(pvals, 1, function(x) {
  sum(!is.na(x))
})

pvals <- pvals[num_non_na > 1, ]
effects <- effects[num_non_na > 1, ]
errors <- errors[num_non_na > 1, ]

# store results
write_feather(pvals, snakemake@output[['pvals']])
write_feather(effects, snakemake@output[['effects']])
write_feather(errors, snakemake@output[['errors']])
