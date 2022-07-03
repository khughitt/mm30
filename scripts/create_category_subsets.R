#!/bin/env Rscript
#
# Creates subsetted versions of the full MM29 p-value / test statistic results which
# include only covariates associated with a particular category (e.g. "survival")
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# load dataset gene- or pathway-level p-values calculcated by fassoc
pvals <- read_feather(snakemake@input[['pvals']])
stats <- read_feather(snakemake@input[['stats']])
coefs <- read_feather(snakemake@input[['coefs']])

# "genes" or "gene sets"
id_field <- colnames(pvals)[1]

# load feature-phenotype association metadata
mdata <- read_feather(snakemake@input[['mdata']])

# get a list of associations of the desired category
pheno_subset <- mdata %>%
  filter(category == snakemake@wildcards$category)

# remove covariates that are not in the specified category
cols_to_keep <- sprintf("%s_%s", pheno_subset$dataset, pheno_subset$phenotype)
cols_to_keep <- c(id_field, cols_to_keep)

mask <- colnames(pvals) %in% cols_to_keep

pvals <- pvals[, mask]
stats <- stats[, mask]
coefs <- coefs[, mask]

# drop any features that no longer have any non-missing values after filtering
num_non_na <- apply(pvals, 1, function(x) {
  sum(!is.na(x))
})

pvals <- pvals[num_non_na > 1, ]
stats <- stats[num_non_na > 1, ]
coefs <- coefs[num_non_na > 1, ]

# store results
write_feather(pvals, snakemake@output[['pvals']])
write_feather(stats, snakemake@output[['stats']])
write_feather(coefs, snakemake@output[['coefs']])
