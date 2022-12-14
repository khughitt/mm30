#!/bin/env Rscript
#
# Creates subsetted versions of the full MM29 p-value / test statistic results which
# include either patient or cell line samples.
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# load dataset gene- or pathway-level p-values calculcated by fassoc
pvals <- read_feather(snakemake@input[['pvals']])
stats <- read_feather(snakemake@input[['stats']])
coefs <- read_feather(snakemake@input[['coefs']])

# "genes" or "gene sets"
id_field <- colnames(pvals)[1]

# load experiment-level metadata
experiment_mdata <- read_tsv("metadata/mm29_experiment_metadata.tsv")

# get relevant experiment ids; excludes single "mixed" experiment
if (snakemake@wildcards$sample_type == "patient") {
  experiment_ids <- experiment_mdata %>%
    filter(sample_type == "Patient") %>%
    pull(accession)

  # add version-less MMRF accession
  experiment_ids <- c(experiment_ids, "MMRF")
} else {
  experiment_ids <- experiment_mdata %>%
    filter(sample_type == "Cell line") %>%
    pull(accession)
}

# load covariate-level metadata
mdata <- read_feather(snakemake@input[['mdata']])

# get a list of entries for the target sample type
pheno_subset <- mdata %>%
  filter(dataset %in% experiment_ids)

# remove covariates that are not of this type
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
