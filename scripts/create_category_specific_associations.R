#!/bin/env Rscript
#
# Creates subsetted versions of the full MM30 p-value, effect size, and standard error tables
# containing only covariates associated with a particular category (e.g. "disease stage")
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

# load dataset gene- or pathway-level p-values calculcated by fassoc
pvals <- read_feather(snek@input[["pvals"]])
effects <- read_feather(snek@input[["effects"]])
errors <- read_feather(snek@input[["errors"]])

# "genes" or "gene sets"
id_field <- colnames(pvals)[1]

# load feature-phenotype association metadata
mdata <- read_feather(snek@input[["mdata"]])

# determine which subset of the data to include based on the requested data subset
category <- snek@wildcards$category

if (category == "disease_stage") {
  pheno_subset <- mdata %>%
    filter(category == "disease_stage")
} else if (category == "survival_os") {
  mask <- grepl("overall", mdata$phenotype)
  pheno_subset <- mdata[mask, ]
} else if (category == "survival_pfs") {
  mask <- grepl("prog_free|pfs|event", mdata$phenotype)
  pheno_subset <- mdata[mask, ]
} else if (category == "treatment_response") {
  mask <- grepl("treatment_response|first_response", mdata$phenotype)
  pheno_subset <- mdata[mask, ]
}

# remove covariates that are not in the specified category
cols_to_keep <- unique(sprintf("%s_%s", pheno_subset$dataset, pheno_subset$phenotype))
cols_to_keep <- c(id_field, cols_to_keep)

mask <- colnames(pvals) %in% cols_to_keep

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
write_feather(pvals, snek@output[["pvals"]])
write_feather(effects, snek@output[["effects"]])
write_feather(errors, snek@output[["errors"]])
