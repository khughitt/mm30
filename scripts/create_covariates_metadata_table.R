#
# create_covariates_metadata_table.R
#
library(arrow)
library(tidyverse)

snek <- snakemake

dataset_mdata <- read_feather(snek@input[[1]])

infile <- file.path(snek@config$fassoc_dir, "metadata/association_metadata.feather")

covariate_mdata <- read_feather(infile) %>%
  filter(feature_level == "genes") %>%
  dplyr::rename(num_genes=num_features) %>%
  dplyr::select(-feature_level, -feature_path, -phenotype_path) %>%
  left_join(dataset_mdata, by="dataset")

# add covariate id + category
covariate_mdata$covariate_id <- sprintf("%s_%s", covariate_mdata$dataset, covariate_mdata$phenotype)
covariate_mdata$category <- factor(covariate_mdata$category)

covariate_mdata %>%
  write_feather(snek@output[[1]])
