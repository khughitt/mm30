#
# create_gene_survival_table.R
#
library(arrow)
library(tidyverse)

effects <- read_feather(snakemake@input[[1]])
errors <- read_feather(snakemake@input[[2]])
pvals <- read_feather(snakemake@input[[3]])
mdata <- read_feather(snakemake@input[[4]])

# get a list of associations of the desired category
mdata <- mdata %>%
  filter(category == "survival")

# remove covariates that are not in the specified category
cols_to_keep <- sprintf("%s_%s", mdata$dataset, mdata$phenotype)
cols_to_keep <- c("symbol", cols_to_keep)

mask <- colnames(effects) %in% cols_to_keep

effects <- effects[, mask]
errors <- errors[, mask]
pvals <- pvals[, mask]

# create long versions
effects <- effects %>% 
  pivot_longer(-symbol, names_to="dataset_covariate", values_to="effect")

errors <- errors %>% 
  pivot_longer(-symbol, names_to="dataset_covariate", values_to="stderror")

pvals <- pvals %>% 
  pivot_longer(-symbol, names_to="dataset_covariate", values_to="pval")

res <- effects

# split dataset & covariate columns
res$dataset <- unlist(lapply(strsplit(pvals$dataset_covariate, "_"), "[", 1))
res$covariate <- sub('[A-Z0-9]+_', '', res$dataset_covariate)

# add error terms & p-values
res$stderror <- errors$stderror
res$pval <- pvals$pval

# drop entries associated with missing values
res <- res[complete.cases(res), ]

res %>%
  select(-dataset_covariate) %>%
  select(symbol, dataset, covariate, effect, stderror, pval) %>%
  write_feather(snakemake@output[[1]])
