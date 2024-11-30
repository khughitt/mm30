#!/bin/env Rscript
#
# Combines expression data from individual experiments making up MM30
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

snek <- snakemake

# create a list of inidividual dataframes
infiles <- c(
  Sys.glob(snek@config$expr_data$geo),
  snek@config$expr_data$mmrf
)

# load data and merge into a single dataframe
dat <- lapply(infiles, read_feather) %>%
  purrr::reduce(full_join, by = "symbol")

# drop any genes which contain a significant number of missing values
num_nas <- apply(dat, 1, function(x) {
  sum(is.na(x))
})

# for now, we will only keep genes with < 50% missing values
mask <- num_nas < .5 * ncol(dat)
dat <- dat[mask, ]

# store combined dataset
write_feather(dat, snek@output[[1]])
