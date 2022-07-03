#!/bin/env Rscript
#
# Combines expression data from individual experiments making up MM29
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# create a list of inidividual dataframes
infiles <- c(
  Sys.glob(snakemake@config$expr_data$geo),
  snakemake@config$expr_data$mmrf
)

# load data and merge into a single dataframe
dat <- lapply(infiles, read_feather) %>%
  purrr::reduce(full_join, by = "symbol")

#dim(dat)
# 102507   5749

# drop any genes which contain a significant number of missing values
num_nas <- apply(dat, 1, function(x) {
  sum(is.na(x))
})

# for now, we will only keep genes with < 50% missing values
mask <- num_nas < .5 * ncol(dat)
dat <- dat[mask, ]

#table(mask)
# mask
# FALSE  TRUE
# 81244 21263

# store combined dataset
write_feather(dat, snakemake@output[[1]])
