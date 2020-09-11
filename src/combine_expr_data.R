#!/bin/env Rscript
#
# Combines expression data from individual experiments making up MM25
#
library(annotables)
library(arrow)
library(tidyverse)

# create a list of inidividual dataframes
infiles <- c(
  Sys.glob(snakemake@config$expr_data$geo),
  snakemake@config$expr_data$mmrf
)

dat_list <- lapply(infiles, read_feather)

# update MMRF symbols from GRCh37 -> 38, where possible
MMRF_IND <- 25

ids <- dat_list[[MMRF_IND]]$symbol
ensgenes <- grch37$ensgene[match(ids, grch37$symbol)]

ind <- !ids %in% grch38$symbol & ids %in% grch37$symbol
ind <- ind & ensgenes %in% grch38$ensgene

dat_list[[MMRF_IND]]$symbol[ind] <- grch38$symbol[match(ensgenes[ind], grch38$ensgene)]

# merge into a single dataframe
dat <- dat_list %>%
  purrr::reduce(full_join, by = "symbol")

#dim(dat)
# [1] 63577  5641

# drop any genes which contain a significant number of missing values; these are
# most likely ones which could not be mapped to GRCh38
num_nas <- apply(dat, 1, function(x) {
  sum(is.na(x))
})

# for now, we will only keep genes with < 50% missing values
mask <- num_nas < .5 * ncol(dat)
dat <- dat[mask, ]

#table(mask)
# mask
# FALSE  TRUE 
# 38767 24810 

# average expression for any duplicated gene symbols; there should not be too many
#table(duplicated(dat$symbol))
# 
# FALSE  TRUE 
# 24590   220 

# split unique and repeated gene entries and average the duplicated ones
dups <- names(table(dat$symbol))[table(dat$symbol) > 1]
dup_mask <- dat$symbol %in% dups

# table(dup_mask)
# dup_mask
# FALSE  TRUE 
# 24375   435 

# dat_dups <- dat[mask, ] %>%
#   group_by(symbol) %>%
#   summarize_all(median, na.rm = TRUE)
dat_dups <- dat[dup_mask, ]

# summarize_all is slow when ncol(dat) is in the thousands or higher, so duplicate
# rows will be averaged manually
collapsed <- NULL

for (gid in unique(dat_dups$symbol)) {
  mask <- dat_dups$symbol == gid

  dat_subset <- dat_dups[mask, ] %>%
    select(-symbol)

  collapsed <- rbind(collapsed, apply(dat_subset, 2, median, na.rm = TRUE))
}

collapsed <- cbind(data.frame(symbol = unique(dat_dups$symbol)), collapsed)

dat <- rbind(dat[!dup_mask, ], collapsed)

# store combined dataset
write_feather(dat, snakemake@output[[1]])
