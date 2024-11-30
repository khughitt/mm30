#
# create gene summary tables
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

means <- read_feather(snakemake@input[[1]])
medians <- read_feather(snakemake@input[[2]])
vars <- read_feather(snakemake@input[[3]])
cvs <- read_feather(snakemake@input[[4]])
ratio_nonzeros <- read_feather(snakemake@input[[5]])

scores_all <- read_feather(snakemake@input[[6]]) %>%
  mutate(rank_all_sumz_wt = dense_rank(sumz_wt_pval),
         rank_all_sumz = dense_rank(sumz_pval)) %>%
  select(symbol, rank_all_sumz_wt, rank_all_sumz)

scores_stage <- read_feather(snakemake@input[[7]]) %>%
  mutate(rank_stage_sumz_wt = dense_rank(sumz_wt_pval),
         rank_stage_sumz = dense_rank(sumz_pval)) %>%
  select(symbol, rank_stage_sumz_wt, rank_stage_sumz)

scores_surv_os <- read_feather(snakemake@input[[8]]) %>%
  mutate(rank_surv_os_sumz_wt = dense_rank(sumz_wt_pval),
         rank_surv_os_sumz = dense_rank(sumz_pval)) %>%
  select(symbol, rank_surv_os_sumz_wt, rank_surv_os_sumz)

scores_surv_pfs <- read_feather(snakemake@input[[9]]) %>%
  mutate(rank_surv_pfs_sumz_wt = dense_rank(sumz_wt_pval),
         rank_surv_pfs_sumz = dense_rank(sumz_pval)) %>%
  select(symbol, rank_surv_pfs_sumz_wt, rank_surv_pfs_sumz)

scores_trmt <- read_feather(snakemake@input[[10]]) %>%
  mutate(rank_trmt_sumz_wt = dense_rank(sumz_wt_pval),
         rank_trmt_sumz = dense_rank(sumz_pval)) %>%
  select(symbol, rank_trmt_sumz_wt, rank_trmt_sumz)

ratio_na <- apply(means[, -1], 1, function(x) {
  sum(is.na(x)) / length(x)
})

df <- data.frame(
  symbol = means$symbol,
  mean = apply(means[, -1], 1, mean, na.rm=TRUE),
  median = apply(medians[, -1], 1, mean, na.rm=TRUE),
  var = apply(vars[, -1], 1, mean, na.rm=TRUE),
  cv = apply(cvs[, -1], 1, mean, na.rm=TRUE),
  ratio_nonzero = apply(ratio_nonzeros[, -1], 1, mean, na.rm=TRUE),
  ratio_missing = ratio_na
)

df <- df %>%
  left_join(scores_all) %>%
  left_join(scores_stage) %>%
  left_join(scores_surv_os) %>%
  left_join(scores_surv_pfs) %>%
  left_join(scores_trmt)

write_feather(df, snakemake@output[[1]])

