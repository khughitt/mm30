#
# create_disease_stage_ranking_tables.R
#
# ranks gene expression magnitude within each sample and computes average rankings across samples
# for each stage.
#
# in addition to comparing individual stages with each other, several aggregate comparisons are also
# performed:
#
# 1. "early" (healthy/mgus/smm) vs. "late" (mm/rrmm)
# 2. "up_to_mgus" (healthy/mgus) vs. "after_mgus" (smm/mm/rrmm)
# 3. relapsed vs. not
#
library(tidyverse)
library(arrow)

snek <- snakemake

acc <- snek@wildcards$stage_dataset
feat_level <- snek@wildcards$feat_level

expr <- read_feather(snek@input[[1]]) %>%
  column_to_rownames(feat_level)

mdat <- read_feather(snek@input[[2]])

mdat <- mdat %>%
  filter(experiment == acc) %>%
  select(sample_id, disease_stage)

stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

early_stages <- c("Healthy", "MGUS", "SMM")
late_stages  <- c("MM", "RRMM")
up_to_mgus_stages <- c("Healthy", "MGUS")
after_mgus_stages <- c("SMM", "MM", "RRMM")

# get expr data for current accession
mask <- colnames(expr) %in% mdat$sample_id
expr <- expr[, mask]

# exclude rows with all missing values
num_nas <- apply(expr, 1, function(x) {
  sum(is.na(x))
})
mask <- num_nas < ncol(expr)

expr <- expr[mask, ]

# convert expression measurements to ranking quantiles within each sample
expr_ranks <- do.call(cbind, lapply(expr, rank, na.last='keep'))

# max ranking may different by sample due to ties so a separate max is used for each column;
# "0" - genes/gene sets with the lowest expr
# "1" - genes/gene sets with the highest expr
max_ranks <- apply(expr_ranks, 2, max, na.rm=TRUE)
expr_rank_quantiles <- sweep(expr_ranks, 2, max_ranks, "/")

# compute average ranking quantiles for each stage
stage_ranks <- list()

for (stage in stages) {
  sample_ids <- mdat %>%
    filter(disease_stage == stage) %>%
    pull(sample_id)
  
  if (length(sample_ids) == 0) {
    next
  }
  stage_ranks[[stage]] <- apply(expr_rank_quantiles[, sample_ids], 1, median)
}

# compute ranking changes across stages
stage_transitions = list()

for (i in 1:(length(stage_ranks) - 1)) {
  ind1 <- i
  ind2 <- i + 1

  stage1 <- names(stage_ranks)[ind1]
  stage2 <- names(stage_ranks)[ind2]

  label <- paste0(c(stage1, stage2), collapse='_')

  stage_transitions[[label]] <- stage_ranks[[stage2]] - stage_ranks[[stage1]]
}

# "early" vs. "late"?
early_samples <- mdat %>%
    filter(disease_stage %in% early_stages) %>%
    pull(sample_id)

late_samples <- mdat %>%
    filter(disease_stage %in% late_stages) %>%
    pull(sample_id)

if (length(early_samples) > 0 && length(late_samples) > 0) {
  early_ranks <- apply(expr_rank_quantiles[, early_samples], 1, median)
  late_ranks <- apply(expr_rank_quantiles[, late_samples], 1, median)

  stage_transitions[["early_vs_late"]] <- late_ranks - early_ranks
}

# Healthy -> MGUS
up_to_mgus_samples <- mdat %>%
    filter(disease_stage %in% up_to_mgus_stages) %>%
    pull(sample_id)

# SMM -> RRMM
after_mgus_samples <- mdat %>%
    filter(disease_stage %in% after_mgus_stages) %>%
    pull(sample_id)

if (length(up_to_mgus_samples) > 0 && length(after_mgus_samples) > 0) {
  up_to_mgus_ranks <- apply(expr_rank_quantiles[, up_to_mgus_samples], 1, median)
  after_mgus_ranks <- apply(expr_rank_quantiles[, after_mgus_samples], 1, median)

  stage_transitions[["before_vs_after_smm"]] <- after_mgus_ranks - up_to_mgus_ranks
}

# relapsed vs. not?
if ("RRMM" %in% mdat$disease_stage) {
  not_relapsed_samples <- mdat %>%
    filter(disease_stage != "RRMM") %>%
    pull(sample_id)

  relapsed_samples <- mdat %>%
    filter(disease_stage == "RRMM") %>%
    pull(sample_id)

  not_relapsed <- apply(expr_rank_quantiles[, not_relapsed_samples], 1, median)
  relapsed <- apply(expr_rank_quantiles[, relapsed_samples], 1, median)

  stage_transitions[["before_vs_after_relapse"]] <- relapsed - not_relapsed
}

df <- data.frame(stage_transitions)
df <- cbind(rownames(expr), df)

write_feather(df, snek@output[[1]])
