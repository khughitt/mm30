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
# 2. "earlier" (healthy/mgus) vs. "less_late" (smm/mm/rrmm)
# 3. relapsed vs. not
#
library(tidyverse)
library(arrow)

snek <- snakemake

expr <- read_feather(snek@input[[1]])
mdat <- read_feather(snek@input[[2]])

acc <- snek@wildcards$stage_dataset
feat_level <- snek@wildcards$feat_level

mdat <- mdat %>%
  filter(experiment == acc) %>%
  select(sample_id, disease_stage)

stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

early_stages <- c("Healthy", "MGUS", "SMM")
earlier_stages <- c("Healthy", "MGUS")
late_stages  <- c("MM", "RRMM")
less_late_stages  <- c("SMM", "MM", "RRMM")
not_relapsed <- c("Healthy", "MGUS", "SMM", "MM")

# extract samples from dataset convert expression values to rankings within the same
mask <- colnames(expr) %in% mdat$sample_id

expr_ranks <- do.call(cbind, lapply(-expr[, mask], rank, na.last='keep'))

# compute average ranks for each stage
stage_ranks <- list()

for (stage in stages) {
  sample_ids <- mdat %>%
    filter(disease_stage == stage) %>%
    pull(sample_id)
  
  if (length(sample_ids) == 0) {
    next
  }
  stage_ranks[[stage]] <- apply(expr_ranks[, sample_ids], 1, median)
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
  early_ranks <- apply(expr_ranks[, early_samples], 1, median)
  late_ranks <- apply(expr_ranks[, late_samples], 1, median)

  stage_transitions[["early_vs_late"]] <- late_ranks - early_ranks
}

# "earlier" vs. "less late"?
earlier_samples <- mdat %>%
    filter(disease_stage %in% earlier_stages) %>%
    pull(sample_id)

less_late_samples <- mdat %>%
    filter(disease_stage %in% less_late_stages) %>%
    pull(sample_id)

if (length(earlier_samples) > 0 && length(less_late_samples) > 0) {
  earlier_ranks <- apply(expr_ranks[, earlier_samples], 1, median)
  less_late_ranks <- apply(expr_ranks[, less_late_samples], 1, median)

  stage_transitions[["earlier_vs_less_late"]] <- less_late_ranks - earlier_ranks
}

# relapsed vs. not?
if ("RRMM" %in% mdat$disease_stage) {
  not_relapsed_samples <- mdat %>%
    filter(disease_stage != "RRMM") %>%
    pull(sample_id)

  relapsed_samples <- mdat %>%
    filter(disease_stage == "RRMM") %>%
    pull(sample_id)

  not_relapsed <- apply(expr_ranks[, not_relapsed_samples], 1, median)
  relapsed <- apply(expr_ranks[, relapsed_samples], 1, median)

  stage_transitions[["before_vs_after_relapse"]] <- relapsed - not_relapsed
}

df <- data.frame(stage_transitions)
df <- cbind(expr[, 1], df)

write_feather(df, snek@output[[1]])
