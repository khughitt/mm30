#
# compute_scaled_disease_stage_expression.R
#
# scales gene/gene set expression within each sample and then computes average scaled expression
# levels for each disease stage.
#
# a similar calculation is performed for several aggregate/composite disease stages as well: 
# "early", "late", and "not relapsed"
#
library(tidyverse)
library(arrow)

snek <- snakemake

expr <- read_feather(snek@input[[1]])
mdat <- read_feather(snek@input[[2]])

acc <- snek@wildcards$stage_dataset
feat_level <- snek@wildcards$feat_level
feat_names <- expr[, 1]

mdat <- mdat %>%
  filter(experiment == acc) %>%
  select(sample_id, disease_stage)

stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

early_stages <- c("Healthy", "MGUS", "SMM")
late_stages  <- c("MM", "RRMM")
not_relapsed <- c("Healthy", "MGUS", "SMM", "MM")

# extract samples from dataset convert expression values to ratios
mask <- colnames(expr) %in% mdat$sample_id

expr <- expr[, mask]
expr_scaled <- sweep(expr, 2, colSums(expr, na.rm = TRUE), "/") * 1E6

# compute average scaled expr for each stage
stage_expr <- list()

for (stage in stages) {
  sample_ids <- mdat %>%
    filter(disease_stage == stage) %>%
    pull(sample_id)
  
  if (length(sample_ids) == 0) {
    next
  }
  stage_expr[[stage]] <- apply(expr_scaled[, sample_ids], 1, median)
}

# add aggregate stages

# "early" vs. "late"?
early_samples <- mdat %>%
    filter(disease_stage %in% early_stages) %>%
    pull(sample_id)

late_samples <- mdat %>%
    filter(disease_stage %in% late_stages) %>%
    pull(sample_id)

if (length(early_samples) > 0) {
  stage_expr[["early"]] <- apply(expr_scaled[, early_samples], 1, median)
}

if (length(late_samples) > 0) {
  late_ranks <- apply(expr_scaled[, late_samples], 1, median)
  stage_expr[["late"]] <- apply(expr_scaled[, late_samples], 1, median)
}

if ("RRMM" %in% mdat$disease_stage) {
  not_relapsed_samples <- mdat %>%
    filter(disease_stage != "RRMM") %>%
    pull(sample_id)

  not_relapsed <- apply(expr_scaled[, not_relapsed_samples], 1, median)
  stage_expr[["pre_relapsed"]] <- not_relapsed
}

df <- data.frame(stage_expr)
df <- cbind(feat_names, df)

write_feather(df, snek@output[[1]])
