#
# create_combined_disease_stage_transition_tables.R
#
library(tidyverse)
library(arrow)

snek <- snakemake

feat_level <- snek@wildcards$feat_level
feat_key <- ifelse(feat_level == "gene", "symbol", "gene_set")

dat_lst <- lapply(snek@input, read_feather)

accessions <- sub(".feather", "", basename(unlist(snek@input)))
names(dat_lst) <- accessions

# iterate over transitions and create a dataframe with all measurements associated with each
res <- list()

for (acc in names(dat_lst)) {
  df <- dat_lst[[acc]]

  for (transition in colnames(df)[-1]) {

    df_subset <- df[, c(feat_key, transition)]

    colnames(df_subset)[2] <- acc

    if (transition %in% names(res)) {
      res[[transition]] <- res[[transition]] %>%
        full_join(df_subset, by=feat_key)

    } else {
      res[[transition]] <- df_subset 
    }
  }
}


# add summary statistics and store resulting dataframes
out_dir <- dirname(snek@output[[1]])

num_non_missing <- function(x) {
  return(sum(!is.na(x)))
}

for (transition in names(res)) {
  res_df <- res[[transition]]

  num_datasets <- apply(res_df[, -1], 1, num_non_missing)
  mean_change <- apply(res_df[, -1], 1, mean, na.rm=TRUE)
  median_change <- apply(res_df[, -1], 1, median, na.rm=TRUE)

  res_df$num_datasets <- num_datasets
  res_df$mean_change <- mean_change
  res_df$median_change <- median_change

  outfile <- file.path(out_dir, sprintf("%s.feather", transition))

  res_df %>%
    write_feather(outfile)
}
