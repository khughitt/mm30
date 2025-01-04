#
# create_combined_disease_stage_transition_table.R
#
library(tidyverse)
library(arrow)

snek <- snakemake

dat_lst <- lapply(snek@input, read_feather)

# get a list of all transitions measured
transitions <- unlist(lapply(dat_lst, colnames))
transitions <- transitions[!transitions %in% c('symbol', 'gene_set')]

# create a dataframe with one column corresponding to the average values for each transition
combined_lst <- list()

for (transition in unique(transitions)) {
  combined_lst[[transition]] <- list()

  for (df in dat_lst) {
    if (transition %in% colnames(df)) {
      combined_lst[[transition]] <- c(combined_lst[[transition]], list(pull(df, transition)))
    }
  }
}

res_lst <- list()

for (transition in names(combined_lst)) {
  df <- do.call(cbind, combined_lst[[transition]])
  res_lst[[transition]] <- apply(df, 1, median, na.rm=TRUE)
}

res_df <- data.frame(res_lst)

feat_names <- dat_lst[[1]][,1]
res_df <- cbind(feat_names, res_df)

res_df %>%
  write_feather(snek@output[[1]])

# store a table with the counts of each transition
transitions %>%
  table() %>%
  enframe(name="transition", value="num") %>%
  write_feather(snek@output[[2]])
