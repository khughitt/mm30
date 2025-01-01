#
# create_average_scaled_disease_stage_expression_table.R
#
library(tidyverse)
library(arrow)

snek <- snakemake

dat_lst <- lapply(snek@input, read_feather)
names(dat_lst) <- stringr::str_extract(snek@input, "GSE[0-9]+")

for (acc in names(dat_lst)) {
  dat_lst[[acc]][["dataset"]] <- acc

  dat_lst[[acc]] <- dat_lst[[acc]] %>% 
    pivot_longer(-c(1, dataset), names_to="stage", values_to="expr")

  dat_lst[[acc]] <- dat_lst[[acc]][complete.cases(dat_lst[[acc]]), ]
}

# combine into a single long df
combined_dat <- do.call(rbind, dat_lst)

combined_dat$dataset <- factor(combined_dat$dataset)
combined_dat$stage <- factor(combined_dat$stage,
                             levels=c('Healthy', 'MGUS', 'SMM', 'MM', 'RRMM', 'early', 'late', 'pre_relapsed'))

combined_dat %>%
  write_feather(snek@output[[1]])
