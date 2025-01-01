#
# create_dataset_metadata.R
#
library(arrow)
library(tidyverse)

snek <- snakemake

dataset_mdata <- read_tsv(snek@input[[1]], show_col_types=FALSE)

geo_mdata <- read_feather(file.path(snek@config$geo_dir, "metadata.feather")) %>%
  select(accession=geo_id, geo_name=name, abstract, overall_design, 
         geo_submission_date=submission_date, geo_last_update_date=last_update_date, urls,
         pubmed_ids, supplementary_files)

dataset_mdata <- dataset_mdata %>%
  left_join(geo_mdata, by='accession') %>%
  select(accession, name=dataset, everything())

dataset_mdata %>%
  write_feather(snek@output[[1]])
