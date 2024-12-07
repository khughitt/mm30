#
# create_dataset_metadata.R
#
library(arrow)
library(tidyverse)

snek <- snakemake

geo_mdata <- read_feather(file.path(snek@config$geo_dir, "metadata.feather"))

dataset_mdata <- rbind(geo_mdata, c("MMRF", "MMRF CoMMpass Study IA22", NA, NA, NA, NA, NA,
                                    "GPL11154", NA, "https://research.themmrf.org/", "", ""))

dataset_mdata <- dataset_mdata %>%
  rename(dataset=geo_id)

dataset_mdata %>%
  write_feather(snek@output[[1]])
