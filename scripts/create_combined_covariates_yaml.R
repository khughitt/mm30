#
# create_combined_covariates_yaml.R
#
library(yaml)

snek <- snakemake

covariates <- list()

for (fname in list.files(snek@config$fassoc_dir)) {
  cfg <- read_yaml(file.path(snek@config$fassoc_dir, fname))
  covariates[[cfg$name]] <- cfg$phenotypes$associations
}

write_yaml(covariates, snek@output[[1]])
