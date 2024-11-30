#
# summarize gene statistics
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

means <- read_feather(snakemake@input[[1]])
medians <- read_feather(snakemake@input[[2]])
vars <- read_feather(snakemake@input[[3]])
cvs <- read_feather(snakemake@input[[4]])
ratio_nonzeros <- read_feather(snakemake@input[[5]])

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

write_feather(df, snakemake@output[[1]])

