#!/bin/env Rscript
#
# Creates a custom filtered gmt file containing the top N highest ranking MM29
# pathways/gene sets.
#
suppressMessages(library(GSEABase))
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)

# largest custom gmt file to generate (# gene sets)
max_gmt_size <- max(snakemake@config$gene_sets$sizes)

# load MM29 pathway rankings
gene_set_ranks <- read_feather(snakemake@input[[1]])
gene_set_ranks$gene_set <- as.character(gene_set_ranks$gene_set)

# grab top N gene sets (already in sorted order)
top_gene_sets <- gene_set_ranks %>%
  head(max_gmt_size) %>%
  pull(gene_set)

# load input gmts
gmt_paths <- Sys.glob(snakemake@config$gene_sets$gmts)

gene_sets <- lapply(gmt_paths, function(x) {
  res <- geneIds(getGmt(x))
  lapply(res, function(gset) { gset[gset != ''] })
})
names(gene_sets) <- tools::file_path_sans_ext(basename(gmt_paths))

# remove gene set :length suffixes, if present
names(gene_sets) <- sub(':\\d+$', '', names(gene_sets))

# remove any leading or trailing whitespace in gene set names;
# encountered at least one instance of this ("HEK 293 T-rex" in NCI-60 gene set
# collection)
for (gset in names(gene_sets)) {
  names(gene_sets[[gset]]) <- trimws(names(gene_sets[[gset]]))
}

# remove gene set weights, if present
# e.g. "ANXA1,1.0" -> "ANXA1"
for (gset in names(gene_sets)) {
  gene_sets[[gset]] <- lapply(gene_sets[[gset]], function(x) { sub(',\\d+\\.\\d+$', '', x) })
}

# create a mapping from <collection, gene_set> pairs to the merged "<collection>_<gene_set>"
# identifiers used in MM29
mapping <- NULL

for (collection in names(gene_sets)) {
  mapping <- rbind(mapping, cbind(collection, gene_set = names(gene_sets[[collection]])))
}
mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)

mapping$combined_id <- paste(mapping$collection, mapping$gene_set, sep = '_')

# collapse nested collection/gene set list into a single list indexed
# by combined <collection>_<gene_set> identifiers to make it easier to find
# entries using the MM29 combined identifiers
combined_gene_sets <- list()

for (collection in names(gene_sets)) {
  collection_list <- gene_sets[[collection]]
  names(collection_list) <- paste(collection, names(collection_list), sep = '_')

  combined_gene_sets <- c(combined_gene_sets, collection_list)
}

# order by MM29 gene set rank and limit to size of largest custom gmt to be created
combined_gene_sets <- combined_gene_sets[top_gene_sets]

max_genes <- max(unlist(lapply(combined_gene_sets, length)))

# dataframe to store gmt output; 2 added to account for collection / gene set
# id columns
res <- matrix("", nrow = length(combined_gene_sets), ncol = max_genes + 2)

# iterate over top gene sets and add to custom GMT file
for (i in 1:length(combined_gene_sets)) {
    combined_id <- names(combined_gene_sets)[i]

    # split combined collection / gene set identifier back into individual components
    ind <- match(combined_id, mapping$combined_id)

    collection <- mapping$collection[ind]
    gene_set <- mapping$gene_set[ind]

    # get genes in gene set
    genes <- combined_gene_sets[[combined_id]]

    res[i, 1:(2 + length(genes))] <- c(gene_set, collection, genes)
}

# iterate over entries and write to file
out_dir <- dirname(snakemake@output[[1]])

for (num_gene_sets in snakemake@config$gene_sets$sizes) {
  outfile <- file.path(out_dir, sprintf('mm29_top_%d_pathways.gmt', num_gene_sets))

  # get the first N entries 
  output <- apply(res[1:num_gene_sets, ], 1, function(x) { 
    x <- x[x != '']
    paste(x, collapse = '\t')
  })

  # save custom gmt
  write(output, file = outfile)
}
