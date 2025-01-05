#
# create_feature_metadata_table.R
# 

library(AnnotationDbi)
library(org.Hs.eg.db)
library(arrow)
library(tidyverse)
library(annotables)

snek <- snakemake

# load metap results
metap_df <- read_feather(snek@input[[1]])

# load cell cycle score data
ccs <- read_tsv("data/cell-cycle-score/cell-cycle-score-signature-lundberg2020.tsv", show_col_types=FALSE) %>%
  dplyr::select(symbol, cell_cycle_phase=phase)

# load DGIdb data
dgidb <- read_tsv("data/dgidb/categories.tsv", show_col_types=FALSE) %>%
  dplyr::select(symbol=name, category) %>%
  distinct() %>%
  group_by(symbol) %>%
  summarise(dgidb_categories=paste(category, collapse = ", "))

# retrieve basic gene metadata from annotables
gene_annot <- grch38 %>%
  dplyr::select(symbol, description, biotype) %>%
  group_by(symbol) %>%
  dplyr::slice(1) %>%
  ungroup()

# exclude genes not present in mm30
to_keep <- metap_df %>%
  pull(symbol)

gene_annot <- gene_annot %>%
  filter(symbol %in% to_keep)

# remove the source info from gene description
gene_annot$description <- str_split(gene_annot$description, " \\[", simplify=TRUE)[, 1]

# add chromosome bands
chr_bands <- AnnotationDbi::select(org.Hs.eg.db,
                                   columns=c("MAP"),
                                   keys=unique(grch38$symbol),
                                   keytype="SYMBOL")

colnames(chr_bands) <- c("symbol", "chr_subband")

gene_annot <- gene_annot %>%
  left_join(chr_bands, by="symbol") %>%
  filter(symbol != "")

# for the small number gene symbols mapped to multiple loci, set values to "NA"
dups <- gene_annot$symbol[duplicated(gene_annot$symbol)]
gene_annot <- gene_annot[!duplicated(gene_annot$symbol), ]
gene_annot[gene_annot$symbol %in% dups, "chr_subband"] <- NA

# add region-level chromosome annotations
gene_annot$chr_region <- sub("\\..+", "", gene_annot$chr_subband)

gene_annot <- gene_annot %>%
  left_join(dgidb, by="symbol") %>%
  left_join(ccs, by="symbol") %>%
  mutate(biotype=as_factor(biotype),
         chr_subband=as_factor(chr_subband),
         chr_region=as_factor(chr_region),
         cell_cycle_phase=fct_relevel(cell_cycle_phase, 'G1', 'G1/S', 'S', 'G2', 'G2/M', 'M'))

gene_annot %>%
  write_feather(snek@output[[1]])
