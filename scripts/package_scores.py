#!/bin/env python
"""
Package MM30 scores
"""
import os
import yaml
import numpy as np
import pandas as pd

#import seaborn as sns
#import umap
#from matplotlib import pyplot as plt
from nodes import BioDataset

snek = snakemake

df = pd.read_feather(snek.input[0])
#df = df.set_index(df.columns[0])

# load package metadata
#  with open("metadata/metadata.yml") as fp:
#      metadata = yaml.load(fp, Loader=yaml.FullLoader)

# package metadata
metadata = {
    "id": snek.params["id"],
    "title": snek.params["title"],
    "version": snek.config["version"],
    "profile": "bio-dataset",
    "urls": ["https://github.com/khughitt/mm30"],
    "row_type": "symbol",
    "datatype": "gene-weights",
    "diseases": ["D009101"],
    "species": [9606],
    "reference": "grch38",
    "contributors": [{
        "title": "V. Keith Hughitt",
        "email": "keith.hughitt@nih.gov",
        "role": "author"
    }]
}

# include gene info from annotables as row metadata
# todo: for pathways, include MSigDB gene set metadata;
if snek.wildcards['feat_level'] == 'gene':
    row_mdata = pd.read_csv("metadata/annotables-grch38.tsv.gz", sep='\t')
    row_mdata = row_mdata.set_index('symbol')

    # create placeholder rows for genes not found in the annotables annotations
    missing_ids = df.symbol[~df.symbol.isin(row_mdata.index)].values

    placeholder_rows = []

    for missing_id in missing_ids:
        placeholder_rows.append([missing_id] + [None] * 8)

    placeholder_df = pd.DataFrame(placeholder_rows, 
                                  columns=['symbol', 'ensgene', 'entrez', 'chr', 
                                           'start', 'end', 'strand', 'biotype', 'description'])
    placeholder_df = placeholder_df.set_index('symbol')
    row_mdata = pd.concat([row_mdata, placeholder_df], axis=0)

    # limit to genes in scores and reorder
    row_mdata = row_mdata.loc[df.symbol]
else:
    row_mdata = None

# create node instance
node = BioDataset(df, row_metadata=row_mdata, **metadata)

# write out package & metadata
node.to_pkg(os.path.dirname(snek.output[0]))
