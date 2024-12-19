#
# create_gene_set_metadata_table.py
#
import json
import pandas as pd

snek = snakemake

with open(snek.config['msigdb_gmt'], 'rt') as fp:
    lines = fp.readlines()
    mm30_gene_sets = [x.split('\t')[0] for x in lines]

with open(snek.config['msigdb_json'], 'rt') as fp:
    msigdb = json.load(fp)

rows = []

for gene_set,info in msigdb.items():
    if gene_set in mm30_gene_sets:
        rows.append({
            'gene_set': gene_set,
            'collection': info['collection'],
            'pmid': info['pmid'],
            'url': info['msigdbURL'],
            'genes': info['geneSymbols']
        })

mdata = pd.DataFrame(rows)
mdata = mdata.astype({'collection': 'category'})

mdata.to_feather(snek.output[0])
