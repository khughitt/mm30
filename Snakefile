"""
MM25 Score generation pipeline

KH 2020.02.02

Combines output from fassoc feature-phenotype association pipeline to generate final
MM25 gene and pathway weights.
"""
import glob
import re
from os.path import join
import pandas as pd

# base input and output data dirs;
# output comes from the separate feature association ("fassoc") pipeline
out_dir = join(config["out_dir"], config["version"])

# load phenotype metadata
phenotypes_infile = join(config['fassoc_dir'], "metadata", "association_metadata.feather")
phenotypes = pd.read_feather(phenotypes_infile)

# phenotype categories
categories = phenotypes.category.unique()

# filepath for one of the custom gmt outputs
custom_gmt = join(out_dir, "gmt", "mm25_top_{}_pathways.gmt".format(config["gene_sets"]["sizes"][0]))

wildcard_constraints:
    category="|".join(categories)

rule all:
    input:
        expand(join(out_dir, "results", "all", "mm25_{feat_level}_scores.feather"), 
               feat_level=["gene", "pathway"]),
        expand(join(out_dir, "results", "categories", "mm25_{feat_level}_{category}_scores.feather"), 
               feat_level=["gene", "pathway"], category=categories),
        expand(join(out_dir, "results", "clusters", "mm25_{feat_level}_{cluster_num}_scores.feather"),
               feat_level=["gene", "pathway"], cluster_num=range(config['clustering']['num_clusters'])),
        join(out_dir, "metadata.tsv"),
        custom_gmt

rule create_custom_gene_sets:
    input:
        join(out_dir, "results", "all", "mm25_pathway_scores.feather") 
    output:
        custom_gmt
    script:
        "src/create_custom_gmt.R"

rule mm25_all:
    input: 
        pvals=join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather")
    output:
        join(out_dir, "results", "all", "mm25_{feat_level}_scores.feather")
    script:
        "src/build_scores.R"

rule mm25_clustered:
    input: 
        pvals=join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_pvals.feather"),
        stats=join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_stats.feather")
    output:
        join(out_dir, "results", "clusters", "mm25_{feat_level}_{cluster_num}_scores.feather")
    script:
        "src/build_scores.R"

rule mm25_categories:
    input: 
        pvals=join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_pvals.feather"),
        stats=join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_stats.feather")
    output:
        join(out_dir, "results", "categories", "mm25_{feat_level}_{category}_scores.feather")
    script:
        "src/build_scores.R"

rule create_mm25_cluster_subsets:
    input:
        pvals=join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather"),
        clusters=join(out_dir, "clusters", "mm25_{feat_level}_covariate_clusters.feather")
    output:
        pvals=join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_pvals.feather"),
        stats=join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_stats.feather")
    script:
        "src/create_cluster_subsets.R"

rule create_mm25_category_subsets:
    input:
        pvals=join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather"),
        mdata=join(config['fassoc_dir'], "metadata", "association_metadata.feather")
    output:
        pvals=join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_pvals.feather"),
        stats=join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_stats.feather")
    script:
        "src/create_category_subsets.R"

rule cluster_covariates:
    input: 
        join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather")
    output:
        join(out_dir, "clusters", "mm25_{feat_level}_covariate_clusters.feather")
    script:
        "src/cluster_covariates.py"

rule create_combined_sample_metadata:
    output: join(out_dir, "metadata.tsv")
    run:
        geo_mdata = glob.glob('/data/human/geo/3.0/*/processed/*_sample_metadata.tsv')
        mmrf_mdata = "/data/human/mmrf-commpass/IA15/clean/clinical_data_tables/MMRF_CoMMpass_IA15_combined_metadata.tsv"

        outfile = "/data/nih/mm25/%s/metadata.tsv" % (config["version"])

        # combine geo/mmrf metadata into a single dataframe with columns for sample id,
        # experiment, platform, and platform_type
        mdat = pd.read_csv(mmrf_mdata, sep='\t')[["public_id", "disease_stage", "cell_type"]]

        mdat.rename(columns={'public_id': 'sample_id'}, inplace=True)

        # mmrf
        mdat['experiment'] = ['MMRF'] * mdat.shape[0]
        mdat['platform_id'] = ['GPL16791'] * mdat.shape[0] # HiSeq 2500
        mdat['platform_type'] = ['RNA-Seq'] * mdat.shape[0]

        # microarray platforms included in MM25 (v3.0)
        mm25_microarray_platforms = ["GPL96", "GPL97", "GPL570", "GPL10558", "GPL5175", "GPL6244", "GPL25401", "GPL4819"]

        # geo
        for infile in geo_mdata:
            # determine GEO identifier
            geo_id = re.findall('GSE[0-9]+', infile)[0]

            # load geo metadata
            geo_mdat = pd.read_csv(infile, sep = '\t')

            # rename sample id column
            geo_mdat.rename(columns={geo_mdat.columns[0]: 'sample_id'}, inplace=True)

            # add experiment column
            geo_mdat['experiment'] = [geo_id] * geo_mdat.shape[0]

            # determine platform type
            if geo_mdat.platform_id.iloc[0] in mm25_microarray_platforms:
                geo_mdat['platform_type'] = ['Microarray'] * geo_mdat.shape[0]
            else:
                geo_mdat['platform_type'] = ['RNA-Seq'] * geo_mdat.shape[0]

            # match column order
            geo_mdat = geo_mdat[mdat.columns]

            # append to combined dataframe
            mdat = pd.concat([mdat, geo_mdat])

        # write combined metadata to disk
        mdat.set_index('sample_id').to_csv(outfile, sep = '\t')

