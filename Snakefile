"""
MM25 Score generation pipeline

Combines output from fassoc feature-phenotype association pipeline to generate final
MM25 gene and pathway weights.
"""
import glob
import os
import re
import pandas as pd

# base input and output data dirs;
# output comes from the separate feature association ("fassoc") pipeline
out_dir = os.path.join(config["out_dir"], config["version"], "scores")

# load phenotype metadata
phenotypes_infile = os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
phenotypes = pd.read_feather(phenotypes_infile)

# phenotype categories
categories = phenotypes.category.unique()

wildcard_constraints:
    category="|".join(categories)

rule all:
    input:
        expand(os.path.join(out_dir, "results", "all", "mm25_{feat_level}_scores.feather"), 
               feat_level=["gene", "pathway"]),
        expand(os.path.join(out_dir, "results", "categories", "mm25_{feat_level}_{category}_scores.feather"), 
               feat_level=["gene", "pathway"], category=categories),
        expand(os.path.join(out_dir, "results", "clusters", "mm25_{feat_level}_{cluster_num}_scores.feather"),
               feat_level=["gene", "pathway"], cluster_num=range(config['clustering']['num_clusters'])),
        expand(os.path.join(out_dir, "results", "categories", "mm25_{feat_level}_survival_stats.feather"),
                feat_level=["gene", "pathway"]),
        os.path.join(out_dir, 'expr', 'mm25_combined_expr_data.feather'),
        os.path.join(out_dir, "summary", "gene_score_cor_mat.feather"),
        os.path.join(out_dir, "metadata.feather")

rule compute_mm25_ranking_correlations:
    input:
        os.path.join(out_dir, "results", "all", "mm25_gene_scores.feather"), 
        expand(os.path.join(out_dir, "results", "categories", "mm25_gene_{category}_scores.feather"), 
               category=categories),
        expand(os.path.join(out_dir, "results", "clusters", "mm25_gene_{cluster_num}_scores.feather"),
               cluster_num=range(config['clustering']['num_clusters']))
    output:
        os.path.join(out_dir, "summary", "gene_score_cor_mat.feather")
    run:
        # create a matrix of alternate gene scores created from different subsets
        # of the MM25 datasets
        gene_scores = pd.read_feather(input[0])
        gene_scores = gene_scores.set_index('symbol')[['sumz_wt_pval']]
        gene_scores.columns = [os.path.basename(input[0])]

        for infile in input[1:]:
            dat = pd.read_feather(infile)
            dat = dat.set_index('symbol')[['sumz_wt_pval']]
            dat.columns = [os.path.basename(infile)]

            gene_scores = gene_scores.join(dat, on='symbol')

        gene_scores.corr().reset_index().rename(columns={"index": "file"}).to_feather(output[0])

rule create_combined_expr:
    output: os.path.join(out_dir, 'expr', 'mm25_combined_expr_data.feather')
    script:
        "scripts/combine_expr_data.R"

rule mm25_surv_stats:
    input: 
        stats=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather"),
        coefs=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_coefs.feather"),
        mdata=os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
    output:
        os.path.join(out_dir, "results", "categories", "mm25_{feat_level}_survival_stats.feather")
    script:
        "scripts/build_survival_stats.R"

rule mm25_all:
    input: 
        pvals=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather"),
        coefs=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_coefs.feather")
    output:
        os.path.join(out_dir, "results", "all", "mm25_{feat_level}_scores.feather")
    script:
        "scripts/build_scores.R"

rule mm25_clustered:
    input: 
        pvals=os.path.join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_pvals.feather"),
        stats=os.path.join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_stats.feather"),
        coefs=os.path.join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_coefs.feather")
    output:
        os.path.join(out_dir, "results", "clusters", "mm25_{feat_level}_{cluster_num}_scores.feather")
    script:
        "scripts/build_scores.R"

rule mm25_categories:
    input: 
        pvals=os.path.join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_pvals.feather"),
        stats=os.path.join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_stats.feather"),
        coefs=os.path.join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_coefs.feather")
    output:
        os.path.join(out_dir, "results", "categories", "mm25_{feat_level}_{category}_scores.feather")
    script:
        "scripts/build_scores.R"

rule create_mm25_cluster_subsets:
    input:
        pvals=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather"),
        coefs=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_coefs.feather"),
        clusters=os.path.join(out_dir, "clusters", "mm25_{feat_level}_covariate_clusters.feather")
    output:
        pvals=os.path.join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_pvals.feather"),
        stats=os.path.join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_stats.feather"),
        coefs=os.path.join(out_dir, "subsets", "clusters", "mm25_{feat_level}_{cluster_num}_coefs.feather")
    script:
        "scripts/create_cluster_subsets.R"

rule create_mm25_category_subsets:
    input:
        pvals=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_stats.feather"),
        coefs=os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_coefs.feather"),
        mdata=os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
    output:
        pvals=os.path.join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_pvals.feather"),
        stats=os.path.join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_stats.feather"),
        coefs=os.path.join(out_dir, "subsets", "categories", "mm25_{feat_level}_{category}_coefs.feather")
    script:
        "scripts/create_category_subsets.R"

rule cluster_covariates:
    input: 
        os.path.join(config['fassoc_dir'], "merged", "mm25_{feat_level}_association_pvals.feather")
    output:
        os.path.join(out_dir, "clusters", "mm25_{feat_level}_covariate_clusters.feather")
    script:
        "scripts/cluster_covariates.py"

rule create_combined_sample_metadata:
    output:
        os.path.join(out_dir, "metadata.feather")
    run:
        geo_mdata = glob.glob(config['metadata']['geo'])
        mmrf_mdata = config['metadata']['mmrf'] 

        # combine geo/mmrf metadata into a single dataframe with columns for sample id,
        # experiment, platform, and platform_type
        mdat = pd.read_feather(mmrf_mdata)[["public_id", "disease_stage"]]

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
            geo_mdat = pd.read_feather(infile)

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
            try:
                geo_mdat = geo_mdat[mdat.columns]
            except:
                print(infile)
                breakpoint()

            # append to combined dataframe
            mdat = pd.concat([mdat, geo_mdat])

        # write combined metadata to disk
        mdat.reset_index(drop=True).to_feather(output[0])
