"""
MM25 Score generation pipeline

KH 2020.02.02

Combines output from fassoc feature-phenotype association pipeline to generate final
MM25 gene and pathway weights.
"""
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

