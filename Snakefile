"""
MM25 Score generation pipeline

KH 2020.02.02

Combines output from fassoc feature-phenotype association pipeline to generate final
MM25 gene and pathway weights.
"""
import os
import pandas as pd

# base input and output data dirs;
# output comes from the separate feature association ("fassoc") pipeline
input_dir = os.path.join(config["fassoc_dir"], "mm25", config["version"])
out_dir = os.path.join(config["out_dir"], config["version"])

# load phenotype metadata
phenotypes_infile = os.path.join(input_dir, "metadata", "association_metadata.feather")
phenotypes = pd.read_feather(phenotypes_infile)

# phenotype categories
categories = phenotypes.category.unique()

# filepath for one of the custom gmt outputs
custom_gmt = os.path.join(out_dir, "gmt", "mm25_top_{}_pathways.gmt".format(config["gene_sets"]["sizes"][0]))

wildcard_constraints:
    category="|".join(categories)

rule all:
    input:
        expand(os.path.join(out_dir, "all", "mm25_{feat_level}_scores.feather"), 
               feat_level=["gene", "pathway"]),
        expand(os.path.join(out_dir, "targeted", "mm25_{feat_level}_{category}_scores.feather"), 
               feat_level=["gene", "pathway"], category=categories),
        custom_gmt,
        expand(os.path.join(out_dir, "clusters", "mm25_{feat_level}_covariate_clusters.feather"),
               feat_level=["gene", "pathway"])

rule create_custom_gene_sets:
    input:
        os.path.join(out_dir, "all", "mm25_pathway_scores.feather") 
    output:
        custom_gmt
    script:
        "src/create_custom_gmt.R"

rule mm25_targeted:
    input: 
        pvals=os.path.join(input_dir, "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=os.path.join(input_dir, "merged", "mm25_{feat_level}_association_stats.feather")
    output:
        os.path.join(out_dir, "targeted", "mm25_{feat_level}_{category}_scores.feather")
    params:
        category="{category}"
    script:
        "src/build_scores.R"

rule mm25_all:
    input: 
        pvals=os.path.join(input_dir, "merged", "mm25_{feat_level}_association_pvals.feather"),
        stats=os.path.join(input_dir, "merged", "mm25_{feat_level}_association_stats.feather")
    output:
        os.path.join(out_dir, "all", "mm25_{feat_level}_scores.feather")
    script:
        "src/build_scores.R"

rule cluster_covariates:
    input: 
        os.path.join(input_dir, "merged", "mm25_{feat_level}_association_pvals.feather")
    output:
        os.path.join(out_dir, "clusters", "mm25_{feat_level}_covariate_clusters.feather")
    script:
        "src/cluster_covariates.py"

