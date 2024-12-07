"""
MM30 Score generation pipeline

Combines output from fassoc feature-phenotype association pipeline to generate final
MM30 gene and pathway weights.
"""
import os
import pandas as pd

configfile: "config/config-v7.2.yml"

# base input and output data dirs;
# output comes from the separate feature association ("fassoc") pipeline
out_dir = os.path.join(config["out_dir"], config["version"], "results")

# rankings to create for specific dataset/covariate categories
categories = ["disease_stage", "survival_os", "survival_pfs", "treatment_response"]

wildcard_constraints:
    category="|".join(categories),

rule all:
    input:
        expand(os.path.join(out_dir, "scores", "all", "mm30_{feat_level}_scores.feather"), 
               feat_level=["gene", "pathway"]),
        expand(os.path.join(out_dir, "scores", "categories", "mm30_{feat_level}_{category}_scores.feather"),
               feat_level=["gene", "pathway"], category=categories),
        os.path.join(out_dir, 'expr', 'mm30_combined_expr_data.feather'),
        os.path.join(out_dir, 'genes', 'summary.feather'),
        os.path.join(out_dir, 'genes', 'survival.feather'),
        os.path.join(out_dir, "scores", "gene_score_cor_mat.feather"),
        os.path.join(out_dir, "metadata", "samples.feather")

rule build_packages:
    input: 
        expand(os.path.join(out_dir, "packages", "scores", "mm30_{feat_level}_scores", "datapackage.yml"),
               feat_level=["gene", "pathway"])

rule package_scores:
    input:
        os.path.join(out_dir, "scores", "all", "mm30_{feat_level}_scores.feather"), 
    output:
        os.path.join(out_dir, "packages", "scores", "mm30_{feat_level}_scores", "datapackage.yml"),
        os.path.join(out_dir, "packages", "scores", "mm30_{feat_level}_scores", "data.csv"), 
    params:
        id=lambda w: f"mm30_{w.feat_level}_all_scores",
        title=lambda w: f"MM30 {w.feat_level.capitalize()} scores"
    script:
        "scripts/package_scores.py"

rule compute_mm30_ranking_correlations:
    input:
        os.path.join(out_dir, "scores/all/mm30_gene_scores.feather"), 
        expand(os.path.join(out_dir, "scores/categories/mm30_gene_{category}_scores.feather"), 
               category=categories)
    output:
        os.path.join(out_dir, "scores", "gene_score_cor_mat.feather")
    run:
        gene_scores = pd.read_feather(input[0])
        gene_scores = gene_scores.set_index('symbol')[['sumz_wt_pval']]
        gene_scores.columns = [os.path.basename(input[0])]

        for infile in input[1:]:
            dat = pd.read_feather(infile)
            dat = dat.set_index('symbol')[['sumz_wt_pval']]
            dat.columns = [os.path.basename(infile)]

            gene_scores = gene_scores.join(dat, on='symbol')

        gene_scores.corr().reset_index().rename(columns={"index": "file"}).to_feather(output[0])

rule create_gene_survival_table:
  input:
    os.path.join(config['fassoc_dir'], "merged", "gene_association_effects.feather"),
    os.path.join(config['fassoc_dir'], "merged", "gene_association_errors.feather"),
    os.path.join(config['fassoc_dir'], "merged", "gene_association_pvals.feather"),
    os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
  output:
    os.path.join(out_dir, 'genes', 'survival.feather')
  script:
    "scripts/create_gene_survival_table.R"

rule create_gene_summary_table:
  input:
      os.path.join(out_dir, 'genes', 'mean.feather'),
      os.path.join(out_dir, 'genes', 'median.feather'),
      os.path.join(out_dir, 'genes', 'var.feather'),
      os.path.join(out_dir, 'genes', 'cv.feather'),
      os.path.join(out_dir, 'genes', 'ratio_nonzero.feather'),
      os.path.join(out_dir, 'scores', 'all', 'mm30_gene_scores.feather'),
      os.path.join(out_dir, 'scores', 'categories', 'mm30_gene_disease_stage_scores.feather'),
      os.path.join(out_dir, 'scores', 'categories', 'mm30_gene_survival_os_scores.feather'),
      os.path.join(out_dir, 'scores', 'categories', 'mm30_gene_survival_pfs_scores.feather'),
      os.path.join(out_dir, 'scores', 'categories', 'mm30_gene_treatment_response_scores.feather')
  output:
      os.path.join(out_dir, 'genes', 'summary.feather')
  script:
    "scripts/create_gene_summary_table.R"

rule compute_gene_stats:
    output: 
        os.path.join(out_dir, 'genes', 'mean.feather'),
        os.path.join(out_dir, 'genes', 'median.feather'),
        os.path.join(out_dir, 'genes', 'var.feather'),
        os.path.join(out_dir, 'genes', 'cv.feather'),
        os.path.join(out_dir, 'genes', 'ratio_nonzero.feather'),
    script:
        "scripts/compute_gene_stats.R"

rule create_combined_expr:
    input:
        os.path.join(out_dir, "scores", "all", "mm30_gene_scores.feather")
    output: 
        os.path.join(out_dir, 'expr', 'mm30_combined_expr_data.feather'),
        os.path.join(out_dir, 'expr', 'mm30_combined_expr_data_scaled.feather')
    script:
        "scripts/create_combined_expr.R"

rule build_scores:
    input: 
        pvals=os.path.join(config['fassoc_dir'], "merged", "{feat_level}_association_pvals.feather"),
        mdata=os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
    output:
        os.path.join(out_dir, "scores", "all", "mm30_{feat_level}_scores.feather")
    script:
        "scripts/build_scores.R"

rule build_category_specific_scores:
    input: 
        pvals=os.path.join(out_dir, "subsets", "mm30_{feat_level}_{category}_pvals.feather"),
        effects=os.path.join(out_dir, "subsets", "mm30_{feat_level}_{category}_effects.feather"),
        errors=os.path.join(out_dir, "subsets", "mm30_{feat_level}_{category}_errors.feather"),
        mdata=os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
    output:
        os.path.join(out_dir, "scores", "categories", "mm30_{feat_level}_{category}_scores.feather")
    script:
        "scripts/build_subset_scores.R"

rule create_mm30_category_subsets:
    input:
        pvals=os.path.join(config['fassoc_dir'], "merged", "{feat_level}_association_pvals.feather"),
        effects=os.path.join(config['fassoc_dir'], "merged", "{feat_level}_association_effects.feather"),
        errors=os.path.join(config['fassoc_dir'], "merged", "{feat_level}_association_errors.feather"),
        mdata=os.path.join(config['fassoc_dir'], "metadata", "association_metadata.feather")
    output:
        pvals=os.path.join(out_dir, "subsets", "mm30_{feat_level}_{category}_pvals.feather"),
        effects=os.path.join(out_dir, "subsets", "mm30_{feat_level}_{category}_effects.feather"),
        errors=os.path.join(out_dir, "subsets", "mm30_{feat_level}_{category}_errors.feather")
    script:
        "scripts/create_category_subsets.R"

rule create_combined_covariates_yaml:
    output:
        os.path.join(out_dir, "metadata", "covariates.yml")
    script:
        "scripts/create_combined_covariates_yaml.R"

rule create_dataset_metadata:
    output:
        os.path.join(out_dir, "metadata", "datasets.feather")
    script:
        "scripts/create_dataset_metadata.R"

rule create_combined_sample_metadata:
    output:
        os.path.join(out_dir, "metadata", "metadata.feather")
    script:
        "scripts/create_combined_sample_metadata.py"

# vi:ft=snakemake
