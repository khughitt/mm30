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

feat_levels = ["gene", "pathway"]

wildcard_constraints:
    category="|".join(categories),

rule all:
    input:
        expand(os.path.join(out_dir, "scores", "metap", "all", "{feat_level}.feather"), 
               feat_level=feat_levels),
        expand(os.path.join(out_dir, "scores", "metap", "{category}", "{feat_level}.feather"),
               feat_level=feat_levels, category=categories),
        expand(os.path.join(out_dir, "scores", "metafor", "{category}", "{feat_level}.feather"),
               feat_level=feat_levels, category=categories),
        os.path.join(out_dir, "expr", "gene", "expr.feather"),
        os.path.join(out_dir, "summary", "gene.feather"),
        os.path.join(out_dir, "gene", "survival.feather"),
        os.path.join(out_dir, "scores", "gene_score_cor_mat.feather"),
        os.path.join(out_dir, "metadata", "samples.feather"),
        os.path.join(out_dir, "metadata", "datasets.feather"),
        os.path.join(out_dir, "metadata", "covariates.yml")

rule build_packages:
    input: 
        expand(os.path.join(out_dir, "packages", "{feat_level}", "datapackage.yml"),
               feat_level=feat_levels)

rule package_scores:
    input:
        os.path.join(out_dir, "scores", "metap", "all", "{feat_level}.feather"), 
    output:
        os.path.join(out_dir, "packages", "{feat_level}", "datapackage.yml"),
        os.path.join(out_dir, "packages", "{feat_level}", "data.csv"), 
    params:
        id=lambda w: f"{w.feat_level}_all_scores",
        title=lambda w: f"MM30 {w.feat_level.capitalize()} scores"
    script:
        "scripts/package_scores.py"

rule compute_mm30_ranking_correlations:
    input:
        os.path.join(out_dir, "scores", "metap", "all", "gene.feather"), 
        expand(os.path.join(out_dir, "scores", "metap", "{category}", "gene.feather"), category=categories)
    output:
        os.path.join(out_dir, "scores", "gene_score_cor_mat.feather")
    run:
        gene_scores = pd.read_feather(input[0])
        gene_scores = gene_scores.set_index("symbol")[["sumz_wt_pval"]]
        gene_scores.columns = ["all"]

        for infile in input[1:]:
            dat = pd.read_feather(infile)
            dat = dat.set_index("symbol")[["sumz_wt_pval"]]
            category = os.path.basename(os.path.dirname(infile))
            dat.columns = [category]

            gene_scores = gene_scores.join(dat, on="symbol")

        gene_scores.corr().reset_index().rename(columns={"index": "category"}).to_feather(output[0])

rule create_gene_survival_table:
  input:
    os.path.join(config["fassoc_dir"], "merged", "gene_association_effects.feather"),
    os.path.join(config["fassoc_dir"], "merged", "gene_association_errors.feather"),
    os.path.join(config["fassoc_dir"], "merged", "gene_association_pvals.feather"),
    os.path.join(config["fassoc_dir"], "metadata", "association_metadata.feather")
  output:
    os.path.join(out_dir, "gene", "survival.feather")
  script:
    "scripts/create_gene_survival_table.R"

rule create_gene_summary_table:
  input:
      os.path.join(out_dir, "expr", "gene", "stats", "mean.feather"),
      os.path.join(out_dir, "expr", "gene", "stats", "median.feather"),
      os.path.join(out_dir, "expr", "gene", "stats", "var.feather"),
      os.path.join(out_dir, "expr", "gene", "stats", "cv.feather"),
      os.path.join(out_dir, "expr", "gene", "stats", "ratio_nonzero.feather"),
      os.path.join(out_dir, "scores", "metap", "all", "gene.feather"),
      os.path.join(out_dir, "scores", "metap", "disease_stage", "gene.feather"),
      os.path.join(out_dir, "scores", "metap", "survival_os", "gene.feather"),
      os.path.join(out_dir, "scores", "metap", "survival_pfs", "gene.feather"),
      os.path.join(out_dir, "scores", "metap", "treatment_response", "gene.feather")
  output:
      os.path.join(out_dir, "summary", "gene.feather")
  script:
    "scripts/create_gene_summary_table.R"

rule compute_gene_stats:
    output: 
        os.path.join(out_dir, "expr", "gene", "stats", "mean.feather"),
        os.path.join(out_dir, "expr", "gene", "stats", "median.feather"),
        os.path.join(out_dir, "expr", "gene", "stats", "var.feather"),
        os.path.join(out_dir, "expr", "gene", "stats", "cv.feather"),
        os.path.join(out_dir, "expr", "gene", "stats", "ratio_nonzero.feather"),
    script:
        "scripts/compute_gene_stats.R"

rule create_combined_expr:
    input:
        os.path.join(out_dir, "scores", "metap", "all", "gene.feather")
    output: 
        os.path.join(out_dir, "expr", "gene", "expr.feather"),
        os.path.join(out_dir, "expr", "gene", "expr_scaled.feather")
    script:
        "scripts/create_combined_expr.R"

rule compute_metap_scores:
    input: 
        pvals=os.path.join(config["fassoc_dir"], "merged", "{feat_level}_association_pvals.feather"),
        mdata=os.path.join(config["fassoc_dir"], "metadata", "association_metadata.feather")
    output:
        os.path.join(out_dir, "scores", "metap", "all", "{feat_level}.feather")
    script:
        "scripts/compute_metap_scores.R"

rule compute_category_specific_metap_scores:
    input: 
        pvals=os.path.join(out_dir, "associations", "{category}", "{feat_level}", "pvals.feather"),
        mdata=os.path.join(config["fassoc_dir"], "metadata", "association_metadata.feather")
    output:
        os.path.join(out_dir, "scores", "metap", "{category}", "{feat_level}.feather")
    script:
        "scripts/compute_metap_scores.R"

rule compute_category_specific_metafor_scores:
    input: 
        effects=os.path.join(out_dir, "associations", "{category}", "{feat_level}", "effects.feather"),
        errors =os.path.join(out_dir, "associations", "{category}", "{feat_level}", "errors.feather"),
        mdata  =os.path.join(config["fassoc_dir"], "metadata", "association_metadata.feather")
    output:
        os.path.join(out_dir, "scores", "metafor", "{category}", "{feat_level}.feather")
    script:
        "scripts/compute_metafor_scores.R"

rule create_category_specific_associations:
    input:
        pvals  =os.path.join(config["fassoc_dir"], "merged", "{feat_level}_association_pvals.feather"),
        effects=os.path.join(config["fassoc_dir"], "merged", "{feat_level}_association_effects.feather"),
        errors =os.path.join(config["fassoc_dir"], "merged", "{feat_level}_association_errors.feather"),
        mdata  =os.path.join(config["fassoc_dir"], "metadata", "association_metadata.feather")
    output:
        pvals  =os.path.join(out_dir, "associations", "{category}", "{feat_level}", "pvals.feather"),
        effects=os.path.join(out_dir, "associations", "{category}", "{feat_level}", "effects.feather"),
        errors =os.path.join(out_dir, "associations", "{category}", "{feat_level}", "errors.feather")
    script:
        "scripts/create_category_specific_associations.R"

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
        os.path.join(out_dir, "metadata", "samples.feather")
    script:
        "scripts/create_combined_sample_metadata.py"

# vi:ft=snakemake
