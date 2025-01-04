"""
MM30 Score generation pipeline

Combines output from fassoc feature-phenotype association pipeline to generate final
MM30 gene and gene set weights.
"""
import os
import pandas as pd

configfile: "config/config-v7.3.yml"

# base input and output data dirs;
# output comes from the separate feature association ("fassoc") pipeline
out_dir = os.path.join(config["out_dir"], config["version"], "results")

# rankings to create for specific dataset/covariate categories
categories = ["disease_stage", "survival_os", "survival_pfs", "treatment_response"]

# expression / co-expression variants
expr_subsets = ["100", "500", "1k", "5k", 
                "scaled-100", "scaled-500", "scaled-1k", "scaled-5k"]
expr_versions = ["full", "scaled-full"] + expr_subsets

feat_levels = ["gene", "gene_set"]

mm_stage_datasets = config['mm_stage_datasets']

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
        expand(os.path.join(out_dir, "scores", "combined", "{category}", "{feat_level}.feather"),
               feat_level=feat_levels, category=categories),
        expand(os.path.join(out_dir, "coex", "{feat_level}", "coex-{expr_version}.feather"),
               feat_level=feat_levels, expr_version=expr_subsets),
        expand(os.path.join(out_dir, "expr", "{feat_level}", "expr-{expr_version}.feather"),
               feat_level=feat_levels, expr_version=expr_versions),
        expand(os.path.join(out_dir, "disease_stage", "scaled", "{feat_level}", "combined.feather"), 
               feat_level=feat_levels),
        expand(os.path.join(out_dir, "disease_stage", "transitions", "{feat_level}", "combined.feather"),
               feat_level=feat_levels),
        os.path.join(out_dir, "gene", "survival.feather"),
        os.path.join(out_dir, "scores", "gene_score_cor_mat.feather"),
        os.path.join(out_dir, "metadata", "covariates.yml"),
        os.path.join(out_dir, "metadata", "covariates.feather"),
        os.path.join(out_dir, "metadata", "datasets.feather"),
        os.path.join(out_dir, "metadata", "genes.feather"),
        os.path.join(out_dir, "metadata", "gene_sets.feather"),
        os.path.join(out_dir, "metadata", "samples.feather")

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

rule create_combined_disease_stage_transition_table:
    input:
      expand(os.path.join(out_dir, "disease_stage", "transitions", "{{feat_level}}", "indiv", "{stage_dataset}.feather"),
             stage_dataset=mm_stage_datasets)
    output:
        os.path.join(out_dir, "disease_stage", "transitions", "{feat_level}", "combined.feather"),
        os.path.join(out_dir, "disease_stage", "transitions", "{feat_level}", "counts.feather")
    script:
        "scripts/create_combined_disease_stage_transition_table.R"

rule create_disease_stage_transition_tables:
    input:
        os.path.join(out_dir, "expr", "{feat_level}", "expr-full.feather"),
        os.path.join(out_dir, "metadata", "samples.feather")
    output:
        os.path.join(out_dir, "disease_stage", "transitions", "{feat_level}", "indiv", "{stage_dataset}.feather")
    script:
        "scripts/create_disease_stage_transition_tables.R"

rule create_combined_disease_stage_table:
    input:
      expand(os.path.join(out_dir, "disease_stage", "scaled", "{{feat_level}}", "indiv", "{stage_dataset}.feather"), 
             stage_dataset=mm_stage_datasets)
    output:
      os.path.join(out_dir, "disease_stage", "scaled", "{feat_level}", "combined.feather")
    script:
        "scripts/create_combined_disease_stage_table.R"

rule compute_scaled_disease_stage_expression:
    input:
        os.path.join(out_dir, "expr", "{feat_level}", "expr-full.feather"),
        os.path.join(out_dir, "metadata", "samples.feather")
    output:
        os.path.join(out_dir, "disease_stage", "scaled", "{feat_level}", "indiv", "{stage_dataset}.feather")
    script:
        "scripts/compute_scaled_disease_stage_expression.R"

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

rule create_combined_score_tables:
  input:
      os.path.join(out_dir, "scores", "metap", "{category}", "{feat_level}.feather"),
      os.path.join(out_dir, "scores", "metafor", "{category}", "{feat_level}.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "stats", "mean.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "stats", "median.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "stats", "var.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "stats", "cv.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "stats", "ratio_nonzero.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "stats", "ratio_missing.feather"),
  output:
      os.path.join(out_dir, "scores", "combined", "{category}", "{feat_level}.feather")
  script:
      "scripts/create_combined_score_tables.R"

rule create_coex_matrices:
    input:
      os.path.join(out_dir, "scores", "metap", "all", "{feat_level}.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "expr-5k.feather"),
      os.path.join(out_dir, "expr", "{feat_level}", "expr-scaled-5k.feather")
    output:
      os.path.join(out_dir, "coex", "{feat_level}", "coex-100.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-500.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-1k.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-5k.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-scaled-100.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-scaled-500.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-scaled-1k.feather"),
      os.path.join(out_dir, "coex", "{feat_level}", "coex-scaled-5k.feather")
    script:
      "scripts/create_coex_matrices.R"

rule create_combined_expr_matrices:
    input:
        os.path.join(out_dir, "scores", "metap", "all", "{feat_level}.feather")
    output: 
        os.path.join(out_dir, "expr", "{feat_level}", "expr-full.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-100.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-500.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-1k.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-5k.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-scaled-full.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-scaled-100.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-scaled-500.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-scaled-1k.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "expr-scaled-5k.feather")
    script:
        "scripts/create_combined_expr_matrices.R"

rule compute_expr_stats:
    output: 
        os.path.join(out_dir, "expr", "{feat_level}", "stats", "mean.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "stats", "median.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "stats", "var.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "stats", "cv.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "stats", "ratio_nonzero.feather"),
        os.path.join(out_dir, "expr", "{feat_level}", "stats", "ratio_missing.feather"),
    script:
        "scripts/compute_expr_stats.R"

rule create_gene_metadata_table:
    input:
        os.path.join(out_dir, "scores", "metap", "all", "gene.feather")
    output:
        os.path.join(out_dir, "metadata", "genes.feather")
    script:
       "scripts/create_gene_metadata_table.R"

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

rule create_covariates_metadata_table:
    input:
        os.path.join(out_dir, "metadata", "datasets.feather")
    output:
        os.path.join(out_dir, "metadata", "covariates.feather")
    script:
        "scripts/create_covariates_metadata_table.R"

rule create_combined_covariates_yaml:
    output:
        os.path.join(out_dir, "metadata", "covariates.yml")
    script:
        "scripts/create_combined_covariates_yaml.R"

rule create_gene_set_metadata_table:
    output:
      os.path.join(out_dir, "metadata", "gene_sets.feather")
    script:
      "scripts/create_gene_set_metadata_table.py"

rule create_dataset_metadata_table:
    input:
        "metadata/mm30_experiment_metadata.tsv"
    output:
        os.path.join(out_dir, "metadata", "datasets.feather")
    script:
        "scripts/create_dataset_metadata.R"

rule create_combined_sample_metadata_table:
    output:
        os.path.join(out_dir, "metadata", "samples.feather")
    script:
        "scripts/create_combined_sample_metadata.py"

# vi:ft=snakemake
