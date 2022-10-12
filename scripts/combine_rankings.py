"""
Combines phenotype-specific and "all" rankings into a single table.
"""
import pandas as pd

# all
all = pd.read_feather(snakemake.input[0])

# gene or gene_set?
feature = "symbol" if snakemake.wildcards[0] == "gene" else "gene_set"

all = all.rename(columns={"sumz_wt_pval": "all_pval"})
all = all[[feature, 'all_pval', 'num_present', 'num_missing']]

# disease stage
stage = pd.read_feather(snakemake.input[1])
stage = stage.rename(columns={"sumz_wt_pval": "disease_stage_pval"})
stage = stage[[feature, 'disease_stage_pval']]


# survival
surv = pd.read_feather(snakemake.input[2])

surv = surv.rename(columns={"sumz_wt_pval": "survival_pval"})
surv = surv[[feature, 'survival_pval']]

# treatment
treatment = pd.read_feather(snakemake.input[3])

treatment = treatment.rename(columns={"sumz_wt_pval": "treatment_pval"})
treatment = treatment[[feature, 'treatment_pval']]

# combine and write result
res = (all.merge(stage, on=feature, how='outer')
          .merge(surv, on=feature, how='outer')
          .merge(treatment, on=feature, how='outer'))

res = res[[feature, 'all_pval', 'disease_stage_pval', 'survival_pval', 'treatment_pval',
           'num_present', 'num_missing']]

res = res.to_feather(snakemake.output[0])
