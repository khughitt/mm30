#!/bin/env python
#
# Clusters MM25 covariates
# V. Keith Hughitt
#
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

np.random.seed(1)

# load p-values
dat = pd.read_feather(snakemake.input[0])

feat_id = dat.columns[0]
dat = dat.drop([feat_id], axis = 1)

# store column names
cnames = dat.columns.tolist()

# exclude genes with missing values and convert to a numpy array
dat = dat.dropna().to_numpy()

# fit gaussian mixture model and retrieve cluster assignments
gmm = GaussianMixture(n_components = 2)
gmm.fit(dat.T)

clusters = gmm.predict(dat.T)

# store result
res = pd.DataFrame({'covariate': cnames, 'cluster': clusters})
res.to_feather(snakemake.output[0])
