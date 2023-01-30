#!/bin/env python
#
# Clusters MM29 covariates
# V. Keith Hughitt
#
import sys
import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

snek = snakemake

NUM_CLUSTERS:int = snek.config['clustering']['num_clusters']

np.random.seed(1)

# load p-values
dat = pd.read_feather(snek.input[0])

feat_id = dat.columns[0]
dat = dat.drop([feat_id], axis = 1)

# store column names
cnames = dat.columns.tolist()

if dat.isnull().sum().sum() > 0:
    num_nas = dat.shape[1] - dat.apply(lambda x: x.count(), axis=1)

    mask = num_nas <= snek.config['clustering']['max_missing']

    dat = dat[mask]

    # impute remaining missing values using k-NN
    imputer = KNNImputer(n_neighbors=4, weights="uniform")
    dat = imputer.fit_transform(dat.T).T
else:
    # if no missing values are present, simply convert to an ndarray
    dat = dat.to_numpy()

# fit gaussian mixture model and retrieve cluster assignments
if snek.config['clustering']['method'] == 'gaussian_mixture':
    from sklearn.mixture import GaussianMixture
    gmm = GaussianMixture(n_components = NUM_CLUSTERS)
    gmm.fit(dat.T)

    clusters = gmm.predict(dat.T)
elif snek.config['clustering']['method'] == 'kmeans':
    from sklearn.cluster import KMeans
    kmeans = KMeans(init='k-means++', n_clusters=NUM_CLUSTERS).fit(dat.T)

    clusters = kmeans.labels_
else:
    sys.exit("Invalid clustering method specified!")

# store result
res = pd.DataFrame({'covariate': cnames, 'cluster': clusters})
res.to_feather(snek.output[0])
