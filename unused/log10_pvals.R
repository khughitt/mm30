
# in addition to computing p-value summary statistics, we will also generate some
# p-values based on -log10-transformed p-values
mean_neg_log10_pvals <- apply(scores_mat, 1, function(x) {
  mean(-log10(pmax(x, 1E-20)), na.rm = TRUE)
})

median_neg_log10_pvals <- apply(scores_mat, 1, function(x) {
  median(-log10(pmax(x, 1E-20)), na.rm = TRUE)
})

max_neg_log10_pvals <- apply(scores_mat, 1, function(x) {
  max(-log10(pmax(x, 1E-20)), na.rm = TRUE)
})

# -log10 p-value summary dataframe
res <- data.frame(
  pull(pvals, id_field),
  mean_neg_log10_pval,
  median_neg_log10_pval,
  max_neg_log10_pval,
  num_present,
  num_missing
)
colnames(res)[1] <- id_field
