# BayesLinReg

This R package implements Bayesian linear regression with independent normal,
spike-and-slab, and beta-prime global-local priors on the regression
coefficients. The global-local family defaults to the Strawderman-Berger prior
and includes the horseshoe as a special case. Multiple regression fits can
combine multiple predictor blocks through a BGLR-style `ETA` interface. Each
block controls its own standardization, prior family, and prior parameters,
while coefficients are always returned on their original scale. Gibbs sampling
is available in R and Rcpp with optional parallel chains.

It was created by OpenAI Codex, supervised by Fabio Morgante. It has not been reviewed and tested carefully.
