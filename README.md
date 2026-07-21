# BayesLinReg

This R package implements Bayesian linear regression with independent normal,
spike-and-slab, and beta-prime global-local priors on the regression
coefficients. The global-local family defaults to the Strawderman-Berger prior
and includes the horseshoe as a special case. Multiple regression fits can
combine multiple predictor blocks through a BGLR-style `ETA` interface. Each
block controls its own standardization, prior family, and prior parameters,
while coefficients are always returned on their original scale. Gibbs sampling
is available in R and Rcpp with optional parallel chains.

`blm()` always receives predictors and coefficient priors through `ETA`. For a
single predictor block, use the single-block shorthand:

```r
fit <- blm(
  y,
  ETA = list(X = X, model = "Normal", var = 10),
  residual_var = 1
)
```

The available models are `"Normal"`, `"SpikeSlab"`, and `"GlobalLocal"`.
For mixed priors, use a named `ETA` list whose blocks specify their own
predictors, model, standardization, and prior parameters.

It was created by OpenAI Codex, supervised by Fabio Morgante. It has not been reviewed and tested carefully.
