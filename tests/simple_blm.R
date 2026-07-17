library(blm)

fit <- simple_blm(
  y = c(3, 5, 7, 9, 11),
  x = 1:5,
  prior_var = 10,
  residual_var = 1
)

stopifnot(
  identical(names(fit), c("posterior_mean", "posterior_var")),
  isTRUE(all.equal(fit$posterior_var, 1 / 10.1)),
  isTRUE(all.equal(fit$posterior_mean, 20 / 10.1))
)

# Centering makes the result invariant to shifts in x and y.
shifted_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11) + 100,
  x = 1:5 + 50,
  prior_var = 10,
  residual_var = 1
)
stopifnot(isTRUE(all.equal(fit, shifted_fit)))

# With an inverse-gamma prior, the residual variance is learned from the data.
learned_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11),
  x = 1:5,
  prior_var = 10,
  residual_shape = 2,
  residual_scale = 1
)
expected_relative_var <- 1 / 10.1
expected_mean <- 20 / 10.1
expected_shape <- 4
expected_scale <- 1 + 0.5 * (
  sum((c(-4, -2, 0, 2, 4) - expected_mean * c(-2, -1, 0, 1, 2))^2) +
    expected_mean^2 / 10
)
stopifnot(
  isTRUE(all.equal(learned_fit$posterior_mean, expected_mean)),
  isTRUE(all.equal(learned_fit$residual_shape, expected_shape)),
  isTRUE(all.equal(learned_fit$residual_scale, expected_scale)),
  isTRUE(all.equal(learned_fit$posterior_df, 2 * expected_shape)),
  isTRUE(all.equal(
    learned_fit$posterior_var,
    expected_scale / (expected_shape - 1) * expected_relative_var
  ))
)

learned_shifted_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11) + 100,
  x = 1:5 + 50,
  prior_var = 10,
  residual_shape = 2,
  residual_scale = 1
)
stopifnot(isTRUE(all.equal(learned_fit, learned_shifted_fit)))

# Invalid variances and incompatible inputs are rejected.
stopifnot(
  inherits(try(simple_blm(1:3, 1:2, 1, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 0, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1, NA_real_), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1), silent = TRUE), "try-error"),
  inherits(
    try(simple_blm(1:3, 1:3, 1, 1, residual_shape = 2), silent = TRUE),
    "try-error"
  )
)
