#' Bayesian simple linear regression
#'
#' Computes the posterior distribution of the slope in a simple linear
#' regression with known residual variance. The intercept is integrated out by
#' centering both the predictor and response. The slope has a normal prior with
#' mean zero and variance `prior_var`.
#'
#' @param y A numeric vector containing the response values.
#' @param x A numeric vector containing the predictor values.
#' @param prior_var A positive numeric scalar giving the prior variance of the
#'   regression coefficient.
#' @param residual_var A positive numeric scalar giving the known residual
#'   variance.
#'
#' @return A named list with `posterior_mean` and `posterior_var`, the mean and
#'   variance of the normal posterior distribution of the regression
#'   coefficient.
#' @export
#'
#' @examples
#' x <- 1:5
#' y <- 1 + 2 * x
#' simple_blm(y, x, prior_var = 10, residual_var = 1)
simple_blm <- function(y, x, prior_var, residual_var) {
  if (!is.numeric(y) || !is.atomic(y) || is.object(y) || !is.null(dim(y))) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(x) || !is.atomic(x) || is.object(x) || !is.null(dim(x))) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (length(y) != length(x)) {
    stop("`y` and `x` must have the same length.", call. = FALSE)
  }
  if (length(y) < 2L) {
    stop("`y` and `x` must contain at least two observations.", call. = FALSE)
  }
  if (anyNA(y) || anyNA(x) || any(!is.finite(y)) || any(!is.finite(x))) {
    stop("`y` and `x` must contain only finite, non-missing values.", call. = FALSE)
  }

  validate_variance <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value <= 0) {
      stop(sprintf("`%s` must be a positive, finite numeric scalar.", name),
           call. = FALSE)
    }
  }

  validate_variance(prior_var, "prior_var")
  validate_variance(residual_var, "residual_var")

  x_centered <- x - mean(x)
  y_centered <- y - mean(y)

  posterior_var <- 1 / (sum(x_centered^2) / residual_var + 1 / prior_var)
  posterior_mean <- posterior_var * sum(x_centered * y_centered) / residual_var

  list(
    posterior_mean = posterior_mean,
    posterior_var = posterior_var
  )
}
