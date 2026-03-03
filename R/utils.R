#' Studentized test statistic
#'
#' Computes the studentized mean statistic T_n = sqrt(n) * mean(R) / sd(R),
#' where sd uses the denominator n (population variance) as in the papers.
#'
#' @param R Numeric vector of residual products.
#' @return Scalar studentized statistic.
#' @keywords internal
studentized_stat <- function(R) {
  n <- length(R)
  m <- mean(R)
  v <- mean(R^2) - m^2
  if (v <= 0) return(0)
  sqrt(n) * m / sqrt(v)
}

#' Permutation p-value
#'
#' Computes permutation p-value with the (1 + ...) / (1 + B) formula.
#'
#' @param T_obs Observed test statistic (scalar).
#' @param T_perm Vector of permutation test statistics.
#' @param alternative One of "two.sided", "greater", "less".
#' @return Permutation p-value.
#' @keywords internal
perm_pvalue <- function(T_obs, T_perm, alternative = "greater") {
  B <- length(T_perm)
  switch(alternative,
    "greater" = (1 + sum(T_perm >= T_obs)) / (1 + B),
    "less" = (1 + sum(T_perm <= T_obs)) / (1 + B),
    "two.sided" = (1 + sum(abs(T_perm) >= abs(T_obs))) / (1 + B)
  )
}

#' Fit a regression model using comets infrastructure
#'
#' Thin wrapper around comets regression functions. Fits a model of y on x
#' using the specified method string (e.g., "rf", "lrm", "lasso").
#'
#' @param y Response vector or matrix.
#' @param x Predictor matrix.
#' @param method Character string specifying the regression method.
#' @param args Optional list of additional arguments passed to the fitting function.
#' @return Fitted model object (comets S3 class).
#' @keywords internal
fit_regression <- function(y, x, method = "rf", args = NULL) {
  fn <- get(method, envir = asNamespace("comets"))
  do.call(fn, c(list(y = y, x = x), args))
}

#' Check and standardize data using comets infrastructure
#'
#' @param x Data to check.
#' @param mode One of "Y", "X", "Z".
#' @param test One of "gcm", "pcm".
#' @return Standardized data.
#' @keywords internal
check_data <- function(x, mode = "Y", test = "gcm") {
  comets:::.check_data(x, mode = mode, test = test)
}

#' Split sample into two folds
#'
#' @param Y Response.
#' @param X Covariates.
#' @param Z Conditioning variables.
#' @param frac Fraction for first fold.
#' @param indices Optional pre-specified indices for first fold.
#' @return List with train and test data.
#' @keywords internal
split_sample <- function(Y, X, Z, frac = 0.5, indices = NULL) {
  comets:::.split_sample(Y, X, Z, frac = frac, indices = indices)
}
