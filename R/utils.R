#' Studentized test statistic
#'
#' Computes the studentized mean statistic T_n = sqrt(n) * mean(R) / sd(R),
#' where sd uses the denominator n (population variance) as in the papers.
#'
#' @param R Numeric vector of residual products.
#' @return Scalar studentized statistic.
#' @keywords internal
#' @noRd
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
#' @noRd
perm_pvalue <- function(T_obs, T_perm, alternative = "greater") {
  B <- length(T_perm)
  switch(alternative,
    "greater" = (1 + sum(T_perm >= T_obs)) / (1 + B),
    "less" = (1 + sum(T_perm <= T_obs)) / (1 + B),
    "two.sided" = (1 + sum(abs(T_perm) >= abs(T_obs))) / (1 + B)
  )
}

#' Score map for adaptive permutation
#'
#' Maps a test statistic to a score based on the alternative hypothesis.
#' For one-sided tests, the score is the statistic itself (or negated).
#' For two-sided tests, the score is the absolute value.
#'
#' @param T_stat Numeric scalar or vector of test statistics.
#' @param alternative One of "two.sided", "greater", "less".
#' @return Numeric score(s).
#' @keywords internal
#' @noRd
score_map <- function(T_stat, alternative) {
  switch(alternative,
    "greater" = T_stat,
    "less" = -T_stat,
    "two.sided" = abs(T_stat)
  )
}

#' Anytime-valid Monte Carlo p-value
#'
#' Computes the anytime-valid p-value from the adaptive permutation procedure.
#' Given the current count of losses K (permutation scores >= observed score)
#' at step t, returns \code{h / (t + h - K)} if \code{K < h}, or
#' \code{h / t_star} if \code{K >= h}.
#'
#' @param t Current step (number of permutations drawn so far).
#' @param K Number of losses (permutation scores >= observed score).
#' @param h Tuning parameter (positive integer).
#' @param t_star Step at which the h-th loss occurred (Inf if not yet).
#' @return Scalar anytime-valid p-value.
#' @keywords internal
#' @noRd
adaptive_pvalue <- function(t, K, h, t_star) {
  if (K < h) {
    h / (t + h - K)
  } else {
    h / t_star
  }
}

#' Fit a random forest nuisance regression
#'
#' Uses the random forest regression implementation from comets. This release
#' supports only the method string \code{"rf"}; additional arguments are passed
#' unchanged to that learner.
#'
#' @param y Response vector or matrix.
#' @param x Predictor matrix.
#' @param method The regression method, which must be \code{"rf"}.
#' @param args Optional list of additional arguments for the random forest.
#' @return A fitted comets random forest model.
#' @keywords internal
#' @noRd
fit_regression <- function(y, x, method = "rf", args = NULL) {
  if (!identical(method, "rf")) {
    stop('This release supports only reg_* = "rf".', call. = FALSE)
  }
  fn <- get("rf", envir = asNamespace("comets"))
  do.call(fn, c(list(y = y, x = x), args))
}

#' Check and standardize data
#'
#' Validates and coerces input data to the expected format.
#' Y is returned as a vector (for PCM, must be univariate).
#' X and Z are returned as matrices with column names.
#'
#' @param x Data to check.
#' @param mode One of "Y", "X", "Z".
#' @param test One of "gcm", "pcm".
#' @return Standardized data.
#' @keywords internal
#' @noRd
check_data <- function(x, mode = "Y", test = "gcm") {
  mode <- match.arg(mode, c("Y", "X", "Z"))
  if (mode == "Y") {
    if (!is.matrix(x) && !is.data.frame(x)) {
      ret <- c(x)
      if (is.factor(ret) && length(levels(ret)) > 2 && test == "pcm")
        stop("Only binary factors are allowed for Y.")
      return(ret)
    } else {
      if (NCOL(x) > 1 && test == "pcm")
        stop("Please provide Y as a vector.")
      return(check_data(x, mode = "X"))
    }
  }
  if (inherits(x, "tbl_df")) x <- as.data.frame(x)
  if (NCOL(x) == 0) x <- matrix(1, nrow = NROW(x), ncol = 1)
  if (NCOL(x) == 1) x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0(mode, seq_len(NCOL(x)))
  x
}
