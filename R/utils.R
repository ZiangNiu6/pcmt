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
  idx <- indices
  if (is.null(idx)) idx <- sample.int(NROW(Y), ceiling(frac * NROW(Y)))
  X <- as.data.frame(X)
  Z <- as.data.frame(Z)
  list(
    Ytr = Y[idx],
    Xtr = as.matrix(X[idx, , drop = FALSE]),
    Ztr = as.matrix(Z[idx, , drop = FALSE]),
    Yte = Y[-idx],
    Xte = as.matrix(X[-idx, , drop = FALSE]),
    Zte = as.matrix(Z[-idx, , drop = FALSE]),
    idx = idx
  )
}

# ---- Vectorized helpers for lrm fast path ----

#' Check whether the lrm fast path applies
#'
#' Returns TRUE when the regression method is \code{"lrm"} and no extra
#' arguments are supplied, so that residuals can be obtained via a single
#' QR-based matrix projection instead of repeated \code{lm()} calls.
#'
#' @param method Regression method string.
#' @param args Optional argument list passed to the regression.
#' @return Logical scalar.
#' @keywords internal
is_lrm_fastpath <- function(method, args) {
  identical(method, "lrm") && (is.null(args) || length(args) == 0L)
}

#' Precompute Q matrix from QR decomposition of design matrix
#'
#' Computes the (thin) Q factor of the QR decomposition of \code{cbind(1, Z)},
#' truncated to the numerical rank.  Done once, cost O(np^2).
#'
#' @param Z Predictor matrix (n x p).
#' @return Q matrix (n x rank), suitable for \code{\link{fast_lrm_residuals}}.
#' @keywords internal
precompute_qr_projection <- function(Z) {
  design <- cbind(1, as.matrix(Z))
  qr_obj <- qr(design)
  qr.Q(qr_obj)[, seq_len(qr_obj$rank), drop = FALSE]
}

#' Compute lrm residuals for multiple response columns via matrix ops
#'
#' Given the Q factor from the design matrix, computes
#' \code{Y_mat - Q \%*\% crossprod(Q, Y_mat)} for an n x K response matrix.
#' Cost is O(npK) instead of K separate O(np^2) QR calls.
#'
#' @param Q Q matrix from \code{\link{precompute_qr_projection}}.
#' @param Y_mat Response matrix (n x K).
#' @return Residual matrix (n x K).
#' @keywords internal
fast_lrm_residuals <- function(Q, Y_mat) {
  Y_mat - Q %*% crossprod(Q, Y_mat)
}

#' Build n x B permuted response matrix via vectorized indexing
#'
#' @param x_vec Numeric vector of length n.
#' @param perm_mat Integer matrix (n x B) of permutation indices.
#' @return Numeric matrix (n x B).
#' @keywords internal
apply_permutations <- function(x_vec, perm_mat) {
  matrix(x_vec[perm_mat], nrow = nrow(perm_mat), ncol = ncol(perm_mat))
}

#' Vectorized studentized statistic across columns
#'
#' Computes \code{sqrt(n) * colMeans(R) / sqrt(colMeans(R^2) - colMeans(R)^2)}
#' for each column of R_mat.
#'
#' @param R_mat Numeric matrix (n x B) of residual products.
#' @return Numeric vector of length B.
#' @keywords internal
vectorized_studentized_stat <- function(R_mat) {
  n <- nrow(R_mat)
  cm <- colMeans(R_mat)
  v <- colMeans(R_mat^2) - cm^2
  v[v <= 0] <- Inf
  sqrt(n) * cm / sqrt(v)
}

#' Expand a data frame to a numeric matrix
#'
#' Converts factor/character columns to dummy variables via
#' \code{model.matrix} and leaves numeric columns as-is.
#' Used internally by \code{perm_icp} to handle mixed-type
#' covariates and environments.
#'
#' @param df A data frame.
#' @return A numeric matrix with no intercept column.
#' @keywords internal
.expand_to_numeric <- function(df) {
  df <- as.data.frame(df)
  if (ncol(df) == 0L) return(matrix(nrow = nrow(df), ncol = 0L))
  has_factor <- any(vapply(df, function(x) is.factor(x) || is.character(x),
                           logical(1)))
  if (!has_factor) return(as.matrix(df))
  # Use model.matrix to expand factors to dummies (no intercept)
  mm <- stats::model.matrix(~ . - 1, data = df)
  unname(mm)
}
