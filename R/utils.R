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
#' Also supports \code{"gam"} for a generalized additive model with negative
#' binomial family via \code{mgcv::gam}, which is handled internally rather
#' than through comets.
#'
#' @param y Response vector or matrix.
#' @param x Predictor matrix.
#' @param method Character string specifying the regression method.
#'   Use \code{"gam"} for \code{mgcv::gam} with \code{nb()} family.
#' @param args Optional list of additional arguments passed to the fitting function.
#'   For \code{"gam"}, supported arguments include \code{family} (default
#'   \code{mgcv::nb()}), \code{k} (basis dimension, default -1), and
#'   \code{formula} (custom GAM formula, e.g.
#'   \code{~ s(Z1) + s(Z2) + ti(Z1, Z2)}).
#' @return Fitted model object (comets S3 class, or \code{"gam_fit"} for GAM).
#' @keywords internal
fit_regression <- function(y, x, method = "rf", args = NULL) {
  if (identical(method, "gam")) {
    return(do.call(fit_gam, c(list(y = y, x = x), args)))
  }
  if (identical(method, "nslrm")) {
    return(do.call(fit_nslrm, c(list(y = y, x = x), args)))
  }
  fn <- get(method, envir = asNamespace("comets"))
  do.call(fn, c(list(y = y, x = x), args))
}

#' Fit a GAM with negative binomial family
#'
#' Fits a generalized additive model using \code{mgcv::gam} with smooth terms
#' for each predictor. Intended for count responses where a flexible
#' mean function is desired (e.g., Y | Z regression in the GCM test).
#'
#' By default every column of \code{x} gets its own \code{s()} smooth term.
#' Supply a custom \code{formula} to override this, e.g. to add tensor
#' interaction smooths.  The formula should use \code{.y} as the response
#' and the column names of \code{x} as predictors.  If \code{x} has no
#' column names, they are set to \code{V1, V2, ...}.
#'
#' @param y Response vector (non-negative integers for NB family).
#' @param x Predictor matrix. Each column gets a smooth term \code{s()}
#'   unless \code{formula} is provided.
#' @param family GAM family (default \code{mgcv::nb()} for negative binomial
#'   with estimated theta).
#' @param k Basis dimension for smooth terms (default -1, automatic).
#'   Ignored when \code{formula} is supplied.
#' @param formula Optional one-sided or two-sided formula for the GAM
#'   right-hand side.  Use column names of \code{x} (or \code{V1, V2, ...})
#'   as predictors and \code{.y} as the response.
#'   Examples: \code{~ s(Z1) + s(Z2) + ti(Z1, Z2)},
#'   \code{.y ~ s(Z1) + s(Z2) + ti(Z1, Z2)}.
#' @return An object of class \code{"gam_fit"} containing the fitted
#'   \code{mgcv::gam} model and training data.
#' @keywords internal
fit_gam <- function(y, x, family = NULL, k = -1, formula = NULL) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for method = \"gam\". ",
         "Install it with install.packages(\"mgcv\").")
  }
  if (is.null(family)) family <- mgcv::nb()
  x <- as.matrix(x)
  p <- ncol(x)
  col_nms <- if (!is.null(colnames(x))) colnames(x) else paste0("V", seq_len(p))
  colnames(x) <- col_nms
  dat <- as.data.frame(x)
  dat$.y <- c(y)
  if (!is.null(formula)) {
    ## Accept one-sided (~ rhs) or two-sided (.y ~ rhs) formulas
    if (length(formula) == 2L) {
      ## one-sided: ~ rhs  ->  .y ~ rhs
      fml <- stats::as.formula(paste(".y ~", deparse(formula[[2]])))
    } else {
      fml <- formula
    }
  } else {
    smooth_terms <- paste0("s(", col_nms, ", k=", k, ")", collapse = " + ")
    fml <- stats::as.formula(paste(".y ~", smooth_terms))
  }
  gam_model <- mgcv::gam(fml, data = dat, family = family)
  structure(
    list(model = gam_model, train_x = x, fitted = fitted(gam_model)),
    class = "gam_fit"
  )
}

#' @export
residuals.gam_fit <- function(object, response = NULL, data = NULL, ...) {
  if (!is.null(data)) {
    data <- as.matrix(data)
    if (!is.null(colnames(object$train_x))) {
      colnames(data) <- colnames(object$train_x)
    }
    newdat <- as.data.frame(data)
    pred <- predict(object$model, newdata = newdat, type = "response")
    return(c(response) - c(pred))
  }
  c(response) - object$fitted
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

# ---- Natural spline + OLS (nslrm) ----

#' Fit a natural-spline expanded OLS model
#'
#' Expands each predictor column into a natural spline basis via
#' \code{splines::ns()}, then fits OLS on the combined expanded design.
#' Because the model is linear in the expanded basis, the existing QR-based
#' fast path (\code{\link{precompute_qr_projection}} /
#' \code{\link{fast_lrm_residuals}}) applies directly.
#'
#' When \code{interactions = TRUE} and there are at least two predictors,
#' all pairwise tensor products of the marginal spline bases are appended
#' to the design matrix, adding \code{df^2} columns per pair.  This allows
#' the model to capture interaction effects (e.g. Z1 * Z2) without
#' adaptive smoothing parameter estimation.
#'
#' @param y Response vector or matrix.
#' @param x Predictor matrix.  Each column is expanded independently.
#' @param df Degrees of freedom for each natural spline basis (default 4).
#' @param interactions Logical. If TRUE, include pairwise tensor product
#'   interaction terms between all predictor spline bases (default FALSE).
#' @return An object of class \code{"nslrm_fit"} containing the fitted
#'   \code{lm} model, the list of \code{ns} basis objects (for prediction on
#'   new data), and the original column names.
#' @keywords internal
fit_nslrm <- function(y, x, df = 4, interactions = FALSE) {
  x <- as.matrix(x)
  p <- ncol(x)
  col_nms <- if (!is.null(colnames(x))) colnames(x) else paste0("V", seq_len(p))
  colnames(x) <- col_nms

  bases <- vector("list", p)
  basis_mats <- vector("list", p)
  for (j in seq_len(p)) {
    bases[[j]] <- splines::ns(x[, j], df = df)
    basis_mats[[j]] <- bases[[j]]
  }
  names(bases) <- col_nms

  X_expanded <- do.call(cbind, basis_mats)

  ## Pairwise tensor product interactions
  if (interactions && p >= 2L) {
    ia_mats <- vector("list", p * (p - 1L) / 2L)
    k <- 1L
    for (i in seq_len(p - 1L)) {
      for (j in (i + 1L):p) {
        ia_mats[[k]] <- tensor_basis(basis_mats[[i]], basis_mats[[j]])
        k <- k + 1L
      }
    }
    X_expanded <- cbind(X_expanded, do.call(cbind, ia_mats))
  }

  design <- cbind(1, X_expanded)
  y_mat <- if (is.matrix(y)) y else as.matrix(y)
  fit <- stats::lm.fit(design, y_mat)

  structure(
    list(model = fit, bases = bases, col_names = col_nms,
         interactions = interactions, fitted = fit$fitted.values),
    class = "nslrm_fit"
  )
}

#' Compute row-wise tensor product of two basis matrices
#'
#' Given two matrices A (n x a) and B (n x b), returns an n x (a*b) matrix
#' where each row is the Kronecker product of the corresponding rows of A and B.
#'
#' @param A Basis matrix (n x a).
#' @param B Basis matrix (n x b).
#' @return Matrix (n x a*b).
#' @keywords internal
tensor_basis <- function(A, B) {
  n <- nrow(A)
  a <- ncol(A)
  b <- ncol(B)
  ## A[,i] * B[,j] for all (i,j) pairs
  A[, rep(seq_len(a), each = b)] * B[, rep(seq_len(b), times = a)]
}

#' Expand new data using stored natural spline bases
#'
#' @param x New predictor matrix (n x p).
#' @param bases List of \code{ns} basis objects from \code{\link{fit_nslrm}}.
#' @param interactions Logical. If TRUE, append pairwise tensor products.
#' @return Expanded matrix.
#' @keywords internal
expand_ns_basis <- function(x, bases, interactions = FALSE) {
  x <- as.matrix(x)
  p <- length(bases)
  mats <- vector("list", p)
  for (j in seq_len(p)) {
    mats[[j]] <- predict(bases[[j]], x[, j])
  }
  out <- do.call(cbind, mats)

  if (interactions && p >= 2L) {
    ia_mats <- vector("list", p * (p - 1L) / 2L)
    k <- 1L
    for (i in seq_len(p - 1L)) {
      for (j in (i + 1L):p) {
        ia_mats[[k]] <- tensor_basis(mats[[i]], mats[[j]])
        k <- k + 1L
      }
    }
    out <- cbind(out, do.call(cbind, ia_mats))
  }
  out
}

#' @export
predict.nslrm_fit <- function(object, data = NULL, ...) {
  if (is.null(data)) return(object$fitted)
  ia <- isTRUE(object$interactions)
  X_expanded <- expand_ns_basis(data, object$bases, interactions = ia)
  cfx <- object$model$coefficients
  cfx[is.na(cfx)] <- 0
  c(cbind(1, X_expanded) %*% cfx)
}

#' @export
residuals.nslrm_fit <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict(object, data = data)
  c(response) - c(preds)
}

#' Check whether the nslrm fast path applies
#'
#' Returns TRUE when the regression method is \code{"nslrm"}.  Unlike
#' \code{is_lrm_fastpath}, extra \code{args} (e.g. \code{df}) are allowed
#' because they only affect the one-time basis construction.
#'
#' @param method Regression method string.
#' @param args Optional argument list (ignored for the check).
#' @return Logical scalar.
#' @keywords internal
is_nslrm_fastpath <- function(method, args) {
  identical(method, "nslrm")
}

#' Precompute QR projection from an nslrm fit
#'
#' Extracts the stored natural spline bases, reconstructs the expanded design
#' matrix (including interaction terms if applicable), and calls
#' \code{\link{precompute_qr_projection}} on it.
#'
#' @param nslrm_fit An \code{"nslrm_fit"} object from \code{\link{fit_nslrm}}.
#' @return Q matrix suitable for \code{\link{fast_lrm_residuals}}.
#' @keywords internal
precompute_nslrm_projection <- function(nslrm_fit) {
  basis_mats <- lapply(nslrm_fit$bases, as.matrix)
  X_expanded <- do.call(cbind, basis_mats)
  if (isTRUE(nslrm_fit$interactions)) {
    p <- length(basis_mats)
    if (p >= 2L) {
      ia_mats <- vector("list", p * (p - 1L) / 2L)
      k <- 1L
      for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
          ia_mats[[k]] <- tensor_basis(basis_mats[[i]], basis_mats[[j]])
          k <- k + 1L
        }
      }
      X_expanded <- cbind(X_expanded, do.call(cbind, ia_mats))
    }
  }
  precompute_qr_projection(X_expanded)
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


