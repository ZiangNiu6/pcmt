## Timing benchmark: current perm_gcm vs naive for-loop
## Setup (a) from Shah & Peters (2020): Z ~ N(0,1), X = fa(Z,2)+0.3*eps,
##                                      Y = fa(Z,2)+0.3*eps
## Regression: xgboost with 10-fold CV over maxdepth (paper's method)
##
## Both implementations use the SAME pre-generated permutation matrix
## so that p-values must coincide exactly.

library(pcmt)
library(xgboost)

fa <- function(z, a) exp(-z^2 / 2) * sin(a * z)

# ---------------------------------------------------------------------------
# xgboost with 10-fold CV over maxdepth  (paper's method)
# ---------------------------------------------------------------------------
xgb_cv_maxdepth <- function(y, x, max_depths = 1:6, nfold = 10L,
                            nrounds = 50L, verbose = 0L, ...) {
  dtrain <- xgb.DMatrix(data = as.matrix(x), label = y)
  best_err <- Inf; best_md <- 1L
  for (md in max_depths) {
    cv <- xgb.cv(data = dtrain,
      params = list(max_depth = md, eta = 0.1, objective = "reg:squarederror"),
      nrounds = nrounds, nfold = nfold, verbose = verbose)
    err <- min(cv$evaluation_log$test_rmse_mean)
    if (err < best_err) { best_err <- err; best_md <- md }
  }
  bst <- xgb.train(data = dtrain,
    params = list(max_depth = best_md, eta = 0.1, objective = "reg:squarederror"),
    nrounds = nrounds, verbose = verbose)
  class(bst) <- c("xgb", class(bst))
  bst
}

# Patch fit_regression so it can find xgb_cv_maxdepth in globalenv
assignInNamespace("fit_regression", function(y, x, method = "rf", args = NULL) {
  fn <- tryCatch(get(method, envir = asNamespace("comets")),
                 error = function(e) get(method, envir = globalenv()))
  do.call(fn, c(list(y = y, x = x), args))
}, ns = "pcmt")

# ---------------------------------------------------------------------------
# Naive permutation loop  (the code before our optimization)
# Takes a pre-generated perm_mat (n x B) so results are reproducible.
# ---------------------------------------------------------------------------
naive_perm_gcm <- function(Y, X, Z, reg_YonZ, reg_XonZ, perm_mat,
                           args_YonZ = NULL, args_XonZ = NULL) {
  fit_reg <- get("fit_regression", envir = asNamespace("pcmt"))
  student <- get("studentized_stat", envir = asNamespace("pcmt"))
  perm_pv <- get("perm_pvalue", envir = asNamespace("pcmt"))

  n <- NROW(Y)
  B <- ncol(perm_mat)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  mY <- fit_reg(y = Y, x = Z, method = reg_YonZ, args = args_YonZ)
  rY <- stats::residuals(mY, response = Y, data = Z)

  mX <- fit_reg(y = X, x = Z, method = reg_XonZ, args = args_XonZ)
  rX <- stats::residuals(mX, response = X, data = Z)

  R_obs <- c(rX) * c(rY)
  T_obs <- student(R_obs)

  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    X_perm <- X[perm_mat[, b], , drop = FALSE]
    mX_perm <- fit_reg(y = X_perm, x = Z, method = reg_XonZ, args = args_XonZ)
    rX_perm <- stats::residuals(mX_perm, response = X_perm, data = Z)
    R_perm <- c(rX_perm) * c(rY)
    T_perm[b] <- student(R_perm)
  }

  list(p.value = perm_pv(T_obs, T_perm, "two.sided"),
       T_obs = T_obs, T_perm = T_perm)
}

# ---------------------------------------------------------------------------
# Current (optimized) implementation, also using a fixed perm_mat.
# For lrm fast-path: we feed the same perm_mat through the vectorized code.
# For xgb: falls through to the sequential loop with the same perm_mat.
# ---------------------------------------------------------------------------
current_perm_gcm <- function(Y, X, Z, reg_YonZ, reg_XonZ, perm_mat,
                             args_YonZ = NULL, args_XonZ = NULL) {
  fit_reg  <- get("fit_regression", envir = asNamespace("pcmt"))
  student  <- get("studentized_stat", envir = asNamespace("pcmt"))
  perm_pv  <- get("perm_pvalue", envir = asNamespace("pcmt"))
  is_fast  <- get("is_lrm_fastpath", envir = asNamespace("pcmt"))
  precomp  <- get("precompute_qr_projection", envir = asNamespace("pcmt"))
  fast_res <- get("fast_lrm_residuals", envir = asNamespace("pcmt"))
  apply_pm <- get("apply_permutations", envir = asNamespace("pcmt"))
  vec_stat <- get("vectorized_studentized_stat", envir = asNamespace("pcmt"))

  n <- NROW(Y)
  B <- ncol(perm_mat)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  mY <- fit_reg(y = Y, x = Z, method = reg_YonZ, args = args_YonZ)
  rY <- stats::residuals(mY, response = Y, data = Z)

  mX <- fit_reg(y = X, x = Z, method = reg_XonZ, args = args_XonZ)
  rX <- stats::residuals(mX, response = X, data = Z)

  R_obs <- c(rX) * c(rY)
  T_obs <- student(R_obs)

  if (is_fast(reg_XonZ, args_XonZ)) {
    ## Vectorized lrm fast path
    Q <- precomp(Z)
    if (NCOL(X) == 1L) {
      X_perm_mat <- apply_pm(c(X), perm_mat)
      rX_perm_mat <- fast_res(Q, X_perm_mat)
      R_perm_mat <- rX_perm_mat * c(rY)
      T_perm <- vec_stat(R_perm_mat)
    } else {
      d <- NCOL(X); nd <- as.numeric(n * d)
      S1 <- numeric(B); S2 <- numeric(B)
      for (j in seq_len(d)) {
        Xj_perm_mat <- apply_pm(X[, j], perm_mat)
        rXj_perm_mat <- fast_res(Q, Xj_perm_mat)
        Rj <- rXj_perm_mat * c(rY)
        S1 <- S1 + colSums(Rj); S2 <- S2 + colSums(Rj^2)
      }
      mean_v <- S1 / nd; var_v <- S2 / nd - mean_v^2
      var_v[var_v <= 0] <- Inf
      T_perm <- sqrt(nd) * mean_v / sqrt(var_v)
    }
  } else {
    ## Sequential loop (same as naive)
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      X_perm <- X[perm_mat[, b], , drop = FALSE]
      mX_perm <- fit_reg(y = X_perm, x = Z, method = reg_XonZ, args = args_XonZ)
      rX_perm <- stats::residuals(mX_perm, response = X_perm, data = Z)
      R_perm <- c(rX_perm) * c(rY)
      T_perm[b] <- student(R_perm)
    }
  }

  list(p.value = perm_pv(T_obs, T_perm, "two.sided"),
       T_obs = T_obs, T_perm = T_perm)
}

# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------
sample_sizes <- c(100L)
B <- 499L

cat("================================================================\n")
cat("Setup (a) variant: Z~N(0,I_5), X=fa(Z1,2)+0.3eps, Y=fa(Z1,2)+0.3eps\n")
cat(sprintf("B = %d permutations (shared permutation matrix)\n", B))
cat("================================================================\n\n")

for (n in sample_sizes) {
  set.seed(42)
  p_Z <- 5L
  Z <- matrix(rnorm(n * p_Z), ncol = p_Z)
  X <- fa(Z[, 1], 2) + 0.3 * rnorm(n)
  Y <- fa(Z[, 1], 2) + 0.3 * rnorm(n)
  perm_mat <- replicate(B, sample.int(n))   # shared across both

  cat(sprintf("n = %d\n", n))
  cat("--------\n")

  ## --- lasso (cv.glmnet) ---
  cat("  reg = lasso (cv.glmnet)\n")

  set.seed(999)
  t3 <- system.time(suppressWarnings(
    r_naive_lasso <- naive_perm_gcm(Y, X, Z, "lasso", "lasso", perm_mat)
  ))[["elapsed"]]

  set.seed(999)
  t4 <- system.time(suppressWarnings(
    r_curr_lasso <- current_perm_gcm(Y, X, Z, "lasso", "lasso", perm_mat)
  ))[["elapsed"]]

  cat(sprintf("    naive loop:   %7.2f s   p = %.4f\n", t3, r_naive_lasso$p.value))
  cat(sprintf("    current:      %7.2f s   p = %.4f\n", t4, r_curr_lasso$p.value))
  cat(sprintf("    p-values match: %s\n",
              identical(r_naive_lasso$p.value, r_curr_lasso$p.value)))
  cat(sprintf("    max |T_perm diff|: %.2e\n",
              max(abs(r_naive_lasso$T_perm - r_curr_lasso$T_perm))))
  cat(sprintf("    speedup: %.2fx\n\n", t3 / t4))

  ## --- lrm ---
  cat("  reg = lrm\n")

  t5 <- system.time(
    r_naive_lrm <- naive_perm_gcm(Y, X, Z, "lrm", "lrm", perm_mat)
  )[["elapsed"]]

  t6 <- system.time(
    r_curr_lrm <- current_perm_gcm(Y, X, Z, "lrm", "lrm", perm_mat)
  )[["elapsed"]]

  cat(sprintf("    naive loop:   %7.3f s   p = %.4f\n", t5, r_naive_lrm$p.value))
  cat(sprintf("    current:      %7.3f s   p = %.4f\n", t6, r_curr_lrm$p.value))
  cat(sprintf("    p-values match: %s\n",
              identical(r_naive_lrm$p.value, r_curr_lrm$p.value)))
  cat(sprintf("    max |T_perm diff|: %.2e\n",
              max(abs(r_naive_lrm$T_perm - r_curr_lrm$T_perm))))
  cat(sprintf("    speedup: %.2fx\n\n", t5 / t6))
}
