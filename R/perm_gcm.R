#' Permutation-Calibrated Generalised Covariance Measure Test
#'
#' Tests conditional independence X _||_ Y | Z using a permutation-calibrated
#' variant of the GCM test. The test is exact under the strong null
#' X _||_ (Y, Z) and asymptotically valid under the weaker null X _||_ Y | Z
#' (under permutation stability conditions).
#'
#' @param Y Response variable. Numeric vector or single-column matrix.
#' @param X Covariate of interest. Numeric vector or matrix.
#' @param Z Conditioning variables. Numeric vector or matrix.
#' @param reg_YonZ Regression method for Y on Z. Character string (default "rf")
#'   or a function. See \code{comets} package for available methods.
#' @param reg_XonZ Regression method for X on Z. Character string (default "rf")
#'   or a function.
#' @param args_YonZ Optional list of additional arguments for the Y on Z regression.
#' @param args_XonZ Optional list of additional arguments for the X on Z regression.
#' @param B Number of permutations (default 499). When \code{h} is specified,
#'   this serves as the deterministic cap \eqn{T_{max}} on the number of
#'   sequential permutations.
#' @param n_cores Number of cores for parallel permutation loop (default 1,
#'   sequential). Only used when the lrm fast path is not applicable and
#'   adaptive permutation is not active. Requires the \code{parallel} package.
#'   Uses \code{parallel::mclapply} (not available on Windows).
#' @param h Tuning parameter for adaptive permutation (default NULL, which
#'   uses the standard fixed-B permutation). A positive integer controlling
#'   early stopping: the procedure stops once \code{h} permutation scores
#'   exceed the observed score (futility) or the anytime-valid p-value drops
#'   below \code{alpha} (rejection). Recommended values are between 10 and 30.
#' @param alpha Significance level for adaptive stopping rules (default 0.05).
#'   Only used when \code{h} is not NULL.
#' @param alternative Type of alternative hypothesis. One of \code{"two.sided"},
#'   \code{"greater"} (default), or \code{"less"}.
#' @param return_fitted_models Logical. If TRUE, return fitted regression models
#'   in the output.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Fit Y on Z to obtain residuals rY (kept fixed).
#'   \item Fit X on Z to obtain residuals rX.
#'   \item Compute residual products and the studentized statistic T_obs.
#'   \item For each of B permutations: permute X, refit X on Z to get new
#'     residuals, recompute the studentized statistic.
#'   \item Compute permutation p-value.
#' }
#'
#' The key insight is that rY is kept fixed across permutations while X is
#' permuted and refitted, ensuring exactness under the strong null by
#' exchangeability of X.
#'
#' @return An object of class \code{c("perm_gcm", "htest")} with components:
#' \describe{
#'   \item{statistic}{The observed studentized test statistic.}
#'   \item{p.value}{Permutation p-value.}
#'   \item{parameter}{Number of permutations B.}
#'   \item{method}{Description of the test method.}
#'   \item{data.name}{Description of the data.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{rY}{Residuals from Y on Z regression.}
#'   \item{rX}{Residuals from X on Z regression (original, unpermuted).}
#'   \item{perm_stats}{Vector of permutation test statistics (length B, or
#'     shorter if adaptive permutation stopped early).}
#'   \item{adaptive}{List with adaptive permutation details (if \code{h} was
#'     specified): \code{h}, \code{alpha}, \code{K} (loss count),
#'     \code{t_star} (step of h-th loss), \code{stopped_early} (logical),
#'     \code{stop_reason} ("rejection", "futility", or "cap"), and
#'     \code{trajectory} (data frame with columns \code{t}, \code{S},
#'     \code{K}, \code{p_AV}).}
#'   \item{models}{Fitted models (if \code{return_fitted_models = TRUE}).}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 200
#' Z <- matrix(rnorm(n), ncol = 1)
#' X <- Z + rnorm(n)
#' Y <- Z + rnorm(n)  # X _||_ Y | Z
#' result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm", B = 199)
#' print(result)
#' }
#'
#' @export
perm_gcm <- function(Y, X, Z,
                     reg_YonZ = "rf", reg_XonZ = "rf",
                     args_YonZ = NULL, args_XonZ = NULL,
                     B = 499L,
                     n_cores = 1L,
                     h = NULL,
                     alpha = 0.05,
                     alternative = c("two.sided", "greater", "less"),
                     return_fitted_models = FALSE) {
  alternative <- match.arg(alternative)
  B <- as.integer(B)
  adaptive <- !is.null(h)
  if (adaptive) {
    h <- as.integer(h)
    stopifnot(h >= 1L, alpha > 0, alpha < 1)
  }

  ## Data validation
  Y <- check_data(Y, "Y")
  X <- check_data(X, "X")
  Z <- check_data(Z, "Z")
  n <- NROW(Y)

  ## Step 1: Fit Y on Z -> residuals rY (FIXED across permutations)
  mY <- fit_regression(y = Y, x = Z, method = reg_YonZ, args = args_YonZ)
  rY <- stats::residuals(mY, response = Y, data = Z)

  ## Step 2: Fit X on Z -> residuals rX (original)
  ## For gols: inject Y|Z fitted values and theta into args
  if (identical(reg_XonZ, "gols")) {
    mu_Y <- c(Y) - c(rY)  # fitted values from Y|Z model
    if (is.null(args_XonZ[["theta"]])) {
      stop("reg_XonZ = \"gols\" requires args_XonZ$theta ",
           "(NB dispersion parameter).")
    }
    args_XonZ[["mu_Y"]] <- mu_Y
  }
  mX <- fit_regression(y = X, x = Z, method = reg_XonZ, args = args_XonZ)
  rX <- stats::residuals(mX, response = X, data = Z)

  ## Step 3: Observed test statistic
  R_obs <- c(rX) * c(rY)
  T_obs <- studentized_stat(R_obs)

  ## Step 4: Permutation loop
  adaptive_info <- NULL

  if (adaptive) {
    ## --- Adaptive sequential loop ---
    S0 <- score_map(T_obs, alternative)
    K <- 0L
    t_star <- Inf
    T_perm <- numeric(B)
    traj_S <- numeric(B)
    traj_K <- integer(B)
    traj_pAV <- numeric(B)
    stop_reason <- "cap"
    actual_B <- B

    for (b in seq_len(B)) {
      perm <- sample.int(n)
      X_perm <- X[perm, , drop = FALSE]
      mX_perm <- fit_regression(y = X_perm, x = Z, method = reg_XonZ,
                                args = args_XonZ)
      rX_perm <- stats::residuals(mX_perm, response = X_perm, data = Z)
      R_perm <- c(rX_perm) * c(rY)
      T_perm[b] <- studentized_stat(R_perm)

      S_b <- score_map(T_perm[b], alternative)
      traj_S[b] <- S_b
      if (S_b >= S0) {
        K <- K + 1L
        if (K == h && is.infinite(t_star)) t_star <- b
      }
      traj_K[b] <- K
      traj_pAV[b] <- adaptive_pvalue(b, K, h, t_star)

      if (traj_pAV[b] <= alpha) {
        actual_B <- b
        stop_reason <- "rejection"
        break
      }
      if (K == h && traj_pAV[b] > alpha) {
        actual_B <- b
        stop_reason <- "futility"
        break
      }
    }

    T_perm <- T_perm[seq_len(actual_B)]
    pval <- traj_pAV[actual_B]
    adaptive_info <- list(
      h = h, alpha = alpha, K = K, t_star = t_star,
      stopped_early = actual_B < B, stop_reason = stop_reason,
      trajectory = data.frame(
        t = seq_len(actual_B),
        S = traj_S[seq_len(actual_B)],
        K = traj_K[seq_len(actual_B)],
        p_AV = traj_pAV[seq_len(actual_B)]
      )
    )
  } else if (is_lrm_fastpath(reg_XonZ, args_XonZ)) {
    ## --- Vectorized lrm fast path ---
    Q <- precompute_qr_projection(Z)
    perm_mat <- replicate(B, sample.int(n))

    if (NCOL(X) == 1L) {
      ## Univariate X: single batch matrix multiply
      X_perm_mat <- apply_permutations(c(X), perm_mat)
      rX_perm_mat <- fast_lrm_residuals(Q, X_perm_mat)
      R_perm_mat <- rX_perm_mat * c(rY)
      T_perm <- vectorized_studentized_stat(R_perm_mat)
    } else {
      ## Multivariate X: accumulate sums across columns
      d <- NCOL(X)
      nd <- as.numeric(n * d)
      S1 <- numeric(B)
      S2 <- numeric(B)
      for (j in seq_len(d)) {
        Xj_perm_mat <- apply_permutations(X[, j], perm_mat)
        rXj_perm_mat <- fast_lrm_residuals(Q, Xj_perm_mat)
        Rj <- rXj_perm_mat * c(rY)
        S1 <- S1 + colSums(Rj)
        S2 <- S2 + colSums(Rj^2)
      }
      mean_vals <- S1 / nd
      var_vals <- S2 / nd - mean_vals^2
      var_vals[var_vals <= 0] <- Inf
      T_perm <- sqrt(nd) * mean_vals / sqrt(var_vals)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else if (is_nslrm_fastpath(reg_XonZ, args_XonZ)) {
    ## --- Vectorized nslrm fast path (natural spline + OLS) ---
    Q <- precompute_nslrm_projection(mX)
    perm_mat <- replicate(B, sample.int(n))

    if (NCOL(X) == 1L) {
      X_perm_mat <- apply_permutations(c(X), perm_mat)
      rX_perm_mat <- fast_lrm_residuals(Q, X_perm_mat)
      R_perm_mat <- rX_perm_mat * c(rY)
      T_perm <- vectorized_studentized_stat(R_perm_mat)
    } else {
      d <- NCOL(X)
      nd <- as.numeric(n * d)
      S1 <- numeric(B)
      S2 <- numeric(B)
      for (j in seq_len(d)) {
        Xj_perm_mat <- apply_permutations(X[, j], perm_mat)
        rXj_perm_mat <- fast_lrm_residuals(Q, Xj_perm_mat)
        Rj <- rXj_perm_mat * c(rY)
        S1 <- S1 + colSums(Rj)
        S2 <- S2 + colSums(Rj^2)
      }
      mean_vals <- S1 / nd
      var_vals <- S2 / nd - mean_vals^2
      var_vals[var_vals <= 0] <- Inf
      T_perm <- sqrt(nd) * mean_vals / sqrt(var_vals)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else if (is_gols_fastpath(reg_XonZ, args_XonZ)) {
    ## --- Vectorized gols fast path (GLS with NB weights) ---
    gols_proj <- precompute_gols_projection(mX)
    perm_mat <- replicate(B, sample.int(n))

    if (NCOL(X) == 1L) {
      X_perm_mat <- apply_permutations(c(X), perm_mat)
      rX_perm_mat <- fast_gols_residuals(gols_proj, X_perm_mat)
      R_perm_mat <- rX_perm_mat * c(rY)
      T_perm <- vectorized_studentized_stat(R_perm_mat)
    } else {
      d <- NCOL(X)
      nd <- as.numeric(n * d)
      S1 <- numeric(B)
      S2 <- numeric(B)
      for (j in seq_len(d)) {
        Xj_perm_mat <- apply_permutations(X[, j], perm_mat)
        rXj_perm_mat <- fast_gols_residuals(gols_proj, Xj_perm_mat)
        Rj <- rXj_perm_mat * c(rY)
        S1 <- S1 + colSums(Rj)
        S2 <- S2 + colSums(Rj^2)
      }
      mean_vals <- S1 / nd
      var_vals <- S2 / nd - mean_vals^2
      var_vals[var_vals <= 0] <- Inf
      T_perm <- sqrt(nd) * mean_vals / sqrt(var_vals)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    ## --- Parallel loop for non-lrm methods ---
    T_perm <- unlist(parallel::mclapply(seq_len(B), function(b) {
      perm <- sample.int(n)
      X_perm <- X[perm, , drop = FALSE]
      mX_perm <- fit_regression(y = X_perm, x = Z, method = reg_XonZ,
                                args = args_XonZ)
      rX_perm <- stats::residuals(mX_perm, response = X_perm, data = Z)
      R_perm <- c(rX_perm) * c(rY)
      studentized_stat(R_perm)
    }, mc.cores = n_cores))
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else {
    ## --- Sequential loop (original) ---
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      perm <- sample.int(n)
      X_perm <- X[perm, , drop = FALSE]
      mX_perm <- fit_regression(y = X_perm, x = Z, method = reg_XonZ,
                                args = args_XonZ)
      rX_perm <- stats::residuals(mX_perm, response = X_perm, data = Z)
      R_perm <- c(rX_perm) * c(rY)
      T_perm[b] <- studentized_stat(R_perm)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  }

  ## Return htest object
  stat <- T_obs
  names(stat) <- "T"
  param <- c("B" = length(T_perm))

  models <- if (return_fitted_models) {
    list(reg_YonZ = mY, reg_XonZ = mX)
  } else {
    NULL
  }

  method_str <- if (adaptive) {
    "Adaptive permutation-calibrated generalised covariance measure test"
  } else {
    "Permutation-calibrated generalised covariance measure test"
  }

  structure(list(
    statistic = stat,
    p.value = pval,
    parameter = param,
    method = method_str,
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    alternative = alternative,
    rY = rY,
    rX = rX,
    perm_stats = T_perm,
    adaptive = adaptive_info,
    models = models
  ), class = c("perm_gcm", "htest"))
}
