#' Permutation-Calibrated Weighted Generalised Covariance Measure Test
#'
#' Tests conditional independence X _||_ Y | Z using a permutation-calibrated
#' variant of the weighted GCM (WGCM) test. The weight function is estimated on
#' a held-out fold to provide power against alternatives where the sign of the
#' conditional covariance varies with Z. The test is exact under the strong null
#' X _||_ (Y, Z) and asymptotically valid under the weaker null X _||_ Y | Z.
#'
#' @param Y Response variable. Numeric vector.
#' @param X Covariate of interest. Numeric vector or matrix.
#' @param Z Conditioning variables. Numeric vector or matrix.
#' @param reg_YonZ Regression method for Y on Z. Character string (default "rf").
#' @param reg_XonZ Regression method for X on Z. Character string (default "rf").
#' @param reg_wfun Regression method for weight function estimation. Character
#'   string (default "rf").
#' @param args_YonZ,args_XonZ,args_wfun Optional lists of additional arguments
#'   for each regression.
#' @param frac Fraction of data for weight estimation fold D1 (default 0.5).
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
#' @param return_fitted_models Logical. If TRUE, return fitted regression models.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Split data into D1 (weight estimation) and D2 (testing).
#'   \item On D1 (fitted once, kept fixed):
#'     \itemize{
#'       \item Fit Y on Z and X on Z, compute residuals.
#'       \item Estimate weight function W(z) = sign of conditional covariance.
#'     }
#'   \item On D2 (testing with permutation calibration):
#'     \itemize{
#'       \item Fit Y on Z to get rY (kept fixed).
#'       \item Fit X on Z to get rX.
#'       \item Apply weights: R = (rX * W) * rY, T_obs = studentized(R).
#'     }
#'   \item Permutation loop on D2: permute X, refit X on Z, apply fixed weights,
#'     compute T_perm.
#' }
#'
#' @return An object of class \code{c("perm_wgcm", "htest")} with components:
#' \describe{
#'   \item{statistic}{The observed studentized test statistic.}
#'   \item{p.value}{Permutation p-value.}
#'   \item{parameter}{Number of permutations B.}
#'   \item{method}{Description of the test method.}
#'   \item{data.name}{Description of the data.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{rY}{Residuals from Y on Z on D2.}
#'   \item{rX}{Residuals from X on Z on D2 (original, unpermuted).}
#'   \item{W}{Estimated weight matrix on D2.}
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
#' result <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
#'   reg_wfun = "lrm", B = 199)
#' print(result)
#' }
#'
#' @export
perm_wgcm <- function(Y, X, Z,
                      reg_YonZ = "rf", reg_XonZ = "rf",
                      reg_wfun = "rf",
                      args_YonZ = NULL, args_XonZ = NULL,
                      args_wfun = NULL,
                      frac = 0.5,
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
  Y <- check_data(Y, "Y", "pcm")  # scalar Y for WGCM
  X <- check_data(X, "X")
  Z <- check_data(Z, "Z")

  ## Sample splitting
  dsp <- split_sample(Y, X, Z, frac = frac)
  Y_D1 <- dsp$Ytr; X_D1 <- dsp$Xtr; Z_D1 <- dsp$Ztr
  Y_D2 <- dsp$Yte; X_D2 <- dsp$Xte; Z_D2 <- dsp$Zte
  n <- length(Y_D2)

  ## ---- D1: Weight estimation (FIXED) ----

  ## Fit Y on Z and X on Z on D1
  wYZ <- fit_regression(y = Y_D1, x = Z_D1, method = reg_YonZ, args = args_YonZ)
  wrY <- stats::residuals(wYZ, response = Y_D1, data = Z_D1)

  ## Handle multivariate X: fit each column separately
  wrX <- apply(X_D1, 2, function(tX) {
    m <- fit_regression(y = tX, x = Z_D1, method = reg_XonZ, args = args_XonZ)
    stats::residuals(m, response = tX, data = Z_D1)
  })
  wrX <- as.matrix(wrX)

  ## Residual products and weight function
  RP <- c(wrY) * wrX
  mW <- apply(as.data.frame(RP), 2, function(tRP) {
    fit_regression(y = tRP, x = Z_D1, method = reg_wfun, args = args_wfun)
  })
  W <- do.call("cbind", lapply(mW, function(mZ) {
    sign(stats::predict(mZ, data = Z_D2))
  }))

  ## ---- D2: Permutation-calibrated GCM with weights ----

  ## Fit Y on Z on D2 -> rY (FIXED across permutations)
  mY <- fit_regression(y = Y_D2, x = Z_D2, method = reg_YonZ, args = args_YonZ)
  rY <- stats::residuals(mY, response = Y_D2, data = Z_D2)

  ## Fit X on Z on D2 -> rX
  rX <- apply(X_D2, 2, function(tX) {
    m <- fit_regression(y = tX, x = Z_D2, method = reg_XonZ, args = args_XonZ)
    stats::residuals(m, response = tX, data = Z_D2)
  })
  rX <- as.matrix(rX)

  ## Apply weights and compute observed statistic
  rX_w <- rX * W
  R_obs <- c(rY) * rX_w
  ## For multivariate X, aggregate across columns
  if (NCOL(R_obs) == 1) {
    T_obs <- studentized_stat(c(R_obs))
  } else {
    ## Use max-type statistic across X columns
    T_obs <- max(apply(R_obs, 2, studentized_stat))
  }

  ## ---- Permutation loop on D2 ----
  adaptive_info <- NULL

  ## Helper: compute one permutation statistic for WGCM
  compute_wgcm_perm_stat <- function(perm_idx) {
    X_perm <- X_D2[perm_idx, , drop = FALSE]
    rX_perm <- apply(X_perm, 2, function(tX) {
      m <- fit_regression(y = tX, x = Z_D2, method = reg_XonZ,
                          args = args_XonZ)
      stats::residuals(m, response = tX, data = Z_D2)
    })
    rX_perm <- as.matrix(rX_perm)
    rX_w_perm <- rX_perm * W
    R_perm <- c(rY) * rX_w_perm
    if (NCOL(R_perm) == 1) {
      studentized_stat(c(R_perm))
    } else {
      max(apply(R_perm, 2, studentized_stat))
    }
  }

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
      T_perm[b] <- compute_wgcm_perm_stat(perm)

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
    Q_D2 <- precompute_qr_projection(Z_D2)
    perm_mat <- replicate(B, sample.int(n))
    d <- NCOL(X_D2)

    if (d == 1L) {
      ## Univariate X
      X_perm_mat <- apply_permutations(c(X_D2), perm_mat)
      rX_perm_mat <- fast_lrm_residuals(Q_D2, X_perm_mat)
      rX_w_perm_mat <- rX_perm_mat * c(W)
      R_perm_mat <- c(rY) * rX_w_perm_mat
      T_perm <- vectorized_studentized_stat(R_perm_mat)
    } else {
      ## Multivariate X: max-type statistic across columns
      Tj_mat <- matrix(NA_real_, d, B)
      for (j in seq_len(d)) {
        Xj_perm_mat <- apply_permutations(X_D2[, j], perm_mat)
        rXj_perm_mat <- fast_lrm_residuals(Q_D2, Xj_perm_mat)
        rXj_w_perm_mat <- rXj_perm_mat * W[, j]
        Rj_perm_mat <- c(rY) * rXj_w_perm_mat
        Tj_mat[j, ] <- vectorized_studentized_stat(Rj_perm_mat)
      }
      T_perm <- apply(Tj_mat, 2, max)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else if (is_nslrm_fastpath(reg_XonZ, args_XonZ)) {
    ## --- Vectorized nslrm fast path (natural spline + OLS) ---
    ns_df <- if (!is.null(args_XonZ$df)) args_XonZ$df else 4L
    tmp_fit <- fit_nslrm(y = X_D2[, 1], x = Z_D2, df = ns_df)
    Q_D2 <- precompute_nslrm_projection(tmp_fit)
    perm_mat <- replicate(B, sample.int(n))
    d <- NCOL(X_D2)

    if (d == 1L) {
      X_perm_mat <- apply_permutations(c(X_D2), perm_mat)
      rX_perm_mat <- fast_lrm_residuals(Q_D2, X_perm_mat)
      rX_w_perm_mat <- rX_perm_mat * c(W)
      R_perm_mat <- c(rY) * rX_w_perm_mat
      T_perm <- vectorized_studentized_stat(R_perm_mat)
    } else {
      Tj_mat <- matrix(NA_real_, d, B)
      for (j in seq_len(d)) {
        Xj_perm_mat <- apply_permutations(X_D2[, j], perm_mat)
        rXj_perm_mat <- fast_lrm_residuals(Q_D2, Xj_perm_mat)
        rXj_w_perm_mat <- rXj_perm_mat * W[, j]
        Rj_perm_mat <- c(rY) * rXj_w_perm_mat
        Tj_mat[j, ] <- vectorized_studentized_stat(Rj_perm_mat)
      }
      T_perm <- apply(Tj_mat, 2, max)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    ## --- Parallel loop for non-lrm methods ---
    T_perm <- unlist(parallel::mclapply(seq_len(B), function(b) {
      compute_wgcm_perm_stat(sample.int(n))
    }, mc.cores = n_cores))
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else {
    ## --- Sequential loop (original) ---
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      T_perm[b] <- compute_wgcm_perm_stat(sample.int(n))
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  }

  ## Return htest object
  stat <- T_obs
  names(stat) <- "T"
  param <- c("B" = length(T_perm))

  models_out <- if (return_fitted_models) {
    list(reg_YonZ_weight = wYZ, reg_wfun = mW, reg_YonZ_test = mY)
  } else {
    NULL
  }

  method_str <- if (adaptive) {
    "Adaptive permutation-calibrated weighted generalised covariance measure test"
  } else {
    "Permutation-calibrated weighted generalised covariance measure test"
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
    W = W,
    perm_stats = T_perm,
    adaptive = adaptive_info,
    models = models_out
  ), class = c("perm_wgcm", "htest"))
}
