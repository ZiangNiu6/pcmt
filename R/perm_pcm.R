#' Permutation-Calibrated Projected Covariance Measure Test
#'
#' Tests the conditional independence hypothesis Y indep X_j given the
#' remaining covariates, using a permutation-calibrated variant of the vanilla
#' PCM (vPCM) test. The test is exact under the strong exchangeable null
#' and asymptotically valid under the weaker conditional null (under
#' permutation stability conditions).
#'
#' @param Y Response variable. Numeric vector.
#' @param X Covariate matrix. The \code{target} column is the variable of interest.
#' @param Z Optional additional conditioning variables. Numeric vector or matrix.
#'   If NULL, tests Y indep X_j given X_rest.
#' @param target Integer index of the target variable in X (default 1).
#' @param reg_YonXZ Regression method for Y on (X, Z), the full model on D2.
#' @param reg_YonZ Regression method for Y on the conditioning set, reduced model on D1.
#' @param reg_YhatonZ Regression method for projecting full-model predictions
#'   onto the conditioning set on D2.
#' @param reg_RonZ Regression method for projecting the direction onto the
#'   conditioning set on D1.
#' @param args_YonXZ,args_YonZ,args_YhatonZ,args_RonZ Optional lists of
#'   additional arguments for each regression.
#' @param frac Fraction of data for D1 (test fold). Default 0.5.
#' @param indices Optional pre-specified indices for D2 (direction estimation fold).
#' @param B Number of permutations (default 499). When \code{h} is specified,
#'   this serves as the deterministic cap \eqn{T_{max}} on the number of
#'   sequential permutations.
#' @param n_cores Number of cores for parallel permutation loop (default 1,
#'   sequential). Only used when no lrm fast path is applicable and
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
#' The algorithm (from the permuted PCM mini-paper, Algorithm 1) proceeds:
#' \enumerate{
#'   \item Split data into D1 (test, size n) and D2 (direction, size n).
#'   \item On D2 (fitted once, kept fixed): fit full model and projection to
#'     estimate the direction function.
#'   \item On D1: fit reduced model for Y residuals, compute direction values,
#'     fit projection of direction, and form the residual product statistic.
#'   \item Permutation loop: permute X_j within D1, recompute direction
#'     values using fixed D2 models, refit projection on D1, compute
#'     permutation statistics.
#' }
#'
#' @return An object of class \code{c("perm_pcm", "htest")} with components:
#' \describe{
#'   \item{statistic}{The observed studentized test statistic.}
#'   \item{p.value}{Permutation p-value.}
#'   \item{parameter}{Number of permutations B.}
#'   \item{method}{Description of the test method.}
#'   \item{data.name}{Description of the data.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{rY}{Residuals of Y from reduced model on D1.}
#'   \item{rF}{Residuals of the direction from projection on D1.}
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
#' n <- 400
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- X[, 2] + X[, 3] + rnorm(n)  # Y _||_ X_1 | X_{-1}
#' result <- perm_pcm(Y, X, target = 1, reg_YonXZ = "lrm",
#'   reg_YonZ = "lrm", reg_YhatonZ = "lrm", reg_RonZ = "lrm", B = 199)
#' print(result)
#' }
#'
#' @export
perm_pcm <- function(Y, X, Z = NULL,
                     target = 1L,
                     reg_YonXZ = "rf", reg_YonZ = "rf",
                     reg_YhatonZ = "rf", reg_RonZ = "rf",
                     args_YonXZ = NULL, args_YonZ = NULL,
                     args_YhatonZ = NULL, args_RonZ = NULL,
                     frac = 0.5, indices = NULL,
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
  Y <- check_data(Y, "Y", "pcm")
  X <- check_data(X, "X", "pcm")
  n_total <- NROW(X)
  target <- as.integer(target)

  ## Separate X into X_j (target) and X_{-j}
  Xj_all <- X[, target, drop = FALSE]
  Xmj_all <- X[, -target, drop = FALSE]

  ## Form the conditioning set: (X_{-j}, Z) or just X_{-j}
  if (!is.null(Z)) {
    Z <- check_data(Z, "Z", "pcm")
    W_all <- cbind(Xmj_all, Z)  # full conditioning set
  } else {
    W_all <- Xmj_all
  }

  ## Full covariate set for full model: (X_j, X_{-j}, Z) = (X_j, W)
  XW_all <- cbind(Xj_all, W_all)

  ## Sample splitting: D2 = direction (indices), D1 = test (complement)
  ## Following the paper: |D1| = |D2| = n
  idx_D2 <- indices
  if (is.null(idx_D2)) {
    idx_D2 <- sample.int(n_total, ceiling(frac * n_total))
  }
  idx_D1 <- setdiff(seq_len(n_total), idx_D2)
  n <- length(idx_D1)

  Y_D1 <- Y[idx_D1]
  Y_D2 <- Y[idx_D2]
  Xj_D1 <- Xj_all[idx_D1, , drop = FALSE]
  Xmj_D1 <- Xmj_all[idx_D1, , drop = FALSE]
  W_D1 <- W_all[idx_D1, , drop = FALSE]
  XW_D1 <- XW_all[idx_D1, , drop = FALSE]
  W_D2 <- W_all[idx_D2, , drop = FALSE]
  XW_D2 <- XW_all[idx_D2, , drop = FALSE]

  ## ---- D2: Direction estimation (fitted ONCE, kept FIXED) ----

  ## Fit m_hat_p: Y ~ (X_j, X_{-j}, Z) on D2
  m_hat_p <- fit_regression(y = Y_D2, x = XW_D2, method = reg_YonXZ, args = args_YonXZ)
  pghat_D2 <- stats::predict(m_hat_p, data = XW_D2)

  ## Fit m_hat_qj: m_hat_p(X) ~ (X_{-j}, Z) on D2
  m_hat_qj <- fit_regression(y = pghat_D2, x = W_D2, method = reg_YhatonZ, args = args_YhatonZ)

  ## Direction function: f_hat(x) = m_hat_p(x) - m_hat_qj(x_{-j})
  ## This function takes the full covariate matrix and the conditioning set
  compute_fhat <- function(XW, W) {
    stats::predict(m_hat_p, data = XW) - stats::predict(m_hat_qj, data = W)
  }

  ## ---- D1: Testing ----

  ## Fit m_hat_rj: Y ~ (X_{-j}, Z) on D1 (FIXED, doesn't depend on X_j)
  m_hat_rj <- fit_regression(y = Y_D1, x = W_D1, method = reg_YonZ, args = args_YonZ)
  rY <- c(Y_D1) - c(stats::predict(m_hat_rj, data = W_D1))

  ## Compute f_hat on D1 (original data)
  fhat_D1 <- compute_fhat(XW_D1, W_D1)

  ## Fit m_hat_pfj: f_hat(X) ~ (X_{-j}, Z) on D1
  m_hat_pfj <- fit_regression(y = fhat_D1, x = W_D1, method = reg_RonZ, args = args_RonZ)
  rF <- c(fhat_D1) - c(stats::predict(m_hat_pfj, data = W_D1))

  ## Observed test statistic
  L_obs <- rY * rF
  T_obs <- studentized_stat(L_obs)

  ## ---- Permutation loop ----
  adaptive_info <- NULL
  lrm_full <- is_lrm_fastpath(reg_YonXZ, args_YonXZ) &&
              is_lrm_fastpath(reg_YhatonZ, args_YhatonZ) &&
              is_lrm_fastpath(reg_RonZ, args_RonZ)
  lrm_proj <- is_lrm_fastpath(reg_RonZ, args_RonZ)
  nslrm_proj <- is_nslrm_fastpath(reg_RonZ, args_RonZ)

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
      Xj_perm <- Xj_D1[perm, , drop = FALSE]
      XW_D1_perm <- cbind(Xj_perm, W_D1)
      fhat_perm <- compute_fhat(XW_D1_perm, W_D1)
      m_hat_pfj_perm <- fit_regression(y = fhat_perm, x = W_D1,
                                       method = reg_RonZ, args = args_RonZ)
      rF_perm <- c(fhat_perm) - c(stats::predict(m_hat_pfj_perm, data = W_D1))
      L_perm <- rY * rF_perm
      T_perm[b] <- studentized_stat(L_perm)

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
  } else if (lrm_full) {
    ## --- Full lrm fast path ---
    ## Decompose fhat_perm = beta_j * Xj_perm + const_fhat
    cfx_p <- stats::coef(m_hat_p)
    cfx_p[is.na(cfx_p)] <- 0
    beta_j <- cfx_p[2]  # coefficient for Xj (1st predictor after intercept)

    ## Constant part: predict(m_hat_p, Xj=0, W) - predict(m_hat_qj, W)
    cfx_W <- cfx_p[-(1:2)]
    if (length(cfx_W) > 0L) {
      pred_W_const <- c(cfx_p[1] + as.matrix(W_D1) %*% cfx_W)
    } else {
      pred_W_const <- rep(cfx_p[1], n)
    }
    qj_pred <- c(stats::predict(m_hat_qj, data = W_D1))
    const_fhat <- pred_W_const - qj_pred

    ## Precompute projection matrix Q_W and fixed residual component
    Q_W <- precompute_qr_projection(W_D1)
    M_W_const <- c(fast_lrm_residuals(Q_W, matrix(const_fhat, ncol = 1)))

    ## All permutations at once
    perm_mat <- replicate(B, sample.int(n))
    Xj_perm_mat <- apply_permutations(c(Xj_D1), perm_mat)

    ## rF_perm = beta_j * M_W %*% Xj_perm + M_W_const
    M_W_Xj_perm <- fast_lrm_residuals(Q_W, Xj_perm_mat)
    rF_perm_mat <- beta_j * M_W_Xj_perm + M_W_const

    L_perm_mat <- c(rY) * rF_perm_mat
    T_perm <- vectorized_studentized_stat(L_perm_mat)
    pval <- perm_pvalue(T_obs, T_perm, alternative)

  } else if (lrm_proj) {
    ## --- Partial fast path: loop over compute_fhat, batch projection ---
    Q_W <- precompute_qr_projection(W_D1)
    perm_mat <- replicate(B, sample.int(n))
    fhat_perm_mat <- matrix(NA_real_, n, B)

    for (b in seq_len(B)) {
      Xj_perm <- Xj_D1[perm_mat[, b], , drop = FALSE]
      XW_D1_perm <- cbind(Xj_perm, W_D1)
      fhat_perm_mat[, b] <- compute_fhat(XW_D1_perm, W_D1)
    }

    rF_perm_mat <- fast_lrm_residuals(Q_W, fhat_perm_mat)
    L_perm_mat <- c(rY) * rF_perm_mat
    T_perm <- vectorized_studentized_stat(L_perm_mat)
    pval <- perm_pvalue(T_obs, T_perm, alternative)

  } else if (nslrm_proj) {
    ## --- Partial nslrm fast path: loop over compute_fhat, batch projection ---
    ns_df <- if (!is.null(args_RonZ$df)) args_RonZ$df else 4L
    tmp_fit <- fit_nslrm(y = fhat_D1, x = W_D1, df = ns_df)
    Q_W <- precompute_nslrm_projection(tmp_fit)
    perm_mat <- replicate(B, sample.int(n))
    fhat_perm_mat <- matrix(NA_real_, n, B)

    for (b in seq_len(B)) {
      Xj_perm <- Xj_D1[perm_mat[, b], , drop = FALSE]
      XW_D1_perm <- cbind(Xj_perm, W_D1)
      fhat_perm_mat[, b] <- compute_fhat(XW_D1_perm, W_D1)
    }

    rF_perm_mat <- fast_lrm_residuals(Q_W, fhat_perm_mat)
    L_perm_mat <- c(rY) * rF_perm_mat
    T_perm <- vectorized_studentized_stat(L_perm_mat)
    pval <- perm_pvalue(T_obs, T_perm, alternative)

  } else if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    ## --- Parallel loop for non-lrm methods ---
    T_perm <- unlist(parallel::mclapply(seq_len(B), function(b) {
      perm <- sample.int(n)
      Xj_perm <- Xj_D1[perm, , drop = FALSE]
      XW_D1_perm <- cbind(Xj_perm, W_D1)
      fhat_perm <- compute_fhat(XW_D1_perm, W_D1)
      m_pfj_perm <- fit_regression(y = fhat_perm, x = W_D1,
                                   method = reg_RonZ, args = args_RonZ)
      rF_perm <- c(fhat_perm) - c(stats::predict(m_pfj_perm, data = W_D1))
      L_perm <- rY * rF_perm
      studentized_stat(L_perm)
    }, mc.cores = n_cores))
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else {
    ## --- Sequential loop (original) ---
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      perm <- sample.int(n)
      Xj_perm <- Xj_D1[perm, , drop = FALSE]
      XW_D1_perm <- cbind(Xj_perm, W_D1)
      fhat_perm <- compute_fhat(XW_D1_perm, W_D1)
      m_hat_pfj_perm <- fit_regression(y = fhat_perm, x = W_D1,
                                       method = reg_RonZ, args = args_RonZ)
      rF_perm <- c(fhat_perm) - c(stats::predict(m_hat_pfj_perm, data = W_D1))
      L_perm <- rY * rF_perm
      T_perm[b] <- studentized_stat(L_perm)
    }
    pval <- perm_pvalue(T_obs, T_perm, alternative)
  }

  ## Return htest object
  stat <- T_obs
  names(stat) <- "T"
  param <- c("B" = length(T_perm))

  models <- if (return_fitted_models) {
    list(m_hat_p = m_hat_p, m_hat_qj = m_hat_qj,
         m_hat_rj = m_hat_rj, m_hat_pfj = m_hat_pfj)
  } else {
    NULL
  }

  method_str <- if (adaptive) {
    "Adaptive permutation-calibrated projected covariance measure test"
  } else {
    "Permutation-calibrated projected covariance measure test"
  }

  structure(list(
    statistic = stat,
    p.value = pval,
    parameter = param,
    method = method_str,
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    alternative = alternative,
    rY = rY,
    rF = rF,
    perm_stats = T_perm,
    D1_indices = idx_D1,
    D2_indices = idx_D2,
    adaptive = adaptive_info,
    models = models
  ), class = c("perm_pcm", "htest"))
}
