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
#' @param B Number of permutations (default 499).
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
#'   \item{perm_stats}{Vector of permutation test statistics (length B).}
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
                     alternative = c("two.sided", "greater", "less"),
                     return_fitted_models = FALSE) {
  alternative <- match.arg(alternative)
  B <- as.integer(B)

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
  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    ## Permute X_j within D1
    perm <- sample.int(n)
    Xj_perm <- Xj_D1[perm, , drop = FALSE]

    ## Form permuted full covariates on D1
    XW_D1_perm <- cbind(Xj_perm, W_D1)

    ## Recompute direction on D1 using FIXED D2 models
    fhat_perm <- compute_fhat(XW_D1_perm, W_D1)

    ## Refit projection: f_hat^pi ~ (X_{-j}, Z) on D1
    m_hat_pfj_perm <- fit_regression(y = fhat_perm, x = W_D1, method = reg_RonZ, args = args_RonZ)
    rF_perm <- c(fhat_perm) - c(stats::predict(m_hat_pfj_perm, data = W_D1))

    ## rY stays FIXED
    L_perm <- rY * rF_perm
    T_perm[b] <- studentized_stat(L_perm)
  }

  ## Permutation p-value
  pval <- perm_pvalue(T_obs, T_perm, alternative)

  ## Return htest object
  stat <- T_obs
  names(stat) <- "T"
  param <- c("B" = B)

  models <- if (return_fitted_models) {
    list(m_hat_p = m_hat_p, m_hat_qj = m_hat_qj,
         m_hat_rj = m_hat_rj, m_hat_pfj = m_hat_pfj)
  } else {
    NULL
  }

  structure(list(
    statistic = stat,
    p.value = pval,
    parameter = param,
    method = "Permutation-calibrated projected covariance measure test",
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    alternative = alternative,
    rY = rY,
    rF = rF,
    perm_stats = T_perm,
    D1_indices = idx_D1,
    D2_indices = idx_D2,
    models = models
  ), class = c("perm_pcm", "htest"))
}
