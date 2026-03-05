#' Permutation-Calibrated TRAM-GCM Invariance Test
#'
#' Tests whether score residuals and residualized environments are uncorrelated,
#' using a permutation-calibrated variant of the TRAM-GCM test from
#' Kook et al. (2024). This is the core invariance test used inside
#' \code{\link{perm_icp}}.
#'
#' @param score_resid Numeric vector of length n: score residuals from a fitted
#'   model (e.g., from \code{\link{compute_score_residuals}}).
#' @param resid_env Numeric vector or matrix: residualized environment(s)
#'   (e.g., from \code{\link{residualize_environment}}).
#' @param B Number of permutations (default 499).
#' @param alternative Type of alternative hypothesis. One of \code{"two.sided"}
#'   (default), \code{"greater"}, or \code{"less"}.
#'
#' @details
#' Implements Algorithm 2 of Kook et al. (2024) with permutation calibration:
#' \enumerate{
#'   \item Compute residual product terms: \eqn{L_i = R_i \cdot \xi_i} where
#'     \eqn{R_i} are score residuals and \eqn{\xi_i} are residualized environments.
#'   \item Estimate the covariance matrix \eqn{\hat{\Sigma}} of \eqn{L}.
#'   \item Compute test statistic: \eqn{T_n = \hat{\Sigma}^{-1/2} n^{-1/2} \sum L_i}.
#'   \item For univariate environments, use a studentized scalar statistic.
#'     For multivariate, use \eqn{\|T_n\|^2}.
#'   \item Permute the indices of \eqn{\xi} to break the association with R,
#'     recompute the statistic B times.
#'   \item Compute a permutation p-value.
#' }
#'
#' @return An object of class \code{c("perm_tram_gcm", "htest")} with components:
#' \describe{
#'   \item{statistic}{The observed test statistic.}
#'   \item{p.value}{Permutation p-value.}
#'   \item{parameter}{Number of permutations B.}
#'   \item{method}{Description of the test.}
#'   \item{perm_stats}{Vector of permutation statistics (length B).}
#' }
#'
#' @references
#' Kook, L., Saengkyongam, S., Lundborg, A. R., Hothorn, T., & Peters, J. (2024).
#' Model-based causal feature selection for general response types.
#' \emph{arXiv:2309.12833v4}.
#'
#' @keywords internal
perm_tram_gcm <- function(score_resid, resid_env,
                          B = 499L, alternative = "two.sided") {
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  B <- as.integer(B)
  n <- length(score_resid)
  R <- c(score_resid)
  xi <- as.matrix(resid_env)
  q <- ncol(xi)  # dimension of environment

  ## Compute residual product terms L_i = R_i * xi_i (n x q matrix)
  L <- R * xi  # broadcasting: n-vector * (n x q) matrix

  if (q == 1L) {
    ## --- Univariate environment: studentized scalar statistic ---
    T_obs <- studentized_stat(c(L))

    ## Permutation: permute xi indices
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      perm <- sample.int(n)
      L_perm <- R * xi[perm, 1L]
      T_perm[b] <- studentized_stat(L_perm)
    }

    pval <- perm_pvalue(T_obs, T_perm, alternative)
  } else {
    ## --- Multivariate environment: chi-squared-type statistic ---
    ## T_n = Sigma_hat^{-1/2} * n^{-1/2} * sum(L_i)
    ## Test statistic = ||T_n||^2
    L_bar <- colMeans(L)
    Sigma_hat <- crossprod(L) / n - tcrossprod(L_bar)

    ## Regularize if needed
    eig <- eigen(Sigma_hat, symmetric = TRUE)
    tol <- max(eig$values) * .Machine$double.eps * n
    keep <- eig$values > tol
    if (sum(keep) == 0L) {
      ## Degenerate case
      return(structure(list(
        statistic = c("T" = 0),
        p.value = 1,
        parameter = c("B" = B),
        method = "Permutation-calibrated TRAM-GCM invariance test",
        perm_stats = rep(0, B)
      ), class = c("perm_tram_gcm", "htest")))
    }

    Sigma_inv_sqrt <- eig$vectors[, keep, drop = FALSE] %*%
      diag(1 / sqrt(eig$values[keep]), nrow = sum(keep)) %*%
      t(eig$vectors[, keep, drop = FALSE])

    T_vec <- c(Sigma_inv_sqrt %*% (sqrt(n) * L_bar))
    T_obs <- sum(T_vec^2)

    ## Permutation loop
    T_perm <- numeric(B)
    for (b in seq_len(B)) {
      perm <- sample.int(n)
      L_p <- R * xi[perm, , drop = FALSE]
      L_bar_p <- colMeans(L_p)
      T_vec_p <- c(Sigma_inv_sqrt %*% (sqrt(n) * L_bar_p))
      T_perm[b] <- sum(T_vec_p^2)
    }

    ## For chi-squared-type statistic, always use "greater" for p-value
    pval <- (1 + sum(T_perm >= T_obs)) / (1 + B)
  }

  stat <- T_obs
  names(stat) <- "T"

  structure(list(
    statistic = stat,
    p.value = pval,
    parameter = c("B" = B),
    method = "Permutation-calibrated TRAM-GCM invariance test",
    data.name = "",
    alternative = alternative,
    perm_stats = T_perm
  ), class = c("perm_tram_gcm", "htest"))
}
