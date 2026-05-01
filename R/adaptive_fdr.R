## Adaptive permutation testing for FDR control
##
## Implements Algorithm 5 of Barry, Niu, Katsevich, and Lin (2025),
## "The permuted score test for robust differential expression analysis"
## (arXiv:2501.03530).

#' Benjamini-Hochberg rejection threshold
#'
#' Returns the p-value cutoff used by the BH procedure at level \code{alpha}:
#' the largest \eqn{k \alpha / m} such that the \eqn{k}-th order statistic of
#' \code{p} is at most \eqn{k \alpha / m}, or 0 if no such \eqn{k} exists.
#'
#' @param p Numeric vector of p-values.
#' @param alpha FDR level.
#' @return Scalar threshold (0 if no rejections).
#' @keywords internal
bh_threshold <- function(p, alpha) {
  m <- length(p)
  if (m == 0L) return(0)
  p_sorted <- sort(p)
  cuts <- seq_len(m) * alpha / m
  ks <- which(p_sorted <= cuts)
  if (length(ks) == 0L) 0 else max(ks) * alpha / m
}

#' Adaptive permutation test for FDR control (Algorithm 5)
#'
#' Generic implementation of the adaptive permutation FDR-control procedure of
#' Fischer, Barry, and Ramdas (2024), as described in Algorithm 5 of Barry,
#' Niu, Katsevich, and Lin (2025).  The procedure processes \eqn{m} hypotheses
#' simultaneously: each round it draws one permuted statistic per active
#' hypothesis, updates a per-hypothesis anytime-valid p-value
#' \eqn{p^i_t = h / (t + h - K^i_t)}, applies the BH rejection threshold to
#' the current p-values, and stops a hypothesis as soon as either (i) its
#' p-value falls at or below the BH threshold (rejection) or (ii) its loss
#' counter \eqn{K^i_t} hits \code{h} (futility).  Rejection and futility are
#' both sticky: once a hypothesis leaves the active set its p-value is frozen
#' at its last value.
#'
#' @param T_orig Numeric vector of length \code{m} of original test statistics.
#' @param permute_stat A function with signature
#'   \code{function(active_idx, t)} returning a numeric vector of length
#'   \code{length(active_idx)} containing the permuted statistics for each
#'   currently-active hypothesis at round \code{t}.  Implementations are
#'   encouraged to share a single permutation across hypotheses for speed.
#' @param alpha FDR level (default 0.05).
#' @param h Tuning parameter (positive integer).  Larger \code{h} gives more
#'   power but takes longer; the paper recommends 10--30.
#' @param B_max Cap on the number of permutation rounds (default 25000).  If
#'   the active set is non-empty at \code{t = B_max} the procedure stops; the
#'   remaining hypotheses are not rejected and their p-values are returned at
#'   their last anytime-valid value.
#' @param alternative One of \code{"greater"} (default), \code{"less"}, or
#'   \code{"two.sided"}.  Statistics are mapped to scores via
#'   \code{\link{score_map}} before comparison.
#'
#' @details
#' The anytime-valid p-value of Fischer and Ramdas (2024) is
#' \deqn{p^i_t = h / (t + h - K^i_t),}
#' where \eqn{K^i_t = \sum_{s=1}^{t} \mathbb{1}\{S^i_s \ge S^i_{\text{orig}}\}}
#' and \eqn{S} is the score returned by \code{score_map}.  When
#' \eqn{K^i_t = h}, the hypothesis is moved to the futility set and its
#' p-value is frozen at \eqn{h / t}.  Active-hypothesis p-values are
#' monotonically decreasing in \eqn{t}, so the BH threshold is monotonically
#' increasing and rejection is sticky.
#'
#' @return A list with components:
#' \describe{
#'   \item{rejected}{Integer vector of indices of rejected hypotheses.}
#'   \item{futility}{Integer vector of indices stopped for futility.}
#'   \item{active}{Integer vector of indices still active at termination
#'     (non-empty only if \code{B_max} was hit).}
#'   \item{p_value}{Numeric vector (length \code{m}) of final anytime-valid
#'     p-values.}
#'   \item{K}{Integer vector (length \code{m}) of final loss counts.}
#'   \item{stop_round}{Integer vector (length \code{m}) of the round at which
#'     each hypothesis was stopped (\code{NA} if still active).}
#'   \item{t_total}{Total number of rounds performed.}
#'   \item{R_threshold}{Final BH threshold.}
#' }
#'
#' @references
#' Fischer, L., Barry, T., and Ramdas, A. (2024). Multiple testing with
#' anytime-valid Monte Carlo p-values. arXiv:2404.15586.
#'
#' Barry, T., Niu, Z., Katsevich, E., and Lin, X. (2025). The permuted score
#' test for robust differential expression analysis.  arXiv:2501.03530.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- 50
#' T_orig <- c(rnorm(45), rnorm(5, mean = 4))   # 5 alternatives
#' permute_stat <- function(active, t) rnorm(length(active))
#' res <- adaptive_fdr_bh(T_orig, permute_stat, alpha = 0.1, h = 15L)
#' res$rejected
#' }
#'
#' @export
adaptive_fdr_bh <- function(T_orig, permute_stat,
                            alpha = 0.05,
                            h = 20L,
                            B_max = 25000L,
                            alternative = c("greater", "less", "two.sided")) {
  alternative <- match.arg(alternative)
  h <- as.integer(h)
  B_max <- as.integer(B_max)
  stopifnot(h >= 1L, alpha > 0, alpha < 1, B_max >= 1L)
  stopifnot(is.function(permute_stat))

  m <- length(T_orig)
  if (m == 0L) {
    return(list(rejected = integer(0), futility = integer(0),
                active = integer(0), p_value = numeric(0), K = integer(0),
                stop_round = integer(0), t_total = 0L, R_threshold = 0))
  }

  S0 <- score_map(T_orig, alternative)

  ## State (length m)
  K          <- integer(m)
  p_val      <- rep(1, m)            # frozen value before any losses
  status     <- rep(0L, m)           # 0 = active, 1 = futility, 2 = rejection
  stop_round <- rep(NA_integer_, m)

  active <- seq_len(m)
  t <- 0L
  R_threshold <- 0

  while (length(active) > 0L && t < B_max) {
    t <- t + 1L

    ## Draw one permuted statistic for each currently-active hypothesis.
    z_perm <- permute_stat(active, t)
    if (length(z_perm) != length(active)) {
      stop("permute_stat must return a vector of length length(active_idx).")
    }
    S_perm <- score_map(z_perm, alternative)

    ## Update K and the AV p-value for active hypotheses.
    is_loss <- S_perm >= S0[active]
    K[active] <- K[active] + as.integer(is_loss)

    ## Anytime-valid p-value for active hypotheses.
    p_val[active] <- h / (t + h - K[active])

    ## Move hypotheses with K == h into futility (frozen at h/t).
    futile_now <- active[K[active] >= h]
    if (length(futile_now) > 0L) {
      status[futile_now] <- 1L
      stop_round[futile_now] <- t
      ## p_val already equals h / (t + h - h) = h / t.  Freeze.
    }

    ## Drop futile hypotheses from the active set.
    active <- active[status[active] == 0L]

    ## BH threshold on the current p-values across all m hypotheses.
    R_threshold <- bh_threshold(p_val, alpha)

    ## Reject any active hypothesis whose p-value is at or below threshold.
    if (R_threshold > 0 && length(active) > 0L) {
      reject_now <- active[p_val[active] <= R_threshold]
      if (length(reject_now) > 0L) {
        status[reject_now] <- 2L
        stop_round[reject_now] <- t
        active <- setdiff(active, reject_now)
      }
    }
  }

  list(
    rejected    = which(status == 2L),
    futility    = which(status == 1L),
    active      = active,
    p_value     = p_val,
    K           = K,
    stop_round  = stop_round,
    t_total     = t,
    R_threshold = R_threshold
  )
}

#' Adaptive FDR-controlled permuted GCM test for many response columns
#'
#' Convenience wrapper that runs \code{\link{adaptive_fdr_bh}} for the common
#' setup in which a treatment vector \code{X} and nuisance covariates \code{Z}
#' are shared across \eqn{m} hypotheses, with one response column per
#' hypothesis.  For each column \eqn{j}, the test statistic is the studentized
#' GCM statistic \eqn{T_j = \sqrt{n}\,\bar{R}_j / \mathrm{sd}(R_j)} where
#' \eqn{R_j = \tilde r_X \cdot \tilde r_{Y_j}} are the residual products from
#' fitting \eqn{X \mid Z} and \eqn{Y_j \mid Z}.
#'
#' Permutations of \code{X} are shared across hypotheses for speed: at each
#' round, one permuted X-residual vector \eqn{\tilde r_X^{(t)}} is computed
#' (via the lrm/nslrm/gols fast path when applicable) and combined with the
#' fixed Y-residuals of all currently-active hypotheses in a single
#' matrix-multiply.
#'
#' @param Y Response matrix of dimension \eqn{n \times m}.  Each column is one
#'   hypothesis's response.
#' @param X Treatment vector or single-column matrix of length \eqn{n}, shared
#'   across hypotheses.
#' @param Z Nuisance covariate matrix \eqn{n \times p}, shared across
#'   hypotheses.
#' @param reg_YonZ,reg_XonZ Regression methods for \eqn{Y_j \mid Z} and
#'   \eqn{X \mid Z}.  See \code{\link{perm_gcm}}.
#' @param args_YonZ,args_XonZ Optional argument lists for the regressions.
#'   For \code{reg_XonZ = "gols"}, \code{theta} (NB dispersion) is required;
#'   the per-column \eqn{\mu_{Y_j}} are computed automatically.
#' @param alpha FDR level (default 0.05).
#' @param h Adaptive tuning parameter (default 20).
#' @param B_max Cap on the number of permutation rounds (default 25000).
#' @param alternative One of \code{"greater"} (default), \code{"less"},
#'   \code{"two.sided"}.
#'
#' @return A list inheriting all components of \code{\link{adaptive_fdr_bh}}
#'   plus:
#'   \describe{
#'     \item{T_orig}{Vector of length \code{m} of original studentized
#'       statistics.}
#'     \item{rY}{Matrix of \eqn{Y_j \mid Z} residuals (\eqn{n \times m}).}
#'     \item{rX}{Vector of \eqn{X \mid Z} residuals (length \eqn{n}).}
#'   }
#'
#' @export
perm_gcm_fdr <- function(Y, X, Z,
                         reg_YonZ = "lrm", reg_XonZ = "lrm",
                         args_YonZ = NULL, args_XonZ = NULL,
                         alpha = 0.05,
                         h = 20L,
                         B_max = 25000L,
                         alternative = c("greater", "less", "two.sided")) {
  alternative <- match.arg(alternative)

  Y <- as.matrix(Y)
  X <- check_data(X, "X")
  Z <- check_data(Z, "Z")
  n <- NROW(Y)
  m <- NCOL(Y)
  if (NROW(X) != n) stop("X and Y must have the same number of rows.")
  if (NROW(Z) != n) stop("Z and Y must have the same number of rows.")
  if (NCOL(X) != 1L) stop("perm_gcm_fdr currently supports a univariate X.")

  ## Per-hypothesis Y_j | Z fits (only the residual matrix is retained).
  ## We also keep mu_Y if reg_XonZ == "gols", as gols needs it per column;
  ## but the X|Z regression is shared across hypotheses, so a single mu_Y
  ## is not directly applicable when Y varies per column.  Instead, we use
  ## the *average* mu_Y across columns (or simply require the user to pass
  ## a fixed theta and a representative mu_Y via args_XonZ if they want
  ## gols).  For correctness, we do not enable per-column gols here.
  rY <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m)) {
    Yj <- Y[, j, drop = FALSE]
    mYj <- fit_regression(y = Yj, x = Z, method = reg_YonZ, args = args_YonZ)
    rY[, j] <- stats::residuals(mYj, response = Yj, data = Z)
  }

  ## X | Z fit (shared across hypotheses).
  args_X <- args_XonZ
  if (identical(reg_XonZ, "gols")) {
    if (is.null(args_X[["theta"]])) {
      stop("reg_XonZ = \"gols\" requires args_XonZ$theta ",
           "(NB dispersion parameter).")
    }
    if (is.null(args_X[["mu_Y"]])) {
      ## Use the row-wise mean of fitted Y values across hypotheses as a
      ## shared mu_Y for the working weights.
      mu_Y_shared <- rowMeans(Y - rY)
      args_X[["mu_Y"]] <- mu_Y_shared
    }
  }
  mX <- fit_regression(y = X, x = Z, method = reg_XonZ, args = args_X)
  rX <- stats::residuals(mX, response = X, data = Z)

  ## Original studentized statistics across hypotheses.
  R_obs <- as.numeric(rX) * rY                # n x m
  T_orig <- vectorized_studentized_stat(R_obs)

  ## Determine the fast path for permuted X residuals.
  fast_kind <- if (is_lrm_fastpath(reg_XonZ, args_X)) {
    "lrm"
  } else if (is_nslrm_fastpath(reg_XonZ, args_X)) {
    "nslrm"
  } else if (is_gols_fastpath(reg_XonZ, args_X)) {
    "gols"
  } else {
    "slow"
  }
  proj <- switch(fast_kind,
    "lrm"   = precompute_qr_projection(Z),
    "nslrm" = precompute_nslrm_projection(mX),
    "gols"  = precompute_gols_projection(mX),
    NULL
  )
  X_vec <- as.numeric(X)

  ## Closure: permute X once, compute permuted X-residuals once, then
  ## combine with rY[, active] in a single matrix multiply.
  permute_stat <- function(active_idx, t) {
    perm <- sample.int(n)
    X_perm <- X_vec[perm]
    if (fast_kind == "lrm" || fast_kind == "nslrm") {
      rX_perm <- as.numeric(fast_lrm_residuals(proj, matrix(X_perm, ncol = 1L)))
    } else if (fast_kind == "gols") {
      rX_perm <- as.numeric(fast_gols_residuals(proj, matrix(X_perm, ncol = 1L)))
    } else {
      Xp <- matrix(X_perm, ncol = 1L)
      args_p <- args_X
      if (identical(reg_XonZ, "gols")) {
        ## mu_Y already injected into args_X above.
      }
      mXp <- fit_regression(y = Xp, x = Z, method = reg_XonZ, args = args_p)
      rX_perm <- as.numeric(stats::residuals(mXp, response = Xp, data = Z))
    }
    R_mat <- rX_perm * rY[, active_idx, drop = FALSE]
    vectorized_studentized_stat(R_mat)
  }

  res <- adaptive_fdr_bh(T_orig, permute_stat,
                         alpha = alpha, h = h, B_max = B_max,
                         alternative = alternative)
  res$T_orig <- T_orig
  res$rY <- rY
  res$rX <- rX
  res
}
