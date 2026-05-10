## perm_pcm_fdr.R
##
## Adaptive FDR-controlled permuted PCM test (and a matching asymptotic
## variant) built on top of `adaptive_fdr_bh`.  The cross-pair speedup
## relies on the test statistic being linear in the permuted X column on
## D1, which holds when reg_YonXZ in {"lrm", "gols"} and reg_YhatonZ /
## reg_RonZ are "lrm".  Under those choices the only X-dependent term in
## fhat is beta_X * X_j, so a single shared sigma (over D1) updates all
## active hypotheses without refitting any regression.

.PCM_LIN_REGS <- c("lrm", "gols")

## Internal: extract a length-(1 + ncol(x)) coefficient vector from a
## fitted lrm or gols object.  Coefficient layout is [intercept, x1, ...].
.pcm_extract_coefs <- function(fit, y) {
  if (inherits(fit, "lrm")) {
    cfx <- stats::coef(fit)
    cfx[is.na(cfx)] <- 0
    return(as.numeric(cfx))
  }
  if (inherits(fit, "gols_fit")) {
    Z_d <- fit$Z_d
    ZtWy <- crossprod(fit$W_hat * Z_d, as.matrix(y))
    coefs <- backsolve(fit$R_chol,
                       forwardsolve(t(fit$R_chol), ZtWy))
    return(as.numeric(coefs))
  }
  stop("Cannot extract coefficients from regression of class: ",
       paste(class(fit), collapse = "/"))
}

#' Precompute a per-pair PCM cache for fast permutation
#'
#' Builds the per-hypothesis cache used by \code{\link{perm_pcm_fdr}} and
#' \code{\link{asym_pcm}}.  The cache stores everything needed to update
#' the PCM test statistic under permutation of \code{X[, target]} on D1
#' without refitting any regression: an \code{lrm}/\code{gols} fit of
#' \eqn{Y \mid X, Z} on D2 is linear in X, so under permutation the only
#' X-dependent term in \eqn{\hat f} is \eqn{\beta_{\!X} \cdot X_j}.
#'
#' Multiple caches can share an \code{idx_D2} so a single permutation of
#' \eqn{\{1, \ldots, n_{D1}\}} updates all of them in lockstep -- this is
#' what \code{perm_pcm_fdr} exploits across response columns.  Callers
#' building per-pair pipelines (where \code{Y} and \code{X} both vary per
#' pair) supply \code{idx_D2} themselves to keep the splits aligned.
#'
#' @param Y Numeric response vector of length \eqn{n}.
#' @param X Matrix or data frame of covariates (\eqn{n \times d}).  The
#'   \code{target} column is the variable of interest.
#' @param Z Optional additional conditioning matrix (\eqn{n \times p}).  If
#'   \code{NULL}, the conditioning set is \eqn{X_{-\text{target}}} alone.
#' @param target Integer index into the columns of \code{X} (default 1).
#' @param reg_YonXZ Regression for \eqn{Y \mid X, Z} on D2.  Must be one of
#'   \code{"lrm"} or \code{"gols"} -- the cache assumes linearity in
#'   \code{X[, target]}.
#' @param reg_YonZ Regression for \eqn{Y \mid Z} on D1 (any method with
#'   proper \code{residuals} dispatch, e.g. \code{"lrm"}, \code{"glrm"}).
#' @param reg_YhatonZ Regression for the Yhat-on-Z projection on D2.  Must
#'   be \code{"lrm"} (we extract closed-form coefficients).
#' @param reg_RonZ Regression for the direction-on-Z projection on D1.
#'   Must be \code{"lrm"} (uses the QR fast path).
#' @param args_YonXZ,args_YonZ,args_YhatonZ,args_RonZ Optional argument
#'   lists.  For \code{reg_YonXZ = "gols"}, \code{args_YonXZ$theta} is
#'   required and \code{args_YonXZ$mu_Y} is auto-extracted from a fresh
#'   \code{reg_YonZ} fit on D2 if not supplied.
#' @param idx_D2 Integer vector of D2 (direction-estimation half) indices.
#'   If \code{NULL}, drawn at random with size \code{ceiling(frac * n)}.
#' @param frac Fraction of \eqn{n} used for D2 (default 0.5).  Ignored
#'   when \code{idx_D2} is supplied.
#'
#' @return A list with components used by \code{pcm_permute_stat_active}:
#' \describe{
#'   \item{idx_D1, idx_D2}{Integer index vectors for the two halves.}
#'   \item{n_D1}{\code{length(idx_D1)}.}
#'   \item{Xj_D1}{Length-\eqn{n_{D1}} numeric vector of \code{X[, target]}
#'     restricted to D1.}
#'   \item{rY_D1}{Numeric vector of \eqn{Y \mid Z} residuals on D1.}
#'   \item{beta_X_p}{Coefficient on \code{X[, target]} from the D2 fit.}
#'   \item{fhat_const}{Length-\eqn{n_{D1}} vector of the X-independent part
#'     of \eqn{\hat f} on D1: \eqn{\hat g_{-X}(W) - \hat q(W)}.}
#'   \item{Q_W}{Q factor from \code{precompute_qr_projection(W_D1)} for the
#'     residualization of \eqn{\hat f^\sigma} on D1.}
#'   \item{T_obs}{Studentized statistic at the unpermuted X.}
#' }
#' @export
precompute_pcm_pair <- function(Y, X, Z = NULL, target = 1L,
                                reg_YonXZ = "lrm", reg_YonZ = "lrm",
                                reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                                args_YonXZ = NULL, args_YonZ = NULL,
                                args_YhatonZ = NULL, args_RonZ = NULL,
                                idx_D2 = NULL, frac = 0.5) {
  if (!reg_YonXZ %in% .PCM_LIN_REGS) {
    stop("perm_pcm_fdr requires reg_YonXZ in {",
         paste(.PCM_LIN_REGS, collapse = ", "), "}; got \"",
         reg_YonXZ, "\".")
  }
  if (!identical(reg_YhatonZ, "lrm")) {
    stop("perm_pcm_fdr requires reg_YhatonZ = \"lrm\"; got \"",
         reg_YhatonZ, "\".")
  }
  if (!identical(reg_RonZ, "lrm")) {
    stop("perm_pcm_fdr requires reg_RonZ = \"lrm\"; got \"",
         reg_RonZ, "\".")
  }

  Y <- check_data(Y, "Y", "pcm")
  X <- check_data(X, "X", "pcm")
  target <- as.integer(target)
  if (target < 1L || target > NCOL(X)) {
    stop("`target` (", target, ") is out of range [1, ", NCOL(X), "].")
  }
  n <- NROW(Y)
  if (NROW(X) != n) stop("X and Y must have the same number of rows.")

  Xj_all  <- X[, target, drop = FALSE]
  Xmj_all <- X[, -target, drop = FALSE]
  if (!is.null(Z)) {
    Z <- check_data(Z, "Z", "pcm")
    if (NROW(Z) != n) stop("Z and Y must have the same number of rows.")
    W_all <- cbind(Xmj_all, Z)
  } else {
    W_all <- Xmj_all
  }
  XW_all <- cbind(Xj_all, W_all)

  if (is.null(idx_D2)) {
    idx_D2 <- sample.int(n, ceiling(frac * n))
  } else {
    idx_D2 <- as.integer(idx_D2)
    if (!all(idx_D2 %in% seq_len(n))) {
      stop("`idx_D2` contains out-of-range indices.")
    }
  }
  idx_D1 <- setdiff(seq_len(n), idx_D2)
  n_D1   <- length(idx_D1)

  Y_D2  <- Y[idx_D2];  Y_D1  <- Y[idx_D1]
  Xj_D2 <- Xj_all[idx_D2, , drop = FALSE]
  Xj_D1 <- Xj_all[idx_D1, , drop = FALSE]
  W_D2  <- W_all[idx_D2, , drop = FALSE]
  W_D1  <- W_all[idx_D1, , drop = FALSE]
  XW_D2 <- XW_all[idx_D2, , drop = FALSE]
  XW_D1 <- XW_all[idx_D1, , drop = FALSE]

  ## ---- D1: m_hat_rj = Y | Z fit (any reg_YonZ) ----
  m_hat_rj <- fit_regression(y = Y_D1, x = W_D1,
                             method = reg_YonZ, args = args_YonZ)
  rY_D1 <- as.numeric(stats::residuals(m_hat_rj,
                                       response = Y_D1, data = W_D1))

  ## ---- gols: auto-extract mu_Y on D2 from a fresh reg_YonZ fit ----
  args_YonXZ_use <- args_YonXZ
  if (identical(reg_YonXZ, "gols")) {
    if (is.null(args_YonXZ_use[["theta"]])) {
      stop("reg_YonXZ = \"gols\" requires args_YonXZ$theta ",
           "(NB dispersion parameter).")
    }
    if (is.null(args_YonXZ_use[["mu_Y"]])) {
      m_hat_rj_D2 <- fit_regression(y = Y_D2, x = W_D2,
                                    method = reg_YonZ, args = args_YonZ)
      r_D2 <- as.numeric(stats::residuals(m_hat_rj_D2,
                                          response = Y_D2, data = W_D2))
      args_YonXZ_use[["mu_Y"]] <- as.numeric(Y_D2) - r_D2
    }
  }

  ## ---- D2: m_hat_p = Y | (X, Z) fit (lrm or gols) ----
  m_hat_p <- fit_regression(y = Y_D2, x = XW_D2,
                            method = reg_YonXZ, args = args_YonXZ_use)
  coefs_p <- .pcm_extract_coefs(m_hat_p, Y_D2)
  ## coefs_p layout: [intercept, beta_target, coefs over W_D2 columns]
  beta_X_p <- coefs_p[2]

  ## ---- D2: m_hat_qj = Yhat | Z fit (lrm) ----
  pghat_D2 <- as.numeric(cbind(1, XW_D2) %*% coefs_p)
  m_hat_qj <- fit_regression(y = pghat_D2, x = W_D2,
                             method = reg_YhatonZ, args = args_YhatonZ)
  coefs_qj <- .pcm_extract_coefs(m_hat_qj, pghat_D2)

  ## ---- D1: split fhat = beta_X * Xj + (X-independent constant) ----
  ## pghat on D1 with the original Xj_D1, then subtract beta_X * Xj_D1 to
  ## isolate the X-independent contribution.
  pghat_D1   <- as.numeric(cbind(1, XW_D1) %*% coefs_p)
  pghat_D1_0 <- pghat_D1 - beta_X_p * as.numeric(Xj_D1[, 1])
  qjhat_D1   <- as.numeric(coefs_qj[1] + W_D1 %*% coefs_qj[-1])
  fhat_const <- pghat_D1_0 - qjhat_D1

  ## ---- D1: precompute Q_W for fast lrm residualization on W_D1 ----
  Q_W <- precompute_qr_projection(W_D1)

  ## ---- T_obs at the unpermuted Xj_D1 ----
  fhat_orig <- beta_X_p * as.numeric(Xj_D1[, 1]) + fhat_const
  rF_orig <- as.numeric(fast_lrm_residuals(Q_W,
                                            matrix(fhat_orig, ncol = 1L)))
  T_obs <- studentized_stat(rY_D1 * rF_orig)

  list(
    idx_D1     = idx_D1,
    idx_D2     = idx_D2,
    n_D1       = n_D1,
    Xj_D1      = as.numeric(Xj_D1[, 1]),
    rY_D1      = rY_D1,
    beta_X_p   = beta_X_p,
    fhat_const = fhat_const,
    Q_W        = Q_W,
    T_obs      = T_obs
  )
}

#' Permutation statistic on a list of pair caches sharing a sigma
#'
#' Given a list of \code{\link{precompute_pcm_pair}} caches that all share
#' the same D1 sample-split size, applies a single permutation \code{sigma}
#' (a permutation of \code{seq_len(n_D1)}) to every cache and returns the
#' updated studentized statistic for each.
#'
#' Callers writing custom multi-pair pipelines (where \code{X} varies per
#' pair, e.g. Perturb-seq pair construction) compose this with
#' \code{\link{adaptive_fdr_bh}} via a closure that draws a fresh sigma per
#' round and indexes into their cache list with \code{active_idx}.
#'
#' @param caches A list of pair caches.
#' @param sigma Integer permutation of \code{seq_len(n_D1)}, where
#'   \code{n_D1} matches every cache.
#' @return Numeric vector of length \code{length(caches)}.
#' @export
pcm_permute_stat_active <- function(caches, sigma) {
  vapply(caches, function(pc) {
    Xj_perm  <- pc$Xj_D1[sigma]
    fhat_perm <- pc$beta_X_p * Xj_perm + pc$fhat_const
    rF_perm  <- as.numeric(fast_lrm_residuals(pc$Q_W,
                                               matrix(fhat_perm, ncol = 1L)))
    studentized_stat(pc$rY_D1 * rF_perm)
  }, numeric(1))
}

#' Adaptive FDR-controlled permuted PCM test for many response columns
#'
#' Convenience wrapper that runs \code{\link{adaptive_fdr_bh}} for the
#' common multi-hypothesis setup in which a covariate matrix \code{X} (with
#' a single \code{target} column) and conditioning matrix \code{Z} are
#' shared across \eqn{m} hypotheses, with one response column per
#' hypothesis.  Per-pair PCM caches are built once with a shared D2
#' sample-split, then a single permutation is drawn per round and applied
#' to all currently-active hypotheses.
#'
#' Regression restrictions: \code{reg_YonXZ} must be \code{"lrm"} or
#' \code{"gols"}, and both projection regressions must be \code{"lrm"}.
#' These enable the linear-in-X structure that makes per-permutation
#' updates O(n_D1) per active hypothesis with no regression refit.
#'
#' @param Y Response matrix of dimension \eqn{n \times m}.
#' @param X Covariate matrix of dimension \eqn{n \times d}.
#' @param Z Optional additional conditioning matrix \eqn{n \times p}.
#' @param target Integer index of the target column in \code{X} (default 1).
#' @param reg_YonXZ,reg_YonZ,reg_YhatonZ,reg_RonZ See
#'   \code{\link{precompute_pcm_pair}}.
#' @param args_YonXZ,args_YonZ,args_YhatonZ,args_RonZ Optional argument
#'   lists.  For \code{reg_YonXZ = "gols"}, \code{args_YonXZ$theta} is
#'   required (shared across \eqn{m}); \code{mu_Y} is auto-extracted per
#'   column.
#' @param frac Fraction of rows for D2 (default 0.5).
#' @param alpha FDR level (default 0.05).
#' @param h Adaptive tuning parameter (default 20).
#' @param B_max Cap on permutation rounds (default 25000).
#' @param alternative One of \code{"greater"} (default), \code{"less"},
#'   \code{"two.sided"}.
#'
#' @return A list inheriting all components of
#'   \code{\link{adaptive_fdr_bh}} plus:
#'   \describe{
#'     \item{T_orig}{Vector of length \eqn{m} of studentized statistics.}
#'     \item{idx_D2}{Integer index vector of the shared D2 split.}
#'   }
#' @export
perm_pcm_fdr <- function(Y, X, Z = NULL, target = 1L,
                         reg_YonXZ = "lrm", reg_YonZ = "lrm",
                         reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                         args_YonXZ = NULL, args_YonZ = NULL,
                         args_YhatonZ = NULL, args_RonZ = NULL,
                         frac = 0.5,
                         alpha = 0.05, h = 20L, B_max = 25000L,
                         alternative = c("greater", "less", "two.sided")) {
  alternative <- match.arg(alternative)
  Y <- as.matrix(Y)
  X <- check_data(X, "X", "pcm")
  if (!is.null(Z)) Z <- check_data(Z, "Z", "pcm")
  n <- NROW(Y); m <- NCOL(Y)
  if (NROW(X) != n) stop("X and Y must have the same number of rows.")
  if (!is.null(Z) && NROW(Z) != n)
    stop("Z and Y must have the same number of rows.")

  ## Shared sample split so a single sigma updates all hypotheses.
  idx_D2 <- sample.int(n, ceiling(frac * n))

  caches <- vector("list", m)
  for (j in seq_len(m)) {
    caches[[j]] <- precompute_pcm_pair(
      Y           = Y[, j],
      X           = X, Z = Z, target = target,
      reg_YonXZ   = reg_YonXZ, args_YonXZ   = args_YonXZ,
      reg_YonZ    = reg_YonZ,  args_YonZ    = args_YonZ,
      reg_YhatonZ = reg_YhatonZ, args_YhatonZ = args_YhatonZ,
      reg_RonZ    = reg_RonZ,  args_RonZ    = args_RonZ,
      idx_D2      = idx_D2
    )
  }
  T_orig <- vapply(caches, function(pc) pc$T_obs, numeric(1))
  n_D1   <- caches[[1L]]$n_D1

  permute_stat <- function(active_idx, t) {
    sigma <- sample.int(n_D1)
    pcm_permute_stat_active(caches[active_idx], sigma)
  }

  res <- adaptive_fdr_bh(T_orig, permute_stat,
                         alpha = alpha, h = h, B_max = B_max,
                         alternative = alternative)
  res$T_orig <- T_orig
  res$idx_D2 <- idx_D2
  res
}

#' Asymptotic projected covariance measure test
#'
#' Single-hypothesis asymptotic PCM test computed from a
#' \code{\link{precompute_pcm_pair}} cache: the studentized statistic
#' \eqn{T} is asymptotically standard normal under the null, so the
#' p-value is \code{2 * pnorm(abs(T), lower.tail = FALSE)} for the
#' two-sided alternative.  Sibling to \code{\link{perm_pcm}}.
#'
#' Same regression restrictions as \code{\link{perm_pcm_fdr}}: this fills
#' the gap that \code{comets::pcm} does not support \code{"gols"} (no
#' \code{predict.gols_fit} method), so for NB-distributed \eqn{Y} the
#' user can run asymptotic PCM with NB working weights via this entry
#' point.
#'
#' @param Y Response vector of length \eqn{n}.
#' @param X Covariate matrix.
#' @param Z Optional conditioning matrix.
#' @param target Integer index of the target column in \code{X}.
#' @param reg_YonXZ,reg_YonZ,reg_YhatonZ,reg_RonZ See
#'   \code{\link{precompute_pcm_pair}}.
#' @param args_YonXZ,args_YonZ,args_YhatonZ,args_RonZ Optional argument
#'   lists.
#' @param frac Fraction of rows for D2 (default 0.5).
#' @param indices Optional pre-specified D2 indices (overrides \code{frac}).
#' @param alternative One of \code{"two.sided"} (default), \code{"greater"},
#'   or \code{"less"}.
#'
#' @return An object of class \code{c("asym_pcm", "htest")} with components
#'   \code{statistic}, \code{p.value}, \code{rY}, \code{D1_indices},
#'   \code{D2_indices}.
#' @export
asym_pcm <- function(Y, X, Z = NULL, target = 1L,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     args_YonXZ = NULL, args_YonZ = NULL,
                     args_YhatonZ = NULL, args_RonZ = NULL,
                     frac = 0.5, indices = NULL,
                     alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  cache <- precompute_pcm_pair(
    Y           = Y, X = X, Z = Z, target = target,
    reg_YonXZ   = reg_YonXZ, args_YonXZ   = args_YonXZ,
    reg_YonZ    = reg_YonZ,  args_YonZ    = args_YonZ,
    reg_YhatonZ = reg_YhatonZ, args_YhatonZ = args_YhatonZ,
    reg_RonZ    = reg_RonZ,  args_RonZ    = args_RonZ,
    idx_D2      = indices, frac = frac
  )
  T_obs <- cache$T_obs
  pval <- switch(alternative,
    "two.sided" = 2 * stats::pnorm(abs(T_obs), lower.tail = FALSE),
    "greater"   = stats::pnorm(T_obs, lower.tail = FALSE),
    "less"      = stats::pnorm(T_obs)
  )
  stat <- T_obs; names(stat) <- "T"
  structure(list(
    statistic   = stat,
    p.value     = pval,
    method      = "Asymptotic projected covariance measure test",
    data.name   = paste0(deparse(match.call()), collapse = "\n"),
    alternative = alternative,
    rY          = cache$rY_D1,
    D1_indices  = cache$idx_D1,
    D2_indices  = cache$idx_D2
  ), class = c("asym_pcm", "htest"))
}
