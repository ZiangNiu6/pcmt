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
#' @noRd
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
#' Implements the adaptive resampling procedure used in the manuscript.
#' The procedure processes \eqn{m} hypotheses
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
#'   \code{"two.sided"}. Statistics are compared directly, negated, or taken
#'   in absolute value, respectively.
#'
#' @details
#' The anytime-valid p-value is
#' \deqn{p^i_t = h / (t + h - K^i_t),}
#' where \eqn{K^i_t = \sum_{s=1}^{t} \mathbb{1}\{S^i_s \ge S^i_{\text{orig}}\}}
#' and \eqn{S} is the corresponding one- or two-sided score. When
#' \eqn{K^i_t = h}, the hypothesis is moved to the futility set and its
#' p-value is frozen at \eqn{h / t}.  Active-hypothesis p-values are
#' monotonically decreasing in \eqn{t}, so the BH threshold is monotonically
#' increasing and rejection is sticky.
#'
#' An FDR interpretation requires valid resampling p-values and dependence
#' conditions under which the Benjamini--Hochberg procedure controls FDR.
#'
#' @references Ai and Niu (2026), \doi{10.48550/arXiv.2611.03255}.
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
  h <- validate_positive_integer(h, "h")
  B_max <- validate_positive_integer(B_max, "B_max")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be one number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(T_orig) || anyNA(T_orig)) {
    stop("T_orig must contain numeric statistics without missing values.", call. = FALSE)
  }
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
    if (!is.numeric(z_perm) || anyNA(z_perm)) {
      stop("permute_stat must return numeric statistics without missing values.", call. = FALSE)
    }
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
