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
#' @param B Number of permutations (default 499).
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
#'   \item{perm_stats}{Vector of permutation test statistics (length B).}
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
                     alternative = c("two.sided", "greater", "less"),
                     return_fitted_models = FALSE) {
  alternative <- match.arg(alternative)
  B <- as.integer(B)

  ## Data validation
  Y <- check_data(Y, "Y")
  X <- check_data(X, "X")
  Z <- check_data(Z, "Z")
  n <- NROW(Y)

  ## Step 1: Fit Y on Z -> residuals rY (FIXED across permutations)
  mY <- fit_regression(y = Y, x = Z, method = reg_YonZ, args = args_YonZ)
  rY <- stats::residuals(mY, response = Y, data = Z)

  ## Step 2: Fit X on Z -> residuals rX (original)
  mX <- fit_regression(y = X, x = Z, method = reg_XonZ, args = args_XonZ)
  rX <- stats::residuals(mX, response = X, data = Z)

  ## Step 3: Observed test statistic
  R_obs <- c(rX) * c(rY)
  T_obs <- studentized_stat(R_obs)

  ## Step 4: Permutation loop
  T_perm <- numeric(B)
  for (b in seq_len(B)) {
    perm <- sample.int(n)
    X_perm <- X[perm, , drop = FALSE]
    mX_perm <- fit_regression(y = X_perm, x = Z, method = reg_XonZ, args = args_XonZ)
    rX_perm <- stats::residuals(mX_perm, response = X_perm, data = Z)
    R_perm <- c(rX_perm) * c(rY)
    T_perm[b] <- studentized_stat(R_perm)
  }

  ## Step 5: Permutation p-value
  pval <- perm_pvalue(T_obs, T_perm, alternative)

  ## Return htest object
  stat <- T_obs
  names(stat) <- "T"
  param <- c("B" = B)

  models <- if (return_fitted_models) {
    list(reg_YonZ = mY, reg_XonZ = mX)
  } else {
    NULL
  }

  structure(list(
    statistic = stat,
    p.value = pval,
    parameter = param,
    method = "Permutation-calibrated generalised covariance measure test",
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    alternative = alternative,
    rY = rY,
    rX = rX,
    perm_stats = T_perm,
    models = models
  ), class = c("perm_gcm", "htest"))
}
