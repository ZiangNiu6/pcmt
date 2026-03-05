#' Permutation-Calibrated Invariant Causal Prediction
#'
#' Identifies plausible direct causes of a response variable using invariant
#' causal prediction (ICP) with permutation-calibrated invariance tests.
#' For each subset S of covariates, fits a model of Y | X^S, computes score
#' residuals, and tests whether these are uncorrelated with (residualized)
#' environments. Returns the intersection of all non-rejected subsets.
#'
#' @param formula A formula specifying the response and candidate predictors,
#'   e.g., \code{Y ~ X1 + X2 + X3} or \code{Surv(time, status) ~ X1 + X2}.
#' @param data A data frame containing all variables.
#' @param env A one-sided formula specifying the environment variable(s),
#'   e.g., \code{~ E} or \code{~ E1 + E2}.
#' @param modFUN A model-fitting function. Must accept \code{formula} and
#'   \code{data} arguments, e.g., \code{survival::coxph}, \code{stats::glm},
#'   \code{stats::lm}. Additional arguments can be passed via \code{...}.
#' @param alpha Significance level (default 0.05).
#' @param B Number of permutations for the invariance test (default 499).
#' @param mandatory A one-sided formula of predictors that must always be
#'   included in the conditioning set, e.g., \code{~ X1 + X2}. This reduces
#'   the search space (only supersets of the mandatory set are tested).
#' @param max_size Maximum subset size to test (default NULL = all subsets).
#'   Useful for computational feasibility when d is large.
#' @param reg_env Method for residualizing the environment on X^S.
#'   Default \code{"rf"} (random forest, consistent with tramicp).
#'   Can also be \code{"lrm"} for linear regression.
#' @param score_type Type of score residual to extract. Default \code{"auto"}.
#' @param verbose Logical. Print progress (default FALSE).
#' @param ... Additional arguments passed to \code{modFUN} (e.g.,
#'   \code{family = binomial} for \code{glm}).
#'
#' @details
#' Implements Algorithm 1 of Peters et al. (2016) with the TRAM-GCM invariance
#' test of Kook et al. (2024), calibrated via permutations instead of
#' asymptotic chi-squared.
#'
#' The algorithm enumerates all subsets \eqn{S \subseteq \{1, \ldots, d\}}
#' (respecting \code{mandatory} and \code{max_size}), fits the model for each,
#' and tests invariance. The output is the intersection of all non-rejected
#' sets, which is guaranteed to be a subset of the true causal parents with
#' probability at least \eqn{1 - \alpha}.
#'
#' Predictor-specific p-values are computed as:
#' \eqn{\hat{p}_j = \max_{S \subseteq [d]: j \notin S} p_{S,n}} if the
#' maximum set-specific p-value is at least \eqn{\alpha}, and 1 otherwise
#' (i.e., when all sets are rejected).
#'
#' @return An object of class \code{"perm_icp"} with components:
#' \describe{
#'   \item{causal_set}{Character vector of identified causal predictors.}
#'   \item{pvalues_set}{Named numeric vector of set-specific p-values.}
#'   \item{pvalues_predictor}{Named numeric vector of predictor-specific p-values.}
#'   \item{alpha}{The significance level used.}
#'   \item{B}{Number of permutations.}
#'   \item{n_sets_tested}{Number of subsets tested.}
#'   \item{call}{The matched call.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 300
#' E <- rbinom(n, 1, 0.5)
#' X1 <- E + rnorm(n)
#' X2 <- rnorm(n)
#' Y <- 2 * X1 + rnorm(n)  # X1 is causal, X2 is not
#' df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)
#' result <- perm_icp(Y ~ X1 + X2, data = df, env = ~ E,
#'                    modFUN = stats::lm, B = 199)
#' print(result)
#' }
#'
#' @export
perm_icp <- function(formula, data, env, modFUN,
                     alpha = 0.05, B = 499L,
                     mandatory = NULL, max_size = NULL,
                     reg_env = "rf", score_type = "auto",
                     verbose = FALSE, ...) {
  cl <- match.call()
  B <- as.integer(B)

  ## Parse response and predictors from formula
  resp_var <- all.vars(formula)[1L]
  pred_names <- attr(stats::terms(formula, data = data), "term.labels")
  d <- length(pred_names)

  ## Parse environment variables
  env_names <- all.vars(env)
  E_mat <- as.matrix(data[, env_names, drop = FALSE])

  ## Parse mandatory predictors
  mandatory_names <- character(0)
  if (!is.null(mandatory)) {
    mandatory_names <- all.vars(mandatory)
    if (!all(mandatory_names %in% pred_names)) {
      stop("All mandatory predictors must appear in the formula.")
    }
  }
  free_names <- setdiff(pred_names, mandatory_names)
  d_free <- length(free_names)

  ## Enumerate subsets of free predictors
  if (is.null(max_size)) max_size <- d_free
  max_size <- min(max_size, d_free)

  subsets <- list()
  for (k in 0:max_size) {
    combos <- utils::combn(d_free, k, simplify = FALSE)
    for (combo in combos) {
      subsets <- c(subsets, list(c(mandatory_names, free_names[combo])))
    }
  }
  n_sets <- length(subsets)

  if (verbose) {
    message(sprintf("Testing %d subsets (d=%d, mandatory=%d, max_size=%d)",
                    n_sets, d, length(mandatory_names), max_size))
  }

  ## Test each subset
  pval_set <- numeric(n_sets)
  set_labels <- character(n_sets)

  for (i in seq_len(n_sets)) {
    S <- subsets[[i]]
    set_labels[i] <- if (length(S) == 0L) "Empty" else paste(S, collapse = "+")

    if (verbose) {
      message(sprintf("  [%d/%d] Testing S = {%s}", i, n_sets, set_labels[i]))
    }

    ## Build formula for this subset
    if (length(S) == 0L) {
      # Intercept-only model
      resp_part <- deparse(formula[[2]])
      sub_formula <- stats::as.formula(paste(resp_part, "~ 1"))
    } else {
      resp_part <- deparse(formula[[2]])
      sub_formula <- stats::as.formula(paste(resp_part, "~",
                                             paste(S, collapse = " + ")))
    }

    ## Fit the model
    fit <- tryCatch(
      modFUN(formula = sub_formula, data = data, ...),
      error = function(e) {
        if (verbose) message("    Model fitting failed: ", e$message)
        NULL
      }
    )

    if (is.null(fit)) {
      pval_set[i] <- NA_real_
      next
    }

    ## Extract score residuals
    score_r <- tryCatch(
      compute_score_residuals(fit, data = data, type = score_type),
      error = function(e) {
        if (verbose) message("    Score residual extraction failed: ", e$message)
        NULL
      }
    )

    if (is.null(score_r)) {
      pval_set[i] <- NA_real_
      next
    }

    ## Residualize environment on X^S
    if (length(S) == 0L) {
      xi <- scale(E_mat, scale = FALSE)
    } else {
      X_S <- as.matrix(data[, S, drop = FALSE])
      xi <- residualize_environment(E_mat, X_S, method = reg_env)
    }

    ## Run permutation invariance test
    test_result <- perm_tram_gcm(score_r, xi, B = B)
    pval_set[i] <- test_result$p.value
  }

  names(pval_set) <- set_labels

  ## Treat NA p-values (failed fits) as rejected (p = 0)
  pval_set_clean <- pval_set
  pval_set_clean[is.na(pval_set_clean)] <- 0

  ## Compute ICP intersection: S_n = intersect of all S with p > alpha
  non_rejected <- which(pval_set_clean > alpha)
  if (length(non_rejected) == 0L) {
    causal_set <- character(0)
  } else {
    causal_set <- Reduce(intersect, subsets[non_rejected])
  }

  ## Compute predictor-specific p-values
  ## p_hat_j = max_{S: j not in S} p_{S,n}
  ## If max_S p_{S,n} < alpha (all rejected), set p_hat_j = 1
  all_rejected <- all(pval_set_clean <= alpha)
  pval_pred <- numeric(d)
  names(pval_pred) <- pred_names

  for (j in seq_len(d)) {
    xname <- pred_names[j]
    # Find subsets that do NOT contain this predictor
    idx_without_j <- which(vapply(subsets, function(S) !(xname %in% S),
                                  logical(1)))
    if (length(idx_without_j) == 0L || all_rejected) {
      pval_pred[j] <- 1
    } else {
      pval_pred[j] <- max(pval_set_clean[idx_without_j])
    }
  }

  structure(list(
    causal_set = causal_set,
    pvalues_set = pval_set,
    pvalues_predictor = pval_pred,
    alpha = alpha,
    B = B,
    n_sets_tested = n_sets,
    call = cl
  ), class = "perm_icp")
}

#' @export
print.perm_icp <- function(x, ...) {
  cat("Permutation-Calibrated Invariant Causal Prediction\n")
  cat("Call: ")
  print(x$call)
  cat("\nPredictor p-values:\n")
  print(round(x$pvalues_predictor, 4))
  cat("\nSet of plausible causal predictors (alpha =", x$alpha, "):\n")
  if (length(x$causal_set) == 0L) {
    cat("  Empty set\n")
  } else {
    cat(" ", paste(x$causal_set, collapse = ", "), "\n")
  }
  cat(sprintf("\n(%d subsets tested, B = %d permutations)\n",
              x$n_sets_tested, x$B))
  invisible(x)
}

#' Extract p-values from perm_icp result
#'
#' @param object A \code{perm_icp} object.
#' @param which Which p-values: \code{"predictor"} (default) or \code{"set"}.
#' @param ... Ignored.
#' @return Named numeric vector of p-values.
#' @export
pvalues <- function(object, ...) UseMethod("pvalues")

#' @rdname pvalues
#' @export
pvalues.perm_icp <- function(object, which = c("predictor", "set"), ...) {
  which <- match.arg(which)
  if (which == "predictor") object$pvalues_predictor else object$pvalues_set
}

# ---- Convenience wrappers ----

#' ICP with Cox Proportional Hazards Model
#'
#' Wrapper around \code{\link{perm_icp}} using \code{survival::coxph} as the
#' model-fitting function.
#'
#' @inheritParams perm_icp
#' @param ... Additional arguments passed to \code{survival::coxph}.
#' @return A \code{perm_icp} object.
#'
#' @examples
#' \donttest{
#' library(survival)
#' # Simulated survival data
#' set.seed(42)
#' n <- 200
#' E <- rbinom(n, 1, 0.5)
#' X1 <- E + rnorm(n)
#' X2 <- rnorm(n)
#' time <- rexp(n, exp(0.5 * X1))
#' status <- rep(1, n)
#' df <- data.frame(time, status, X1, X2, E)
#' result <- coxphICP(Surv(time, status) ~ X1 + X2, data = df,
#'                    env = ~ E, B = 199)
#' print(result)
#' }
#'
#' @export
coxphICP <- function(formula, data, env, alpha = 0.05, B = 499L,
                     mandatory = NULL, max_size = NULL,
                     reg_env = "rf", verbose = FALSE, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for coxphICP.")
  }
  perm_icp(formula = formula, data = data, env = env,
           modFUN = survival::coxph, alpha = alpha, B = B,
           mandatory = mandatory, max_size = max_size,
           reg_env = reg_env, verbose = verbose, ...)
}

#' ICP with Generalized Linear Model
#'
#' Wrapper around \code{\link{perm_icp}} using \code{stats::glm} as the
#' model-fitting function.
#'
#' @inheritParams perm_icp
#' @param family A \code{family} object for the GLM (default \code{gaussian()}).
#' @param ... Additional arguments passed to \code{stats::glm}.
#' @return A \code{perm_icp} object.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 300
#' E <- rbinom(n, 1, 0.5)
#' X1 <- -E + rnorm(n)
#' X2 <- rnorm(n)
#' Y <- rbinom(n, 1, plogis(0.5 * X1))
#' df <- data.frame(Y, X1, X2, E)
#' result <- glmICP(Y ~ X1 + X2, data = df, env = ~ E,
#'                  family = binomial, B = 199)
#' print(result)
#' }
#'
#' @export
glmICP <- function(formula, data, env, family = gaussian(),
                   alpha = 0.05, B = 499L,
                   mandatory = NULL, max_size = NULL,
                   reg_env = "lrm", verbose = FALSE, ...) {
  perm_icp(formula = formula, data = data, env = env,
           modFUN = stats::glm, alpha = alpha, B = B,
           mandatory = mandatory, max_size = max_size,
           reg_env = reg_env, verbose = verbose,
           family = family, ...)
}

#' ICP with Linear Model
#'
#' Wrapper around \code{\link{perm_icp}} using \code{stats::lm} as the
#' model-fitting function.
#'
#' @inheritParams perm_icp
#' @param ... Additional arguments passed to \code{stats::lm}.
#' @return A \code{perm_icp} object.
#'
#' @export
lmICP <- function(formula, data, env, alpha = 0.05, B = 499L,
                  mandatory = NULL, max_size = NULL,
                  reg_env = "lrm", verbose = FALSE, ...) {
  perm_icp(formula = formula, data = data, env = env,
           modFUN = stats::lm, alpha = alpha, B = B,
           mandatory = mandatory, max_size = max_size,
           reg_env = reg_env, verbose = verbose, ...)
}
