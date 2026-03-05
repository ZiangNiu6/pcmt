#' Extract Score Residuals from a Fitted Model
#'
#' Computes score residuals from a fitted regression model. Score residuals are
#' the derivative of the log-likelihood with respect to a scalar perturbation
#' parameter, evaluated at zero. They generalize ordinary residuals to
#' non-linear and non-additive noise models (e.g., GLMs, Cox PH, TRAMs).
#'
#' @param model A fitted model object (e.g., from \code{coxph}, \code{glm},
#'   \code{lm}, or a \code{tram} model).
#' @param data The data used to fit the model (needed for some model types).
#' @param type Character string specifying the type of score residual.
#'   Default \code{"auto"} selects based on model class.
#' @return A numeric vector of length n (univariate score residuals) or an
#'   n x p matrix (multivariate, e.g., from Cox PH with p covariates).
#'   For the ICP application we typically need a single column, so multivariate
#'   score residuals are summed row-wise when \code{type = "auto"}.
#'
#' @details
#' These are TRAM score residuals: the derivative of the log-likelihood with
#' respect to a scalar perturbation of the transformation function h, evaluated
#' at zero (Definition 6 in Kook et al., 2024). They are NOT the same as
#' regression parameter score contributions.
#'
#' Dispatches on model class:
#' \describe{
#'   \item{\code{coxph}}{Martingale residuals:
#'     \code{residuals(model, type = "martingale")}. These correspond to the
#'     TRAM score residuals for the Cox/Weibull model with minimum extreme
#'     value error distribution.}
#'   \item{\code{glm}}{Response residuals \eqn{y_i - \hat{\mu}_i}:
#'     \code{residuals(model, type = "response")}. For canonical link GLMs,
#'     these equal the TRAM score residuals.}
#'   \item{\code{lm}}{Ordinary residuals \eqn{y_i - \hat{y}_i}:
#'     \code{residuals(model)}. Proportional to TRAM score residuals for
#'     the normal linear model.}
#'   \item{\code{mlt}/\code{tram}}{Uses \code{residuals(model, type = "score")}
#'     which directly returns TRAM score residuals.}
#' }
#'
#' @keywords internal
compute_score_residuals <- function(model, data = NULL, type = "auto") {
  cls <- class(model)

  if ("coxph" %in% cls) {
    ## Martingale residuals = TRAM score residuals for Cox PH
    return(c(residuals(model, type = "martingale")))
  }

  if ("glm" %in% cls) {
    ## Default (deviance) residuals, consistent with tramicp implementation.
    ## For GLMs, stats::residuals() returns deviance residuals, which are
    ## signed square roots of the deviance contributions. These work well
    ## as score residuals for the TRAM-GCM invariance test.
    return(c(residuals(model)))
  }

  if ("lm" %in% cls) {
    ## Ordinary residuals (proportional to TRAM score residuals)
    return(c(residuals(model)))
  }

  # tram / mlt models
  if (any(c("tram", "mlt") %in% cls)) {
    sr <- residuals(model, type = "score")
    if (is.matrix(sr)) return(rowSums(sr))
    return(c(sr))
  }

  # Generic fallback: try residuals()
  r <- tryCatch(c(residuals(model)), error = function(e) NULL)
  if (!is.null(r)) return(r)

  stop("Cannot extract score residuals from model of class: ",
       paste(cls, collapse = ", "))
}

#' Residualize Environment on Covariates
#'
#' Computes residualized environments \eqn{\hat{\xi} = E - \hat{\mu}(X^S)}
#' where \eqn{\hat{\mu}} estimates \eqn{E[E | X^S]}.
#'
#' @param E Numeric vector or matrix of environment variable(s).
#' @param X_S Numeric matrix of conditioning covariates.
#' @param method Character string: regression method for E on X^S.
#'   \code{"lrm"} (default) for linear regression, \code{"rf"} for random forest.
#' @return Numeric vector (if E is univariate) or matrix of residualized
#'   environments.
#'
#' @keywords internal
residualize_environment <- function(E, X_S, method = "lrm") {
  E <- as.matrix(E)
  X_S <- as.matrix(X_S)
  n <- nrow(E)

  if (ncol(X_S) == 0L || all(apply(X_S, 2, stats::var) == 0)) {
    # No conditioning covariates or constant: return centered E
    return(scale(E, scale = FALSE))
  }

  resid_mat <- matrix(0, nrow = n, ncol = ncol(E))

  for (j in seq_len(ncol(E))) {
    e_j <- E[, j]

    if (identical(method, "lrm")) {
      # Fast QR path
      design <- cbind(1, X_S)
      qr_obj <- qr(design)
      Q <- qr.Q(qr_obj)[, seq_len(qr_obj$rank), drop = FALSE]
      resid_mat[, j] <- e_j - Q %*% crossprod(Q, e_j)
    } else if (identical(method, "rf")) {
      # Use ranger directly to avoid comets::rf() converting binary to factor
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("Package 'ranger' is required for reg_env = \"rf\".")
      }
      X_named <- as.data.frame(X_S)
      if (is.null(names(X_named))) {
        names(X_named) <- paste0("V", seq_len(ncol(X_S)))
      }
      X_named[["__env__"]] <- e_j
      rf <- ranger::ranger(`__env__` ~ ., data = X_named)
      resid_mat[, j] <- e_j - rf$predictions
    } else {
      # Use comets regression infrastructure for other methods
      X_named <- X_S
      if (is.null(colnames(X_named))) {
        colnames(X_named) <- paste0("V", seq_len(ncol(X_named)))
      }
      m <- fit_regression(y = e_j, x = X_named, method = method)
      fitted_vals <- if (!is.null(m$predictions)) {
        c(m$predictions)
      } else {
        stats::fitted(m)
      }
      resid_mat[, j] <- e_j - fitted_vals
    }
  }

  if (ncol(E) == 1L) return(c(resid_mat))
  resid_mat
}
