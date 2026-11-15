# Small, seeded RF fixtures and independent full-refit references.
# The reference uses comets directly and never calls a pcmt fit/statistic helper.
paper_data <- function(n = 40L) {
  set.seed(615L)
  Z <- matrix(stats::rnorm(n * 2L), ncol = 2L,
              dimnames = list(NULL, c("Z1", "Z2")))
  X <- matrix(0.4 * Z[, 1L] + stats::rnorm(n), ncol = 1L,
              dimnames = list(NULL, "X1"))
  Y <- sin(Z[, 1L]) + 0.3 * X[, 1L] + stats::rnorm(n)
  X_pcm <- cbind(W1 = stats::rnorm(n), target = X[, 1L],
                 W2 = stats::runif(n))
  list(Y = Y, X = X, Z = Z, X_pcm = X_pcm)
}

paper_rf_args <- function(seed) {
  list(seed = as.integer(seed), num.threads = 1L, num.trees = 12L)
}

paper_gcm_args <- function(dat = paper_data()) {
  list(Y = dat$Y, X = dat$X, Z = dat$Z,
       args_YonZ = paper_rf_args(101L),
       args_XonZ = paper_rf_args(102L), n_cores = 1L)
}

paper_pcm_args <- function(dat = paper_data()) {
  list(Y = dat$Y, X = dat$X_pcm, Z = dat$Z, target = 2L,
       indices = rev(seq.int(2L, length(dat$Y), by = 2L)),
       args_YonXZ = paper_rf_args(201L),
       args_YonZ = paper_rf_args(202L),
       args_YhatonZ = paper_rf_args(203L),
       args_RonZ = paper_rf_args(204L), n_cores = 1L)
}

paper_rf <- function(y, x, args) {
  do.call(get("rf", envir = asNamespace("comets")),
          c(list(y = y, x = x), args))
}

paper_stat <- function(products) {
  centered_ss <- sum((products - mean(products))^2)
  if (centered_ss <= 0) return(0)
  sum(products) / sqrt(centered_ss)
}

paper_permutation_p <- function(observed, permuted, alternative) {
  exceed <- switch(alternative,
    two.sided = abs(permuted) >= abs(observed),
    greater = permuted >= observed,
    less = permuted <= observed)
  (1 + sum(exceed)) / (1 + length(permuted))
}

paper_gcm_reference <- function(args, B, alternative) {
  Y <- args$Y; X <- args$X; Z <- args$Z
  fit_y <- paper_rf(Y, Z, args$args_YonZ)
  rY <- Y - as.numeric(stats::predict(fit_y, data = Z))
  fit_x <- paper_rf(X, Z, args$args_XonZ)
  rX <- as.numeric(X) - as.numeric(stats::predict(fit_x, data = Z))
  observed <- paper_stat(rY * rX)
  permuted <- vapply(seq_len(B), function(b) {
    Xp <- X[sample.int(NROW(X)), , drop = FALSE]
    fit <- paper_rf(Xp, Z, args$args_XonZ)
    residual <- as.numeric(Xp) - as.numeric(stats::predict(fit, data = Z))
    paper_stat(rY * residual)
  }, numeric(1))
  list(statistic = observed, perm_stats = permuted,
       p.value = paper_permutation_p(observed, permuted, alternative))
}

paper_pcm_reference <- function(args, B, alternative) {
  Y <- args$Y; X <- args$X; Z <- args$Z; target <- args$target
  D2 <- args$indices; D1 <- setdiff(seq_along(Y), D2)
  Xj <- X[, target, drop = FALSE]
  W <- cbind(X[, -target, drop = FALSE], Z)
  XW <- cbind(Xj, W)
  fit_p <- paper_rf(Y[D2], XW[D2, , drop = FALSE], args$args_YonXZ)
  fit_q <- paper_rf(stats::predict(fit_p, data = XW[D2, , drop = FALSE]),
                   W[D2, , drop = FALSE], args$args_YhatonZ)
  direction <- function(xw) {
    as.numeric(stats::predict(fit_p, data = xw)) -
      as.numeric(stats::predict(fit_q, data = W[D1, , drop = FALSE]))
  }
  fit_y <- paper_rf(Y[D1], W[D1, , drop = FALSE], args$args_YonZ)
  rY <- Y[D1] - as.numeric(stats::predict(fit_y, data = W[D1, , drop = FALSE]))
  fhat <- direction(XW[D1, , drop = FALSE])
  fit_f <- paper_rf(fhat, W[D1, , drop = FALSE], args$args_RonZ)
  rF <- fhat - as.numeric(stats::predict(fit_f, data = W[D1, , drop = FALSE]))
  observed <- paper_stat(rY * rF)
  permuted <- vapply(seq_len(B), function(b) {
    Xp <- Xj[D1, , drop = FALSE][sample.int(length(D1)), , drop = FALSE]
    fp <- direction(cbind(Xp, W[D1, , drop = FALSE]))
    fit <- paper_rf(fp, W[D1, , drop = FALSE], args$args_RonZ)
    residual <- fp - as.numeric(stats::predict(fit, data = W[D1, , drop = FALSE]))
    paper_stat(rY * residual)
  }, numeric(1))
  list(statistic = observed, perm_stats = permuted,
       p.value = paper_permutation_p(observed, permuted, alternative))
}

expect_adaptive_contract <- function(fit, fixed, h, alpha, B) {
  actual <- length(fit$perm_stats)
  expect_equal(unname(fit$statistic), unname(fixed$statistic), tolerance = 1e-12)
  expect_equal(fit$perm_stats, fixed$perm_stats[seq_len(actual)], tolerance = 1e-12)
  expect_equal(unname(fit$parameter[["B"]]), actual)
  tr <- fit$adaptive$trajectory
  scores <- switch(fit$alternative, two.sided = abs(fit$perm_stats),
                   greater = fit$perm_stats, less = -fit$perm_stats)
  observed <- switch(fit$alternative, two.sided = abs(fit$statistic),
                     greater = fit$statistic, less = -fit$statistic)
  losses <- cumsum(scores >= observed)
  expect_equal(tr$t, seq_len(actual))
  expect_equal(tr$K, losses)
  expect_equal(tr$S, scores)
  expected_p <- h / (seq_len(actual) + h - losses)
  expect_equal(tr$p_AV, expected_p, tolerance = 1e-12)
  expect_equal(fit$p.value, utils::tail(expected_p, 1L), tolerance = 1e-12)
  expect_equal(fit$adaptive$h, h)
  expect_equal(fit$adaptive$alpha, alpha)
  expect_equal(fit$adaptive$K, utils::tail(losses, 1L))
  expect_identical(fit$adaptive$stopped_early, actual < B)
  expected_reason <- if (fit$p.value <= alpha) "rejection" else
    if (utils::tail(losses, 1L) == h) "futility" else "cap"
  expect_identical(fit$adaptive$stop_reason, expected_reason)
}
