## Validity tests for adaptive permutation calibration
##
## Under the null hypothesis, the rejection probability should be <= alpha.
## We run many independent replications and check that the empirical rejection
## rate is consistent with the nominal level using a one-sided binomial test.

test_that("adaptive perm_gcm controls type-I error under the null (lrm)", {
  skip_on_cran()

  alpha <- 0.05
  h <- 15L
  B_max <- 499L
  n_reps <- 200L
  n <- 100L

  rejections <- 0L
  for (r in seq_len(n_reps)) {
    set.seed(r)
    Z <- matrix(rnorm(n), ncol = 1)
    X <- rnorm(n)           # X independent of (Y, Z)
    Y <- Z[, 1] + rnorm(n)

    res <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                    B = B_max, h = h, alpha = alpha,
                    alternative = "two.sided")
    if (res$p.value <= alpha) rejections <- rejections + 1L
  }

  ## One-sided binomial test: H0: rejection rate <= alpha
  ## If the adaptive p-value is valid, we should NOT reject this binomial test.
  bt <- binom.test(rejections, n_reps, p = alpha, alternative = "greater")
  expect_gt(bt$p.value, 0.01,
    label = sprintf("Adaptive GCM rejection rate = %d/%d = %.3f (binomial p = %.4f)",
                    rejections, n_reps, rejections / n_reps, bt$p.value))
})

test_that("adaptive perm_gcm controls type-I error under the null (strong null, rf)", {
  skip_on_cran()

  alpha <- 0.05
  h <- 20L
  B_max <- 499L
  n_reps <- 200L
  n <- 100L

  rejections <- 0L
  for (r in seq_len(n_reps)) {
    set.seed(1000 + r)
    Z <- matrix(rnorm(n * 2), ncol = 2)
    X <- rnorm(n)           # X independent of (Y, Z)
    Y <- sin(Z[, 1]) + Z[, 2]^2 + rnorm(n)

    res <- perm_gcm(Y, X, Z, reg_YonZ = "rf", reg_XonZ = "rf",
                    B = B_max, h = h, alpha = alpha,
                    alternative = "two.sided")
    if (res$p.value <= alpha) rejections <- rejections + 1L
  }

  bt <- binom.test(rejections, n_reps, p = alpha, alternative = "greater")
  expect_gt(bt$p.value, 0.01,
    label = sprintf("Adaptive GCM (rf) rejection rate = %d/%d = %.3f",
                    rejections, n_reps, rejections / n_reps))
})

test_that("adaptive perm_pcm controls type-I error under the null (lrm)", {
  skip_on_cran()

  alpha <- 0.05
  h <- 15L
  B_max <- 499L
  n_reps <- 200L
  n <- 200L  # need 2n total for sample split

  rejections <- 0L
  for (r in seq_len(n_reps)) {
    set.seed(2000 + r)
    X <- matrix(rnorm(n * 3), ncol = 3)
    Y <- X[, 2] + X[, 3] + rnorm(n)  # Y _||_ X_1 | X_{-1}

    res <- perm_pcm(Y, X, target = 1,
                    reg_YonXZ = "lrm", reg_YonZ = "lrm",
                    reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                    B = B_max, h = h, alpha = alpha,
                    alternative = "two.sided")
    if (res$p.value <= alpha) rejections <- rejections + 1L
  }

  bt <- binom.test(rejections, n_reps, p = alpha, alternative = "greater")
  expect_gt(bt$p.value, 0.01,
    label = sprintf("Adaptive PCM rejection rate = %d/%d = %.3f",
                    rejections, n_reps, rejections / n_reps))
})

test_that("adaptive perm_wgcm controls type-I error under the null (lrm)", {
  skip_on_cran()

  alpha <- 0.05
  h <- 15L
  B_max <- 499L
  n_reps <- 200L
  n <- 150L

  rejections <- 0L
  for (r in seq_len(n_reps)) {
    set.seed(3000 + r)
    Z <- matrix(rnorm(n), ncol = 1)
    X <- rnorm(n)           # X independent of (Y, Z)
    Y <- Z[, 1] + rnorm(n)

    res <- perm_wgcm(Y, X, Z,
                     reg_YonZ = "lrm", reg_XonZ = "lrm",
                     reg_wfun = "lrm",
                     B = B_max, h = h, alpha = alpha,
                     alternative = "two.sided")
    if (res$p.value <= alpha) rejections <- rejections + 1L
  }

  bt <- binom.test(rejections, n_reps, p = alpha, alternative = "greater")
  expect_gt(bt$p.value, 0.01,
    label = sprintf("Adaptive WGCM rejection rate = %d/%d = %.3f",
                    rejections, n_reps, rejections / n_reps))
})

test_that("adaptive perm_gcm controls type-I error at different h values", {
  skip_on_cran()

  alpha <- 0.05
  B_max <- 499L
  n_reps <- 200L
  n <- 100L

  for (h_val in c(5L, 10L, 30L)) {
    rejections <- 0L
    for (r in seq_len(n_reps)) {
      set.seed(4000 + h_val * 1000 + r)
      Z <- matrix(rnorm(n), ncol = 1)
      X <- rnorm(n)
      Y <- Z[, 1] + rnorm(n)

      res <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      B = B_max, h = h_val, alpha = alpha,
                      alternative = "two.sided")
      if (res$p.value <= alpha) rejections <- rejections + 1L
    }

    bt <- binom.test(rejections, n_reps, p = alpha, alternative = "greater")
    expect_gt(bt$p.value, 0.01,
      label = sprintf("h=%d: rejection rate = %d/%d = %.3f",
                      h_val, rejections, n_reps, rejections / n_reps))
  }
})

test_that("adaptive perm_gcm controls type-I error with one-sided alternatives", {
  skip_on_cran()

  alpha <- 0.05
  h <- 15L
  B_max <- 499L
  n_reps <- 200L
  n <- 100L

  for (alt in c("greater", "less")) {
    rejections <- 0L
    for (r in seq_len(n_reps)) {
      set.seed(5000 + match(alt, c("greater", "less")) * 1000 + r)
      Z <- matrix(rnorm(n), ncol = 1)
      X <- rnorm(n)
      Y <- Z[, 1] + rnorm(n)

      res <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      B = B_max, h = h, alpha = alpha,
                      alternative = alt)
      if (res$p.value <= alpha) rejections <- rejections + 1L
    }

    bt <- binom.test(rejections, n_reps, p = alpha, alternative = "greater")
    expect_gt(bt$p.value, 0.01,
      label = sprintf("alt=%s: rejection rate = %d/%d = %.3f",
                      alt, rejections, n_reps, rejections / n_reps))
  }
})
