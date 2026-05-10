test_that("precompute_pcm_pair errors on disallowed regression types", {
  set.seed(1)
  n <- 100
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- X[, 2] + rnorm(n)

  expect_error(
    precompute_pcm_pair(Y, X, target = 1, reg_YonXZ = "rf",
                         reg_YonZ = "lrm", reg_YhatonZ = "lrm",
                         reg_RonZ = "lrm"),
    "reg_YonXZ"
  )
  expect_error(
    precompute_pcm_pair(Y, X, target = 1, reg_YonXZ = "lrm",
                         reg_YonZ = "lrm", reg_YhatonZ = "rf",
                         reg_RonZ = "lrm"),
    "reg_YhatonZ"
  )
  expect_error(
    precompute_pcm_pair(Y, X, target = 1, reg_YonXZ = "lrm",
                         reg_YonZ = "lrm", reg_YhatonZ = "lrm",
                         reg_RonZ = "rf"),
    "reg_RonZ"
  )
})

test_that("precompute_pcm_pair returns the expected cache fields", {
  set.seed(2)
  n <- 200
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)
  cache <- precompute_pcm_pair(Y, X, target = 1,
                                reg_YonXZ = "lrm", reg_YonZ = "lrm",
                                reg_YhatonZ = "lrm", reg_RonZ = "lrm")

  expect_named(cache,
    c("idx_D1", "idx_D2", "n_D1", "Xj_D1", "rY_D1",
      "beta_X_p", "fhat_const", "Q_W", "T_obs"))
  expect_equal(length(cache$Xj_D1),     cache$n_D1)
  expect_equal(length(cache$rY_D1),     cache$n_D1)
  expect_equal(length(cache$fhat_const), cache$n_D1)
  expect_equal(length(intersect(cache$idx_D1, cache$idx_D2)), 0L)
  expect_equal(sort(c(cache$idx_D1, cache$idx_D2)), seq_len(n))
  expect_true(is.finite(cache$T_obs))
})

test_that("asym_pcm returns a well-formed htest object", {
  set.seed(3)
  n <- 200
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)

  res <- asym_pcm(Y, X, target = 1,
                   reg_YonXZ = "lrm", reg_YonZ = "lrm",
                   reg_YhatonZ = "lrm", reg_RonZ = "lrm")

  expect_s3_class(res, "htest")
  expect_s3_class(res, "asym_pcm")
  expect_true(is.finite(res$statistic))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_equal(length(res$rY), length(res$D1_indices))
  expect_equal(length(intersect(res$D1_indices, res$D2_indices)), 0L)
})

test_that("cache T_obs matches perm_pcm()$statistic on the same split", {
  set.seed(4)
  n <- 300
  X <- matrix(rnorm(n * 4), n, 4)
  Y <- X[, 2] + X[, 4] + rnorm(n)

  idx_D2 <- sample.int(n, ceiling(0.5 * n))

  cache <- precompute_pcm_pair(Y, X, target = 1,
                                reg_YonXZ = "lrm", reg_YonZ = "lrm",
                                reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                                idx_D2 = idx_D2)
  ref <- perm_pcm(Y, X, target = 1,
                  reg_YonXZ = "lrm", reg_YonZ = "lrm",
                  reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                  indices = idx_D2, B = 1L)

  expect_equal(unname(ref$statistic), cache$T_obs, tolerance = 1e-10)
})

test_that("perm_pcm_fdr controls FDR under the global null", {
  ## Global null: every Y_j is independent of (X, Z).
  skip_on_cran()
  set.seed(5)
  n <- 200; m <- 40
  X <- matrix(rnorm(n * 3), n, 3)
  Z <- matrix(rnorm(n * 2), n, 2)
  Y <- matrix(rnorm(n * m), n, m)
  q <- 0.2

  res <- perm_pcm_fdr(Y, X, Z, target = 1,
                      reg_YonXZ = "lrm", reg_YonZ = "lrm",
                      reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                      alpha = q, h = 10L, B_max = 500L,
                      alternative = "two.sided")
  expect_length(res$T_orig, m)
  expect_true(length(res$rejected) <= ceiling(q * m) + 5L)
})

test_that("perm_pcm_fdr has power against a true alternative", {
  skip_on_cran()
  set.seed(6)
  n <- 300; m <- 20
  X <- matrix(rnorm(n * 2), n, 2)
  Z <- matrix(rnorm(n * 2), n, 2)
  Y <- matrix(rnorm(n * m), n, m)
  ## Inject 5 true alternatives
  alt_idx <- 1:5
  Y[, alt_idx] <- Y[, alt_idx] + matrix(rep(0.6 * X[, 1], length(alt_idx)),
                                         n, length(alt_idx))

  res <- perm_pcm_fdr(Y, X, Z, target = 1,
                      reg_YonXZ = "lrm", reg_YonZ = "lrm",
                      reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                      alpha = 0.1, h = 10L, B_max = 500L,
                      alternative = "two.sided")
  ## At least 2 of 5 alternatives recovered (loose, just to catch regressions)
  expect_true(sum(res$rejected %in% alt_idx) >= 2L)
})

test_that("gols path produces a calibrated p-value on simulated NB data", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(7)
  n <- 500
  Z <- matrix(rnorm(n * 2), n, 2)
  X <- matrix(rnorm(n * 1), n, 1)
  theta <- 2
  mu <- exp(0.3 + Z %*% c(0.1, -0.05))
  Y <- as.numeric(MASS::rnegbin(n, mu = mu, theta = theta))

  res <- asym_pcm(Y, X, Z = Z, target = 1,
                   reg_YonXZ = "gols",
                   args_YonXZ = list(theta = theta),
                   reg_YonZ = "glrm",
                   args_YonZ = list(family = MASS::neg.bin(theta)),
                   reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                   alternative = "two.sided")
  expect_s3_class(res, "asym_pcm")
  expect_true(res$p.value > 0.01)
})
