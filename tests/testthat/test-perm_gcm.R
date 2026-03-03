test_that("perm_gcm returns correct htest structure", {
  set.seed(123)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm", B = 49L)

  expect_s3_class(result, "htest")
  expect_s3_class(result, "perm_gcm")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_equal(result$parameter[["B"]], 49L)
  expect_length(result$perm_stats, 49L)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_equal(result$alternative, "two.sided")
  expect_length(result$rY, n)
  expect_length(result$rX, n)
})

test_that("perm_gcm respects alternative argument", {
  set.seed(456)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  res_g <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 49L, alternative = "greater")
  res_l <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 49L, alternative = "less")
  res_t <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 49L, alternative = "two.sided")

  expect_equal(res_g$alternative, "greater")
  expect_equal(res_l$alternative, "less")
  expect_equal(res_t$alternative, "two.sided")
})

test_that("perm_gcm does not reject under strong null (linear)", {
  # Under strong null X _||_ (Y, Z), p-values should be roughly uniform.
  # We check that at least most runs do not reject at alpha = 0.05.
  set.seed(789)
  n <- 150
  n_rep <- 20
  pvals <- numeric(n_rep)

  for (i in seq_len(n_rep)) {
    Z <- matrix(rnorm(n), ncol = 1)
    X <- rnorm(n)  # X independent of (Y, Z)
    Y <- 2 * Z + rnorm(n)
    pvals[i] <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                          B = 99L, alternative = "greater")$p.value
  }

  # Under exactness, rejection rate should be near 0.05
  reject_rate <- mean(pvals <= 0.05)
  expect_true(reject_rate <= 0.3)  # conservative check
})

test_that("perm_gcm has power under alternative", {
  set.seed(101)
  n <- 200
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + 0.5 * X + rnorm(n)  # Y depends on X given Z

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 199L, alternative = "greater")
  expect_true(result$p.value < 0.1)
})

test_that("perm_gcm works with multivariate X", {
  set.seed(202)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- cbind(Z + rnorm(n), rnorm(n))
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm", B = 49L)
  expect_s3_class(result, "perm_gcm")
  expect_true(!is.null(result$p.value))
})

test_that("perm_gcm return_fitted_models works", {
  set.seed(303)
  n <- 80
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 19L, return_fitted_models = TRUE)
  expect_true(!is.null(result$models))
  expect_true("reg_YonZ" %in% names(result$models))
  expect_true("reg_XonZ" %in% names(result$models))
})
