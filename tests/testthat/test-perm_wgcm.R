test_that("perm_wgcm returns correct htest structure", {
  set.seed(123)
  n <- 150
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      reg_wfun = "lrm", B = 49L)

  expect_s3_class(result, "htest")
  expect_s3_class(result, "perm_wgcm")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_equal(result$parameter[["B"]], 49L)
  expect_length(result$perm_stats, 49L)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(!is.null(result$rY))
  expect_true(!is.null(result$rX))
  expect_true(!is.null(result$W))
})

test_that("perm_wgcm does not reject under strong null", {
  set.seed(456)
  n <- 150
  n_rep <- 15
  pvals <- numeric(n_rep)

  for (i in seq_len(n_rep)) {
    Z <- matrix(rnorm(n), ncol = 1)
    X <- rnorm(n)  # X independent of (Y, Z)
    Y <- 2 * Z + rnorm(n)
    pvals[i] <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                           reg_wfun = "lrm", B = 99L,
                           alternative = "greater")$p.value
  }

  reject_rate <- mean(pvals <= 0.05)
  expect_true(reject_rate <= 0.3)
})

test_that("perm_wgcm has power under alternative", {
  set.seed(789)
  n <- 200
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + 0.5 * X + rnorm(n)

  result <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      reg_wfun = "lrm", B = 199L, alternative = "greater")
  expect_true(result$p.value < 0.15)
})

test_that("perm_wgcm return_fitted_models works", {
  set.seed(101)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      reg_wfun = "lrm", B = 19L, return_fitted_models = TRUE)
  expect_true(!is.null(result$models))
})

test_that("perm_wgcm alternative options work", {
  set.seed(202)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- rnorm(n)
  Y <- rnorm(n)

  res_g <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      reg_wfun = "lrm", B = 19L, alternative = "greater")
  res_t <- perm_wgcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                      reg_wfun = "lrm", B = 19L, alternative = "two.sided")

  expect_equal(res_g$alternative, "greater")
  expect_equal(res_t$alternative, "two.sided")
})
