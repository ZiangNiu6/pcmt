# ---- nslrm fit and S3 methods ----

test_that("fit_nslrm returns correct structure", {
  set.seed(42)
  n <- 100
  x <- matrix(rnorm(n * 2), ncol = 2)
  y <- x[, 1] + rnorm(n)

  fit <- pcmt:::fit_nslrm(y, x, df = 4)
  expect_s3_class(fit, "nslrm_fit")
  expect_length(fit$bases, 2)
  expect_length(fit$fitted, n)
})

test_that("predict.nslrm_fit works on new data", {
  set.seed(43)
  n <- 100
  x <- matrix(rnorm(n * 2), ncol = 2)
  y <- x[, 1] + rnorm(n)

  fit <- pcmt:::fit_nslrm(y, x, df = 4)
  preds <- predict(fit, data = x)
  expect_length(preds, n)
  expect_equal(preds, fit$fitted, tolerance = 1e-10)
})

test_that("residuals.nslrm_fit returns correct residuals", {
  set.seed(44)
  n <- 100
  x <- matrix(rnorm(n * 2), ncol = 2)
  y <- x[, 1] + rnorm(n)

  fit <- pcmt:::fit_nslrm(y, x, df = 4)
  resids <- residuals(fit, response = y, data = x)
  expect_length(resids, n)
  expect_equal(resids, y - fit$fitted, tolerance = 1e-10)
})

# ---- perm_gcm with nslrm ----

test_that("perm_gcm nslrm fast path returns valid htest", {
  set.seed(100)
  n <- 150
  Z <- matrix(rnorm(n * 2), ncol = 2)
  X <- sin(Z[, 1]) + rnorm(n)
  Y <- cos(Z[, 2]) + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "nslrm", reg_XonZ = "nslrm", B = 49L)
  expect_s3_class(result, "perm_gcm")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_length(result$perm_stats, 49L)
})

test_that("perm_gcm nslrm does not reject under strong null", {
  set.seed(200)
  n <- 150
  n_rep <- 20
  pvals <- numeric(n_rep)

  for (i in seq_len(n_rep)) {
    Z <- matrix(rnorm(n * 2), ncol = 2)
    X <- rnorm(n)
    Y <- sin(Z[, 1]) + rnorm(n)
    pvals[i] <- perm_gcm(Y, X, Z, reg_YonZ = "nslrm", reg_XonZ = "nslrm",
                          B = 99L, alternative = "greater")$p.value
  }
  reject_rate <- mean(pvals <= 0.05)
  expect_true(reject_rate <= 0.3)
})

test_that("perm_gcm nslrm has power under alternative", {
  set.seed(300)
  n <- 200
  Z <- matrix(rnorm(n), ncol = 1)
  X <- sin(Z) + rnorm(n, sd = 0.5)
  Y <- Z + 0.5 * X + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "nslrm", reg_XonZ = "nslrm",
                     B = 199L, alternative = "greater")
  expect_true(result$p.value < 0.1)
})

test_that("perm_gcm nslrm works with multivariate X", {
  set.seed(400)
  n <- 120
  Z <- matrix(rnorm(n * 2), ncol = 2)
  X <- cbind(rnorm(n), rnorm(n))
  Y <- sin(Z[, 1]) + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "nslrm", reg_XonZ = "nslrm", B = 49L)
  expect_s3_class(result, "perm_gcm")
  expect_true(result$p.value > 0.01)
})

test_that("perm_gcm nslrm respects custom df", {
  set.seed(500)
  n <- 150
  Z <- matrix(rnorm(n), ncol = 1)
  X <- rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "nslrm", reg_XonZ = "nslrm",
                     args_XonZ = list(df = 6), B = 49L)
  expect_s3_class(result, "perm_gcm")
})

# ---- perm_wgcm with nslrm ----

test_that("perm_wgcm nslrm fast path returns valid htest", {
  set.seed(600)
  n <- 200
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_wgcm(Y, X, Z, reg_YonZ = "nslrm", reg_XonZ = "nslrm",
                      reg_wfun = "lrm", B = 49L)
  expect_s3_class(result, "perm_wgcm")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

# ---- perm_pcm with nslrm ----

test_that("perm_pcm nslrm partial fast path returns valid htest", {
  set.seed(700)
  n <- 400
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)

  result <- perm_pcm(Y, X, target = 1, reg_YonXZ = "lrm",
                     reg_YonZ = "lrm", reg_YhatonZ = "lrm",
                     reg_RonZ = "nslrm", B = 49L)
  expect_s3_class(result, "perm_pcm")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})
