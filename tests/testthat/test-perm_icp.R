test_that("lmICP identifies correct causal parent", {
  set.seed(42)
  n <- 300
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  Y <- 2 * X1 + rnorm(n)
  df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)

  result <- lmICP(Y ~ X1 + X2, data = df, env = ~ E, B = 199)

  expect_s3_class(result, "perm_icp")
  expect_true("X1" %in% result$causal_set)
  expect_true(result$pvalues_predictor["X1"] < 0.05)
  expect_true(result$pvalues_predictor["X2"] > 0.05)
})

test_that("glmICP works with binomial family", {
  set.seed(42)
  n <- 400
  E <- rbinom(n, 1, 0.5)
  X1 <- -E + rnorm(n)
  X2 <- rnorm(n)
  Y <- as.integer(0.5 * X1 > rnorm(n))
  df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)

  result <- glmICP(Y ~ X1 + X2, data = df, env = ~ E,
                   family = binomial, B = 199)

  expect_s3_class(result, "perm_icp")
  expect_equal(result$n_sets_tested, 4L)
  # X1 should be identified as causal
  expect_true(result$pvalues_predictor["X1"] < 0.05)
})

test_that("coxphICP identifies causal parent in survival data", {
  skip_if_not_installed("survival")
  set.seed(42)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  time <- rexp(n, exp(0.5 * X1))
  status <- rep(1L, n)
  df <- data.frame(time = time, status = status, X1 = X1, X2 = X2, E = E)

  result <- coxphICP(survival::Surv(time, status) ~ X1 + X2,
                     data = df, env = ~ E, B = 199)

  expect_s3_class(result, "perm_icp")
  expect_true("X1" %in% result$causal_set)
})

test_that("mandatory predictors are always included", {
  set.seed(42)
  n <- 300
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- 2 * X1 + X2 + rnorm(n)
  df <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E)

  result <- lmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E,
                  mandatory = ~ X1, B = 199)

  # With X1 mandatory, only supersets of {X1} are tested
  # Should test: {X1}, {X1,X2}, {X1,X3}, {X1,X2,X3} = 4 subsets
  expect_equal(result$n_sets_tested, 4L)
})

test_that("max_size limits subset enumeration", {
  set.seed(42)
  n <- 200
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- X1 + rnorm(n)
  df <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E)

  result <- lmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E,
                  max_size = 1, B = 99)

  # max_size=1: Empty, {X1}, {X2}, {X3} = 4 subsets
  expect_equal(result$n_sets_tested, 4L)
})

test_that("pvalues method works", {
  set.seed(42)
  n <- 200
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  Y <- X1 + rnorm(n)
  df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)

  result <- lmICP(Y ~ X1 + X2, data = df, env = ~ E, B = 99)

  p_pred <- pvalues(result, which = "predictor")
  p_set <- pvalues(result, which = "set")

  expect_named(p_pred, c("X1", "X2"))
  expect_length(p_set, 4L)
})

test_that("print.perm_icp works without error", {
  set.seed(42)
  n <- 200
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  Y <- X1 + rnorm(n)
  df <- data.frame(Y = Y, X1 = X1, E = E)

  result <- lmICP(Y ~ X1, data = df, env = ~ E, B = 99)
  expect_output(print(result), "Permutation-Calibrated")
})

test_that("perm_tram_gcm returns correct structure", {
  set.seed(42)
  n <- 100
  R <- rnorm(n)
  xi <- rnorm(n)

  result <- pcmt:::perm_tram_gcm(R, xi, B = 99)

  expect_s3_class(result, "htest")
  expect_true("p.value" %in% names(result))
  expect_length(result$perm_stats, 99L)
})

test_that("perm_tram_gcm with multivariate environment", {
  set.seed(42)
  n <- 100
  R <- rnorm(n)
  xi <- matrix(rnorm(n * 2), ncol = 2)

  result <- pcmt:::perm_tram_gcm(R, xi, B = 99)

  expect_s3_class(result, "htest")
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("compute_score_residuals works for lm", {
  m <- lm(mpg ~ wt + hp, data = mtcars)
  sr <- pcmt:::compute_score_residuals(m)
  expect_length(sr, nrow(mtcars))
})

test_that("compute_score_residuals works for glm", {
  m <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  sr <- pcmt:::compute_score_residuals(m)
  expect_length(sr, nrow(mtcars))
})

test_that("residualize_environment with lrm method", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  E <- X %*% c(1, 2) + rnorm(n)

  xi <- pcmt:::residualize_environment(E, X, method = "lrm")

  # Residuals should be uncorrelated with X
  expect_true(abs(cor(xi, X[, 1])) < 0.1)
  expect_true(abs(cor(xi, X[, 2])) < 0.1)
})
