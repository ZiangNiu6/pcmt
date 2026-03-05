# Validation tests: compare pcmt::perm_icp against tramicp
#
# These tests verify that our permutation-calibrated ICP implementation
# produces results consistent with the tramicp package. We check:
# 1. Score residuals match between implementations
# 2. Both identify the correct causal predictor on clear-signal problems
# 3. Both reject the same non-invariant sets
# 4. Correct type I error control

skip_if_not_installed("tramicp")
skip_if_not_installed("survival")

# ============================================================================
# Test 1: GLM score residuals match tramicp (deviance residuals)
# ============================================================================
test_that("GLM score residuals match tramicp (deviance residuals)", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- rbinom(n, 1, plogis(X1))
  df <- data.frame(Y = Y, X1 = X1, X2 = X2)

  m <- glm(Y ~ X1 + X2, data = df, family = binomial)

  sr_pcmt <- pcmt:::compute_score_residuals(m)
  sr_default <- residuals(m)  # deviance by default, same as tramicp

  expect_equal(sr_pcmt, c(sr_default), tolerance = 1e-10)
})

# ============================================================================
# Test 2: Cox PH score residuals = martingale residuals
# ============================================================================
test_that("Cox PH score residuals are martingale residuals", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  time <- rexp(n, exp(0.5 * X1))
  status <- rep(1L, n)
  df <- data.frame(time = time, status = status, X1 = X1)

  m <- survival::coxph(survival::Surv(time, status) ~ X1, data = df)

  sr_pcmt <- pcmt:::compute_score_residuals(m)
  sr_martingale <- residuals(m, type = "martingale")

  expect_equal(sr_pcmt, c(sr_martingale), tolerance = 1e-10)
})

# ============================================================================
# Test 3: lm score residuals = ordinary residuals
# ============================================================================
test_that("lm score residuals are ordinary residuals", {
  m <- lm(mpg ~ wt + hp, data = mtcars)

  sr_pcmt <- pcmt:::compute_score_residuals(m)
  sr_lm <- residuals(m)

  expect_equal(sr_pcmt, c(sr_lm), tolerance = 1e-10)
})

# ============================================================================
# Test 4: lmICP agrees with tramicp on strong-signal linear model
# ============================================================================
test_that("lmICP: pcmt and tramicp agree on causal set (linear model)", {
  set.seed(456)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- 2 * X1 + rnorm(n)

  df <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E)

  # tramicp
  res_tramicp <- tramicp::lmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E,
                                 type = "residual", test = "gcm.test",
                                 verbose = FALSE)
  set_tramicp <- res_tramicp$candidate_causal_predictors
  p_tramicp <- tramicp::pvalues(res_tramicp, which = "predictor")

  # pcmt — use B=999 for stability
  set.seed(1)
  res_pcmt <- pcmt::lmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E, B = 999)
  set_pcmt <- res_pcmt$causal_set
  p_pcmt <- pcmt::pvalues(res_pcmt, which = "predictor")

  # Both should identify X1 as causal
  expect_true("X1" %in% set_tramicp)
  expect_true("X1" %in% set_pcmt)

  # X1 should have small predictor p-value in both
  expect_true(p_tramicp["X1"] < 0.05)
  expect_true(p_pcmt["X1"] < 0.05)
})

# ============================================================================
# Test 5: glmICP agrees with tramicp (binary logistic, strong signal)
# ============================================================================
test_that("glmICP: pcmt and tramicp agree on causal predictor", {
  set.seed(42)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- -E + rnorm(n)
  X2 <- rnorm(n)
  Y <- as.integer(0.5 * X1 > rnorm(n))
  X2_obs <- Y + 0.8 * E + rnorm(n)

  df <- data.frame(Y = Y, X1 = X1, X2 = X2_obs, E = E)

  # tramicp
  res_tramicp <- tramicp::glmICP(Y ~ X1 + X2, data = df, env = ~ E,
                                  family = binomial, type = "residual",
                                  test = "gcm.test", verbose = FALSE)
  p_tramicp <- tramicp::pvalues(res_tramicp, which = "predictor")

  # pcmt
  set.seed(1)
  res_pcmt <- pcmt::glmICP(Y ~ X1 + X2, data = df, env = ~ E,
                            family = binomial, B = 999)
  p_pcmt <- pcmt::pvalues(res_pcmt, which = "predictor")

  # X1 should have small predictor p-value in both
  expect_true(p_tramicp["X1"] < 0.05)
  expect_true(p_pcmt["X1"] < 0.05)

  # Both should identify X1 as causal
  expect_true("X1" %in% res_tramicp$candidate_causal_predictors)
  expect_true("X1" %in% res_pcmt$causal_set)
})

# ============================================================================
# Test 6: coxphICP correctly identifies causal predictor
# (tramicp::coxphICP has issues with Surv formula parsing, so we only
# validate pcmt's coxphICP against known ground truth)
# ============================================================================
test_that("coxphICP: pcmt correctly identifies causal predictor", {
  set.seed(789)
  n <- 800
  E <- rbinom(n, 1, 0.5)
  X1 <- 1.5 * E + rnorm(n)
  X2 <- rnorm(n)
  time <- rexp(n, exp(0.5 * X1))
  status <- rep(1L, n)

  df <- data.frame(time = time, status = status, X1 = X1, X2 = X2, E = E)

  set.seed(1)
  res_pcmt <- pcmt::coxphICP(
    survival::Surv(time, status) ~ X1 + X2,
    data = df, env = ~ E, B = 999)
  p_pcmt <- pcmt::pvalues(res_pcmt, which = "predictor")

  # X1 should be identified as causal
  expect_true(p_pcmt["X1"] < 0.05)
  expect_true("X1" %in% res_pcmt$causal_set)
})

# ============================================================================
# Test 7: Null case - both maintain correct level
# ============================================================================
test_that("Null case: full invariant set is not rejected", {
  set.seed(321)
  n <- 300
  E <- rbinom(n, 1, 0.5)
  X1 <- rnorm(n)  # independent of E
  X2 <- rnorm(n)
  Y <- X1 + X2 + rnorm(n)

  df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)

  # Both should NOT reject the full set {X1, X2}
  set.seed(1)
  res_pcmt <- pcmt::lmICP(Y ~ X1 + X2, data = df, env = ~ E, B = 499)
  full_set_pval <- res_pcmt$pvalues_set["X1+X2"]
  expect_true(full_set_pval > 0.05)

  res_tramicp <- tramicp::lmICP(Y ~ X1 + X2, data = df, env = ~ E,
                                 type = "residual", test = "gcm.test",
                                 verbose = FALSE)
  p_t_set <- tramicp::pvalues(res_tramicp, which = "set")
  expect_true(p_t_set["X1+X2"] > 0.05)
})

# ============================================================================
# Test 8: Predictor p-value ordering is consistent
# ============================================================================
test_that("Predictor p-values: causal < non-causal in both", {
  set.seed(111)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  Y <- 2 * X1 + rnorm(n)

  df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)

  # tramicp
  res_tramicp <- tramicp::lmICP(Y ~ X1 + X2, data = df, env = ~ E,
                                 type = "residual", test = "gcm.test",
                                 verbose = FALSE)
  p_tramicp <- tramicp::pvalues(res_tramicp, which = "predictor")

  # pcmt
  set.seed(1)
  res_pcmt <- pcmt::lmICP(Y ~ X1 + X2, data = df, env = ~ E, B = 999)
  p_pcmt <- pcmt::pvalues(res_pcmt, which = "predictor")

  # X1 should have smaller p-value than X2 in both
  expect_true(p_tramicp["X1"] < p_tramicp["X2"])
  expect_true(p_pcmt["X1"] < p_pcmt["X2"])
})

# ============================================================================
# Test 9: Set-specific p-values agree on invariant vs non-invariant
# ============================================================================
test_that("Set p-values: invariant set not rejected, non-invariant rejected", {
  set.seed(222)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  Y <- 2 * X1 + rnorm(n)

  df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)

  # tramicp
  res_tramicp <- tramicp::lmICP(Y ~ X1 + X2, data = df, env = ~ E,
                                 type = "residual", test = "gcm.test",
                                 verbose = FALSE)
  p_set_t <- tramicp::pvalues(res_tramicp, which = "set")

  # pcmt
  set.seed(1)
  res_pcmt <- pcmt::lmICP(Y ~ X1 + X2, data = df, env = ~ E, B = 999)
  p_set_p <- pcmt::pvalues(res_pcmt, which = "set")

  # Both: {X1} should not be rejected
  expect_true(p_set_t["X1"] > 0.05)
  expect_true(p_set_p["X1"] > 0.05)

  # Both: Empty and {X2} should be rejected
  expect_true(p_set_t["Empty"] < 0.05)
  expect_true(p_set_p["Empty"] < 0.05)
  expect_true(p_set_t["X2"] < 0.05)
  expect_true(p_set_p["X2"] < 0.05)
})

# ============================================================================
# Test 10: Residualized environments via lrm match lm residuals
# ============================================================================
test_that("Residualized environments: lrm matches lm residuals", {
  set.seed(444)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  E <- X %*% c(1, 0.5, -1) + rnorm(n)

  xi_pcmt <- pcmt:::residualize_environment(E, X, method = "lrm")
  xi_lm <- residuals(lm(c(E) ~ X))

  expect_equal(c(xi_pcmt), unname(c(xi_lm)), tolerance = 1e-10)
})

# ============================================================================
# Test 11: Complex structure - descendant not identified as causal
# ============================================================================
test_that("Descendant X3 not in causal set for both methods", {
  set.seed(555)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- -E + rnorm(n)
  X2 <- rnorm(n)
  Y <- as.integer(X1 + 0.5 * X2 > rnorm(n))
  X3 <- Y + 0.5 * E + rnorm(n)  # descendant

  df <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E)

  # tramicp
  res_t <- tramicp::glmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E,
                            family = binomial, type = "residual",
                            test = "gcm.test", verbose = FALSE)

  # pcmt
  set.seed(1)
  res_p <- pcmt::glmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E,
                         family = binomial, B = 999)

  # X3 should NOT be in the causal set for either
  expect_false("X3" %in% res_t$candidate_causal_predictors)
  expect_false("X3" %in% res_p$causal_set)
})

# ============================================================================
# Test 12: mandatory argument works consistently
# ============================================================================
test_that("mandatory argument produces consistent results", {
  set.seed(333)
  n <- 500
  E <- rbinom(n, 1, 0.5)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  Y <- 2 * X1 + X2 + rnorm(n)

  df <- data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E)

  # pcmt with mandatory X1
  set.seed(1)
  res_pcmt <- pcmt::lmICP(Y ~ X1 + X2 + X3, data = df, env = ~ E,
                           mandatory = ~ X1, B = 999)

  # X1 should always be in the causal set (it's mandatory and causal)
  expect_true("X1" %in% res_pcmt$causal_set)
  expect_equal(res_pcmt$n_sets_tested, 4L)
})
