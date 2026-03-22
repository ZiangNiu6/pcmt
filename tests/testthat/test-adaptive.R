## Tests for adaptive permutation calibration

# ---- perm_gcm adaptive ----

test_that("adaptive perm_gcm returns correct structure", {
  set.seed(42)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 199L, h = 15L, alpha = 0.05)

  expect_s3_class(result, "htest")
  expect_s3_class(result, "perm_gcm")
  expect_true(!is.null(result$adaptive))
  expect_equal(result$adaptive$h, 15L)
  expect_equal(result$adaptive$alpha, 0.05)
  expect_true(result$adaptive$K <= 15L)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(result$adaptive$stop_reason %in% c("rejection", "futility", "cap"))

  ## Trajectory is a data frame with correct columns

  traj <- result$adaptive$trajectory
  expect_true(is.data.frame(traj))
  expect_true(all(c("t", "S", "K", "p_AV") %in% names(traj)))
  expect_equal(nrow(traj), length(result$perm_stats))

  ## B parameter matches actual number of permutations
  expect_equal(result$parameter[["B"]], length(result$perm_stats))
})

test_that("adaptive perm_gcm stops for futility under the null", {
  set.seed(100)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)  # X _||_ Y | Z

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 999L, h = 10L, alpha = 0.05)

  ## Under null, likely stops for futility with fewer than B permutations
  if (result$adaptive$stop_reason == "futility") {
    expect_true(result$adaptive$stopped_early)
    expect_true(length(result$perm_stats) < 999L)
    expect_equal(result$adaptive$K, 10L)
    expect_true(result$p.value > 0.05)
  }
})

test_that("adaptive perm_gcm can reject under alternative", {
  set.seed(200)
  n <- 200
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n, sd = 0.5)  # strong signal: Y depends on X given Z

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 999L, h = 15L, alpha = 0.05,
                     alternative = "greater")

  ## With strong alternative, should reject early
  expect_true(result$p.value <= 0.05)
  if (result$adaptive$stop_reason == "rejection") {
    expect_true(result$adaptive$stopped_early)
  }
})

test_that("adaptive perm_gcm with h=NULL behaves like standard", {
  set.seed(300)
  n <- 80
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 49L, h = NULL)

  expect_null(result$adaptive)
  expect_equal(result$parameter[["B"]], 49L)
  expect_length(result$perm_stats, 49L)
})

test_that("adaptive perm_gcm works with two.sided alternative", {
  set.seed(400)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 199L, h = 10L, alpha = 0.05,
                     alternative = "two.sided")

  expect_equal(result$alternative, "two.sided")
  expect_true(!is.null(result$adaptive))
  ## Scores should be non-negative (absolute values)
  expect_true(all(result$adaptive$trajectory$S >= 0))
})

# ---- perm_pcm adaptive ----

test_that("adaptive perm_pcm returns correct structure", {
  set.seed(42)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 199L, h = 15L, alpha = 0.05)

  expect_s3_class(result, "perm_pcm")
  expect_true(!is.null(result$adaptive))
  expect_equal(result$adaptive$h, 15L)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(result$adaptive$stop_reason %in% c("rejection", "futility", "cap"))
  expect_equal(result$parameter[["B"]], length(result$perm_stats))
})

test_that("adaptive perm_pcm stops for futility under the null", {
  set.seed(101)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)  # Y _||_ X_1 | X_{-1}

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 999L, h = 10L, alpha = 0.05)

  if (result$adaptive$stop_reason == "futility") {
    expect_true(result$adaptive$stopped_early)
    expect_true(result$p.value > 0.05)
  }
})

# ---- perm_wgcm adaptive ----

test_that("adaptive perm_wgcm returns correct structure", {
  set.seed(42)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_wgcm(Y, X, Z,
                      reg_YonZ = "lrm", reg_XonZ = "lrm",
                      reg_wfun = "lrm",
                      B = 199L, h = 15L, alpha = 0.05)

  expect_s3_class(result, "perm_wgcm")
  expect_true(!is.null(result$adaptive))
  expect_equal(result$adaptive$h, 15L)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(result$adaptive$stop_reason %in% c("rejection", "futility", "cap"))
  expect_equal(result$parameter[["B"]], length(result$perm_stats))
})

# ---- Utility function tests ----

test_that("score_map computes correctly", {
  expect_equal(score_map(3.5, "greater"), 3.5)
  expect_equal(score_map(3.5, "less"), -3.5)
  expect_equal(score_map(-2.0, "two.sided"), 2.0)
  expect_equal(score_map(0, "two.sided"), 0)
})

test_that("adaptive_pvalue formula is correct", {
  ## Before h-th loss: p = h / (t + h - K)
  expect_equal(adaptive_pvalue(t = 5, K = 2, h = 10, t_star = Inf),
               10 / (5 + 10 - 2))
  ## After h-th loss: p = h / t_star
  expect_equal(adaptive_pvalue(t = 20, K = 10, h = 10, t_star = 15),
               10 / 15)
  ## At boundary: K = h for the first time
  expect_equal(adaptive_pvalue(t = 15, K = 10, h = 10, t_star = 15),
               10 / 15)
})

test_that("adaptive perm_gcm validates h parameter", {
  n <- 50
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  expect_error(
    perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
             B = 49L, h = 0L),
    "h >= 1"
  )
  expect_error(
    perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
             B = 49L, h = 10L, alpha = 0),
    "alpha > 0"
  )
})

test_that("adaptive trajectory K is monotonically non-decreasing", {
  set.seed(500)
  n <- 100
  Z <- matrix(rnorm(n), ncol = 1)
  X <- Z + rnorm(n)
  Y <- Z + rnorm(n)

  result <- perm_gcm(Y, X, Z, reg_YonZ = "lrm", reg_XonZ = "lrm",
                     B = 199L, h = 20L, alpha = 0.05)

  traj <- result$adaptive$trajectory
  expect_true(all(diff(traj$K) >= 0))
})
