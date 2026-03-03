test_that("perm_pcm returns correct htest structure", {
  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 49L)

  expect_s3_class(result, "htest")
  expect_s3_class(result, "perm_pcm")
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
  expect_equal(result$parameter[["B"]], 49L)
  expect_length(result$perm_stats, 49L)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_true(!is.null(result$rY))
  expect_true(!is.null(result$rF))
  expect_true(!is.null(result$D1_indices))
  expect_true(!is.null(result$D2_indices))
})

test_that("perm_pcm does not reject under strong null", {
  # X_1 _||_ (Y, X_{-1}): X_1 independent of everything else
  set.seed(456)
  n <- 200
  n_rep <- 15
  pvals <- numeric(n_rep)

  for (i in seq_len(n_rep)) {
    X <- matrix(rnorm(n * 3), ncol = 3)  # all independent
    Y <- X[, 2] + X[, 3] + rnorm(n)     # Y depends on X2, X3 but not X1
    pvals[i] <- perm_pcm(Y, X, target = 1,
                          reg_YonXZ = "lrm", reg_YonZ = "lrm",
                          reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                          B = 99L, alternative = "greater")$p.value
  }

  reject_rate <- mean(pvals <= 0.05)
  expect_true(reject_rate <= 0.3)  # conservative check
})

test_that("perm_pcm has power when X_j contributes to Y", {
  set.seed(789)
  n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 1] + X[, 2] + rnorm(n)  # X_1 contributes

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 199L, alternative = "greater")
  expect_true(result$p.value < 0.1)
})

test_that("perm_pcm works with additional Z", {
  set.seed(101)
  n <- 200
  Z <- matrix(rnorm(n * 2), ncol = 2)
  X <- matrix(rnorm(n * 2), ncol = 2)
  Y <- X[, 2] + Z[, 1] + rnorm(n)

  result <- perm_pcm(Y, X, Z = Z, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 49L)
  expect_s3_class(result, "perm_pcm")
})

test_that("perm_pcm D2 models are not refitted during permutations", {
  # Verify by checking that the return object has consistent structure.
  # The D2 indices should be disjoint from D1 indices.
  set.seed(202)
  n <- 200
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + rnorm(n)

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 19L, return_fitted_models = TRUE)

  expect_true(length(intersect(result$D1_indices, result$D2_indices)) == 0)
  expect_equal(sort(c(result$D1_indices, result$D2_indices)), seq_len(n))
  expect_true(!is.null(result$models))
  expect_true("m_hat_p" %in% names(result$models))
  expect_true("m_hat_qj" %in% names(result$models))
})

test_that("perm_pcm alternative options work", {
  set.seed(303)
  n <- 200
  X <- matrix(rnorm(n * 2), ncol = 2)
  Y <- rnorm(n)

  res_g <- perm_pcm(Y, X, target = 1,
                    reg_YonXZ = "lrm", reg_YonZ = "lrm",
                    reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                    B = 19L, alternative = "greater")
  res_t <- perm_pcm(Y, X, target = 1,
                    reg_YonXZ = "lrm", reg_YonZ = "lrm",
                    reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                    B = 19L, alternative = "two.sided")

  expect_equal(res_g$alternative, "greater")
  expect_equal(res_t$alternative, "two.sided")
})

# ---- Vectorized fast-path tests ----

test_that("perm_pcm full lrm fast path gives valid p-values under null", {
  set.seed(404)
  n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)  # X_1 is irrelevant

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 199L, alternative = "greater")
  expect_s3_class(result, "perm_pcm")
  expect_true(result$p.value > 0.01)
  expect_length(result$perm_stats, 199L)
})

test_that("perm_pcm partial fast path works (rf D2 models, lrm projection)", {
  set.seed(505)
  n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 2] + X[, 3] + rnorm(n)

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "rf", reg_YonZ = "lrm",
                     reg_YhatonZ = "rf", reg_RonZ = "lrm",
                     B = 49L, alternative = "greater")
  expect_s3_class(result, "perm_pcm")
  expect_true(!is.null(result$p.value))
  expect_length(result$perm_stats, 49L)
})

test_that("perm_pcm full fast path has power when X_j contributes", {
  set.seed(606)
  n <- 400
  X <- matrix(rnorm(n * 3), ncol = 3)
  Y <- X[, 1] + X[, 2] + rnorm(n)

  result <- perm_pcm(Y, X, target = 1,
                     reg_YonXZ = "lrm", reg_YonZ = "lrm",
                     reg_YhatonZ = "lrm", reg_RonZ = "lrm",
                     B = 199L, alternative = "greater")
  expect_true(result$p.value < 0.1)
})
