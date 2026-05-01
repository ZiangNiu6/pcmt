## Tests for the adaptive FDR-control procedure (Algorithm 5)

# ---- bh_threshold ----

test_that("bh_threshold matches stats::p.adjust", {
  set.seed(1)
  for (m in c(5L, 20L, 100L)) {
    p <- runif(m)
    alpha <- 0.1
    R <- pcmt:::bh_threshold(p, alpha)
    ## p.adjust gives BH-adjusted p-values; compare rejection sets.
    rej_pcmt <- which(p <= R)
    rej_padj <- which(stats::p.adjust(p, method = "BH") <= alpha)
    expect_setequal(rej_pcmt, rej_padj)
  }
})

test_that("bh_threshold handles edge cases", {
  expect_equal(pcmt:::bh_threshold(numeric(0), 0.05), 0)
  expect_equal(pcmt:::bh_threshold(c(0.9, 0.95, 0.99), 0.05), 0)
  ## All zeros => everything rejected at threshold m * alpha / m = alpha.
  expect_equal(pcmt:::bh_threshold(rep(0, 10), 0.1), 0.1)
})

# ---- adaptive_fdr_bh: structural properties ----

test_that("adaptive_fdr_bh returns expected structure", {
  set.seed(42)
  m <- 30L
  T_orig <- c(rnorm(25), rnorm(5, mean = 4))   # 5 alternatives
  permute_stat <- function(active, t) rnorm(length(active))

  res <- adaptive_fdr_bh(T_orig, permute_stat,
                         alpha = 0.1, h = 10L, B_max = 5000L)

  expect_named(res,
               c("rejected", "futility", "active", "p_value", "K",
                 "stop_round", "t_total", "R_threshold"),
               ignore.order = TRUE)
  expect_length(res$p_value, m)
  expect_length(res$K, m)
  expect_length(res$stop_round, m)
  expect_true(all(res$p_value >= 0 & res$p_value <= 1))
  expect_true(all(res$K >= 0L & res$K <= 10L))
  ## Every hypothesis must end up rejected, futile, or still active.
  expect_setequal(c(res$rejected, res$futility, res$active), seq_len(m))
})

test_that("adaptive_fdr_bh: futility hypotheses have p = h / stop_round", {
  set.seed(7)
  m <- 40L
  T_orig <- rnorm(m)                            # all nulls
  permute_stat <- function(active, t) rnorm(length(active))

  h <- 8L
  res <- adaptive_fdr_bh(T_orig, permute_stat,
                         alpha = 0.05, h = h, B_max = 10000L)

  expect_true(length(res$futility) > 0L)
  for (i in res$futility) {
    expect_equal(res$K[i], h)
    expect_equal(res$p_value[i], h / res$stop_round[i],
                 tolerance = 1e-12,
                 info = sprintf("hypothesis %d", i))
  }
})

test_that("adaptive_fdr_bh: active p-value formula holds at termination", {
  set.seed(99)
  m <- 200L
  T_orig <- rnorm(m)
  permute_stat <- function(active, t) rnorm(length(active))

  h <- 50L
  ## Force early termination via B_max so some hypotheses are still active.
  res <- adaptive_fdr_bh(T_orig, permute_stat,
                         alpha = 0.001, h = h, B_max = 30L)

  expect_gt(length(res$active), 0L)
  for (i in res$active) {
    expect_equal(res$p_value[i], h / (res$t_total + h - res$K[i]),
                 tolerance = 1e-12)
    expect_true(res$K[i] < h)
    expect_true(is.na(res$stop_round[i]))
  }
})

test_that("adaptive_fdr_bh rejects strong signals", {
  set.seed(123)
  m <- 50L
  ## 5 strong alternatives, 45 nulls.
  T_orig <- c(rnorm(45), rep(8, 5))
  permute_stat <- function(active, t) rnorm(length(active))

  res <- adaptive_fdr_bh(T_orig, permute_stat,
                         alpha = 0.1, h = 15L, B_max = 5000L)

  ## All five strong signals should be rejected.
  expect_true(all((m - 4):m %in% res$rejected))
})

test_that("adaptive_fdr_bh: handles two.sided alternative", {
  set.seed(8)
  m <- 30L
  T_orig <- c(rnorm(25), rep(-6, 5))            # left-tail signals
  permute_stat <- function(active, t) rnorm(length(active))

  res_two <- adaptive_fdr_bh(T_orig, permute_stat,
                             alpha = 0.1, h = 12L, B_max = 5000L,
                             alternative = "two.sided")
  ## Strong negative signals should be picked up by two-sided test.
  expect_true(all((m - 4):m %in% res_two$rejected))

  res_less <- adaptive_fdr_bh(T_orig, permute_stat,
                              alpha = 0.1, h = 12L, B_max = 5000L,
                              alternative = "less")
  expect_true(all((m - 4):m %in% res_less$rejected))
})

# ---- FDR control: synthetic Gaussian z-test simulation ----

test_that("adaptive_fdr_bh empirically controls FDR on synthetic data", {
  skip_on_cran()
  set.seed(2025)
  alpha <- 0.1
  m <- 100L
  m1 <- 20L                                     # alternatives
  mu_alt <- 4
  reps <- 80L

  fdp <- numeric(reps)
  power <- numeric(reps)
  for (r in seq_len(reps)) {
    T_orig <- c(rnorm(m - m1), rnorm(m1, mean = mu_alt))
    truth <- c(rep(FALSE, m - m1), rep(TRUE, m1))
    permute_stat <- function(active, t) rnorm(length(active))
    res <- adaptive_fdr_bh(T_orig, permute_stat,
                           alpha = alpha, h = 15L, B_max = 5000L)
    R <- length(res$rejected)
    V <- sum(!truth[res$rejected])
    fdp[r] <- if (R == 0L) 0 else V / R
    power[r] <- sum(truth[res$rejected]) / m1
  }

  ## Empirical FDR should respect the level (with margin for Monte Carlo).
  expect_lt(mean(fdp), alpha + 0.05)
  ## And we should have meaningful power on the strong alternatives.
  expect_gt(mean(power), 0.5)
})

# ---- perm_gcm_fdr wrapper ----

test_that("perm_gcm_fdr runs and recovers signals (lrm)", {
  skip_on_cran()
  set.seed(11)
  n <- 200L
  m <- 30L
  Z <- matrix(rnorm(n * 2L), ncol = 2L)
  X <- Z[, 1] + rnorm(n)
  Y <- matrix(0, nrow = n, ncol = m)
  for (j in seq_len(m - 5L)) Y[, j] <- Z[, 1] + rnorm(n)              # null
  for (j in (m - 4L):m)        Y[, j] <- Z[, 1] + 0.6 * X + rnorm(n)   # alt

  res <- perm_gcm_fdr(Y, X, Z,
                      reg_YonZ = "lrm", reg_XonZ = "lrm",
                      alpha = 0.1, h = 15L, B_max = 2000L)

  expect_length(res$T_orig, m)
  expect_equal(dim(res$rY), c(n, m))
  expect_length(res$rX, n)
  ## Most of the strong signals should be rejected.
  expect_gte(sum((m - 4L):m %in% res$rejected), 4L)
  ## Empirical FDP on this single replicate should be modest.
  R <- length(res$rejected)
  V <- sum(res$rejected <= m - 5L)
  fdp <- if (R == 0L) 0 else V / R
  expect_lt(fdp, 0.3)
})
