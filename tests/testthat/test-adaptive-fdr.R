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
