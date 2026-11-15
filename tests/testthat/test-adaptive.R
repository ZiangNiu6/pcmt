test_that("score maps and anytime p-values retain their mathematical formulas", {
  expect_equal(pcmt:::score_map(c(-2, 3), "greater"), c(-2, 3))
  expect_equal(pcmt:::score_map(c(-2, 3), "less"), c(2, -3))
  expect_equal(pcmt:::score_map(c(-2, 3), "two.sided"), c(2, 3))
  expect_equal(pcmt:::adaptive_pvalue(5, 2, 10, Inf), 10 / 13)
  expect_equal(pcmt:::adaptive_pvalue(20, 10, 10, 15), 10 / 15)
  products <- c(-0.4, 0.5, 2, 1.3)
  expect_equal(pcmt:::studentized_stat(products), paper_stat(products))
  expect_equal(pcmt:::studentized_stat(rep(1, 5)), 0)
})

test_that("plus-one fixed permutation p-values count ties conservatively", {
  permuted <- c(-2, -1, 0, 1, 2)
  expect_equal(pcmt:::perm_pvalue(2, permuted, "two.sided"), 3 / 6)
  expect_equal(pcmt:::perm_pvalue(2, permuted, "greater"), 2 / 6)
  expect_equal(pcmt:::perm_pvalue(2, permuted, "less"), 1)
})

test_that("adaptive FDR callback receives only active hypotheses", {
  seen <- list()
  callback <- function(active_idx, t) {
    seen[[t]] <<- active_idx
    ifelse(active_idx == 1L, 100, -100)
  }
  fit <- adaptive_fdr_bh(c(0, 10), callback,
                        alpha = 0.1, h = 2L, B_max = 50L)
  expect_identical(seen[[1L]], c(1L, 2L))
  expect_equal(fit$futility, 1L)
  expect_equal(fit$rejected, 2L)
  expect_true(all(vapply(seen[-seq_len(2L)], function(active) {
    !1L %in% active
  }, logical(1))))
  expect_length(fit$active, 0L)
  expect_equal(fit$p_value[1L], 1)
  expect_equal(fit$K[1L], 2L)
  expect_equal(fit$stop_round[1L], 2L)
})
