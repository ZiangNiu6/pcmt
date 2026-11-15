test_that("RF PCM matches full refits on the supplied ordered split", {
  args <- paper_pcm_args()
  for (alternative in c("two.sided", "greater", "less")) {
    set.seed(801L)
    actual <- do.call(perm_pcm, c(args, list(B = 7L, alternative = alternative)))
    set.seed(801L)
    expected <- paper_pcm_reference(args, B = 7L, alternative = alternative)
    expect_s3_class(actual, "perm_pcm")
    expect_s3_class(actual, "htest")
    expect_equal(unname(actual$statistic), expected$statistic, tolerance = 1e-10)
    expect_equal(actual$perm_stats, expected$perm_stats, tolerance = 1e-10)
    expect_equal(actual$p.value, expected$p.value)
    expect_identical(actual$D2_indices, args$indices)
    expect_identical(actual$D1_indices, setdiff(seq_along(args$Y), args$indices))
    expect_identical(actual$alternative, alternative)
    expect_equal(unname(actual$parameter[["B"]]), 7L)
    expect_null(actual$adaptive)
  }
})

test_that("RF PCM adaptive calibration uses the same ordered refits", {
  args <- paper_pcm_args()
  for (alternative in c("two.sided", "greater")) {
    set.seed(802L)
    fixed <- do.call(perm_pcm, c(args, list(B = 11L, alternative = alternative)))
    set.seed(802L)
    adaptive <- do.call(perm_pcm, c(args, list(B = 11L, h = 2L,
                        alpha = 0.2, alternative = alternative)))
    expect_adaptive_contract(adaptive, fixed, h = 2L, alpha = 0.2, B = 11L)
    expect_identical(adaptive$D2_indices, args$indices)
    expect_equal(adaptive$rY, fixed$rY)
    expect_equal(adaptive$rF, fixed$rF)
  }
})

test_that("PCM frac controls D2 and explicit indices override it", {
  args <- paper_pcm_args()
  args$indices <- NULL
  set.seed(803L)
  expected <- sample.int(length(args$Y), ceiling(0.4 * length(args$Y)))
  set.seed(803L)
  fit <- do.call(perm_pcm, c(args, list(B = 2L, frac = 0.4)))
  expect_identical(fit$D2_indices, expected)
  expect_length(fit$D1_indices, length(args$Y) - length(expected))

  args$indices <- rev(expected)
  fit <- do.call(perm_pcm, c(args, list(B = 2L, frac = 0.75)))
  expect_identical(fit$D2_indices, rev(expected))
})

test_that("RF PCM retains the four fitted nuisance models", {
  fit <- do.call(perm_pcm, c(paper_pcm_args(),
                 list(B = 2L, return_fitted_models = TRUE)))
  expect_named(fit$models, c("m_hat_p", "m_hat_qj", "m_hat_rj", "m_hat_pfj"))
  expect_true(all(vapply(fit$models, inherits, logical(1), what = "rf")))
})
