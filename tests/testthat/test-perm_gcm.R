test_that("RF GCM matches genuine nuisance refitting for each alternative", {
  args <- paper_gcm_args()
  for (alternative in c("two.sided", "greater", "less")) {
    set.seed(701L)
    actual <- do.call(perm_gcm, c(args, list(B = 7L, alternative = alternative)))
    set.seed(701L)
    expected <- paper_gcm_reference(args, B = 7L, alternative = alternative)
    expect_s3_class(actual, "perm_gcm")
    expect_s3_class(actual, "htest")
    expect_equal(unname(actual$statistic), expected$statistic, tolerance = 1e-10)
    expect_equal(actual$perm_stats, expected$perm_stats, tolerance = 1e-10)
    expect_equal(actual$p.value, expected$p.value)
    expect_identical(actual$alternative, alternative)
    expect_equal(unname(actual$parameter[["B"]]), 7L)
    expect_null(actual$adaptive)
  }
})

test_that("RF GCM adaptive calibration uses the same ordered refits", {
  args <- paper_gcm_args()
  for (alternative in c("two.sided", "greater")) {
    set.seed(702L)
    fixed <- do.call(perm_gcm, c(args, list(B = 11L, alternative = alternative)))
    set.seed(702L)
    adaptive <- do.call(perm_gcm, c(args, list(B = 11L, h = 2L,
                        alpha = 0.2, alternative = alternative)))
    expect_adaptive_contract(adaptive, fixed, h = 2L, alpha = 0.2, B = 11L)
    expect_equal(adaptive$rY, fixed$rY)
    expect_equal(adaptive$rX, fixed$rX)
  }
})

test_that("RF GCM retains optional fitted models", {
  fit <- do.call(perm_gcm, c(paper_gcm_args(),
                 list(B = 2L, return_fitted_models = TRUE)))
  expect_named(fit$models, c("reg_YonZ", "reg_XonZ"))
  expect_s3_class(fit$models$reg_YonZ, "rf")
  expect_s3_class(fit$models$reg_XonZ, "rf")
})
