test_that("the publication API contains only the three used functions", {
  expect_setequal(getNamespaceExports("pcmt"),
                  c("perm_gcm", "perm_pcm", "adaptive_fdr_bh"))
  expect_true(all(c("args_YonZ", "args_XonZ", "B", "n_cores", "h", "alpha",
                    "alternative") %in% names(formals(perm_gcm))))
  expect_true(all(c("indices", "frac", "target", "args_YonXZ", "args_RonZ") %in%
                    names(formals(perm_pcm))))
})

test_that("permutation budgets, stopping rules and cores reject invalid inputs", {
  cases <- list(gcm = list(fun = perm_gcm, args = paper_gcm_args()),
                pcm = list(fun = perm_pcm, args = paper_pcm_args()))
  for (case in cases) {
    for (invalid in list(0, -1, 1.5, NA_real_, Inf, numeric(0), c(1, 2))) {
      args <- case$args; args$B <- invalid
      expect_error(do.call(case$fun, args), "B|positive integer")
      args <- case$args; args$B <- 2L; args$h <- invalid
      expect_error(do.call(case$fun, args), "h|positive integer")
    }
    for (invalid in list(0, 2, 1.5, NA_real_, numeric(0), c(1, 2))) {
      args <- case$args; args$B <- 2L; args$n_cores <- invalid
      expect_error(do.call(case$fun, args), "n_cores")
    }
    for (invalid in list(0, 1, -0.1, NA_real_, Inf)) {
      args <- case$args; args$B <- 2L; args$h <- 2L; args$alpha <- invalid
      expect_error(do.call(case$fun, args), "alpha")
    }
  }
})

test_that("removed regression learners are explicitly rejected", {
  for (method in c("lrm", "nslrm", "gols", "gam", "krr", "lasso", "glrm")) {
    for (field in c("reg_YonZ", "reg_XonZ")) {
      args <- paper_gcm_args(); args$B <- 2L; args[[field]] <- method
      expect_error(do.call(perm_gcm, args), "rf|regression")
    }
    for (field in c("reg_YonXZ", "reg_YonZ", "reg_YhatonZ", "reg_RonZ")) {
      args <- paper_pcm_args(); args$B <- 2L; args[[field]] <- method
      expect_error(do.call(perm_pcm, args), "rf|regression")
    }
  }
})

test_that("GCM and PCM reject inconsistent or nonfinite data", {
  for (fun in list(perm_gcm, perm_pcm)) {
    dat <- paper_data()
    args <- if (identical(fun, perm_gcm)) paper_gcm_args(dat) else paper_pcm_args(dat)
    args$B <- 2L
    bad <- args; bad$Y <- bad$Y[-1L]
    expect_error(do.call(fun, bad), "Y|rows|length")
    bad <- args; bad$Z <- bad$Z[-1L, , drop = FALSE]
    expect_error(do.call(fun, bad), "Z|rows|length")
    bad <- args; bad$Y[1L] <- NA_real_
    expect_error(do.call(fun, bad), "Y|finite")
    bad <- args; bad$X[1L, 1L] <- Inf
    expect_error(do.call(fun, bad), "X|finite")
    bad <- args; bad$Y <- cbind(bad$Y, bad$Y)
    expect_error(do.call(fun, bad), "Y|vector|column")
  }
  args <- paper_gcm_args(); args$B <- 2L
  args$X <- cbind(args$X, args$X)
  expect_error(do.call(perm_gcm, args), "X|column|univariate")
})

test_that("PCM validates target, random split fraction, and explicit indices", {
  for (target in list(0, 4, 1.5, NA_real_, c(1, 2))) {
    args <- paper_pcm_args(); args$B <- 2L; args$target <- target
    expect_error(do.call(perm_pcm, args), "target")
  }
  for (frac in list(0, 1, -0.1, 1.1, NA_real_, Inf, c(0.4, 0.6))) {
    args <- paper_pcm_args(); args$B <- 2L; args$indices <- NULL; args$frac <- frac
    expect_error(do.call(perm_pcm, args), "frac|fold")
  }
  for (indices in list(c(1L, 1L, 2L), c(0L, 2L, 3L), c(1L, 2L, 41L),
                        c(1, 2, 3.5), c(1, 2, NA_real_), integer(0), 1L,
                        seq_len(39L))) {
    args <- paper_pcm_args(); args$B <- 2L; args$indices <- indices
    expect_error(do.call(perm_pcm, args), "indices|fold|split")
  }
})

test_that("adaptive FDR validates budgets and the scalar FDR level", {
  callback <- function(active_idx, t) rep(0, length(active_idx))
  for (field in c("h", "B_max")) {
    for (invalid in list(0, -1, 1.5, NA_real_, Inf, numeric(0), c(1, 2))) {
      args <- list(T_orig = c(0, 1), permute_stat = callback,
                   alpha = 0.1, h = 2L, B_max = 3L)
      args[[field]] <- invalid
      expect_error(do.call(adaptive_fdr_bh, args),
                   paste0(field, "|positive integer"))
    }
  }
  for (invalid in list(0, 1, -0.1, 1.1, NA_real_, Inf,
                       numeric(0), c(0.05, 0.1), "0.1")) {
    expect_error(adaptive_fdr_bh(c(0, 1), callback, alpha = invalid,
                                 h = 2L, B_max = 3L), "alpha")
  }
})

test_that("adaptive FDR rejects nonnumeric or missing statistics", {
  callback <- function(active_idx, t) rep(0, length(active_idx))
  for (invalid in list(c("a", "b"), list(0, 1), c(0, NA_real_), c(0, NaN))) {
    expect_error(adaptive_fdr_bh(invalid, callback, h = 2L, B_max = 3L),
                 "T_orig|numeric|NA|missing")
  }
  for (invalid in list(c("a", "b"), list(0, 1), c(0, NA_real_), c(0, NaN))) {
    callback <- function(active_idx, t) invalid
    expect_error(adaptive_fdr_bh(c(0, 1), callback, h = 2L, B_max = 3L),
                 "permute_stat|numeric|NA|missing")
  }
})

test_that("adaptive FDR accepts infinite statistics and an empty family", {
  res <- adaptive_fdr_bh(c(Inf, -Inf),
                         function(active_idx, t) rep(0, length(active_idx)),
                         alpha = 0.1, h = 2L, B_max = 3L,
                         alternative = "two.sided")
  expect_equal(res$p_value, rep(2 / 5, 2L))
  expect_equal(res$active, c(1L, 2L))
  expect_length(res$rejected, 0L)

  res <- adaptive_fdr_bh(c(1, -1),
                         function(active_idx, t) rep(Inf, length(active_idx)),
                         alpha = 0.1, h = 2L, B_max = 3L,
                         alternative = "two.sided")
  expect_equal(res$p_value, c(1, 1))
  expect_equal(res$futility, c(1L, 2L))
  expect_length(res$active, 0L)

  empty <- adaptive_fdr_bh(numeric(0), function(active_idx, t) {
    stop("An empty hypothesis family must not invoke the callback.")
  })
  expect_true(all(lengths(empty[c("rejected", "futility", "active", "p_value",
                                  "K", "stop_round")]) == 0L))
  expect_identical(empty$t_total, 0L)
  expect_identical(empty$R_threshold, 0)
})
