## Benchmark: adaptive vs fixed-B permutation
##
## Compares wall-clock time and number of permutations used under:
##   (a) the null hypothesis (adaptive stops for futility)
##   (b) the alternative hypothesis (adaptive stops for rejection)

library(pcmt)

cat("=== Adaptive permutation benchmark ===\n\n")

run_benchmark <- function(label, n, B, h, n_reps, gen_data, test_fn, ...) {
  cat(sprintf("--- %s (n=%d, B=%d, h=%d, reps=%d) ---\n", label, n, B, h, n_reps))

  times_fixed <- numeric(n_reps)
  times_adaptive <- numeric(n_reps)
  perms_adaptive <- numeric(n_reps)
  stop_reasons <- character(n_reps)

  for (r in seq_len(n_reps)) {
    set.seed(r)
    dat <- gen_data(n)

    extra_args <- list(...)

    ## Fixed-B
    set.seed(r + 10000L)
    t0 <- proc.time()["elapsed"]
    res_fixed <- do.call(test_fn, c(dat, list(B = B), extra_args))
    times_fixed[r] <- proc.time()["elapsed"] - t0

    ## Adaptive
    set.seed(r + 10000L)
    t0 <- proc.time()["elapsed"]
    res_adapt <- do.call(test_fn, c(dat, list(B = B, h = h, alpha = 0.05), extra_args))
    times_adaptive[r] <- proc.time()["elapsed"] - t0

    perms_adaptive[r] <- length(res_adapt$perm_stats)
    stop_reasons[r] <- res_adapt$adaptive$stop_reason
  }

  speedup <- times_fixed / times_adaptive
  perm_ratio <- perms_adaptive / B

  cat(sprintf("  Fixed-B time:     mean=%.3fs, median=%.3fs\n",
              mean(times_fixed), median(times_fixed)))
  cat(sprintf("  Adaptive time:    mean=%.3fs, median=%.3fs\n",
              mean(times_adaptive), median(times_adaptive)))
  cat(sprintf("  Speedup:          mean=%.1fx, median=%.1fx\n",
              mean(speedup), median(speedup)))
  cat(sprintf("  Perms used:       mean=%.0f/%d (%.0f%%), median=%.0f/%d (%.0f%%)\n",
              mean(perms_adaptive), B, 100 * mean(perm_ratio),
              median(perms_adaptive), B, 100 * median(perm_ratio)))
  cat(sprintf("  Stop reasons:     rejection=%d, futility=%d, cap=%d\n",
              sum(stop_reasons == "rejection"),
              sum(stop_reasons == "futility"),
              sum(stop_reasons == "cap")))
  cat("\n")

  invisible(list(
    times_fixed = times_fixed, times_adaptive = times_adaptive,
    speedup = speedup, perms_adaptive = perms_adaptive,
    stop_reasons = stop_reasons
  ))
}

## ============================================================
## 1. perm_gcm with random forest (expensive refit)
## ============================================================

cat("\n========== perm_gcm (rf) ==========\n\n")

## Null: X independent of (Y, Z)
run_benchmark(
  label = "GCM-rf NULL",
  n = 100, B = 499L, h = 15L, n_reps = 20L,
  gen_data = function(n) {
    Z <- matrix(rnorm(n), ncol = 1)
    list(Y = Z[,1] + rnorm(n), X = rnorm(n), Z = Z)
  },
  test_fn = perm_gcm,
  reg_YonZ = "rf", reg_XonZ = "rf", alternative = "two.sided"
)

## Alternative: Y depends on X given Z
run_benchmark(
  label = "GCM-rf ALTERNATIVE",
  n = 100, B = 499L, h = 15L, n_reps = 20L,
  gen_data = function(n) {
    Z <- matrix(rnorm(n), ncol = 1)
    X <- Z[,1] + rnorm(n)
    Y <- X + Z[,1] + rnorm(n, sd = 0.5)
    list(Y = Y, X = X, Z = Z)
  },
  test_fn = perm_gcm,
  reg_YonZ = "rf", reg_XonZ = "rf", alternative = "two.sided"
)

## ============================================================
## 2. perm_gcm with lrm (cheap refit -- less room for speedup)
## ============================================================

cat("\n========== perm_gcm (lrm, sequential) ==========\n\n")

## For lrm, adaptive forces sequential mode; fixed-B uses vectorized fast path.
## This isolates the permutation-count savings from the fast-path advantage.

## Null
run_benchmark(
  label = "GCM-lrm NULL",
  n = 200, B = 999L, h = 15L, n_reps = 50L,
  gen_data = function(n) {
    Z <- matrix(rnorm(n * 3), ncol = 3)
    list(Y = Z[,1] + rnorm(n), X = rnorm(n), Z = Z)
  },
  test_fn = perm_gcm,
  reg_YonZ = "lrm", reg_XonZ = "lrm", alternative = "two.sided"
)

## Alternative
run_benchmark(
  label = "GCM-lrm ALTERNATIVE",
  n = 200, B = 999L, h = 15L, n_reps = 50L,
  gen_data = function(n) {
    Z <- matrix(rnorm(n * 3), ncol = 3)
    X <- Z[,1] + rnorm(n)
    Y <- X + Z[,2] + rnorm(n, sd = 0.5)
    list(Y = Y, X = X, Z = Z)
  },
  test_fn = perm_gcm,
  reg_YonZ = "lrm", reg_XonZ = "lrm", alternative = "two.sided"
)

## ============================================================
## 3. perm_pcm with lrm
## ============================================================

cat("\n========== perm_pcm (lrm) ==========\n\n")

## Null: Y _||_ X_1 | X_{-1}
run_benchmark(
  label = "PCM-lrm NULL",
  n = 200, B = 499L, h = 15L, n_reps = 50L,
  gen_data = function(n) {
    X <- matrix(rnorm(n * 3), ncol = 3)
    Y <- X[,2] + X[,3] + rnorm(n)
    list(Y = Y, X = X, target = 1L)
  },
  test_fn = perm_pcm,
  reg_YonXZ = "lrm", reg_YonZ = "lrm",
  reg_YhatonZ = "lrm", reg_RonZ = "lrm",
  alternative = "two.sided"
)

## Alternative: Y depends on X_1
run_benchmark(
  label = "PCM-lrm ALTERNATIVE",
  n = 200, B = 499L, h = 15L, n_reps = 50L,
  gen_data = function(n) {
    X <- matrix(rnorm(n * 3), ncol = 3)
    Y <- X[,1] + X[,2] + X[,3] + rnorm(n, sd = 0.5)
    list(Y = Y, X = X, target = 1L)
  },
  test_fn = perm_pcm,
  reg_YonXZ = "lrm", reg_YonZ = "lrm",
  reg_YhatonZ = "lrm", reg_RonZ = "lrm",
  alternative = "two.sided"
)

## ============================================================
## 4. perm_pcm with rf (expensive)
## ============================================================

cat("\n========== perm_pcm (rf) ==========\n\n")

## Null
run_benchmark(
  label = "PCM-rf NULL",
  n = 200, B = 199L, h = 15L, n_reps = 10L,
  gen_data = function(n) {
    X <- matrix(rnorm(n * 3), ncol = 3)
    Y <- X[,2] + X[,3] + rnorm(n)
    list(Y = Y, X = X, target = 1L)
  },
  test_fn = perm_pcm,
  reg_YonXZ = "rf", reg_YonZ = "rf",
  reg_YhatonZ = "rf", reg_RonZ = "rf",
  alternative = "two.sided"
)

## Alternative
run_benchmark(
  label = "PCM-rf ALTERNATIVE",
  n = 200, B = 199L, h = 15L, n_reps = 10L,
  gen_data = function(n) {
    X <- matrix(rnorm(n * 3), ncol = 3)
    Y <- X[,1] + X[,2] + X[,3] + rnorm(n, sd = 0.5)
    list(Y = Y, X = X, target = 1L)
  },
  test_fn = perm_pcm,
  reg_YonXZ = "rf", reg_YonZ = "rf",
  reg_YhatonZ = "rf", reg_RonZ = "rf",
  alternative = "two.sided"
)

## ============================================================
## 5. perm_wgcm with lrm
## ============================================================

cat("\n========== perm_wgcm (lrm) ==========\n\n")

## Null
run_benchmark(
  label = "WGCM-lrm NULL",
  n = 150, B = 499L, h = 15L, n_reps = 50L,
  gen_data = function(n) {
    Z <- matrix(rnorm(n), ncol = 1)
    list(Y = Z[,1] + rnorm(n), X = rnorm(n), Z = Z)
  },
  test_fn = perm_wgcm,
  reg_YonZ = "lrm", reg_XonZ = "lrm", reg_wfun = "lrm",
  alternative = "two.sided"
)

## Alternative
run_benchmark(
  label = "WGCM-lrm ALTERNATIVE",
  n = 150, B = 499L, h = 15L, n_reps = 50L,
  gen_data = function(n) {
    Z <- matrix(rnorm(n), ncol = 1)
    X <- Z[,1] + rnorm(n)
    Y <- X + Z[,1] + rnorm(n, sd = 0.5)
    list(Y = Y, X = X, Z = Z)
  },
  test_fn = perm_wgcm,
  reg_YonZ = "lrm", reg_XonZ = "lrm", reg_wfun = "lrm",
  alternative = "two.sided"
)

cat("\n=== Benchmark complete ===\n")
