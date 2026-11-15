# pcmt

A minimal R package supporting the experiments in
[Conditional independence testing through marginal permutation](https://arxiv.org/abs/2611.03255).
The public interface contains three functions:

- `perm_gcm()`: permutation calibration of the generalised covariance measure.
- `perm_pcm()`: permutation calibration of the projected covariance measure.
- `adaptive_fdr_bh()`: adaptive resampling with a Benjamini--Hochberg threshold.

Nuisance regressions use random forests through `comets`. Resampling is
sequential; `n_cores = 1` is retained for compatibility with the reproduction
scripts. Finite-sample validity requires independence of the exposure from the
outcome and covariates. Under conditional independence, validity is asymptotic
and requires the nuisance-estimation and permutation-stability conditions in
the paper.

## Installation

After the GitHub repository is public, install with:

```r
install.packages("remotes")
remotes::install_github("ZiangNiu6/pcmt")
```

## GCM and PCM

These examples use small forests and resampling budgets to keep execution
short. Use the settings in the paper's reproduction code for its experiments.

```r
library(pcmt)
set.seed(1)
n <- 80L
Z <- matrix(rnorm(n * 2L), ncol = 2L)
X <- Z[, 1] + rnorm(n)
Y <- sin(Z[, 1]) + Z[, 2] + rnorm(n)
rf_args <- list(num.trees = 50L, num.threads = 1L)

gcm <- perm_gcm(
  Y, X, Z, B = 49L,
  args_YonZ = rf_args, args_XonZ = rf_args
)
gcm$p.value

pcm <- perm_pcm(
  Y, cbind(X, Z), target = 1L, B = 49L, alternative = "greater",
  args_YonXZ = rf_args, args_YonZ = rf_args,
  args_YhatonZ = rf_args, args_RonZ = rf_args
)
pcm$p.value
```

`perm_gcm()` and `perm_pcm()` use a two-sided alternative by default. The
paper's PCM simulations set `alternative = "greater"`, as in the example. Supplying
`h` enables adaptive stopping, with `B` as the maximum number of draws and
`alpha` as the stopping threshold. For PCM, `frac` is the fraction assigned to
the direction-estimation fold D2; `indices` supplies that fold directly.

## Adaptive multiple testing

`adaptive_fdr_bh()` accepts observed statistics and a callback returning one
null-resampled statistic per active hypothesis. The callback must implement
the intended null resampling law. This example permutes a common exposure
against independent response columns:

```r
set.seed(2)
X_test <- rnorm(n)
Y_test <- matrix(rnorm(n * 4L), nrow = n)
T_orig <- colSums(Y_test * X_test)
permute_stat <- function(active_idx, t) {
  colSums(Y_test[, active_idx, drop = FALSE] * sample(X_test))
}

multiple <- adaptive_fdr_bh(
  T_orig, permute_stat, alpha = 0.10, h = 15L,
  B_max = 199L, alternative = "two.sided"
)
multiple$rejected
```

See the function help pages for returned components and resampling details.
Licensed under GPL (>= 3).
