# Internal input guards. These do not fit models or consume random numbers.
validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be one positive integer.", call. = FALSE)
  }
  as.integer(x)
}

validate_permutation_controls <- function(B, n_cores, h, alpha) {
  B <- validate_positive_integer(B, "B")
  n_cores <- validate_positive_integer(n_cores, "n_cores")
  if (n_cores != 1L) {
    stop("n_cores must be 1 in this paper-specific package.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be one number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(h)) h <- validate_positive_integer(h, "h")
  list(B = B, h = h)
}

validate_rf_methods <- function(methods) {
  for (name in names(methods)) {
    if (!identical(methods[[name]], "rf")) {
      stop(name, " must be \"rf\" in this paper-specific package.", call. = FALSE)
    }
  }
  invisible(NULL)
}

validate_numeric_data <- function(x, name) {
  numeric_columns <- if (is.data.frame(x)) {
    all(vapply(x, is.numeric, logical(1)))
  } else {
    is.numeric(x)
  }
  if (!numeric_columns || NROW(x) == 0L || NCOL(x) == 0L ||
      any(!is.finite(as.matrix(x)))) {
    stop(name, " must contain finite numeric data.", call. = FALSE)
  }
  invisible(NULL)
}

validate_fraction <- function(frac) {
  if (!is.numeric(frac) || length(frac) != 1L || !is.finite(frac) ||
      frac <= 0 || frac >= 1) {
    stop("frac must be one number strictly between 0 and 1.", call. = FALSE)
  }
  invisible(NULL)
}

validate_split_indices <- function(indices, n) {
  if (!is.numeric(indices) || any(!is.finite(indices)) ||
      any(indices != floor(indices)) || anyDuplicated(indices) ||
      any(indices < 1L | indices > n) || length(indices) < 2L ||
      n - length(indices) < 2L) {
    stop("indices must be distinct valid integer D2 indices leaving at least two observations in each fold.",
         call. = FALSE)
  }
  invisible(NULL)
}
