.validate_graph_constraint <- function(x, p, name, default = FALSE) {
  ## score_shd accepts logical and 0/1-style constraints; convert once and let
  ## the C++ aggregator enforce acyclicity/conflict rules.
  if (is.null(x)) {
    x <- matrix(default, p, p)
  }
  if (!is.matrix(x) || !identical(dim(x), c(p, p)) || anyNA(x)) {
    stop(name, " must be a p by p matrix with no NA values")
  }
  x <- x != 0
  dimnames(x) <- NULL
  x
}

.prepare_score_constraints <- function(p, whiteList, blackList) {
  whiteList <- .validate_graph_constraint(whiteList, p, "whiteList")
  blackList <- .validate_graph_constraint(blackList, p, "blackList")
  diag(blackList) <- TRUE
  list(whiteList = whiteList, blackList = blackList)
}

score_shd <- function(boot.adj, alpha = 1, freqCutoff = 0.5, whiteList = NULL,
                      blackList = NULL, maxStep = NULL, verbose = FALSE) {
  ## C++ computes edge frequencies from the bootstrap array and then performs
  ## the deterministic greedy generalized-SHD aggregation.
  ## maxStep is retained for backward compatibility with v1 call sites but has
  ## no effect; the C++ aggregator processes all eligible candidates in one pass.
  if (!is.null(maxStep)) {
    warning("maxStep is deprecated and has no effect", call. = FALSE)
  }
  if (!is.numeric(freqCutoff) || length(freqCutoff) != 1L ||
      !is.finite(freqCutoff) || freqCutoff < 0 || freqCutoff > 1) {
    stop("freqCutoff must be a finite scalar in [0, 1]")
  }
  if (is.null(dim(boot.adj)) || length(dim(boot.adj)) != 3L ||
      dim(boot.adj)[1] != dim(boot.adj)[2]) {
    stop("boot.adj must be a p by p by B array")
  }
  p <- dim(boot.adj)[1]
  constraints <- .prepare_score_constraints(p, whiteList, blackList)
  score_shd_cpp(boot.adj, alpha, freqCutoff, constraints$whiteList,
                constraints$blackList, verbose)
}

score_shd_freq <- function(freq, alpha = 1, freqCutoff = 0.5, whiteList = NULL,
                           blackList = NULL, maxStep = NULL, verbose = FALSE) {
  ## Frequency-only aggregation pairs with hc_boot(..., return = "freq") so
  ## large bootstrap runs do not need to retain all individual adjacency arrays.
  ## maxStep is retained for backward compatibility with v1 call sites but has
  ## no effect; the C++ aggregator processes all eligible candidates in one pass.
  if (!is.null(maxStep)) {
    warning("maxStep is deprecated and has no effect", call. = FALSE)
  }
  if (!is.numeric(freqCutoff) || length(freqCutoff) != 1L ||
      !is.finite(freqCutoff) || freqCutoff < 0 || freqCutoff > 1) {
    stop("freqCutoff must be a finite scalar in [0, 1]")
  }
  if (!is.matrix(freq) || nrow(freq) != ncol(freq)) {
    stop("freq must be a square matrix")
  }
  freq <- as.matrix(freq)
  storage.mode(freq) <- "double"
  p <- nrow(freq)
  constraints <- .prepare_score_constraints(p, whiteList, blackList)
  score_shd_freq_cpp(freq, alpha, freqCutoff, constraints$whiteList,
                     constraints$blackList, verbose)
}
