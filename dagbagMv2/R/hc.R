.is_dag_matrix <- function(adj) {
  ## Kahn-style DAG check used for validating whitelist inputs before C++.
  p <- nrow(adj)
  indeg <- colSums(adj != 0)
  queue <- which(indeg == 0)
  seen <- 0L

  while (length(queue) > 0L) {
    node <- queue[1L]
    queue <- queue[-1L]
    seen <- seen + 1L
    children <- which(adj[node, ] != 0)
    for (child in children) {
      indeg[child] <- indeg[child] - 1L
      if (indeg[child] == 0L) {
        queue <- c(queue, child)
      }
    }
  }

  seen == p
}

.validate_controls <- function(tol, maxStep, restart, seed) {
  ## Keep public integer-like controls signed and checked before they cross the
  ## R/C++ boundary. The C++ entry points repeat these checks for direct use.
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("tol must be a finite nonnegative scalar")
  }
  if (!is.numeric(maxStep) || length(maxStep) != 1L || is.na(maxStep) || maxStep < 0) {
    stop("maxStep must be a nonnegative integer")
  }
  if (!is.numeric(restart) || length(restart) != 1L || is.na(restart) || restart < 1) {
    stop("restart must be a positive integer")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0) {
    stop("seed must be a nonnegative integer")
  }

  list(
    maxStep = as.integer(maxStep),
    restart = as.integer(restart),
    seed = as.integer(seed)
  )
}

.validate_hc_inputs <- function(Y, nodeType, whiteList, blackList,
                               tol, maxStep, restart, seed) {
  ## Validate and normalise all inputs shared by hc() and hc_boot().
  ## Returns raw (unscaled) Y; standardization is the caller's responsibility:
  ##   hc()      -- scales the full dataset once, before the single HC run.
  ##   hc_boot() -- scales each bootstrap sample independently inside
  ##                .fit_boot_one(), after row-resampling.
  ## The adjacency convention is matrix[from, to] == TRUE for from -> to.
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("Y must be a numeric matrix or data frame")
  }
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (!is.numeric(Y) || anyNA(Y) || any(!is.finite(Y))) {
    stop("Y must contain only finite numeric values")
  }

  p <- ncol(Y)
  if (p < 1L || nrow(Y) < 1L) {
    stop("Y must have positive numbers of rows and columns")
  }

  if (is.null(nodeType)) {
    nodeType <- rep("c", p)
  }
  if (!is.character(nodeType) || length(nodeType) != p || anyNA(nodeType) ||
      !all(nodeType %in% c("c", "b"))) {
    stop("nodeType must be a character vector of length ncol(Y) with entries 'c' or 'b'")
  }
  for (i in which(nodeType == "b")) {
    if (!all(Y[, i] %in% c(0, 1))) {
      stop("binary nodes must contain only 0/1 values")
    }
  }

  if (is.null(whiteList)) {
    whiteList <- matrix(FALSE, p, p)
  }
  if (!is.matrix(whiteList) || !is.logical(whiteList) ||
      !identical(dim(whiteList), c(p, p)) || anyNA(whiteList)) {
    stop("whiteList must be a p by p logical matrix with no NA values")
  }
  if (any(diag(whiteList))) {
    stop("whiteList diagonal must be FALSE")
  }
  if (any(whiteList & t(whiteList))) {
    stop("whiteList cannot contain both directions of an edge")
  }
  if (!.is_dag_matrix(whiteList)) {
    stop("whiteList must be acyclic")
  }

  if (is.null(blackList)) {
    blackList <- matrix(FALSE, p, p)
    diag(blackList) <- TRUE
  }
  if (!is.matrix(blackList) || !is.logical(blackList) ||
      !identical(dim(blackList), c(p, p)) || anyNA(blackList)) {
    stop("blackList must be a p by p logical matrix with no NA values")
  }
  diag(blackList) <- TRUE
  if (any(whiteList & blackList)) {
    stop("whiteList and blackList conflict")
  }

  controls <- .validate_controls(tol, maxStep, restart, seed)

  list(
    Y = Y,
    nodeType = nodeType,
    whiteList = whiteList,
    blackList = blackList,
    maxStep = controls$maxStep,
    restart = controls$restart,
    seed = controls$seed
  )
}

hc <- function(Y, nodeType = NULL, whiteList = NULL, blackList = NULL,
               standardize = TRUE, tol = 1e-6, maxStep = 2000L,
               restart = 1L, seed = 1L, verbose = FALSE, debug = FALSE,
               addDeleteOnly = FALSE) {
  ## standardize = TRUE: continuous columns of Y are scaled to mean 0 / sd 1
  ## on the full dataset before the HC run.
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }
  if (!is.logical(addDeleteOnly) || length(addDeleteOnly) != 1L || is.na(addDeleteOnly)) {
    stop("addDeleteOnly must be TRUE or FALSE")
  }
  args <- .validate_hc_inputs(Y, nodeType, whiteList, blackList,
                              tol, maxStep, restart, seed)
  if (standardize) {
    for (i in which(args$nodeType == "c")) {
      s <- stats::sd(args$Y[, i])
      if (!is.finite(s) || s <= 0)
        stop("continuous nodes must have positive finite standard deviation when standardize = TRUE")
      args$Y[, i] <- (args$Y[, i] - mean(args$Y[, i])) / s
    }
  }
  hc_(args$Y, args$nodeType, args$whiteList, args$blackList, tol,
      args$maxStep, args$restart, args$seed, verbose, debug, addDeleteOnly)
}

.fit_boot_one <- function(b, Y, nodeType, whiteList, blackList, tol, maxStep,
                          restart, hc_seed, boot_index, node_perm, standardize,
                          verbose, debug, addDeleteOnly) {
  ## Fit one bootstrap HC replicate.
  ## Y is the raw (unscaled) data matrix; boot_index and node_perm are
  ## pre-generated upfront so the result is independent of the backend.
  ##
  ## standardize = TRUE: each bootstrap sample is standardized to mean 0 / sd 1
  ## using its OWN rows, after row-resampling but before the HC run.
  ## This ensures every replicate sees exactly unit-scale data rather than
  ## approximately unit-scale data (which would result from scaling the full
  ## dataset once and then resampling rows).
  p <- ncol(Y)
  Y.B <- Y[boot_index, , drop = FALSE]  # bootstrap resample (raw scale)

  if (standardize) {
    ## Scale each continuous column using the bootstrap sample's own mean/sd.
    ## The is.finite(s) && s > 0 guard handles the rare case where a bootstrap
    ## column is constant (all duplicated rows); such columns are left unscaled
    ## and the C++ BIC score returns +Inf, which hc_() treats as a non-improving
    ## candidate (effectively excluding those edges).
    for (i in which(nodeType == "c")) {
      m <- mean(Y.B[, i])
      s <- stats::sd(Y.B[, i])
      if (is.finite(s) && s > 0) {
        Y.B[, i] <- (Y.B[, i] - m) / s
      }
    }
  }

  if (!is.null(node_perm)) {
    ## Permute node order before fitting and restore it afterward.
    ## HC is a greedy search and its result can depend on the order in which
    ## nodes are visited (especially at score ties).  Randomizing node order
    ## each bootstrap replicate reduces this ordering bias in the aggregated
    ## edge frequencies.
    node.index <- order(node_perm)  # inverse permutation for restoring original order
    Y.B <- Y.B[, node_perm, drop = FALSE]
    node.B <- nodeType[node_perm]
    whiteList.B <- whiteList[node_perm, node_perm, drop = FALSE]
    blackList.B <- blackList[node_perm, node_perm, drop = FALSE]
  } else {
    node.index <- seq_len(p)
    node.B <- nodeType
    whiteList.B <- whiteList
    blackList.B <- blackList
  }

  cur_res <- hc_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep,
                 restart, hc_seed, verbose, debug, addDeleteOnly)
  ## Restore the original node ordering before returning.
  cur_res$adjacency[node.index, node.index, drop = FALSE]
}

.format_boot_result <- function(result, p, n.boot, output_type) {
  ## Collate bootstrap adjacency matrices into the requested output form.
  ## "array"  : raw p x p x B adjacency array (largest memory, keeps all info).
  ## "freq"   : p x p edge-selection frequency matrix, freq[i,j] = fraction of
  ##            bootstrap DAGs containing edge i->j; diagonal forced to zero.
  ## "both"   : both the array and the freq matrix.
  if (output_type == "array") {
    return(array(as.numeric(unlist(result, use.names = FALSE)), dim = c(p, p, n.boot)))
  }

  ## Accumulate frequencies incrementally rather than building the full 3D array
  ## first, saving memory when only frequencies are needed.
  freq <- matrix(0, p, p)
  for (i in seq_len(n.boot)) {
    freq <- freq + result[[i]]
  }
  freq <- freq / n.boot
  diag(freq) <- 0

  if (output_type == "freq") {
    return(freq)
  }
  list(
    adjacency = array(as.numeric(unlist(result, use.names = FALSE)), dim = c(p, p, n.boot)),
    freq = freq
  )
}

.resolve_future_workers <- function(workers, n.boot, verbose) {
  ## Cap parallel workers to (available cores - 1) and to n.boot.
  ## Returns 1 if the machine reports no available cores.
  available <- max(1L, as.integer(future::availableCores()) - 1L)
  requested <- if (is.null(workers)) available else as.integer(workers)
  capped <- max(1L, min(requested, available, n.boot))
  if (verbose && capped < requested) {
    message("Capping workers from ", requested, " to ", capped,
            " based on available cores and n.boot.")
  }
  capped
}

hc_boot <- function(Y, n.boot = 100L, nodeType = NULL, whiteList = NULL,
                    blackList = NULL, standardize = TRUE, tol = 1e-6,
                    maxStep = 2000L, restart = 1L, seed = 1L,
                    nodeShuffle = TRUE,
                    backend = c("sequential", "future"),
                    workers = NULL,
                    verbose = FALSE, debug = FALSE,
                    addDeleteOnly = FALSE,
                    output_type = c("array", "freq", "both")) {
  ## Fit bootstrap HC replicates and aggregate edge frequencies.
  ##
  ## standardize = TRUE: each bootstrap sample is standardized independently
  ## to mean 0 / sd 1 (using that sample's own rows) inside .fit_boot_one(),
  ## immediately after row-resampling and before the HC run.  This is the
  ## correct bootstrap behavior: every replicate sees exactly unit-scale data.
  ## Contrast with hc(), where the full dataset is scaled once before the
  ## single HC run.
  ##
  ## backend = "sequential": replicates run one by one in the current session.
  ## backend = "future":     replicates are dispatched as parallel futures
  ##                         (future::multisession).
  ##
  ## All randomness (bootstrap row indices, HC tie-breaking seeds, node
  ## permutations) is generated upfront from a single set.seed(seed) call so
  ## that sequential and future backends produce bit-identical results.
  backend <- match.arg(backend)
  output_type <- match.arg(output_type)

  if (!is.numeric(n.boot) || length(n.boot) != 1L || is.na(n.boot) || n.boot < 1) {
    stop("n.boot must be a positive integer")
  }
  n.boot <- as.integer(n.boot)
  if (!is.logical(nodeShuffle) || length(nodeShuffle) != 1L || is.na(nodeShuffle)) {
    stop("nodeShuffle must be TRUE or FALSE")
  }
  if (!is.null(workers) && (!is.numeric(workers) || length(workers) != 1L ||
      is.na(workers) || workers < 1)) {
    stop("workers must be NULL or a positive integer")
  }
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }
  if (!is.logical(addDeleteOnly) || length(addDeleteOnly) != 1L || is.na(addDeleteOnly)) {
    stop("addDeleteOnly must be TRUE or FALSE")
  }
  args <- .validate_hc_inputs(Y, nodeType, whiteList, blackList,
                              tol, maxStep, restart, seed)
  n <- nrow(args$Y)
  p <- ncol(args$Y)

  ## Validate constant columns on the full data upfront so the user gets a clear
  ## error before any bootstrap work begins.  The actual transformation is deferred
  ## to .fit_boot_one so each replicate uses its own bootstrap-sample statistics.
  if (standardize) {
    for (i in which(args$nodeType == "c")) {
      s <- stats::sd(args$Y[, i])
      if (!is.finite(s) || s <= 0) {
        stop("continuous nodes must have positive finite standard deviation when standardize = TRUE")
      }
    }
  }

  ## Save and restore the caller's global RNG state so set.seed() below does
  ## not have side effects outside this function.
  old_rng <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (is.null(old_rng)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_rng, envir = .GlobalEnv)
    }
  }, add = TRUE)

  ## Generate all bootstrap indices, HC seeds, and node permutations from one
  ## seed.  Drawing them upfront (before any parallel dispatch) ensures that
  ## sequential and future backends visit the same bootstrap problems.
  set.seed(args$seed)
  boot_indices <- replicate(n.boot, sample.int(n, n, replace = TRUE),
                            simplify = FALSE)
  ## HC tie-breaking seed for bootstrap b = seed + b, giving each replicate a
  ## distinct seed without risk of collisions for typical n.boot values.
  if (args$seed > .Machine$integer.max - n.boot) {
    stop("seed is too large relative to n.boot; reduce seed so that seed + n.boot does not overflow")
  }
  hc_seeds <- args$seed + seq_len(n.boot)
  ## NULL signals the no-shuffle fast path in .fit_boot_one (avoids copying
  ## the permutation matrices for every bootstrap when nodeShuffle = FALSE).
  node_permutations <- if (nodeShuffle) {
    replicate(n.boot, sample.int(p), simplify = FALSE)
  } else {
    NULL
  }

  run_one <- function(b) {
    if (verbose) message("Processing bootstrap sample ", b)
    node_perm <- if (!is.null(node_permutations)) node_permutations[[b]] else NULL
    .fit_boot_one(b, args$Y, args$nodeType, args$whiteList, args$blackList,
                  tol, args$maxStep, args$restart, hc_seeds[b],
                  boot_indices[[b]], node_perm, standardize,
                  verbose, debug, addDeleteOnly)
  }

  if (backend == "sequential") {
    result <- vector("list", n.boot)
    for (b in seq_len(n.boot)) {
      result[[b]] <- run_one(b)
    }
  } else {
    ## future backend: dispatch all replicates as independent async futures and
    ## collect results.  seed = TRUE in future::future() suppresses the package's
    ## L'Ecuyer-CMRG RNG warning; actual reproducibility is controlled by the
    ## upfront set.seed() above, not by future's own RNG management.
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("The future package is required for backend = 'future'.")
    }
    nw <- .resolve_future_workers(workers, n.boot, verbose)
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = nw)
    futures <- lapply(seq_len(n.boot), function(b) {
      future::future(run_one(b), seed = TRUE)
    })
    result <- lapply(futures, future::value)
  }

  .format_boot_result(result, p, n.boot, output_type)
}

## ---------------------------------------------------------------------------
## Original parallel bootstrap implementation (preserved for comparison).
## Differences from the current hc_boot():
##   (1) Standardization: scales the full dataset ONCE before resampling,
##       so each bootstrap sample inherits approximately (not exactly) unit scale.
##       hc_boot() standardizes each bootstrap sample independently after
##       row-resampling, giving every replicate exactly unit-scale data.
##   (2) RNG: each worker calls set.seed(i * 1001L + seed) independently,
##       mixing the global seed with the replicate index.  hc_boot() draws all
##       randomness upfront from a single set.seed(seed) call, which guarantees
##       sequential and future backends are bit-identical.
##   (3) HC seed: i * 11L + seed (here) vs seed + i (hc_boot()).
## ---------------------------------------------------------------------------

.fit_boot_one_parallel <- function(i, Y, nodeType, whiteList, blackList, tol,
                                   maxStep, restart, seed, nodeShuffle,
                                   verbose, debug, addDeleteOnly) {
  ## Original .fit_boot_one logic: generates bootstrap indices and the node
  ## permutation from the same set.seed call inside the worker.
  p <- ncol(Y)
  n <- nrow(Y)
  set.seed(i * 1001L + seed)
  s.pick <- sample.int(n, n, replace = TRUE)

  if (nodeShuffle) {
    node.rand  <- sample.int(p, p, replace = FALSE)
    node.index <- order(node.rand)
  } else {
    node.rand  <- seq_len(p)
    node.index <- seq_len(p)
  }

  Y.B        <- Y[s.pick, node.rand, drop = FALSE]
  node.B     <- nodeType[node.rand]
  whiteList.B <- whiteList[node.rand, node.rand, drop = FALSE]
  blackList.B <- blackList[node.rand, node.rand, drop = FALSE]
  curRes <- hc_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep,
                restart, i * 11L + seed, verbose, debug, addDeleteOnly)
  curRes$adjacency[node.index, node.index, drop = FALSE]
}

hc_boot_parallel <- function(Y, n.boot = 100L, nodeType = NULL, whiteList = NULL,
                             blackList = NULL, standardize = TRUE, tol = 1e-6,
                             maxStep = 2000L, restart = 1L, seed = 1L,
                             nodeShuffle = TRUE, numThread = 2L,
                             verbose = FALSE, debug = FALSE,
                             addDeleteOnly = FALSE,
                             output_type = c("array", "freq", "both")) {
  ## Original parallel bootstrap HC using foreach + doFuture.
  ## Standardizes the full dataset once before resampling (full-data scaling).
  ## Use hc_boot(backend = "future") for per-sample standardization and
  ## upfront reproducible randomness.
  output_type <- match.arg(output_type)

  if (!requireNamespace("foreach",  quietly = TRUE) ||
      !requireNamespace("future",   quietly = TRUE) ||
      !requireNamespace("doFuture", quietly = TRUE)) {
    stop("hc_boot_parallel requires foreach, future, and doFuture")
  }
  if (!is.numeric(n.boot) || length(n.boot) != 1L || is.na(n.boot) || n.boot < 1) {
    stop("n.boot must be a positive integer")
  }
  if (!is.numeric(numThread) || length(numThread) != 1L ||
      is.na(numThread) || numThread < 1) {
    stop("numThread must be a positive integer")
  }
  if (!is.logical(nodeShuffle) || length(nodeShuffle) != 1L || is.na(nodeShuffle)) {
    stop("nodeShuffle must be TRUE or FALSE")
  }
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }
  if (!is.logical(addDeleteOnly) || length(addDeleteOnly) != 1L || is.na(addDeleteOnly)) {
    stop("addDeleteOnly must be TRUE or FALSE")
  }
  n.boot    <- as.integer(n.boot)
  numThread <- as.integer(numThread)

  ## Validate inputs; then apply full-dataset standardization so each worker
  ## receives pre-scaled data (original behavior).
  args <- .validate_hc_inputs(Y, nodeType, whiteList, blackList,
                              tol, maxStep, restart, seed)
  if (standardize) {
    for (i in which(args$nodeType == "c")) {
      s <- stats::sd(args$Y[, i])
      if (!is.finite(s) || s <= 0)
        stop("continuous nodes must have positive finite standard deviation when standardize = TRUE")
      args$Y[, i] <- (args$Y[, i] - mean(args$Y[, i])) / s
    }
  }
  p <- ncol(args$Y)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  cores   <- parallel::detectCores()
  if (is.na(cores) || cores < 1L) cores <- 1L
  workers <- min(numThread, max(1L, cores - 1L))
  future::plan(future::multisession, workers = workers)

  `%dofuture%` <- doFuture::`%dofuture%`
  i <- NULL
  result <- foreach::foreach(
    i = seq_len(n.boot),
    .errorhandling = "stop",
    .options.future = list(seed = TRUE, packages = "dagbagMv2")
  ) %dofuture% {
    .fit_boot_one_parallel(i, args$Y, args$nodeType, args$whiteList,
                           args$blackList, tol, args$maxStep,
                           args$restart, args$seed, nodeShuffle,
                           verbose, debug, addDeleteOnly)
  }

  .format_boot_result(result, p, n.boot, output_type)
}
