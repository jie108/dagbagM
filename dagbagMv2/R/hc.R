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

.prepare_hc_inputs <- function(Y, nodeType, whiteList, blackList, standardize,
                               tol, maxStep, restart, seed) {
  ## Shared input preparation for direct HC and bootstrap HC.
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

  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }
  if (standardize) {
    for (i in which(nodeType == "c")) {
      ## Match the original package behavior for continuous nodes, but fail
      ## early on constant columns instead of creating NaN values.
      s <- stats::sd(Y[, i])
      if (!is.finite(s) || s <= 0) {
        stop("continuous nodes must have positive finite standard deviation when standardize = TRUE")
      }
      Y[, i] <- (Y[, i] - mean(Y[, i])) / s
    }
  }

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
  if (!is.logical(addDeleteOnly) || length(addDeleteOnly) != 1L || is.na(addDeleteOnly)) {
    stop("addDeleteOnly must be TRUE or FALSE")
  }
  args <- .prepare_hc_inputs(Y, nodeType, whiteList, blackList, standardize,
                             tol, maxStep, restart, seed)
  hc_(args$Y, args$nodeType, args$whiteList, args$blackList, tol,
      args$maxStep, args$restart, args$seed, verbose, debug, addDeleteOnly)
}

.fit_boot_one <- function(i, Y, nodeType, whiteList, blackList, tol, maxStep,
                          restart, seed, nodeShuffle, verbose, debug,
                          addDeleteOnly) {
  ## Fit one bootstrap sample. If nodeShuffle is TRUE, fit in permuted node
  ## order and then map the learned adjacency back to the original order.
  p <- ncol(Y)
  n <- nrow(Y)
  set.seed(i * 1001L + seed)
  s.pick <- sample.int(n, n, replace = TRUE)

  if (nodeShuffle) {
    node.rand <- sample.int(p, p, replace = FALSE)
    node.index <- order(node.rand)
  } else {
    node.rand <- seq_len(p)
    node.index <- seq_len(p)
  }

  Y.B <- Y[s.pick, node.rand, drop = FALSE]
  node.B <- nodeType[node.rand]
  whiteList.B <- whiteList[node.rand, node.rand, drop = FALSE]
  blackList.B <- blackList[node.rand, node.rand, drop = FALSE]
  curRes <- hc_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep,
                restart, i * 11L + seed, verbose, debug, addDeleteOnly)
  curRes$adjacency[node.index, node.index, drop = FALSE]
}

.format_boot_result <- function(result, p, n.boot, return) {
  ## Large bootstrap runs can avoid storing the full p x p x B array by
  ## returning only edge selection frequencies.
  if (return == "array") {
    return(array(as.numeric(unlist(result, use.names = FALSE)), dim = c(p, p, n.boot)))
  }

  freq <- Reduce(`+`, lapply(result, function(x) x + 0)) / n.boot
  diag(freq) <- 0
  if (return == "freq") {
    return(freq)
  }
  list(
    adjacency = array(as.numeric(unlist(result, use.names = FALSE)), dim = c(p, p, n.boot)),
    freq = freq
  )
}

hc_boot <- function(Y, n.boot = 1L, nodeType = NULL, whiteList = NULL,
                    blackList = NULL, standardize = TRUE, tol = 1e-6,
                    maxStep = 2000L, restart = 1L, seed = 1L,
                    nodeShuffle = TRUE, verbose = FALSE, debug = FALSE,
                    addDeleteOnly = FALSE,
                    return = c("array", "freq", "both")) {
  return <- match.arg(return)
  if (!is.numeric(n.boot) || length(n.boot) != 1L || is.na(n.boot) || n.boot < 1) {
    stop("n.boot must be a positive integer")
  }
  n.boot <- as.integer(n.boot)
  if (!is.logical(nodeShuffle) || length(nodeShuffle) != 1L || is.na(nodeShuffle)) {
    stop("nodeShuffle must be TRUE or FALSE")
  }
  if (!is.logical(addDeleteOnly) || length(addDeleteOnly) != 1L || is.na(addDeleteOnly)) {
    stop("addDeleteOnly must be TRUE or FALSE")
  }

  args <- .prepare_hc_inputs(Y, nodeType, whiteList, blackList, standardize,
                             tol, maxStep, restart, seed)
  p <- ncol(args$Y)
  result <- vector("list", n.boot)
  for (i in seq_len(n.boot)) {
    if (verbose) {
      message("fit bootstrap sample ", i, "...")
    }
    result[[i]] <- .fit_boot_one(i, args$Y, args$nodeType, args$whiteList,
                                 args$blackList, tol, args$maxStep,
                                 args$restart, args$seed, nodeShuffle,
                                 verbose, debug, addDeleteOnly)
  }
  .format_boot_result(result, p, n.boot, return)
}

hc_boot_parallel <- function(Y, n.boot = 1L, nodeType = NULL, whiteList = NULL,
                             blackList = NULL, standardize = TRUE, tol = 1e-6,
                             maxStep = 2000L, restart = 1L, seed = 1L,
                             nodeShuffle = TRUE, numThread = 2L,
                             verbose = FALSE, debug = FALSE,
                             addDeleteOnly = FALSE,
                             return = c("array", "freq", "both")) {
  return <- match.arg(return)
  if (!requireNamespace("foreach", quietly = TRUE) ||
      !requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("doFuture", quietly = TRUE)) {
    stop("hc_boot_parallel requires foreach, future, and doFuture")
  }
  if (!is.numeric(n.boot) || length(n.boot) != 1L || is.na(n.boot) || n.boot < 1) {
    stop("n.boot must be a positive integer")
  }
  if (!is.numeric(numThread) || length(numThread) != 1L || is.na(numThread) || numThread < 1) {
    stop("numThread must be a positive integer")
  }
  if (!is.logical(nodeShuffle) || length(nodeShuffle) != 1L || is.na(nodeShuffle)) {
    stop("nodeShuffle must be TRUE or FALSE")
  }
  if (!is.logical(addDeleteOnly) || length(addDeleteOnly) != 1L || is.na(addDeleteOnly)) {
    stop("addDeleteOnly must be TRUE or FALSE")
  }
  n.boot <- as.integer(n.boot)
  numThread <- as.integer(numThread)

  args <- .prepare_hc_inputs(Y, nodeType, whiteList, blackList, standardize,
                             tol, maxStep, restart, seed)
  p <- ncol(args$Y)

  old_plan <- future::plan()
  ## Restore the user's previous future plan even if a worker errors.
  on.exit(future::plan(old_plan), add = TRUE)
  cores <- parallel::detectCores()
  if (is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  workers <- min(numThread, max(1L, cores - 1L))
  future::plan(future::multisession, workers = workers)

  `%dofuture%` <- doFuture::`%dofuture%`
  i <- NULL
  result <- foreach::foreach(
    i = seq_len(n.boot),
    .errorhandling = "stop",
    .options.future = list(seed = TRUE, packages = "dagbagMv2")
  ) %dofuture% {
    .fit_boot_one(i, args$Y, args$nodeType, args$whiteList,
                  args$blackList, tol, args$maxStep,
                  args$restart, args$seed, nodeShuffle,
                  verbose, debug, addDeleteOnly)
  }

  .format_boot_result(result, p, n.boot, return)
}
