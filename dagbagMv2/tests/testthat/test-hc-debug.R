is_dag_matrix <- function(adj) {
  p <- nrow(adj)
  indeg <- colSums(adj != 0)
  queue <- which(indeg == 0)
  seen <- 0L

  while (length(queue) > 0L) {
    node <- queue[1L]
    queue <- queue[-1L]
    seen <- seen + 1L
    for (child in which(adj[node, ] != 0)) {
      indeg[child] <- indeg[child] - 1L
      if (indeg[child] == 0L) {
        queue <- c(queue, child)
      }
    }
  }

  seen == p
}

make_mixed_data <- function(n = 100L) {
  set.seed(11)
  x1 <- rnorm(n)
  x2 <- 0.8 * x1 + rnorm(n, 0, 0.6)
  x3 <- -0.4 * x1 + 0.7 * x2 + rnorm(n, 0, 0.7)
  x4 <- rbinom(n, 1, plogis(0.9 * x1 - 0.6 * x2))
  cbind(x1, x2, x3, x4)
}

test_that("HC debug mode returns valid DAGs", {
  Y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")

  fit <- hc(Y, node_type, maxStep = 40, restart = 3, seed = 1, debug = TRUE)
  expect_true(is_dag_matrix(fit$adjacency))
  expect_true(all(is.finite(fit$score)))
  expect_equal(fit$steps, length(fit$operations))
})

test_that("addDeleteOnly skips reversal operations", {
  Y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")

  fit <- hc(Y, node_type, maxStep = 40, restart = 3, seed = 1,
            debug = TRUE, addDeleteOnly = TRUE)
  ops <- if (length(fit$operations)) {
    do.call(rbind, fit$operations)
  } else {
    matrix(integer(), 0, 3)
  }

  expect_true(is_dag_matrix(fit$adjacency))
  expect_false(any(ops[, 3] == 3))
})

test_that("bootstrap and C++ aggregation work with debug HC", {
  Y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")

  arr <- hc_boot(Y, n.boot = 4, nodeType = node_type, maxStep = 20,
                 restart = 2, seed = 2, debug = TRUE, return = "array")
  freq <- hc_boot(Y, n.boot = 4, nodeType = node_type, maxStep = 20,
                  restart = 2, seed = 2, debug = TRUE,
                  addDeleteOnly = TRUE, return = "freq")
  agg_array <- score_shd(arr)
  agg_freq <- score_shd_freq(freq)

  expect_equal(dim(arr), c(4L, 4L, 4L))
  expect_equal(dim(freq), c(4L, 4L))
  expect_true(is_dag_matrix(agg_array))
  expect_true(is_dag_matrix(agg_freq))
})

test_that("validation catches invalid graphs and frequencies", {
  Y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  whitelist <- matrix(FALSE, 4, 4)
  whitelist[1, 2] <- TRUE
  whitelist[2, 1] <- TRUE

  expect_error(hc(Y, node_type, whiteList = whitelist, debug = TRUE),
               "both directions|acyclic")
  expect_error(score_shd_freq(matrix(c(0, 2, 0, 0), 2, 2)),
               "\\[0, 1\\]")
})
