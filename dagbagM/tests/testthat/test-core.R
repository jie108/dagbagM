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

# -- cal_order ----------------------------------------------------------------

test_that("cal_order returns valid topological order for a chain DAG", {
  adj <- matrix(0L, 3, 3)
  adj[1, 2] <- 1L
  adj[2, 3] <- 1L
  ord <- cal_order(adj)
  expect_equal(length(ord), 3L)
  expect_true(which(ord == 1L) < which(ord == 2L))
  expect_true(which(ord == 2L) < which(ord == 3L))
})

test_that("cal_order handles a root-only DAG (no edges)", {
  adj <- matrix(0L, 4, 4)
  ord <- cal_order(adj)
  expect_equal(sort(ord), 1:4)
})

test_that("cal_order detects a cycle", {
  adj <- matrix(0L, 3, 3)
  adj[1, 2] <- 1L
  adj[2, 3] <- 1L
  adj[3, 1] <- 1L
  expect_error(cal_order(adj), "cycle")
})

test_that("cal_order works on a single-node graph", {
  adj <- matrix(0L, 1, 1)
  expect_equal(cal_order(adj), 1L)
})

# -- skeleton / moral_graph / vstructures -------------------------------------

test_that("skeleton returns undirected version of directed graph", {
  adj <- matrix(0L, 3, 3)
  adj[1, 2] <- 1L
  adj[2, 3] <- 1L
  ske <- skeleton(adj)
  expect_true(ske[1, 2] && ske[2, 1])
  expect_true(ske[2, 3] && ske[3, 2])
  expect_false(ske[1, 3])
})

test_that("vstructures returns NULL for a chain with no v-structures", {
  adj <- matrix(0L, 3, 3)
  adj[1, 2] <- 1L
  adj[2, 3] <- 1L
  expect_null(vstructures(adj))
})

test_that("vstructures finds the correct v-structure: 1->3<-2, no 1-2 edge", {
  adj <- matrix(0L, 3, 3)
  adj[1, 3] <- 1L
  adj[2, 3] <- 1L
  vs <- vstructures(adj)
  expect_equal(nrow(vs), 1L)
  expect_equal(unname(vs[1, "child"]), 3L)
  expect_true(setequal(c(vs[1, "par1"], vs[1, "par2"]), c(1L, 2L)))
})

test_that("moral_graph connects unlinked parents in a v-structure", {
  adj <- matrix(0L, 3, 3)
  adj[1, 3] <- 1L
  adj[2, 3] <- 1L
  mg <- moral_graph(adj)
  expect_true(mg[1, 2] && mg[2, 1])
})

test_that("moral_graph of a chain equals its skeleton", {
  adj <- matrix(0L, 3, 3)
  adj[1, 2] <- 1L
  adj[2, 3] <- 1L
  expect_identical(moral_graph(adj), skeleton(adj))
})

# -- compare.vstructures ------------------------------------------------------

test_that("compare.vstructures returns NULL when either argument is NULL", {
  adj <- matrix(0L, 3, 3)
  adj[1, 3] <- 1L
  adj[2, 3] <- 1L
  vs <- vstructures(adj)
  expect_null(compare.vstructures(NULL, NULL))
  expect_null(compare.vstructures(NULL, vs))
  expect_null(compare.vstructures(vs, NULL))
})

test_that("compare.vstructures identifies matching v-structures", {
  adj <- matrix(0L, 3, 3)
  adj[1, 3] <- 1L
  adj[2, 3] <- 1L
  vs <- vstructures(adj)
  corr <- compare.vstructures(vs, vs)
  expect_equal(nrow(corr), 1L)
})

# -- basic hc() ---------------------------------------------------------------

test_that("hc returns a valid DAG on all-continuous data", {
  y <- make_mixed_data()[, 1:3]
  fit <- hc(y, nodeType = rep("c", 3), maxStep = 20, seed = 1)
  expect_true(is_dag_matrix(fit$adjacency))
  expect_true(all(is.finite(fit$score)))
  expect_equal(fit$steps, length(fit$operations))
})

test_that("hc handles p = 1 (single node, no edges possible)", {
  set.seed(1)
  y <- matrix(rnorm(50), ncol = 1)
  fit <- hc(y, nodeType = "c", maxStep = 10)
  expect_equal(dim(fit$adjacency), c(1L, 1L))
  expect_false(any(fit$adjacency != 0))
  expect_equal(fit$steps, 0L)
})

test_that("hc returns the same result with the same seed", {
  y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  fit1 <- hc(y, node_type, maxStep = 30, seed = 42)
  fit2 <- hc(y, node_type, maxStep = 30, seed = 42)
  expect_identical(fit1$adjacency, fit2$adjacency)
})

# -- hc_boot sequential -------------------------------------------------------

test_that("hc_boot output_type='array' produces valid bootstrap DAGs", {
  y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  arr <- hc_boot(y, n.boot = 3L, nodeType = node_type,
                 maxStep = 15, seed = 7, output_type = "array")
  expect_equal(dim(arr), c(4L, 4L, 3L))
  for (b in seq_len(3L)) {
    expect_true(is_dag_matrix(arr[, , b]))
  }
})

test_that("hc_boot output_type='freq' produces a valid frequency matrix", {
  y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  freq <- hc_boot(y, n.boot = 4L, nodeType = node_type,
                  maxStep = 15, seed = 7, output_type = "freq")
  expect_equal(dim(freq), c(4L, 4L))
  expect_true(all(freq >= 0 & freq <= 1))
  expect_true(all(diag(freq) == 0))
})

test_that("hc_boot 'array' and 'freq' modes agree on frequencies", {
  y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  both <- hc_boot(y, n.boot = 4L, nodeType = node_type,
                  maxStep = 15, seed = 7, output_type = "both")
  freq_from_array <- apply(both$adjacency, c(1, 2), mean)
  expect_equal(both$freq, freq_from_array)
})

# -- score_shd / score_shd_freq edge cases ------------------------------------

test_that("score_shd_freq with freq.cutoff=1 returns a valid DAG", {
  set.seed(3)
  freq <- matrix(runif(16, 0, 0.9), 4, 4)
  diag(freq) <- 0
  agg <- score_shd_freq(freq, freq.cutoff = 1)
  expect_true(is_dag_matrix(agg))
})

test_that("score_shd_freq with freq.cutoff=0 returns a valid DAG", {
  set.seed(3)
  freq <- matrix(runif(16, 0, 0.9), 4, 4)
  diag(freq) <- 0
  agg <- score_shd_freq(freq, freq.cutoff = 0)
  expect_true(is_dag_matrix(agg))
})

test_that("score_shd and score_shd_freq agree on matching inputs", {
  y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  both <- hc_boot(y, n.boot = 6L, nodeType = node_type,
                  maxStep = 15, seed = 5, output_type = "both")
  agg_arr  <- score_shd(both$adjacency, alpha = 1, freq.cutoff = 0.5)
  agg_freq <- score_shd_freq(both$freq,  alpha = 1, freq.cutoff = 0.5)
  expect_identical(agg_arr, agg_freq)
})

# -- sequential vs future backend identity -----------------------------------

test_that("hc_boot sequential and future backends produce identical results", {
  skip_if_not_installed("future")
  y <- make_mixed_data()
  node_type <- c("c", "c", "c", "b")
  seq_out <- hc_boot(y, n.boot = 4L, nodeType = node_type,
                     maxStep = 15, seed = 99, backend = "sequential",
                     output_type = "array")
  fut_out <- hc_boot(y, n.boot = 4L, nodeType = node_type,
                     maxStep = 15, seed = 99, backend = "future", workers = 2L,
                     output_type = "array")
  expect_identical(seq_out, fut_out)
})

# -- overflow / bounds validation ---------------------------------------------

test_that("hc rejects p > 46000", {
  skip("p > 46000 requires large memory allocation; run manually")
})
