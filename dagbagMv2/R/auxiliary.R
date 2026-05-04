cal_order <- function(adj_matrix) {
  if (nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("adjacency matrix is not square matrix!")
  }
  p <- nrow(adj_matrix)
  in_degree <- colSums(adj_matrix != 0)
  queue <- which(in_degree == 0)
  node_order <- integer(0)

  while (length(queue) > 0) {
    node <- queue[1L]
    queue <- queue[-1L]
    node_order <- c(node_order, node)
    for (child in which(adj_matrix[node, ] != 0)) {
      in_degree[child] <- in_degree[child] - 1L
      if (in_degree[child] == 0L) {
        queue <- c(queue, child)
      }
    }
  }

  if (length(node_order) != p) {
    stop("adj_matrix contains a cycle; topological order does not exist")
  }
  node_order
}

## Moralization: for each node i with >= 2 parents, add an undirected "spouse"
## edge between every pair of those parents (marry the parents).  The resulting
## moral graph is undirected: it contains the skeleton of the DAG plus these
## added spouse edges.  The final symmetrization (A + t(A)) > 0 ensures the
## matrix is symmetric even if the DAG itself has directed edges between parents.
moral_graph <- function(adj.matrix) {
  p <- nrow(adj.matrix)
  moral.adj.matrix <- adj.matrix

  n.parents <- colSums(adj.matrix)
  for (i in seq_len(p)) {
    if (n.parents[i] >= 2) {
      pars <- which(adj.matrix[, i] > 0)
      ## Add spouse edge for every unordered pair of parents (j < k ensures
      ## each pair is visited once; symmetrization below makes both directions).
      for (j in seq_len(n.parents[i] - 1)) {
        for (k in (j + 1):n.parents[i]) {
          moral.adj.matrix[pars[j], pars[k]] <- 1
        }
      }
    }
  }

  ## Symmetrize: (A + t(A)) > 0 captures both the original directed edges
  ## (their undirected versions) and the newly added spouse edges.
  (moral.adj.matrix + t(moral.adj.matrix)) > 0
}

## Enumerate all v-structures (unshielded colliders) in the DAG.
## A v-structure is a triple par1 -> child <- par2 where par1 and par2 are
## NOT adjacent (no edge in either direction between them).
## Each triple is reported exactly once with par1 < par2 (canonical ordering)
## because the inner loops enumerate unordered pairs (j < k).
## Returns an Nx3 matrix with columns (par1, child, par2), or NULL if none exist.
vstructures <- function(adj.matrix) {
  p <- ncol(adj.matrix)
  res <- NULL
  for (i in seq_len(p)) {
    parent.node <- which(adj.matrix[, i] != 0)
    n.par <- length(parent.node)

    if (n.par > 1) {
      ## j < k ensures par1 < par2 in result; unordered pair avoids duplicates.
      for (j in seq_len(n.par - 1)) {
        for (k in (j + 1):n.par) {
          ## No edge in either direction between the two parents = unshielded.
          if (adj.matrix[parent.node[j], parent.node[k]] == 0 &&
              adj.matrix[parent.node[k], parent.node[j]] == 0) {
            res <- rbind(res, c(parent.node[j], i, parent.node[k]))
          }
        }
      }
    }
  }

  if (!is.null(res)) {
    colnames(res) <- c("par1", "child", "par2")
  }
  res
}

skeleton <- function(adj.matrix) {
  (adj.matrix + t(adj.matrix)) > 0
}

## Find which v-structures in target also appear in true (precision-style lookup).
## Semantics are asymmetric: the outer loop iterates over target; each row in
## target is looked up in true.  The return value is the subset of TRUE rows that
## were matched: at most nrow(target) rows and at most nrow(true) distinct
## rows.  Pass (estimated, truth) to compute true-positive v-structures; swap
## arguments to get false positives from the other direction.
compare.vstructures <- function(target.vstructures, true.vstructures) {
  corr.v <- NULL
  if (!is.null(target.vstructures) && !is.null(true.vstructures)) {
    target.vstructures <- matrix(target.vstructures, ncol = 3)
    true.vstructures <- matrix(true.vstructures, ncol = 3)
    target.l <- nrow(target.vstructures)
    for (i in seq_len(target.l)) {
      res <- apply(true.vstructures, 1, function(x) all(x == target.vstructures[i, ]))
      if (any(res)) {
        corr.v <- rbind(corr.v, true.vstructures[which(res), , drop = FALSE])
      }
    }
  }
  corr.v
}
