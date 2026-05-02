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

moral_graph <- function(adj.matrix) {
  p <- nrow(adj.matrix)
  moral.adj.matrix <- adj.matrix

  n.parents <- colSums(adj.matrix)
  for (i in seq_len(p)) {
    if (n.parents[i] >= 2) {
      pars <- which(adj.matrix[, i] > 0)
      for (j in 1:(n.parents[i] - 1)) {
        for (k in (j + 1):n.parents[i]) {
          moral.adj.matrix[pars[j], pars[k]] <- 1
        }
      }
    }
  }

  (moral.adj.matrix + t(moral.adj.matrix)) > 0
}

vstructures <- function(adj.matrix) {
  p <- ncol(adj.matrix)
  res <- NULL
  for (i in seq_len(p)) {
    parent.node <- which(adj.matrix[, i] != 0)
    n.par <- length(parent.node)

    if (n.par > 1) {
      for (j in 1:(n.par - 1)) {
        for (k in (j + 1):n.par) {
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
