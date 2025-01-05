UNITEST_01 <- function() {
  data(example)
  Y.n <- example$Y # data matrix
  true.dir <- example$true.dir # adjacency matrix of the data generating DAG
  true.ske <- example$true.ske # skeleton graph of the data generating DAG


  ## (i) DAG learning  by hill climbing: no aggregation

  # learn DAG using "BIC" score
  temp <- score(Y = Y.n, n.boot = 0, score.type = "BIC")
  adj <- temp$adj.matrix

  # Find the total number of skeleton edges
  # and the number of correct skeleton edges of the estimated DAG
  tt <- adj + t(adj) ## symmetrization
  correct.c <- sum((tt > 0) & (true.ske > 0)) / 2
  total.c <- sum(tt > 0) / 2
  total.true <- sum(true.ske > 0) / 2

  return(
    paste(c(total.c, correct.c, total.true))
  )
}

test_that("dag pass", {
  expect_equal(UNITEST_01(), c("326", "96", "109"))
})
