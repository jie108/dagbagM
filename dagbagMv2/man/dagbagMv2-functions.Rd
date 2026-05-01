\name{dagbagMv2-functions}
\alias{hc}
\alias{hc_boot}
\alias{hc_boot_parallel}
\alias{score_shd}
\alias{score_shd_freq}
\alias{cal_order}
\alias{moral_graph}
\alias{vstructures}
\alias{skeleton}
\alias{compare.vstructures}
\title{dagbagMv2 user functions}
\description{
User-facing helpers for direct hill-climbing DAG learning, bootstrap DAG
learning, bootstrap aggregation, and graph summaries.
}
\usage{
hc(Y, nodeType = NULL, whiteList = NULL, blackList = NULL,
   standardize = TRUE, tol = 1e-6, maxStep = 2000L,
   restart = 1L, seed = 1L, verbose = FALSE, debug = FALSE,
   addDeleteOnly = FALSE)

hc_boot(Y, n.boot = 1L, nodeType = NULL, whiteList = NULL,
        blackList = NULL, standardize = TRUE, tol = 1e-6,
        maxStep = 2000L, restart = 1L, seed = 1L,
        nodeShuffle = TRUE, verbose = FALSE, debug = FALSE,
        addDeleteOnly = FALSE,
        return = c("array", "freq", "both"))

hc_boot_parallel(Y, n.boot = 1L, nodeType = NULL, whiteList = NULL,
                 blackList = NULL, standardize = TRUE, tol = 1e-6,
                 maxStep = 2000L, restart = 1L, seed = 1L,
                 nodeShuffle = TRUE, numThread = 2L,
                 verbose = FALSE, debug = FALSE,
                 addDeleteOnly = FALSE,
                 return = c("array", "freq", "both"))

score_shd(boot.adj, alpha = 1, threshold = 0, whitelist = NULL,
          blacklist = NULL, max.step = NULL, verbose = FALSE)

score_shd_freq(freq, alpha = 1, threshold = 0, whitelist = NULL,
               blacklist = NULL, max.step = NULL, verbose = FALSE)
}
\arguments{
\item{Y}{Numeric data matrix.}
\item{nodeType}{Character vector with entries \code{"c"} for continuous and
\code{"b"} for binary nodes.}
\item{whiteList, whitelist}{Logical or 0/1 matrix of required directed edges.}
\item{blackList, blacklist}{Logical or 0/1 matrix of forbidden directed edges.}
\item{standardize}{Whether to standardize continuous nodes before fitting.}
\item{tol}{Minimum BIC improvement threshold.}
\item{maxStep, max.step}{Maximum hill-climbing steps; \code{max.step} is a
legacy aggregation argument.}
\item{restart}{Number of hill-climbing restarts.}
\item{seed}{Seed for bootstrap sampling, shuffling, and tie-breaking.}
\item{verbose}{Whether to print progress.}
\item{debug}{Whether to run additional HC validity checks.}
\item{addDeleteOnly}{Whether to skip reversal candidates in HC.}
\item{n.boot}{Number of bootstrap resamples.}
\item{nodeShuffle}{Whether to shuffle node order within each bootstrap fit.}
\item{numThread}{Number of workers for parallel bootstrap fitting.}
\item{return}{Bootstrap output mode: full array, edge frequencies, or both.}
\item{boot.adj}{A \code{p by p by B} array of bootstrap DAGs.}
\item{freq}{A \code{p by p} edge selection frequency matrix.}
\item{alpha}{Generalized SHD reversal weight parameter.}
\item{threshold}{Frequency cutoff adjustment.}
}
\value{
\code{hc} returns a list with adjacency, node scores, operations, deltas, and
step count. Bootstrap functions return an adjacency array, a frequency matrix,
or both. Aggregation functions return an adjacency matrix.
}
\keyword{models}
