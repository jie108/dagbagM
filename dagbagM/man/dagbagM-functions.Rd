\name{dagbagM-functions}
\alias{hc}
\alias{hc_boot}
\alias{score_shd}
\alias{score_shd_freq}
\alias{cal_order}
\alias{moral_graph}
\alias{vstructures}
\alias{skeleton}
\alias{compare.vstructures}
\title{dagbagM user functions}
\description{
User-facing helpers for direct hill-climbing DAG learning, bootstrap DAG
learning, bootstrap aggregation, and graph summaries.
}
\usage{
hc(Y, nodeType = NULL, whiteList = NULL, blackList = NULL,
   standardize = TRUE, tol = 1e-6, maxStep = 2000L,
   restart = 1L, seed = 1L, verbose = FALSE, debug = FALSE,
   addDeleteOnly = FALSE)

hc_boot(Y, n.boot = 100L, nodeType = NULL, whiteList = NULL,
        blackList = NULL, standardize = TRUE, tol = 1e-6,
        maxStep = 2000L, restart = 1L, seed = 1L,
        nodeShuffle = TRUE,
        backend = c("sequential", "future"),
        workers = NULL,
        verbose = FALSE, debug = FALSE,
        addDeleteOnly = FALSE,
        output_type = c("array", "freq", "both"))

score_shd(boot.adj, alpha = 1, freq.cutoff = 0.5, whiteList = NULL,
          blackList = NULL, maxStep = NULL, verbose = FALSE)

score_shd_freq(freq, alpha = 1, freq.cutoff = 0.5, whiteList = NULL,
               blackList = NULL, maxStep = NULL, verbose = FALSE)
}
\arguments{
\item{Y}{Numeric data matrix.}
\item{nodeType}{Character vector with entries \code{"c"} for continuous and
\code{"b"} for binary nodes.}
\item{whiteList}{Logical or 0/1 matrix of required directed edges.}
\item{blackList}{Logical or 0/1 matrix of forbidden directed edges.}
\item{standardize}{Whether to standardize continuous nodes before fitting.}
\item{tol}{Minimum BIC improvement threshold.}
\item{maxStep}{Maximum hill-climbing steps. In \code{score_shd} and
\code{score_shd_freq} this is a deprecated legacy argument with no effect.}
\item{restart}{Number of hill-climbing restarts.}
\item{seed}{Seed for bootstrap sampling, node shuffling, and HC tie-breaking.
All randomness is derived from this one value so sequential and future
backends produce identical results.}
\item{verbose}{Whether to print progress.}
\item{debug}{Whether to run additional HC validity checks.}
\item{addDeleteOnly}{Whether to skip reversal candidates in HC.}
\item{n.boot}{Number of bootstrap resamples.}
\item{nodeShuffle}{Whether to randomize node visit order within each bootstrap
fit to reduce greedy-search ordering bias.}
\item{backend}{Execution backend: \code{"sequential"} (default) runs replicates
one by one; \code{"future"} dispatches them in parallel using
\code{future::multisession}.}
\item{workers}{Number of parallel workers for \code{backend = "future"}.
\code{NULL} (default) uses all available cores minus one, capped at
\code{n.boot}.}
\item{output_type}{Bootstrap output mode: \code{"array"} returns the full
\code{p x p x B} adjacency array; \code{"freq"} returns only the edge
selection frequency matrix (saves memory for large \code{B});
\code{"both"} returns a list with both.}
\item{boot.adj}{A \code{p by p by B} array of bootstrap DAGs.}
\item{freq}{A \code{p by p} edge selection frequency matrix.}
\item{alpha}{Generalized SHD reversal weight parameter.}
\item{freq.cutoff}{Edge inclusion threshold in [0, 1]: an edge is a candidate
iff its generalized selection frequency exceeds this value.
Default 0.5 corresponds to majority-vote skeleton selection.}
}
\value{
\code{hc} returns a list with adjacency, node scores, operations, deltas, and
step count. \code{hc_boot} returns an adjacency array, a frequency matrix, or
a list with both, depending on \code{output_type}. Aggregation functions return
an integer adjacency matrix.
}
\examples{
## Generate small mixed data: 3 continuous + 1 binary node
set.seed(42)
n <- 80
x1 <- rnorm(n)
x2 <- 0.8 * x1 + rnorm(n, 0, 0.6)
x3 <- -0.4 * x1 + 0.7 * x2 + rnorm(n, 0, 0.7)
x4 <- rbinom(n, 1, plogis(0.9 * x1 - 0.6 * x2))
Y <- cbind(x1, x2, x3, x4)
nodeType <- c("c", "c", "c", "b")

## Single hill-climbing fit
fit <- hc(Y, nodeType)
fit$adjacency

## Bootstrap aggregation (small B for speed)
boot_freq <- hc_boot(Y, n.boot = 5, nodeType = nodeType, output_type = "freq")
agg <- score_shd_freq(boot_freq)
agg

## Graph summaries
cal_order(agg)
skeleton(agg)
moral_graph(agg)
vstructures(agg)
}
\keyword{models}
