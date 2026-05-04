//4/30/2026: bug fixes, type checks, and efficiency improvements
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"

#include <RcppEigen.h>
#include <RcppNumerical.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Numer;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef std::vector<std::vector<int> > AdjList;

namespace {

// Adjacency convention throughout this file:
// graph[from, to] = 1 means the directed edge from -> to is present.
// Parent lists use the same convention: parents[to] contains all from nodes.
const double kInf = std::numeric_limits<double>::infinity();
const double kNA = std::numeric_limits<double>::quiet_NaN();

inline int mat_index(int from, int to, int p) {
  return from + to * p;
}

inline bool has_edge(const std::vector<unsigned char>& graph, int from, int to, int p) {
  return graph[mat_index(from, to, p)] != 0;
}

inline void set_edge(std::vector<unsigned char>& graph, int from, int to, int p, bool value) {
  graph[mat_index(from, to, p)] = value ? 1 : 0;
}

void require_index(int node, int p, const char* name) {
  if (node < 0 || node >= p) {
    Rcpp::stop("%s is outside valid node index range", name);
  }
}

AdjList graph_to_parents(const std::vector<unsigned char>& graph, int p) {
  AdjList parents(p);
  for (int to = 0; to < p; ++to) {
    for (int from = 0; from < p; ++from) {
      if (has_edge(graph, from, to, p)) {
        parents[to].push_back(from);
      }
    }
  }
  return parents;
}

AdjList graph_to_children(const std::vector<unsigned char>& graph, int p) {
  AdjList children(p);
  for (int from = 0; from < p; ++from) {
    for (int to = 0; to < p; ++to) {
      if (has_edge(graph, from, to, p)) {
        children[from].push_back(to);
      }
    }
  }
  return children;
}

// Check whether adding fromNode -> toNode would create a cycle. This follows
// the original dagbagM logic: trace backward from fromNode through parents and
// ask whether toNode is already an ancestor of fromNode.
bool edge_on_loop_cpp(int fromNode, int toNode, const AdjList& parents) {
  const int p = static_cast<int>(parents.size());
  require_index(fromNode, p, "fromNode");
  require_index(toNode, p, "toNode");

  std::vector<unsigned char> seen(p, 0);
  std::vector<int> path;
  path.reserve(p);
  path.push_back(fromNode);
  seen[fromNode] = 1;

  for (std::size_t pos = 0; pos < path.size(); ++pos) {
    const int cur = path[pos];
    if (cur == toNode) {
      return true;
    }
    const std::vector<int>& curParents = parents[cur];
    for (std::size_t k = 0; k < curParents.size(); ++k) {
      const int parent = curParents[k];
      if (!seen[parent]) {
        seen[parent] = 1;
        path.push_back(parent);
      }
    }
  }

  return false;
}

// Validate that graph has no directed cycles. Checks each existing edge (from,to)
// by asking whether to is already an ancestor of from (i.e., there is a directed
// path to->...->from). If so, the new edge from->to would close a cycle.
// Complexity: O(E*(p+E)) — acceptable because this is only called in
// validate_hc_inputs() (once per run) and in debug mode (never in production).
// For a DAG validity check in a hot path, use a Kahn topological sort instead.
bool is_dag(const std::vector<unsigned char>& graph, int p) {
  AdjList parents = graph_to_parents(graph, p);
  for (int from = 0; from < p; ++from) {
    for (int to = 0; to < p; ++to) {
      if (has_edge(graph, from, to, p) && edge_on_loop_cpp(from, to, parents)) {
        return false;
      }
    }
  }
  return true;
}

// Convert an R logical adjacency matrix into compact two-state C++ storage.
// NA_LOGICAL is rejected at the boundary so the inner HC loop never has to
// reason about R's three-valued logical type.
std::vector<unsigned char> logical_matrix_to_graph(const Rcpp::LogicalMatrix& x,
                                                   const char* name) {
  const int p = x.nrow();
  if (x.ncol() != p) {
    Rcpp::stop("%s must be square", name);
  }

  std::vector<unsigned char> out(p * p, 0);
  for (int to = 0; to < p; ++to) {
    for (int from = 0; from < p; ++from) {
      const int value = x(from, to);
      if (value == NA_LOGICAL) {
        Rcpp::stop("%s must not contain NA values", name);
      }
      out[mat_index(from, to, p)] = value ? 1 : 0;
    }
  }
  return out;
}

// Validate the C++ entry-point inputs even though the R wrappers also validate.
// Keeping both checks makes direct .Call/Rcpp use safer and prevents NA/nonfinite
// values from reaching score calculations or graph state.
void validate_hc_inputs(const Rcpp::NumericMatrix& Y,
                        const Rcpp::CharacterVector& nodeType,
                        const Rcpp::LogicalMatrix& whiteList,
                        const Rcpp::LogicalMatrix& blackList,
                        double tol,
                        int maxStep,
                        int restart,
                        int seed) {
  const int n = Y.nrow();
  const int p = Y.ncol();
  if (n <= 0 || p <= 0) {
    Rcpp::stop("Y must have positive numbers of rows and columns");
  }
  if (p > 46000) {
    Rcpp::stop("p exceeds safe index limit (46000)");
  }
  if (nodeType.size() != p) {
    Rcpp::stop("nodeType must have length equal to ncol(Y)");
  }
  if (whiteList.nrow() != p || whiteList.ncol() != p ||
      blackList.nrow() != p || blackList.ncol() != p) {
    Rcpp::stop("whiteList and blackList must be p by p logical matrices");
  }
  if (!std::isfinite(tol) || tol < 0.0) {
    Rcpp::stop("tol must be a finite nonnegative scalar");
  }
  if (maxStep < 0) {
    Rcpp::stop("maxStep must be nonnegative");
  }
  if (restart < 1) {
    Rcpp::stop("restart must be at least 1");
  }
  if (seed < 0) {
    Rcpp::stop("seed must be nonnegative");
  }
  for (int col = 0; col < p; ++col) {
    const std::string type = Rcpp::as<std::string>(nodeType[col]);
    if (type != "c" && type != "b") {
      Rcpp::stop("nodeType entries must be either 'c' or 'b'");
    }
    for (int row = 0; row < n; ++row) {
      const double value = Y(row, col);
      if (!std::isfinite(value)) {
        Rcpp::stop("Y must not contain NA, NaN, or infinite values");
      }
      if (type == "b" && value != 0.0 && value != 1.0) {
        Rcpp::stop("binary nodes must contain only 0/1 values");
      }
    }
  }

  std::vector<unsigned char> white = logical_matrix_to_graph(whiteList, "whiteList");
  std::vector<unsigned char> black = logical_matrix_to_graph(blackList, "blackList");
  for (int i = 0; i < p; ++i) {
    if (has_edge(white, i, i, p)) {
      Rcpp::stop("whiteList diagonal must be FALSE");
    }
    for (int j = 0; j < p; ++j) {
      if (has_edge(white, i, j, p) && has_edge(black, i, j, p)) {
        Rcpp::stop("whiteList and blackList conflict");
      }
      if (i != j && has_edge(white, i, j, p) && has_edge(white, j, i, p)) {
        Rcpp::stop("whiteList cannot contain both directions of an edge");
      }
    }
  }
  if (!is_dag(white, p)) {
    Rcpp::stop("whiteList must be acyclic");
  }
}

// Fill a pre-allocated workspace with the design matrix for a candidate parent
// set and return q = total columns written.
//
// Column layout (left to right):
//   [extraParent column, if extraParent >= 0]   ← candidate new parent (add op)
//   [existing parents, skipping dropParent]      ← current parents minus dropped (delete op)
//   [ones column]                                ← intercept, always last
//
// q counts every column including the intercept.  The caller passes
// workspace.leftCols(q) to bic_score, which then uses (q - 1) as the BIC
// penalty term — q-1 equals the number of real parent predictors because the
// intercept is not counted as a free structural parameter.
//
// The workspace must be pre-allocated to at least n rows and (p+1) columns.
// This avoids per-call heap allocation inside the hot HC loop.
int fill_design(MatrixXd& workspace, const Rcpp::NumericMatrix& Y,
                const std::vector<int>& parents,
                int extraParent = -1,
                int dropParent = -1) {
  const int n = Y.nrow();
  int col = 0;
  if (extraParent >= 0) {
    if (std::find(parents.begin(), parents.end(), extraParent) != parents.end()) {
      Rcpp::stop("fill_design: extraParent is already in the parent set (internal error)");
    }
    workspace.col(col) = Eigen::Map<const VectorXd>(&Y(0, extraParent), n);
    ++col;
  }
  for (std::size_t k = 0; k < parents.size(); ++k) {
    if (parents[k] == dropParent) {
      continue;
    }
    workspace.col(col) = Eigen::Map<const VectorXd>(&Y(0, parents[k]), n);
    ++col;
  }
  workspace.col(col).setOnes();
  ++col;
  return col;
}

typedef Eigen::Ref<const MatrixXd> ConstMatRef;
typedef Eigen::Ref<const VectorXd> ConstVecRef;

MatrixXd crossprod_self(const ConstMatRef& A) {
  const int n = A.cols();
  return MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint());
}

// Gaussian MLE log-likelihood for the regression Y ~ X*beta + e, e ~ N(0,sigma^2*I).
// beta_hat = (X'X)^{-1} X'Y  (solved via Cholesky of X'X).
// sigma2_hat = ||Y - X*beta_hat||^2 / n  (MLE, divides by n not n-1).
// log L = -n/2 * (log(2*pi) + log(sigma2_hat) + 1).
// Returns -Inf if X'X is singular, if beta is non-finite, or if sigma2 <= 0.
double continuous_loglik(const ConstMatRef& X, const ConstVecRef& Y) {
  const int n = X.rows();
  LLT<MatrixXd> llt(crossprod_self(X));
  if (llt.info() != Eigen::Success) {
    return -kInf;
  }
  const VectorXd beta = llt.solve(X.adjoint() * Y);
  if (!beta.allFinite()) {
    return -kInf;
  }
  const VectorXd resid = Y - X * beta;
  const double s2 = resid.squaredNorm() / static_cast<double>(n);
  if (!std::isfinite(s2) || s2 <= 0.0) {
    return -kInf;
  }
  return -static_cast<double>(n) / 2.0 * (std::log(2.0 * M_PI) + std::log(s2) + 1.0);
}

// Logistic regression objective using Eigen types directly so the design
// matrix can be a workspace block view without R heap allocation.
class LogisticReg: public MFuncGrad {
private:
  const ConstMatRef X;
  const ConstVecRef Y;
  const int n;
  Eigen::VectorXd xbeta;
  Eigen::VectorXd prob;

public:
  LogisticReg(const ConstMatRef& x_, const ConstVecRef& y_) :
    X(x_), Y(y_), n(X.rows()), xbeta(n), prob(n) {}

  double f_grad(Constvec& beta, Refvec grad) {
    xbeta.noalias() = X * beta;
    const double yxbeta = Y.dot(xbeta);
    for (int i = 0; i < n; ++i) {
      prob[i] = R::log1pexp(xbeta[i]);
    }
    const double f = prob.sum() - yxbeta;
    prob = (xbeta - prob).array().exp();
    grad.noalias() = X.transpose() * (prob - Y);
    return f;
  }
};

double logistic_loglik(const ConstMatRef& X, const ConstVecRef& Y) {
  const int q = X.cols();
  LogisticReg nll(X, Y);
  VectorXd beta = VectorXd::Zero(q);
  double fopt = kInf;
  const int status = optim_lbfgs(nll, beta, fopt, 300, 1e-8, 1e-5);
  if (status < 0 || !std::isfinite(fopt)) {
    return -kInf;
  }
  return -fopt;
}

// BIC = -2*logL + log(n) * (q - 1).
// X has q columns: (q-1) parent predictors + 1 intercept.  The BIC penalty
// counts only the (q-1) structural parameters, not the intercept, which is
// always present regardless of parent set size.
// Smaller BIC is better; returns +Inf when loglik is non-finite.
double bic_score(const ConstMatRef& X,
                 const ConstVecRef& y,
                 const std::string& nodeType) {
  const int n = X.rows();
  const int q = X.cols();
  const double loglik = (nodeType == "c") ? continuous_loglik(X, y) : logistic_loglik(X, y);
  if (!std::isfinite(loglik)) {
    return kInf;
  }
  // q - 1 = number of parent predictors (excludes the intercept column).
  const double score = -2.0 * loglik + std::log(static_cast<double>(n)) * (q - 1);
  return std::isfinite(score) ? score : kInf;
}

// Score one node conditional on its current or candidate parent set.
// Smaller BIC is better, so HC accepts negative score changes.
double node_score(MatrixXd& workspace,
                  const Rcpp::NumericMatrix& Y,
                  const std::vector<std::string>& nodeType,
                  const AdjList& parents,
                  int node,
                  int extraParent = -1,
                  int dropParent = -1) {
  const int n = Y.nrow();
  const int q = fill_design(workspace, Y, parents[node], extraParent, dropParent);
  Eigen::Map<const VectorXd> y(&Y(0, node), n);
  return bic_score(workspace.leftCols(q), y, nodeType[node]);
}

// Cache entry for one directed edge candidate.  Fields serve different roles
// depending on the operation being cached:
//
//   Add from->to:
//     scoreA   = BIC score of node `to` after adding `from` as a parent
//     versionA = version[to] at the time scoreA was computed
//     (value = scoreA - curScore[to]; scoreB/versionB unused)
//
//   Delete from->to:
//     scoreA   = BIC score of node `to` after removing `from` from its parents
//     versionA = version[to] at the time scoreA was computed
//     (value = scoreA - curScore[to]; scoreB/versionB unused)
//
//   Reverse from->to (delete from->to, then add to->from):
//     scoreA   = BIC score of `to` after deletion   (= deleteCache[idx].scoreA)
//     scoreB   = BIC score of `from` after adding `to` as a parent
//     versionA = version[to] stamp for scoreA
//     versionB = version[from] stamp for scoreB
//     value    = (scoreA - curScore[to]) + (scoreB - curScore[from])
//
// Cache entries default to versionA=versionB=-1 so they are always stale
// on first access (version[node] initializes to 0, 0 != -1).
struct OneCache {
  double value;
  double scoreA;
  double scoreB;
  int versionA;
  int versionB;
  OneCache() : value(kNA), scoreA(kNA), scoreB(kNA), versionA(-1), versionB(-1) {}
};

struct EdgeCandidate {
  int from;
  int to;
  double sf;
  double gsf;
};

bool better_delta(double candidate,
                  double currentBest,
                  double tol,
                  bool flip,
                  std::mt19937& rng) {
  if (!std::isfinite(candidate)) {
    return false;
  }
  if (candidate < currentBest - tol) {
    return true;
  }
  if (candidate > currentBest + tol) {
    return false;
  }
  if (flip) {
    std::uniform_int_distribution<int> coin(0, 1);
    return coin(rng) == 1;
  }
  return false;
}

// Debug helper: compare cached candidate scores/deltas with full recomputation
// from the current parent sets. This is intentionally only called when
// debug=TRUE because it removes the speed benefit of caching.
bool same_score(double cached, double recomputed) {
  if (!std::isfinite(cached) || !std::isfinite(recomputed)) {
    return !std::isfinite(cached) && !std::isfinite(recomputed);
  }
  const double scale = std::max(1.0, std::max(std::fabs(cached), std::fabs(recomputed)));
  return std::fabs(cached - recomputed) <= 1e-8 * scale;
}

void check_cached_delta(const char* operation,
                        int from,
                        int to,
                        double cachedScoreA,
                        double recomputedScoreA,
                        double cachedScoreB,
                        double recomputedScoreB,
                        double cachedDelta,
                        double recomputedDelta) {
  if (!same_score(cachedScoreA, recomputedScoreA) ||
      !same_score(cachedScoreB, recomputedScoreB) ||
      !same_score(cachedDelta, recomputedDelta)) {
    Rcpp::stop("debug cache check failed for %s %d -> %d", operation, from, to);
  }
}

bool acyclic_add(int from, int to, const AdjList& parents) {
  return !edge_on_loop_cpp(from, to, parents);
}

// Reverse from -> to is checked by temporarily deleting from from parents[to]
// and asking whether adding to -> from would create a cycle.
bool acyclic_reverse(int from, int to, AdjList parents) {
  std::vector<int>& toParents = parents[to];
  std::vector<int>::iterator it = std::find(toParents.begin(), toParents.end(), from);
  if (it == toParents.end()) {
    return false;
  }
  toParents.erase(it);
  return !edge_on_loop_cpp(to, from, parents);
}

// BFS ancestor indicators for node i, including i itself.
std::vector<char> compute_ancestors(int i, const AdjList& parents) {
  const int p = static_cast<int>(parents.size());
  std::vector<char> out(p, 0);
  std::vector<int> queue;
  queue.reserve(p);
  out[i] = 1;
  queue.push_back(i);
  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const std::vector<int>& par = parents[queue[pos]];
    for (std::size_t k = 0; k < par.size(); ++k) {
      if (!out[par[k]]) {
        out[par[k]] = 1;
        queue.push_back(par[k]);
      }
    }
  }
  return out;
}

// BFS descendant indicators for node i, including i itself.
std::vector<char> compute_descendants(int i, const AdjList& children) {
  const int p = static_cast<int>(children.size());
  std::vector<char> out(p, 0);
  std::vector<int> queue;
  queue.reserve(p);
  out[i] = 1;
  queue.push_back(i);
  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const std::vector<int>& chd = children[queue[pos]];
    for (std::size_t k = 0; k < chd.size(); ++k) {
      if (!out[chd[k]]) {
        out[chd[k]] = 1;
        queue.push_back(chd[k]);
      }
    }
  }
  return out;
}

// Acyclicity check for reversing from->to without copying the full AdjList.
// Uses erase+insert to preserve element order so that fill_design column
// ordering is unaffected by the temporary removal.
bool acyclic_reverse_inplace(int from, int to, AdjList& parents) {
  std::vector<int>& toParents = parents[to];
  const std::vector<int>::iterator it = std::find(toParents.begin(), toParents.end(), from);
  if (it == toParents.end()) {
    return false;
  }
  const std::ptrdiff_t pos = it - toParents.begin();
  toParents.erase(it);
  const bool result = !edge_on_loop_cpp(to, from, parents);
  toParents.insert(toParents.begin() + pos, from);
  return result;
}

// State of the last accepted HC operation, including BFS ancestor/descendant
// indicators computed on the post-operation graph. Default-constructed to
// type==-1 meaning no operation has been accepted yet.
struct LastOpState {
  int from = -1, to = -1, type = -1;
  std::vector<char> fromAn, fromDe, toAn, toDe;
};

// Incrementally update acyStatus for one candidate operation given the last
// accepted operation.
//
// Dual-slot convention (one matrix, two semantics):
//   acyStatus(from, to) = true  iff adding from->to is currently acyclic
//   acyStatus(to, from) = true  iff reversing from->to is currently acyclic
//   NA_LOGICAL           = not yet computed (lazy initialization)
//
// operType: 1 = add candidate, 3 = reverse candidate.
// Only the (operFrom, operTo) slot is touched; all other entries are unchanged.
//
// The function propagates the effect of the last accepted operation:
//   type 1 (add):    new path last.from->last.to may close cycles → set false.
//   type 2 (delete): a blocking path was removed → re-check previously false entries.
//   type 3 (reverse): both effects apply; re-check when affected by the deleted
//                     half; set false when affected only by the new edge.
// If neither ancestor/descendant condition is met the entry is left unchanged.
void acyclic_cache_update(const LastOpState& last,
                          int operFrom, int operTo, int operType,
                          AdjList& parents,
                          Rcpp::LogicalMatrix& acyStatus) {
  const int i = operFrom;
  const int j = operTo;

  if (last.type == 1) {
    // Last op was add last.from->last.to. Paths through the new edge can now
    // create cycles for some previously-acyclic candidates.
    if (operType == 1 && acyStatus(i,j) && last.toDe[i] && last.fromAn[j]) {
      acyStatus(i,j) = false;
    }
    if (operType == 3 && acyStatus(j,i) && last.toDe[j] && last.fromAn[i]) {
      if (i != last.from || j != last.to) {
        acyStatus(j,i) = false;
      }
    }
  } else if (last.type == 2) {
    // Last op was delete last.from->last.to. Some previously-cyclic candidates
    // may now be acyclic; re-check those.
    if (operType == 1 && !acyStatus(i,j) && last.toDe[i] && last.fromAn[j]) {
      acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
    }
    if (operType == 3 && !acyStatus(j,i) && last.toDe[j] && last.fromAn[i]) {
      acyStatus(j,i) = acyclic_reverse_inplace(i, j, parents) ? true : false;
    }
  } else {
    // Last op was reverse last.from->last.to (= delete last.from->last.to +
    // add last.to->last.from). The reverse both creates new cycle paths and
    // removes old ones; handle both effects per candidate.
    if (operType == 1) {
      if (acyStatus(i,j) && last.fromDe[i] && last.toAn[j]) {
        // New edge last.to->last.from may create a cycle for this add.
        // If the deleted edge also affected this pair, BFS is authoritative.
        if (last.toDe[i] && last.fromAn[j])
          acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
        else
          acyStatus(i,j) = false;
      } else if (!acyStatus(i,j) && last.toDe[i] && last.fromAn[j]) {
        // Deleted edge may have freed this add; re-check.
        acyStatus(i,j) = acyclic_add(i, j, parents) ? true : false;
      }
    }
    if (operType == 3) {
      const bool isReversedEdge = (i == last.to && j == last.from);
      if (!isReversedEdge && acyStatus(j,i) && last.fromDe[j] && last.toAn[i]) {
        if (last.toDe[j] && last.fromAn[i])
          acyStatus(j,i) = acyclic_reverse_inplace(i, j, parents) ? true : false;
        else
          acyStatus(j,i) = false;
      } else if (!acyStatus(j,i) && last.toDe[j] && last.fromAn[i]) {
        acyStatus(j,i) = acyclic_reverse_inplace(i, j, parents) ? true : false;
      }
    }
  }
}

void erase_value(std::vector<int>& x, int value, const char* context) {
  std::vector<int>::iterator it = std::find(x.begin(), x.end(), value);
  if (it == x.end()) {
    Rcpp::stop("internal graph inconsistency while updating %s", context);
  }
  x.erase(it);
}

Rcpp::LogicalMatrix wrap_graph(const std::vector<unsigned char>& graph, int p) {
  Rcpp::LogicalMatrix out(p, p);
  for (int to = 0; to < p; ++to) {
    for (int from = 0; from < p; ++from) {
      out(from, to) = has_edge(graph, from, to, p);
    }
  }
  return out;
}

Rcpp::IntegerMatrix wrap_int_graph(const std::vector<unsigned char>& graph, int p) {
  Rcpp::IntegerMatrix out(p, p);
  for (int to = 0; to < p; ++to) {
    for (int from = 0; from < p; ++from) {
      out(from, to) = has_edge(graph, from, to, p) ? 1 : 0;
    }
  }
  return out;
}

// Aggregate bootstrap DAGs from an edge-frequency matrix using generalized SHD.
//
// Generalized score function (GSF) for edge from->to:
//   gsf(from, to) = sf(from, to) + (1 - alpha/2) * sf(to, from)
// where sf(i, j) = fraction of bootstrap DAGs containing edge i->j.
//
// With the default alpha=1:  gsf = sf + 0.5 * sf_reverse.
// Both orientations of an undirected skeleton edge enter the candidate list with
// the same gsf value, so the greedy acyclic pass orients the edge toward
// whichever direction the forward sf favors (higher sf wins the sort tie-break).
// freqCutoff acts as a skeleton-inclusion threshold: an undirected edge appears
// iff max(gsf(i,j), gsf(j,i)) > freqCutoff.
//
// Candidates are sorted by (gsf desc, sf desc, from asc, to asc) for
// deterministic output across platforms.
Rcpp::IntegerMatrix aggregate_freq_cpp(const Rcpp::NumericMatrix& seleFreq,
                                       double alpha,
                                       double freqCutoff,
                                       const Rcpp::LogicalMatrix& whiteList,
                                       const Rcpp::LogicalMatrix& blackList,
                                       bool verbose) {
  const int p = seleFreq.nrow();
  if (seleFreq.ncol() != p) {
    Rcpp::stop("freq must be a square matrix");
  }
  if (p > 46000) {
    Rcpp::stop("p exceeds safe index limit (46000)");
  }
  if (whiteList.nrow() != p || whiteList.ncol() != p ||
      blackList.nrow() != p || blackList.ncol() != p) {
    Rcpp::stop("whiteList and blackList must be p by p matrices");
  }
  if (!std::isfinite(alpha) || alpha <= 0.0) {
    Rcpp::stop("alpha must be a positive finite scalar");
  }
  if (!std::isfinite(freqCutoff) || freqCutoff < 0.0 || freqCutoff > 1.0) {
    Rcpp::stop("freqCutoff must be a finite scalar in [0, 1]");
  }

  std::vector<unsigned char> white = logical_matrix_to_graph(whiteList, "whiteList");
  std::vector<unsigned char> black = logical_matrix_to_graph(blackList, "blackList");
  for (int i = 0; i < p; ++i) {
    set_edge(black, i, i, p, true);
    if (has_edge(white, i, i, p)) {
      Rcpp::stop("whiteList diagonal must be FALSE");
    }
    for (int j = 0; j < p; ++j) {
      if (has_edge(white, i, j, p) && has_edge(black, i, j, p)) {
        Rcpp::stop("whiteList and blackList conflict");
      }
      if (i != j && has_edge(white, i, j, p) && has_edge(white, j, i, p)) {
        Rcpp::stop("whiteList cannot contain both directions of an edge");
      }
      const double value = seleFreq(i, j);
      if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        Rcpp::stop("selection frequencies must be finite values in [0, 1]");
      }
    }
  }
  if (!is_dag(white, p)) {
    Rcpp::stop("whitelist must be acyclic");
  }

  std::vector<EdgeCandidate> candidates;
  candidates.reserve(p * p);
  for (int to = 0; to < p; ++to) {
    for (int from = 0; from < p; ++from) {
      if (from == to) {
        continue;
      }
      const double sf = seleFreq(from, to);
      const double gsf = sf + (1.0 - alpha / 2.0) * seleFreq(to, from);
      if (gsf > freqCutoff) {
        EdgeCandidate candidate;
        candidate.from = from;
        candidate.to = to;
        candidate.sf = sf;
        candidate.gsf = gsf;
        candidates.push_back(candidate);
      }
    }
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const EdgeCandidate& a, const EdgeCandidate& b) {
              if (a.gsf != b.gsf) return a.gsf > b.gsf;
              if (a.sf != b.sf) return a.sf > b.sf;
              if (a.from != b.from) return a.from < b.from;
              return a.to < b.to;
            });

  std::vector<unsigned char> result = white;
  AdjList parents = graph_to_parents(result, p);
  for (std::size_t k = 0; k < candidates.size(); ++k) {
    const int from = candidates[k].from;
    const int to = candidates[k].to;
    if (verbose) {
      Rcpp::Rcout << (k + 1) << "th operation: " << (from + 1) << " -> " << (to + 1) << "\n";
    }
    if (has_edge(black, from, to, p) || has_edge(white, from, to, p) ||
        has_edge(result, to, from, p)) {
      continue;
    }
    if (!edge_on_loop_cpp(from, to, parents)) {
      set_edge(result, from, to, p, true);
      parents[to].push_back(from);
    } else if (verbose) {
      Rcpp::Rcout << "not pass acyclic check\n";
    }
  }

  return wrap_int_graph(result, p);
}

} // namespace

bool edgeOnLoop(int fromNode, int toNode, const Rcpp::List& parSet) {
  const int p = parSet.size();
  AdjList parents(p);
  for (int i = 0; i < p; ++i) {
    Rcpp::IntegerVector cur = parSet[i];
    parents[i].reserve(cur.size());
    for (int k = 0; k < cur.size(); ++k) {
      const int parent = cur[k];
      require_index(parent, p, "parent");
      parents[i].push_back(parent);
    }
  }
  return edge_on_loop_cpp(fromNode, toNode, parents);
}

// Run one hill-climbing search from the whitelist-initialized graph.
//
// Operation types in the returned trace:
//   1: add from -> to
//   2: delete from -> to
//   3: reverse from -> to into to -> from
//
// Cache invalidation via version stamps:
//   version[node] is initialized to 0 and incremented every time that node's
//   parent set changes (i.e., after an accepted add, delete, or either leg of
//   a reverse).  Each OneCache entry stores the version[node] at the time its
//   score was computed.  A hit requires the stored stamp to equal the current
//   version; any mismatch triggers a recomputation.
//   Default stamp -1 guarantees a miss on first access (version starts at 0).
//
// Per-step loop structure:
//   Outer loop: to (the node whose BIC score changes for add/delete/reverse-to).
//   Inner loop: from.
//   Branch 1 — edge (from,to) exists:    evaluate delete and (if eligible) reverse.
//   Branch 2 — edge (from,to) absent and edge (to,from) absent: evaluate add.
//
// reverseCache.scoreA reuse:
//   For a reverse of from->to the "new score of to after deletion" is identical
//   to deleteCache.scoreA.  When deleteCache is fresh (same version[to]), that
//   value is copied directly into reverseCache.scoreA without recomputation.
Rcpp::List hc1(const Rcpp::NumericMatrix& Y,
               const Rcpp::CharacterVector& nodeType,
               const Rcpp::LogicalMatrix& whiteList,
               const Rcpp::LogicalMatrix& blackList,
               double tol = 1e-6,
               int maxStep = 500,
               int seed = 1,
               bool flip = true,
               bool verbose = false,
               bool debug = false,
               bool addDeleteOnly = false) {
  validate_hc_inputs(Y, nodeType, whiteList, blackList, tol, maxStep, 1, seed);

  const int n = Y.nrow();
  const int p = Y.ncol();
  std::vector<std::string> types(p);
  for (int i = 0; i < p; ++i) {
    types[i] = Rcpp::as<std::string>(nodeType[i]);
  }

  std::vector<unsigned char> graph = logical_matrix_to_graph(whiteList, "whiteList");
  const std::vector<unsigned char> white = logical_matrix_to_graph(whiteList, "whiteList");
  const std::vector<unsigned char> black = logical_matrix_to_graph(blackList, "blackList");
  AdjList parents = graph_to_parents(graph, p);
  AdjList children = graph_to_children(graph, p);

  MatrixXd workspace(n, p + 1);
  std::vector<double> curScore(p);
  std::vector<int> version(p, 0);
  for (int i = 0; i < p; ++i) {
    curScore[i] = node_score(workspace, Y, types, parents, i);
    if (!std::isfinite(curScore[i])) {
      Rcpp::stop("initial node score is non-finite; check data and whitelist");
    }
  }

  std::vector<OneCache> addCache(p * p);
  std::vector<OneCache> deleteCache(p * p);
  std::vector<OneCache> reverseCache(p * p);
  Rcpp::List stepOper;
  Rcpp::NumericVector stepDelta;
  std::mt19937 rng(static_cast<unsigned int>(seed));

  // Acyclicity status cache. NA_LOGICAL = not yet computed.
  //   acyStatus(from, to) = true iff adding from->to is acyclic
  //   acyStatus(to, from) = true iff reversing from->to is acyclic
  Rcpp::LogicalMatrix acyStatus(p, p);
  std::fill(acyStatus.begin(), acyStatus.end(), NA_LOGICAL);

  // State for incremental acyStatus updates. Ancestor/descendant sets of the
  // last accepted operation's from/to nodes, computed on the post-op graph.
  LastOpState lastOp;

  int acceptedSteps = 0;
  while (acceptedSteps < maxStep) {
    // Each step scans all currently eligible add/delete/reverse candidates and
    // accepts the best score decrease smaller than -tol.
    double bestDelta = 0.0;
    double bestScoreA = kNA;
    double bestScoreB = kNA;
    int bestFrom = -1;
    int bestTo = -1;
    int bestType = -1;

    for (int to = 0; to < p; ++to) {
      for (int from = 0; from < p; ++from) {
        if (from == to) {
          continue;
        }
        const int idx = mat_index(from, to, p);

        if (has_edge(graph, from, to, p)) {
          // Existing edges are acyclic by invariant; set once on first encounter.
          if (acyStatus(from, to) != true) acyStatus(from, to) = true;
          if (has_edge(white, from, to, p)) {
            continue;
          }

          // Deleting from -> to changes only the parent set and local score of
          // the target node to.
          OneCache& del = deleteCache[idx];
          if (del.versionA != version[to]) {
            del.scoreA = node_score(workspace, Y, types, parents, to, -1, from);
            del.value = del.scoreA - curScore[to];
            del.versionA = version[to];
          }
          if (debug) {
            const double recomputedScore = node_score(workspace, Y, types, parents, to, -1, from);
            check_cached_delta("delete", from, to,
                               del.scoreA, recomputedScore,
                               0.0, 0.0,
                               del.value, recomputedScore - curScore[to]);
          }
          if (verbose) {
            Rcpp::Rcout << "delete " << from << "->" << to << ": delta=" << del.value << "\n";
          }
          if (better_delta(del.value, bestDelta, tol, flip, rng)) {
            bestDelta = del.value;
            bestScoreA = del.scoreA;
            bestFrom = from;
            bestTo = to;
            bestType = 2;
          }

          // Reversing from -> to changes two local scores: to after deletion
          // and from after adding to as a parent. addDeleteOnly skips this block.
          if (!addDeleteOnly && !has_edge(black, to, from, p)) {
            // Resolve reversal acyclicity from cache.
            if (acyStatus(to, from) == NA_LOGICAL) {
              acyStatus(to, from) = acyclic_reverse_inplace(from, to, parents) ? true : false;
            } else {
              acyclic_cache_update(lastOp, from, to, 3, parents, acyStatus);
              if (acyStatus(to, from) == NA_LOGICAL) {
                acyStatus(to, from) = acyclic_reverse_inplace(from, to, parents) ? true : false;
              }
            }
            if (debug) {
              const bool bfsResult = acyclic_reverse(from, to, parents);
              if (bfsResult != static_cast<bool>(acyStatus(to, from))) {
                Rcpp::stop("debug: acyStatus mismatch for reverse %d->%d: "
                           "cache=%d actual=%d",
                           from, to,
                           static_cast<int>(static_cast<bool>(acyStatus(to, from))),
                           static_cast<int>(bfsResult));
              }
            }
            if (acyStatus(to, from)) {
              OneCache& rev = reverseCache[idx];
              if (rev.versionA != version[to] || rev.versionB != version[from]) {
                // Reuse del.scoreA: same computation, already fresh this iteration.
                rev.scoreA = del.scoreA;
                rev.scoreB = node_score(workspace, Y, types, parents, from, to, -1);
                rev.value = (rev.scoreA - curScore[to]) + (rev.scoreB - curScore[from]);
                rev.versionA = version[to];
                rev.versionB = version[from];
              }
              if (debug) {
                const double recomputedScoreA = node_score(workspace, Y, types, parents, to, -1, from);
                const double recomputedScoreB = node_score(workspace, Y, types, parents, from, to, -1);
                const double recomputedDelta =
                  (recomputedScoreA - curScore[to]) + (recomputedScoreB - curScore[from]);
                check_cached_delta("reverse", from, to,
                                   rev.scoreA, recomputedScoreA,
                                   rev.scoreB, recomputedScoreB,
                                   rev.value, recomputedDelta);
              }
              if (verbose) {
                Rcpp::Rcout << "reverse " << from << "->" << to << ": delta=" << rev.value << "\n";
              }
              if (better_delta(rev.value, bestDelta, tol, flip, rng)) {
                bestDelta = rev.value;
                bestScoreA = rev.scoreA;
                bestScoreB = rev.scoreB;
                bestFrom = from;
                bestTo = to;
                bestType = 3;
              }
            }
          }
        } else if (!has_edge(graph, to, from, p) && !has_edge(black, from, to, p)) {
          // Resolve add acyclicity from cache.
          if (acyStatus(from, to) == NA_LOGICAL) {
            acyStatus(from, to) = acyclic_add(from, to, parents) ? true : false;
          } else {
            acyclic_cache_update(lastOp, from, to, 1, parents, acyStatus);
            if (acyStatus(from, to) == NA_LOGICAL) {
              acyStatus(from, to) = acyclic_add(from, to, parents) ? true : false;
            }
          }
          if (debug) {
            const bool bfsResult = acyclic_add(from, to, parents);
            if (bfsResult != static_cast<bool>(acyStatus(from, to))) {
              Rcpp::stop("debug: acyStatus mismatch for add %d->%d: "
                         "cache=%d actual=%d",
                         from, to,
                         static_cast<int>(static_cast<bool>(acyStatus(from, to))),
                         static_cast<int>(bfsResult));
            }
          }
          if (!acyStatus(from, to)) {
            continue;
          }

          // Adding from -> to changes only the parent set and local score of to.
          OneCache& add = addCache[idx];
          if (add.versionA != version[to]) {
            add.scoreA = node_score(workspace, Y, types, parents, to, from, -1);
            add.value = add.scoreA - curScore[to];
            add.versionA = version[to];
          }
          if (debug) {
            const double recomputedScore = node_score(workspace, Y, types, parents, to, from, -1);
            check_cached_delta("add", from, to,
                               add.scoreA, recomputedScore,
                               0.0, 0.0,
                               add.value, recomputedScore - curScore[to]);
          }
          if (verbose) {
            Rcpp::Rcout << "add " << from << "->" << to << ": delta=" << add.value << "\n";
          }
          if (better_delta(add.value, bestDelta, tol, flip, rng)) {
            bestDelta = add.value;
            bestScoreA = add.scoreA;
            bestFrom = from;
            bestTo = to;
            bestType = 1;
          }
        }
      }
    }

    if (bestDelta >= -tol || bestType < 0) {
      break;
    }

    if (debug) {
      if (bestType == 1 && (!acyclic_add(bestFrom, bestTo, parents) ||
          has_edge(graph, bestFrom, bestTo, p) || has_edge(graph, bestTo, bestFrom, p))) {
        Rcpp::stop("debug check failed for selected add operation");
      }
      if (bestType == 3 && !acyclic_reverse(bestFrom, bestTo, parents)) {
        Rcpp::stop("debug check failed for selected reverse operation");
      }
    }

    // Apply the selected operation and bump version stamps for every node whose
    // parent set changed. This invalidates stale cache entries lazily.
    if (bestType == 1) {
      set_edge(graph, bestFrom, bestTo, p, true);
      parents[bestTo].push_back(bestFrom);
      children[bestFrom].push_back(bestTo);
      curScore[bestTo] = bestScoreA;
      ++version[bestTo];
    } else if (bestType == 2) {
      set_edge(graph, bestFrom, bestTo, p, false);
      erase_value(parents[bestTo], bestFrom, "parent list");
      erase_value(children[bestFrom], bestTo, "child list");
      curScore[bestTo] = bestScoreA;
      ++version[bestTo];
    } else if (bestType == 3) {
      set_edge(graph, bestFrom, bestTo, p, false);
      erase_value(parents[bestTo], bestFrom, "parent list");
      erase_value(children[bestFrom], bestTo, "child list");
      curScore[bestTo] = bestScoreA;
      ++version[bestTo];

      set_edge(graph, bestTo, bestFrom, p, true);
      parents[bestFrom].push_back(bestTo);
      children[bestTo].push_back(bestFrom);
      curScore[bestFrom] = bestScoreB;
      ++version[bestFrom];
    }

    // Update the last-operation state for incremental acyStatus propagation.
    // Ancestor/descendant sets are computed on the post-operation graph so that
    // acyclic_cache_update in the next step reasons about the updated topology.
    lastOp.from  = bestFrom;
    lastOp.to    = bestTo;
    lastOp.type  = bestType;
    lastOp.fromAn = compute_ancestors(bestFrom, parents);
    lastOp.fromDe = compute_descendants(bestFrom, children);
    lastOp.toAn   = compute_ancestors(bestTo, parents);
    lastOp.toDe   = compute_descendants(bestTo, children);

    stepOper.push_back(Rcpp::IntegerVector::create(bestFrom, bestTo, bestType));
    stepDelta.push_back(bestDelta);
    ++acceptedSteps;
  }

  if (debug && !is_dag(graph, p)) {
    Rcpp::stop("debug check failed: final graph is cyclic");
  }

  return Rcpp::List::create(
    Rcpp::Named("adjacency") = wrap_graph(graph, p),
    Rcpp::Named("score") = Rcpp::wrap(curScore),
    Rcpp::Named("operations") = stepOper,
    Rcpp::Named("deltaMin") = stepDelta,
    Rcpp::Named("steps") = acceptedSteps);
}

// Restart HC with different tie-breaking seeds and keep the fit with the
// smallest total node-wise BIC.
// [[Rcpp::export]]
Rcpp::List hc_(const Rcpp::NumericMatrix& Y,
               const Rcpp::CharacterVector& nodeType,
               const Rcpp::LogicalMatrix& whiteList,
               const Rcpp::LogicalMatrix& blackList,
               double tol = 1e-6,
               int maxStep = 500,
               int restart = 1,
               int seed = 1,
               bool verbose = false,
               bool debug = false,
               bool addDeleteOnly = false) {
  validate_hc_inputs(Y, nodeType, whiteList, blackList, tol, maxStep, restart, seed);

  Rcpp::List bestRes;
  double bestScore = kInf;
  const bool flip = restart > 1;
  for (int i = 0; i < restart; ++i) {
    // Space restart seeds by 101 to ensure independence across restarts.
    Rcpp::List curRes = hc1(Y, nodeType, whiteList, blackList, tol, maxStep,
                            seed + i * 101, flip, verbose, debug, addDeleteOnly);
    Rcpp::NumericVector curScoreVector = curRes["score"];
    const double curScore = Rcpp::sum(curScoreVector);
    if (!std::isfinite(curScore)) {
      continue;
    }
    if (curScore < bestScore) {
      bestRes = curRes;
      bestScore = curScore;
    }
  }
  if (bestScore == kInf) {
    Rcpp::stop("all HC restarts produced non-finite scores");
  }
  return bestRes;
}

// Aggregate directly from an edge-frequency matrix. This is used by the new
// frequency-only bootstrap path.
// [[Rcpp::export]]
Rcpp::IntegerMatrix score_shd_freq_cpp(const Rcpp::NumericMatrix& freq,
                                       double alpha,
                                       double freqCutoff,
                                       const Rcpp::LogicalMatrix& whiteList,
                                       const Rcpp::LogicalMatrix& blackList,
                                       bool verbose = false) {
  Rcpp::NumericMatrix cleanFreq = Rcpp::clone(freq);
  const int p = cleanFreq.nrow();
  if (cleanFreq.ncol() != p) {
    Rcpp::stop("freq must be a square matrix");
  }
  for (int i = 0; i < p; ++i) {
    cleanFreq(i, i) = 0.0;
  }
  return aggregate_freq_cpp(cleanFreq, alpha, freqCutoff, whiteList, blackList, verbose);
}

// Aggregate from a p x p x B bootstrap adjacency array by first computing edge
// selection frequencies in C++ and then using the same greedy aggregator.
// [[Rcpp::export]]
Rcpp::IntegerMatrix score_shd_cpp(const Rcpp::NumericVector& bootAdj,
                                  double alpha,
                                  double freqCutoff,
                                  const Rcpp::LogicalMatrix& whiteList,
                                  const Rcpp::LogicalMatrix& blackList,
                                  bool verbose = false) {
  Rcpp::IntegerVector dims = bootAdj.attr("dim");
  if (dims.size() != 3) {
    Rcpp::stop("boot.adj must be a p by p by B array");
  }
  const int p = dims[0];
  const int p2 = dims[1];
  const int nb = dims[2];
  if (p <= 0 || p2 != p || nb <= 0) {
    Rcpp::stop("boot.adj must be a p by p by B array with B >= 1");
  }

  Rcpp::NumericMatrix freq(p, p);
  const int sliceSize = p * p;
  for (int b = 0; b < nb; ++b) {
    for (int to = 0; to < p; ++to) {
      for (int from = 0; from < p; ++from) {
        const double value = bootAdj[from + to * p + b * sliceSize];
        if (!std::isfinite(value) || (value != 0.0 && value != 1.0)) {
          Rcpp::stop("boot.adj entries must be finite 0/1 values");
        }
        freq(from, to) += value / static_cast<double>(nb);
      }
    }
  }
  for (int i = 0; i < p; ++i) {
    freq(i, i) = 0.0;
  }
  return aggregate_freq_cpp(freq, alpha, freqCutoff, whiteList, blackList, verbose);
}

#pragma GCC diagnostic pop
