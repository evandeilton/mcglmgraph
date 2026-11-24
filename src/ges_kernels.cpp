// ges_kernels.cpp: greedy equivalence search (RcppArmadillo).
// Part of the mcglmgraph package. License: GPL-3.
//
// The stepwise search of the undirected layer moves over supports; this one
// moves over equivalence classes, which is the space the second-order
// separations identify. Because the working criterion is score equivalent, two
// members of a class score alike and a search over classes is well posed.
//
// The two phases are those of Chickering (2002). The forward phase starts from
// the empty class and applies, at each step, the insertion that raises the
// score most; the backward phase starts from what the forward phase left and
// applies the best deletion. Each operator acts on a CPDAG and its result is
// mapped back to a CPDAG, so the search never leaves the space.
//
//   Insert(X, Y, T)  X and Y nonadjacent, T a subset of the neighbors of Y not
//                    adjacent to X. Valid when NA(Y,X) union T is a clique and
//                    every semi-directed path from Y to X meets it. Writes
//                    X -> Y and T' -> Y for T' in T.
//   Delete(X, Y, H)  X and Y adjacent, H a subset of NA(Y,X). Valid when
//                    NA(Y,X) minus H is a clique. Removes the edge and writes
//                    Y -> H' and X -> H'.
//
// The score is the Gaussian working criterion of the residual layer with a
// penalty on the number of parameters. It is decomposable: the contribution of
// a node given its parents is read off the residual variance of its regression
// on them, which is why a single insertion or deletion is scored by two local
// terms and not by refitting the whole layer.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "equilibrated_solve.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace Rcpp;

namespace {

typedef std::vector<arma::uword> uvec_t;

inline bool has_edge(const arma::umat& A, arma::uword a, arma::uword b) {
  return A(a, b) == 1 || A(b, a) == 1;
}
inline bool is_undir(const arma::umat& A, arma::uword a, arma::uword b) {
  return A(a, b) == 1 && A(b, a) == 1;
}
inline bool is_dir(const arma::umat& A, arma::uword a, arma::uword b) {
  return A(a, b) == 1 && A(b, a) == 0;
}

// Local score of node j given a parent set: the Gaussian working criterion
// contribution, penalized. Larger is better.
double local_score(const arma::mat& S, double n, arma::uword j,
                   const uvec_t& pa, double penalty) {
  double resid = S(j, j);
  if (!pa.empty()) {
    arma::uvec P(pa);
    arma::uvec jj(1); jj(0) = j;
    arma::mat A = S.submat(P, P);
    A = 0.5 * (A + A.t());
    arma::vec b = S.submat(P, jj);
    // A near-singular parent block, judged on the equilibrated block so that
    // the test does not depend on scale, has no identified regression: the
    // parent set is refused.
    const double rcond_tol = std::pow(arma::datum::eps, 2.0 / 3.0);
    arma::vec g;
    bool ill = false, ridged = false;
    if (!mcglmgraph::equilibrated_solve(A, b, rcond_tol, 0.0, false, g, ill,
                                        ridged) || ill)
      return -std::numeric_limits<double>::infinity();
    resid = S(j, j) - arma::dot(b, g);
  }
  if (!(resid > 0.0)) return -std::numeric_limits<double>::infinity();
  // -n/2 log(sigma^2) minus the penalty on the parameters of this node
  return -0.5 * n * std::log(resid) - 0.5 * penalty * (double)(pa.size() + 1);
}

// Neighbors of y joined by an undirected edge.
uvec_t neighbors(const arma::umat& A, arma::uword y) {
  uvec_t out;
  for (arma::uword v = 0; v < A.n_rows; ++v)
    if (v != y && is_undir(A, y, v)) out.push_back(v);
  return out;
}

// Parents of y: arcs into y.
uvec_t parents(const arma::umat& A, arma::uword y) {
  uvec_t out;
  for (arma::uword v = 0; v < A.n_rows; ++v)
    if (v != y && is_dir(A, v, y)) out.push_back(v);
  return out;
}

bool is_clique(const arma::umat& A, const uvec_t& S) {
  for (std::size_t a = 0; a + 1 < S.size(); ++a)
    for (std::size_t b = a + 1; b < S.size(); ++b)
      if (!has_edge(A, S[a], S[b])) return false;
  return true;
}

// Is there a semi-directed path from y to x avoiding every node of block?
// A semi-directed path follows undirected edges or arcs in the direction of
// travel.
bool semi_directed_path(const arma::umat& A, arma::uword y, arma::uword x,
                        const std::vector<char>& block) {
  const arma::uword R = A.n_rows;
  std::vector<char> seen(R, 0);
  std::queue<arma::uword> q;
  seen[y] = 1; q.push(y);
  while (!q.empty()) {
    const arma::uword v = q.front(); q.pop();
    for (arma::uword w = 0; w < R; ++w) {
      if (seen[w] || block[w]) continue;
      if (A(v, w) != 1) continue;              // v -> w or v - w
      if (w == x) return true;
      seen[w] = 1; q.push(w);
    }
  }
  return false;
}

// Enumerate the subsets of a set, as index masks.
void subsets_of(std::size_t m, std::vector<unsigned long long>& out) {
  const unsigned long long total = 1ULL << m;
  out.clear();
  out.reserve(total);
  for (unsigned long long k = 0; k < total; ++k) out.push_back(k);
}

}  // namespace

//' Local score of a node given a parent set (C++ kernel)
//'
//' Returns -(n/2) log(r) - (penalty/2) (|pa| + 1), where r = S[j, j] -
//' S[j, pa] S[pa, pa]^-1 S[pa, j] is the residual variance of the regression of
//' node j on its parents. It is -Inf when the parent block is numerically
//' singular or r is not positive.
//'
//' @param S covariance of the responses.
//' @param n number of units.
//' @param j node, 1-based.
//' @param pa parents, 1-based, possibly empty.
//' @param penalty multiplier on the parameter count; log(n) gives the BIC.
//' @return the penalized Gaussian working contribution of the node.
//' @keywords internal
// [[Rcpp::export]]
double cpp_ges_local_score(const arma::mat& S, double n, arma::uword j,
                           const arma::uvec& pa, double penalty) {
  uvec_t p(pa.n_elem);
  for (arma::uword t = 0; t < pa.n_elem; ++t) p[t] = pa(t) - 1;
  return local_score(S, n, j - 1, p, penalty);
}

//' Best insertion available on a CPDAG (C++ kernel)
//'
//' Scans the operators Insert(X, Y, T) of the greedy equivalence search: X and Y
//' nonadjacent, T a subset of the neighbors of Y not adjacent to X, valid when
//' NA(Y,X) union T is a clique and every semi-directed path from Y to X meets
//' it. The operator writes X -> Y and T' -> Y for T' in T, and its gain is the
//' difference of the local scores of Y with and without X. The operator with the
//' largest gain above 1e-10 is returned.
//'
//' @param S covariance of the responses.
//' @param n number of units.
//' @param amat the current CPDAG.
//' @param penalty multiplier on the parameter count.
//' @param max_set largest T enumerated, as a guard against the exponential.
//' @return list with the gain, the nodes x and y, the set T, and a flag saying
//'   whether any valid insertion was found.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ges_best_insert(const arma::mat& S, double n,
                               const arma::umat& amat, double penalty,
                               int max_set = 12) {
  const arma::uword R = amat.n_rows;
  double best = 0.0;
  bool found = false;
  arma::uword bx = 0, by = 0;
  uvec_t bT;
  std::vector<unsigned long long> masks;

  for (arma::uword y = 0; y < R; ++y) {
    const uvec_t na_y = neighbors(amat, y);
    const uvec_t pa_y = parents(amat, y);
    for (arma::uword x = 0; x < R; ++x) {
      if (x == y || has_edge(amat, x, y)) continue;
      // T ranges over the neighbors of y that are not adjacent to x
      uvec_t cand, na_yx;
      for (std::size_t t = 0; t < na_y.size(); ++t) {
        if (has_edge(amat, na_y[t], x)) na_yx.push_back(na_y[t]);
        else cand.push_back(na_y[t]);
      }
      if ((int)cand.size() > max_set) continue;
      subsets_of(cand.size(), masks);
      for (std::size_t k = 0; k < masks.size(); ++k) {
        uvec_t T;
        for (std::size_t t = 0; t < cand.size(); ++t)
          if (masks[k] & (1ULL << t)) T.push_back(cand[t]);

        uvec_t cond = na_yx;
        cond.insert(cond.end(), T.begin(), T.end());
        if (!is_clique(amat, cond)) continue;
        std::vector<char> block(R, 0);
        for (std::size_t t = 0; t < cond.size(); ++t) block[cond[t]] = 1;
        block[y] = 1;                       // paths start at y
        if (semi_directed_path(amat, y, x, block)) continue;

        uvec_t base = pa_y;
        base.insert(base.end(), cond.begin(), cond.end());
        std::sort(base.begin(), base.end());
        base.erase(std::unique(base.begin(), base.end()), base.end());
        uvec_t with = base;
        with.push_back(x);
        const double gain = local_score(S, n, y, with, penalty) -
          local_score(S, n, y, base, penalty);
        if (gain > best + 1e-10) {
          best = gain; found = true; bx = x; by = y; bT = T;
        }
      }
    }
  }
  arma::uvec Tout(bT.size());
  for (std::size_t t = 0; t < bT.size(); ++t) Tout(t) = bT[t] + 1;
  return Rcpp::List::create(_["found"] = found, _["gain"] = best,
                            _["x"] = (double)(bx + 1), _["y"] = (double)(by + 1),
                            _["set"] = Tout);
}

//' Best deletion available on a CPDAG (C++ kernel)
//'
//' Scans the operators Delete(X, Y, H): X and Y adjacent by X - Y or X -> Y, H a
//' subset of NA(Y,X), valid when NA(Y,X) minus H is a clique. The operator
//' removes the edge and writes Y -> H' and X -> H' for H' in H. The operator
//' with the largest gain above 1e-10 is returned.
//'
//' @inheritParams cpp_ges_best_insert
//' @return list with the gain, the nodes x and y, the set H, and a flag.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ges_best_delete(const arma::mat& S, double n,
                               const arma::umat& amat, double penalty,
                               int max_set = 12) {
  const arma::uword R = amat.n_rows;
  double best = 0.0;
  bool found = false;
  arma::uword bx = 0, by = 0;
  uvec_t bH;
  std::vector<unsigned long long> masks;

  for (arma::uword y = 0; y < R; ++y) {
    const uvec_t pa_y = parents(amat, y);
    for (arma::uword x = 0; x < R; ++x) {
      if (x == y) continue;
      // the edge must be x - y or x -> y
      if (!(is_undir(amat, x, y) || is_dir(amat, x, y))) continue;
      uvec_t na_yx;
      const uvec_t na_y = neighbors(amat, y);
      for (std::size_t t = 0; t < na_y.size(); ++t)
        if (has_edge(amat, na_y[t], x)) na_yx.push_back(na_y[t]);
      if ((int)na_yx.size() > max_set) continue;
      subsets_of(na_yx.size(), masks);
      for (std::size_t k = 0; k < masks.size(); ++k) {
        uvec_t H, keep;
        for (std::size_t t = 0; t < na_yx.size(); ++t) {
          if (masks[k] & (1ULL << t)) H.push_back(na_yx[t]);
          else keep.push_back(na_yx[t]);
        }
        if (!is_clique(amat, keep)) continue;

        uvec_t base = pa_y;
        base.insert(base.end(), na_yx.begin(), na_yx.end());
        std::sort(base.begin(), base.end());
        base.erase(std::unique(base.begin(), base.end()), base.end());
        // remove H from the base, and x as well for the reduced set
        std::vector<char> drop(R, 0);
        for (std::size_t t = 0; t < H.size(); ++t) drop[H[t]] = 1;
        uvec_t before, after;
        for (std::size_t t = 0; t < base.size(); ++t) {
          if (drop[base[t]]) continue;
          before.push_back(base[t]);
          if (base[t] != x) after.push_back(base[t]);
        }
        if (before.size() == after.size()) before.push_back(x);
        const double gain = local_score(S, n, y, after, penalty) -
          local_score(S, n, y, before, penalty);
        if (gain > best + 1e-10) {
          best = gain; found = true; bx = x; by = y; bH = H;
        }
      }
    }
  }
  arma::uvec Hout(bH.size());
  for (std::size_t t = 0; t < bH.size(); ++t) Hout(t) = bH[t] + 1;
  return Rcpp::List::create(_["found"] = found, _["gain"] = best,
                            _["x"] = (double)(bx + 1), _["y"] = (double)(by + 1),
                            _["set"] = Hout);
}

//' A DAG consistent with a partially directed graph (C++ kernel)
//'
//' The algorithm of Dor and Tarsi (1992): repeatedly take a node with no
//' outgoing arc whose undirected neighbors are adjacent to all its other
//' neighbors, orient every undirected edge at it inwards, and remove it.
//'
//' @param amat the partially directed graph.
//' @return list with the DAG and a flag saying whether an extension exists.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_pdag_to_dag(const arma::umat& amat) {
  const arma::uword R = amat.n_rows;
  arma::umat work = amat;
  arma::umat out = amat;
  std::vector<char> alive(R, 1);
  arma::uword left = R;

  while (left > 0) {
    bool moved = false;
    for (arma::uword x = 0; x < R && !moved; ++x) {
      if (!alive[x]) continue;
      bool sink = true;
      for (arma::uword w = 0; w < R; ++w)
        if (alive[w] && w != x && is_dir(work, x, w)) { sink = false; break; }
      if (!sink) continue;
      uvec_t und, nb;
      for (arma::uword w = 0; w < R; ++w) {
        if (!alive[w] || w == x) continue;
        if (is_undir(work, x, w)) und.push_back(w);
        if (has_edge(work, x, w)) nb.push_back(w);
      }
      bool ok = true;
      for (std::size_t a = 0; a < und.size() && ok; ++a)
        for (std::size_t b = 0; b < nb.size(); ++b) {
          if (nb[b] == und[a]) continue;
          if (!has_edge(work, und[a], nb[b])) { ok = false; break; }
        }
      if (!ok) continue;
      for (std::size_t a = 0; a < und.size(); ++a) {
        out(und[a], x) = 1; out(x, und[a]) = 0;      // orient into x
      }
      alive[x] = 0; --left; moved = true;
      for (arma::uword w = 0; w < R; ++w) { work(x, w) = 0; work(w, x) = 0; }
    }
    if (!moved) return Rcpp::List::create(_["dag"] = out, _["ok"] = false);
  }
  return Rcpp::List::create(_["dag"] = out, _["ok"] = true);
}
