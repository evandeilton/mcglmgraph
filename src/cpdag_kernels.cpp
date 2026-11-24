// cpdag_kernels.cpp: exhaustive enumeration of separating sets for the
// order-independent orientation rules (RcppArmadillo).
// Part of the mcglmgraph package. License: GPL-3.
//
// The skeleton search of mc_pc_skeleton() stops at the FIRST conditioning set
// that separates a pair and records only that set. That is enough to delete the
// edge and it is what the standard v-structure rule uses, but it is not enough
// for the conservative rule of Ramsey, Zhang and Spirtes (2006) or the majority
// rule of Colombo and Maathuis (2014): both ask whether the center of an
// unshielded triple belongs to none, to all, or to some of the sets that
// separate the endpoints, and the answer needs the whole collection.
//
// The collection is defined on the FINAL skeleton, not on the neighborhoods
// that happened to be live when the edge fell, so it is rebuilt here from
// scratch: for the pair (i, j), every subset of adj(i) \ {j} and every subset of
// adj(j) \ {i} of size at most max_ord is tested, and the separating ones are
// returned, deduplicated across the two neighborhoods.
//
// The conditional-association test is the same second-order one used by the
// skeleton: the partial correlation read off the inverse of the (i, j, S) block
// with Fisher's z transform, under the degrees-of-freedom contract
// n - |S| - 3 for raw units, or nu_e - |S| - 2 when the residual degrees of
// freedom of the first stage are supplied.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <set>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Rcpp;

namespace {

const double RCOND_TOL = std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0);
const double PCOR_CLAMP = 0.999999;

struct CiOut {
  double pcor;
  double p_value;
  double df;
  double ridge;
  bool   regularized;
  bool   usable;          // false when the block gives no valid test
};

// Partial correlation of (i, j) given S, with Fisher's z p-value, the ridge
// fallback for a near-singular block and the clamp. The single implementation.
CiOut ci_pcor_block(const arma::mat& Sigma, arma::uword i, arma::uword j,
                    const std::vector<arma::uword>& S,
                    double n, double df_residual, bool use_dfr) {
  CiOut out;
  out.ridge = 0.0;
  out.regularized = false;
  out.usable = true;
  out.pcor = NA_REAL;
  out.p_value = NA_REAL;

  const arma::uword m = S.size();
  out.df = use_dfr ? (df_residual - double(m) - 2.0)
                   : (n - double(m) - 3.0);
  if (!std::isfinite(out.df) || out.df <= 0.0) { out.usable = false; return out; }

  arma::uvec idx(m + 2);
  idx(0) = i; idx(1) = j;
  for (arma::uword s = 0; s < m; ++s) idx(s + 2) = S[s];
  arma::mat sub = Sigma.submat(idx, idx);
  sub = 0.5 * (sub + sub.t());

  arma::mat P;
  bool ok = false;
  double rc = NA_REAL;
  try { rc = arma::rcond(sub); } catch (...) { rc = NA_REAL; }
  if (std::isfinite(rc) && rc > RCOND_TOL) ok = arma::inv_sympd(P, sub);
  if (!ok || !P.is_finite()) {
    const double scale = std::max(std::max(arma::abs(sub.diag()).max(),
                                           arma::abs(sub).max()), 1.0);
    out.ridge = 1e-8 * scale;
    ok = arma::inv_sympd(P, sub + out.ridge * arma::eye(sub.n_rows, sub.n_rows));
    out.regularized = true;
  }
  if (!ok || !P.is_finite()) { out.usable = false; return out; }

  const double den = P(0, 0) * P(1, 1);
  if (!std::isfinite(den) || den <= 0.0) { out.usable = false; return out; }

  double pc = -P(0, 1) / std::sqrt(den);
  pc = std::max(std::min(pc, PCOR_CLAMP), -PCOR_CLAMP);
  const double z = 0.5 * std::log((1.0 + pc) / (1.0 - pc)) * std::sqrt(out.df);
  out.pcor = pc;
  out.p_value = 2.0 * R::pnorm(-std::fabs(z), 0.0, 1.0, 1, 0);
  return out;
}

// Enumerate the subsets of nb of size 0..top in lexicographic order of indices.
template <typename F>
void for_each_subset(const std::vector<arma::uword>& nb, int top, F&& f) {
  const int m = static_cast<int>(nb.size());
  const int hi = std::min(top, m);
  std::vector<arma::uword> S;
  for (int s = 0; s <= hi; ++s) {
    std::vector<int> c(static_cast<std::size_t>(s));
    for (int t = 0; t < s; ++t) c[t] = t;
    while (true) {
      S.clear();
      S.reserve(static_cast<std::size_t>(s));
      for (int t = 0; t < s; ++t) S.push_back(nb[c[t]]);
      f(S);
      int t = s - 1;
      while (t >= 0 && c[t] == m - s + t) --t;
      if (t < 0) break;
      ++c[t];
      for (int u = t + 1; u < s; ++u) c[u] = c[u - 1] + 1;
    }
  }
}

}  // namespace

//' Second-order conditional-association test of one pair (C++ kernel)
//'
//' The partial correlation of (i, j) given S, read off the inverse of the
//' (i, j, S) block, with the p-value of Fisher's z, z = atanh(pcor) sqrt(df),
//' where df = n - |S| - 3 or df_residual - |S| - 2. The partial correlation is
//' clamped to [-0.999999, 0.999999], and a nearly singular block is inverted
//' after a ridge, which is reported. It is the test the skeleton and the
//' enumeration of separating sets both apply.
//'
//' @param Sigma R x R covariance matrix of the responses.
//' @param i,j the pair, 1-based.
//' @param S conditioning set, 1-based, possibly empty.
//' @param n number of experimental units.
//' @param df_residual residual degrees of freedom of the first stage; a value
//'   not greater than zero selects the raw contract n - |S| - 3.
//' @return list with pcor, p_value, df, ridge, regularized and usable, the
//'   last one false when the block gives no valid test.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ci_pcor(const arma::mat& Sigma, arma::uword i, arma::uword j,
                       const arma::uvec& S, double n,
                       double df_residual = 0.0) {
  const arma::uword R = Sigma.n_rows;
  if (Sigma.n_cols != R) stop("Sigma must be square.");
  if (i < 1 || j < 1 || i > R || j > R || i == j)
    stop("the pair holds an invalid node index.");
  std::vector<arma::uword> cond(S.n_elem);
  for (arma::uword t = 0; t < S.n_elem; ++t) {
    if (S(t) < 1 || S(t) > R || S(t) == i || S(t) == j)
      stop("the conditioning set holds an invalid node index.");
    cond[t] = S(t) - 1;
  }
  const CiOut r = ci_pcor_block(Sigma, i - 1, j - 1, cond, n, df_residual,
                                df_residual > 0.0);
  return Rcpp::List::create(_["pcor"] = r.pcor, _["p_value"] = r.p_value,
                            _["df"] = r.df, _["ridge"] = r.ridge,
                            _["regularized"] = r.regularized,
                            _["usable"] = r.usable);
}

//' All separating sets of given pairs of nodes (C++ kernel)
//'
//' For each requested pair (i, j), enumerates every subset of the neighbors of
//' i other than j, and every subset of the neighbors of j other than i, of size
//' at most max_ord, and returns those that make the second-order partial
//' correlation non-significant at level alpha. Sets found in both neighborhoods
//' are returned once.
//'
//' @param Sigma R x R covariance matrix of the responses.
//' @param n number of experimental units.
//' @param adj integer R x R symmetric 0/1 adjacency matrix of the skeleton.
//' @param pairs integer q x 2 matrix of pairs (i, j), 1-based, assumed
//'   non-adjacent in adj.
//' @param alpha significance level of the conditional-association test.
//' @param max_ord largest conditioning set size to enumerate.
//' @param df_residual residual degrees of freedom of the first stage; a value
//'   not greater than zero selects the raw contract n - |S| - 3.
//' @return list with set_pair (1-based row of pairs), set_members (n_sets x R,
//'   0/1), set_size, set_pcor, set_p_value, and the diagnostics n_tests,
//'   n_regularized, max_ridge, n_unusable and truncated (one flag per pair).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_sepsets_pairs(const arma::mat& Sigma, double n,
                             const arma::umat& adj, const arma::umat& pairs,
                             double alpha, int max_ord,
                             double df_residual = 0.0) {
  const arma::uword R = Sigma.n_rows;
  if (Sigma.n_cols != R) stop("Sigma must be square.");
  if (adj.n_rows != R || adj.n_cols != R)
    stop("adj must be square with one row per response.");
  if (R > 64) stop("cpp_sepsets_pairs supports at most 64 responses.");
  if (pairs.n_cols != 2 && pairs.n_rows > 0)
    stop("pairs must have exactly two columns.");
  const bool use_dfr = (df_residual > 0.0);
  const arma::uword q = pairs.n_rows;

  std::vector<arma::uword> out_pair;
  std::vector<double> out_size, out_pcor, out_p;
  std::vector<std::vector<arma::uword> > out_sets;
  arma::uvec truncated(std::max<arma::uword>(q, 1), arma::fill::zeros);

  long long n_tests = 0, n_regularized = 0, n_unusable = 0;
  double max_ridge = 0.0;

  for (arma::uword e = 0; e < q; ++e) {
    const arma::uword i = pairs(e, 0) - 1;
    const arma::uword j = pairs(e, 1) - 1;
    if (i >= R || j >= R || i == j) stop("pairs holds an invalid node index.");

    std::vector<arma::uword> nbi, nbj;
    for (arma::uword v = 0; v < R; ++v) {
      if (v != i && v != j && adj(i, v) != 0) nbi.push_back(v);
      if (v != i && v != j && adj(j, v) != 0) nbj.push_back(v);
    }

    std::set<unsigned long long> seen;
    bool trunc = false;

    auto visit = [&](const std::vector<arma::uword>& S) {
      unsigned long long mask = 0ULL;
      for (std::size_t t = 0; t < S.size(); ++t)
        mask |= (1ULL << S[t]);
      if (!seen.insert(mask).second) return;          // already tested
      CiOut r = ci_pcor_block(Sigma, i, j, S, n, df_residual, use_dfr);
      ++n_tests;
      if (!r.usable) { ++n_unusable; trunc = true; return; }
      if (r.regularized) ++n_regularized;
      if (r.ridge > max_ridge) max_ridge = r.ridge;
      if (r.p_value >= alpha) {
        out_pair.push_back(e + 1);
        out_sets.push_back(S);
        out_size.push_back(static_cast<double>(S.size()));
        out_pcor.push_back(r.pcor);
        out_p.push_back(r.p_value);
      }
    };

    for_each_subset(nbi, max_ord, visit);
    for_each_subset(nbj, max_ord, visit);
    if (static_cast<int>(nbi.size()) > max_ord ||
        static_cast<int>(nbj.size()) > max_ord) trunc = true;
    if (q > 0) truncated(e) = trunc ? 1 : 0;
  }

  const arma::uword n_sets = out_sets.size();
  arma::umat members(n_sets, R, arma::fill::zeros);
  for (arma::uword s = 0; s < n_sets; ++s)
    for (std::size_t t = 0; t < out_sets[s].size(); ++t)
      members(s, out_sets[s][t]) = 1;

  return Rcpp::List::create(
    _["set_pair"]      = arma::uvec(out_pair),
    _["set_members"]   = members,
    _["set_size"]      = arma::vec(out_size),
    _["set_pcor"]      = arma::vec(out_pcor),
    _["set_p_value"]   = arma::vec(out_p),
    _["n_tests"]       = static_cast<double>(n_tests),
    _["n_regularized"] = static_cast<double>(n_regularized),
    _["n_unusable"]    = static_cast<double>(n_unusable),
    _["max_ridge"]     = max_ridge,
    _["truncated"]     = (q > 0 ? truncated : arma::uvec())
  );
}

// ---------------------------------------------------------------------------
// Meek's orientation rules.
//
// R1 to R3 complete a graph whose only oriented edges come from v-structures:
// applied to exhaustion they return the CPDAG of the equivalence class. R4 is
// idle in that case and becomes necessary as soon as some edges are oriented
// from outside, by background knowledge, because then the graph need not be the
// CPDAG of any class and the first three rules are no longer complete.
//
//   R1  c -> a,  a - b,  c and b nonadjacent            =>  a -> b
//   R2  a -> c -> b,  a - b                             =>  a -> b
//   R3  a - c,  a - d,  c -> b,  d -> b,  c,d nonadj.   =>  a -> b
//   R4  a - c,  a - d,  c -> d,  d -> b,  c,b nonadj.   =>  a -> b
// In every rule a - b is the undirected edge being oriented.
//
// R1, R3 and R4 all exist to avoid creating a collider at the tail node, so
// each is withheld when the unshielded triple carrying that collider is one a
// conservative or majority rule left ambiguous. R2 only forbids a cycle and is
// applied unconditionally.
//
// Two rules can require opposite directions on the same edge when the input is
// not consistent. The requests of a sweep are therefore collected first,
// contested edges are left undirected, and the clash is reported, rather than
// letting the later write overturn the earlier one.

namespace {

inline bool amb_has(const std::vector<unsigned long long>& keys,
                    arma::uword R, arma::uword a, arma::uword b,
                    arma::uword c) {
  if (keys.empty()) return false;
  const unsigned long long k =
    (unsigned long long)a * R * R + (unsigned long long)b * R + c;
  return std::binary_search(keys.begin(), keys.end(), k);
}

}  // namespace

//' Apply Meek's four orientation rules (C++ kernel)
//'
//' @param amat integer R x R adjacency in the package convention: amat(a,b) and
//'   amat(b,a) both one is an undirected edge, only amat(a,b) one is a -> b.
//'   Each of the four rules is applied in sweeps on the graph as it stood at the
//'   start of the sweep, and an edge two rules orient in opposite directions is
//'   left undirected and flagged in conflicts.
//' @param amb_keys keys a*R*R + b*R + c of the ambiguous triples (a, b, c with
//'   center b), 0-based node indices, both endpoint orders present. Rules R1, R3
//'   and R4 are withheld where they would settle such a triple; R2 is not.
//' @param max_iter cap on the number of sweeps.
//' @return list with amat, conflicts (R x R), iterations and blocked.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_meek_rules(arma::umat amat, const arma::uvec& amb_keys,
                          int max_iter = 1000) {
  const arma::uword R = amat.n_rows;
  if (amat.n_cols != R) stop("amat must be square.");
  std::vector<unsigned long long> keys(amb_keys.n_elem);
  for (arma::uword t = 0; t < amb_keys.n_elem; ++t)
    keys[t] = (unsigned long long)amb_keys(t);
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

  arma::umat conflicts(R, R, arma::fill::zeros);
  long long blocked = 0;
  int iterations = 0;
  bool changed = true;

  while (changed && iterations < max_iter) {
    ++iterations;
    changed = false;
    const arma::umat frozen = amat;

    auto has_edge = [&](arma::uword a, arma::uword b) {
      return frozen(a, b) == 1 || frozen(b, a) == 1;
    };
    auto is_dir = [&](arma::uword a, arma::uword b) {
      return frozen(a, b) == 1 && frozen(b, a) == 0;
    };
    auto is_undir = [&](arma::uword a, arma::uword b) {
      return frozen(a, b) == 1 && frozen(b, a) == 1;
    };

    arma::umat req(R, R, arma::fill::zeros);
    for (arma::uword a = 0; a + 1 < R; ++a) {
      for (arma::uword b = a + 1; b < R; ++b) {
        if (!is_undir(a, b)) continue;
        for (int side = 0; side < 2; ++side) {
          const arma::uword from = (side == 0) ? a : b;
          const arma::uword to   = (side == 0) ? b : a;

          bool fire = false;

          // R1
          for (arma::uword c = 0; c < R && !fire; ++c) {
            if (c == from || c == to) continue;
            if (!(is_dir(c, from) && !has_edge(c, to))) continue;
            if (amb_has(keys, R, c, from, to)) { ++blocked; continue; }
            fire = true;
          }
          // R2
          for (arma::uword c = 0; c < R && !fire; ++c) {
            if (c == from || c == to) continue;
            if (is_dir(from, c) && is_dir(c, to)) fire = true;
          }
          // R3
          for (arma::uword c = 0; c < R && !fire; ++c) {
            if (c == from || c == to) continue;
            if (!(is_undir(from, c) && is_dir(c, to))) continue;
            for (arma::uword d = c + 1; d < R && !fire; ++d) {
              if (d == from || d == to) continue;
              if (!(is_undir(from, d) && is_dir(d, to))) continue;
              if (has_edge(c, d)) continue;
              if (amb_has(keys, R, c, from, d)) { ++blocked; continue; }
              fire = true;
            }
          }
          // R4
          for (arma::uword c = 0; c < R && !fire; ++c) {
            if (c == from || c == to) continue;
            if (!is_undir(from, c) || has_edge(c, to)) continue;
            for (arma::uword d = 0; d < R && !fire; ++d) {
              if (d == from || d == to || d == c) continue;
              if (!is_undir(from, d)) continue;
              if (!(is_dir(c, d) && is_dir(d, to))) continue;
              if (amb_has(keys, R, c, from, to)) { ++blocked; continue; }
              fire = true;
            }
          }

          if (fire) req(from, to) = 1;
        }
      }
    }

    arma::umat clash = req % req.t();
    conflicts = arma::max(conflicts, clash);
    for (arma::uword f = 0; f < R; ++f) {
      for (arma::uword t = 0; t < R; ++t) {
        if (req(f, t) == 1 && clash(f, t) == 0 &&
            amat(f, t) == 1 && amat(t, f) == 1) {
          amat(t, f) = 0;
          changed = true;
        }
      }
    }
  }

  return Rcpp::List::create(_["amat"] = amat, _["conflicts"] = conflicts,
                            _["iterations"] = iterations,
                            _["blocked"] = (double)blocked);
}

// ---------------------------------------------------------------------------
// Structure comparison.
//
// A pair of responses is in one of four states: absent, undirected, or directed
// either way. Coding them once and comparing the codes gives the structural
// Hamming distance of the partially directed graph, its skeleton counts, and
// the v-structures, without a second pass over the adjacency.
// ---------------------------------------------------------------------------

//' State of every pair of nodes (C++ kernel)
//'
//' @param amat integer adjacency in the package convention.
//' @return integer vector over the pairs in upper-triangular order: 1 absent,
//'   2 undirected, 3 i to j, 4 j to i.
//' @keywords internal
// [[Rcpp::export]]
arma::uvec cpp_pair_states(const arma::umat& amat) {
  const arma::uword R = amat.n_rows;
  if (amat.n_cols != R) stop("amat must be square.");
  arma::uvec out(R * (R - 1) / 2);
  arma::uword e = 0;
  for (arma::uword i = 0; i + 1 < R; ++i)
    for (arma::uword j = i + 1; j < R; ++j) {
      const bool a = amat(i, j) == 1, b = amat(j, i) == 1;
      out(e++) = (!a && !b) ? 1 : (a && b) ? 2 : a ? 3 : 4;
    }
  return out;
}

//' V-structures of a partially directed graph (C++ kernel)
//'
//' Triples i -> k <- j with i < j and i, j nonadjacent, each arc directed.
//'
//' @param amat integer adjacency in the package convention.
//' @return integer matrix with columns (i, k, j), 1-based.
//' @keywords internal
// [[Rcpp::export]]
arma::umat cpp_vstructure_set(const arma::umat& amat) {
  const arma::uword R = amat.n_rows;
  std::vector<arma::uword> vi, vk, vj;
  for (arma::uword k = 0; k < R; ++k) {
    std::vector<arma::uword> pa;
    for (arma::uword p = 0; p < R; ++p)
      if (p != k && amat(p, k) == 1 && amat(k, p) == 0) pa.push_back(p);
    for (std::size_t a = 0; a + 1 < pa.size(); ++a)
      for (std::size_t b = a + 1; b < pa.size(); ++b) {
        const arma::uword i = pa[a], j = pa[b];
        if (amat(i, j) == 0 && amat(j, i) == 0) {
          vi.push_back(std::min(i, j) + 1);
          vk.push_back(k + 1);
          vj.push_back(std::max(i, j) + 1);
        }
      }
  }
  arma::umat out(vi.size(), 3);
  for (arma::uword t = 0; t < vi.size(); ++t) {
    out(t, 0) = vi[t]; out(t, 1) = vk[t]; out(t, 2) = vj[t];
  }
  return out;
}

//' Confusion counts of two edge indicators (C++ kernel)
//'
//' @param est logical vector of estimated edges.
//' @param tru logical vector of true edges, same length.
//' @return integer vector (tp, fp, fn, tn) over the positions where both are
//'   not missing.
//' @keywords internal
// [[Rcpp::export]]
arma::uvec cpp_confusion(const arma::uvec& est, const arma::uvec& tru) {
  if (est.n_elem != tru.n_elem)
    stop("the two indicators must have the same length.");
  arma::uvec out(4, arma::fill::zeros);
  for (arma::uword k = 0; k < est.n_elem; ++k) {
    const bool a = est(k) == 1, b = tru(k) == 1;
    if (a && b) ++out(0); else if (a && !b) ++out(1);
    else if (!a && b) ++out(2); else ++out(3);
  }
  return out;
}

// ---------------------------------------------------------------------------
// The two validity gates of a partially directed output.
//
// Both are pure combinatorics on the adjacency, with no floating point in
// them, and both are evaluated once per resample by the structure-frequency
// routines, which is what makes them worth moving out of R.
// ---------------------------------------------------------------------------

//' Does the directed part contain a cycle? (C++ kernel)
//'
//' Kahn's algorithm on the arcs alone; undirected edges are ignored.
//'
//' @param amat integer adjacency in the package convention.
//' @return TRUE when some arc lies on a directed cycle.
//' @keywords internal
// [[Rcpp::export]]
bool cpp_has_directed_cycle(const arma::umat& amat) {
  const arma::uword R = amat.n_rows;
  std::vector<char> removed(R, 0);
  arma::uword left = R;
  while (left > 0) {
    bool moved = false;
    for (arma::uword v = 0; v < R; ++v) {
      if (removed[v]) continue;
      bool has_parent = false;
      for (arma::uword p = 0; p < R; ++p)
        if (!removed[p] && p != v && amat(p, v) == 1 && amat(v, p) == 0) {
          has_parent = true; break;
        }
      if (has_parent) continue;
      removed[v] = 1; --left; moved = true;
    }
    if (!moved) return true;          // every survivor has a parent
  }
  return false;
}

//' Does a consistent DAG extension exist? (C++ kernel)
//'
//' The criterion of Dor and Tarsi (1992): repeatedly remove a node with no
//' outgoing arc whose undirected neighbors are adjacent to all its other
//' neighbors. A graph that empties out admits an extension preserving the arcs
//' and the v-structures.
//'
//' @param amat integer adjacency in the package convention.
//' @return TRUE when an extension exists.
//' @keywords internal
// [[Rcpp::export]]
bool cpp_has_consistent_extension(const arma::umat& amat) {
  const arma::uword R = amat.n_rows;
  std::vector<char> alive(R, 1);
  for (arma::uword step = 0; step < R; ++step) {
    bool found = false;
    for (arma::uword x = 0; x < R && !found; ++x) {
      if (!alive[x]) continue;
      bool sink = true;
      for (arma::uword w = 0; w < R; ++w)
        if (alive[w] && w != x && amat(x, w) == 1 && amat(w, x) == 0) {
          sink = false; break;
        }
      if (!sink) continue;
      std::vector<arma::uword> und, nb;
      for (arma::uword w = 0; w < R; ++w) {
        if (!alive[w] || w == x) continue;
        const bool edge = amat(x, w) == 1 || amat(w, x) == 1;
        if (edge) nb.push_back(w);
        if (amat(x, w) == 1 && amat(w, x) == 1) und.push_back(w);
      }
      bool ok = true;
      for (std::size_t a = 0; a < und.size() && ok; ++a)
        for (std::size_t b = 0; b < nb.size(); ++b) {
          if (nb[b] == und[a]) continue;
          if (amat(und[a], nb[b]) == 0 && amat(nb[b], und[a]) == 0) {
            ok = false; break;
          }
        }
      if (!ok) continue;
      alive[x] = 0; found = true;
    }
    if (!found) break;
  }
  for (arma::uword v = 0; v < R; ++v) if (alive[v]) return false;
  return true;
}
