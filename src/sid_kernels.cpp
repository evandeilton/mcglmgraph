// sid_kernels.cpp: structural intervention distance (RcppArmadillo).
// Part of the mcglmgraph package. License: GPL-3.
//
// The structural Hamming distance counts the edges two graphs disagree about.
// It says nothing about what the disagreement costs, and two graphs at the same
// Hamming distance from the truth can be worth very different things: one may
// still give the right adjustment set for every pair of nodes, the other none.
//
// The structural intervention distance of Peters and Buehlmann (2015) counts
// instead the ordered pairs (i, j) for which the parent set of i in the
// estimated graph is not a valid adjustment set for the effect of i on j in the
// true graph.
//
// Validity is decided by the generalized adjustment criterion rather than by
// the back-door one. The back-door criterion is sufficient and not necessary:
// it refuses any set containing a descendant of i, while a descendant that lies
// on no causal path from i to j does no harm. Using it would count pairs the
// adjustment in fact recovers. Writing cn(i, j) for the nodes on proper causal
// paths from i to j, and forb(i, j) for their descendants together with i, the
// set S is valid when
//   (a) S does not meet forb(i, j), and
//   (b) S d-separates i from j in the graph with the first edges of the proper
//       causal paths, that is the arcs from i into cn(i, j), removed.
// With no causal path from i to j the second condition asks for d-separation in
// the whole graph, which is what recovers a null effect.
//
// One case needs care. When j itself belongs to the adjustment set the formula
// collapses to the marginal of j, since conditioning on j and averaging over it
// returns it unchanged. That is right exactly when the effect is null, and
// wrong otherwise.
//
// D-separation is decided by the ancestral moralization: restrict to the
// ancestors of the three sets, marry the parents of every node, delete the
// conditioning set, and look for a path.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <queue>
using namespace Rcpp;

namespace {

// Nodes reachable from the seeds by following arcs backwards (ancestors,
// inclusive of the seeds themselves).
std::vector<char> ancestors_of(const arma::umat& G,
                               const std::vector<arma::uword>& seeds) {
  const arma::uword R = G.n_rows;
  std::vector<char> in(R, 0);
  std::queue<arma::uword> q;
  for (std::size_t t = 0; t < seeds.size(); ++t)
    if (!in[seeds[t]]) { in[seeds[t]] = 1; q.push(seeds[t]); }
  while (!q.empty()) {
    const arma::uword v = q.front(); q.pop();
    for (arma::uword p = 0; p < R; ++p)
      if (!in[p] && G(p, v) == 1) { in[p] = 1; q.push(p); }
  }
  return in;
}

// Nodes reachable from i by following arcs forwards, i included.
std::vector<char> descendants_of(const arma::umat& G, arma::uword i) {
  const arma::uword R = G.n_rows;
  std::vector<char> in(R, 0);
  std::queue<arma::uword> q;
  in[i] = 1; q.push(i);
  while (!q.empty()) {
    const arma::uword v = q.front(); q.pop();
    for (arma::uword c = 0; c < R; ++c)
      if (!in[c] && G(v, c) == 1) { in[c] = 1; q.push(c); }
  }
  return in;
}

// Are i and j d-separated by S in the DAG G?
bool d_separated(const arma::umat& G, arma::uword i, arma::uword j,
                 const std::vector<arma::uword>& S) {
  const arma::uword R = G.n_rows;
  std::vector<arma::uword> seeds = S;
  seeds.push_back(i);
  seeds.push_back(j);
  const std::vector<char> anc = ancestors_of(G, seeds);

  // moral graph of the ancestral subgraph
  std::vector<std::vector<char> > M(R, std::vector<char>(R, 0));
  for (arma::uword v = 0; v < R; ++v) {
    if (!anc[v]) continue;
    std::vector<arma::uword> pa;
    for (arma::uword p = 0; p < R; ++p)
      if (anc[p] && G(p, v) == 1) pa.push_back(p);
    for (std::size_t a = 0; a < pa.size(); ++a) {
      M[pa[a]][v] = M[v][pa[a]] = 1;
      for (std::size_t b = a + 1; b < pa.size(); ++b)
        M[pa[a]][pa[b]] = M[pa[b]][pa[a]] = 1;
    }
  }
  // delete the conditioning set
  std::vector<char> gone(R, 0);
  for (std::size_t t = 0; t < S.size(); ++t) gone[S[t]] = 1;
  if (gone[i] || gone[j]) return true;      // conditioned on an endpoint

  std::vector<char> seen(R, 0);
  std::queue<arma::uword> q;
  seen[i] = 1; q.push(i);
  while (!q.empty()) {
    const arma::uword v = q.front(); q.pop();
    if (v == j) return false;
    for (arma::uword w = 0; w < R; ++w)
      if (!seen[w] && !gone[w] && anc[w] && M[v][w] == 1) {
        seen[w] = 1; q.push(w);
      }
  }
  return true;
}

}  // namespace

//' Structural intervention distance between two DAGs (C++ kernel)
//'
//' Counts the ordered pairs (i, j), i not j, for which the parents of i in the
//' estimated graph are not a valid adjustment set for the effect of i on j in
//' the true graph.
//'
//' @param Gtrue integer R x R adjacency of the true DAG, one in (a, b) meaning
//'   a -> b.
//' @param Hest integer R x R adjacency of the estimated DAG, same convention.
//' @return list with the count and the R x R matrix of offending pairs.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_sid(const arma::umat& Gtrue, const arma::umat& Hest) {
  const arma::uword R = Gtrue.n_rows;
  if (Gtrue.n_cols != R || Hest.n_rows != R || Hest.n_cols != R)
    stop("both graphs must be square and of the same size.");
  arma::umat bad(R, R, arma::fill::zeros);
  long long count = 0;

  for (arma::uword i = 0; i < R; ++i) {
    // parents of i in the estimate: the adjustment set it prescribes
    std::vector<arma::uword> pa;
    for (arma::uword p = 0; p < R; ++p)
      if (p != i && Hest(p, i) == 1) pa.push_back(p);

    const std::vector<char> desc_i = descendants_of(Gtrue, i);

    for (arma::uword j = 0; j < R; ++j) {
      if (j == i) continue;

      // nodes on proper causal paths from i to j: descendants of i other than
      // i that are ancestors of j or j itself
      std::vector<arma::uword> seed_j(1, j);
      const std::vector<char> anc_j = ancestors_of(Gtrue, seed_j);
      std::vector<arma::uword> cn;
      for (arma::uword w = 0; w < R; ++w)
        if (w != i && desc_i[w] && anc_j[w]) cn.push_back(w);

      // the forbidden set: descendants of those nodes, plus i
      std::vector<char> forb(R, 0);
      forb[i] = 1;
      if (!cn.empty()) {
        std::vector<char> seen(R, 0);
        std::queue<arma::uword> q;
        for (std::size_t t = 0; t < cn.size(); ++t)
          if (!seen[cn[t]]) { seen[cn[t]] = 1; q.push(cn[t]); }
        while (!q.empty()) {
          const arma::uword v = q.front(); q.pop();
          forb[v] = 1;
          for (arma::uword c = 0; c < R; ++c)
            if (!seen[c] && Gtrue(v, c) == 1) { seen[c] = 1; q.push(c); }
        }
      }

      bool offending = false;
      for (std::size_t t = 0; t < pa.size(); ++t)
        if (forb[pa[t]]) { offending = true; break; }

      if (!offending && Hest(j, i) == 1) {
        // j is in the adjustment set: the formula returns the marginal of j,
        // right when the effect is null and wrong when it is not
        offending = (desc_i[j] == 1);
      } else if (!offending) {
        // remove the first edges of the proper causal paths and ask for
        // d-separation there
        arma::umat Gcut = Gtrue;
        for (std::size_t t = 0; t < cn.size(); ++t) Gcut(i, cn[t]) = 0;
        offending = !d_separated(Gcut, i, j, pa);
      }
      if (offending) { bad(i, j) = 1; ++count; }
    }
  }
  return Rcpp::List::create(_["distance"] = (double)count, _["pairs"] = bad);
}
