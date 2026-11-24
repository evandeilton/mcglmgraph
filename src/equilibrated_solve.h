// equilibrated_solve.h: scale-free solve of a symmetric parent block.
// Part of the mcglmgraph package. License: GPL-3.
//
// The normal equations A g = b of a nodal regression are solved on the
// equilibrated block As = D^-1 A D^-1, D = diag(sqrt(diag(A))), which has a
// unit diagonal. Its reciprocal condition number does not change when a
// response is rescaled, so the near-singularity test and the fallback ridge are
// relative to the scale of the data.

#ifndef MCGLMGRAPH_EQUILIBRATED_SOLVE_H
#define MCGLMGRAPH_EQUILIBRATED_SOLVE_H

#include <RcppArmadillo.h>
#include <cmath>

namespace mcglmgraph {

// Returns true when a solution g was found. `ill` is set when the block is
// near singular (rcond(As) <= rcond_tol, or a diagonal entry that is not
// finite and positive). With allow_ridge a block that cannot be solved exactly
// is solved with ridge_rel added to the unit diagonal of As, and `ridged` is
// set; without it no solution is returned.
inline bool equilibrated_solve(const arma::mat& A, const arma::vec& b,
                               double rcond_tol, double ridge_rel,
                               bool allow_ridge, arma::vec& g, bool& ill,
                               bool& ridged) {
  ill = false;
  ridged = false;
  if (!A.is_finite() || !b.is_finite()) { ill = true; return false; }
  const arma::uword m = A.n_rows;
  arma::vec d = A.diag();
  for (arma::uword k = 0; k < m; ++k) {
    if (std::isfinite(d(k)) && d(k) > 0.0) {
      d(k) = std::sqrt(d(k));
    } else {
      // no scale to equilibrate by: flag the block and leave that row as is
      ill = true;
      d(k) = 1.0;
    }
  }
  arma::mat As = A.each_col() / d;
  As.each_row() /= d.t();
  As = 0.5 * (As + As.t());
  const arma::vec bs = b / d;

  arma::vec gs;
  bool ok = false;
  if (!ill) {
    const double rc = arma::rcond(As);
    ill = !(std::isfinite(rc) && rc > rcond_tol);
  }
  if (!ill) ok = arma::solve(gs, As, bs, arma::solve_opts::no_approx);
  if (!ok || !gs.is_finite()) {
    ok = false;
    if (!allow_ridge) return false;
    ridged = true;
    ok = arma::solve(gs, As + ridge_rel * arma::eye(m, m), bs,
                     arma::solve_opts::no_approx);
  }
  if (!ok || !gs.is_finite()) return false;
  g = gs / d;
  return true;
}

}  // namespace mcglmgraph

#endif
