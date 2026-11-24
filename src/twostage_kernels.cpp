// twostage_kernels.cpp: propagation of the first-stage uncertainty into the
// graph layer (RcppArmadillo).
// Part of the mcglmgraph package. License: GPL-3.
//
// Stacking the first-stage and the graph estimating functions gives a block
// triangular sensitivity, and block inversion gives the influence function of
// the second stage,
//   IF_i = H^-1 ( psi_i - B A^-1 psi_alpha_i ),
// with H the observed curvature of the graph block, A the first-stage
// sensitivity and B = -E( d psi_i / d alpha' ) the cross sensitivity.
//
// When every response carries a single dispersion parameter that multiplies
// its variance function, the within-response structure is the identity and
// the prior weights are constant, the first-stage Pearson equations fix the
// average m_r of the squared residuals of each response, and the first-stage
// term has, to first order and when the first stage acts on the residual
// covariance as a diagonal rescaling, the reduced form
//   B A^-1 psi_alpha_i = Q ( r_i o r_i - m ),
//   Q(d, r) = (1/2) [ K Ddot_d K S_e ](r, r) / m_r,
// with Ddot_d = dSigma_b / dtheta_d and S_e the sample covariance of the
// Pearson residuals. The factor 1/m_r is the influence of the dispersion: the
// residuals of response r scale as tau_r^-1/2 and the first stage fixes their
// average square at m_r. This kernel forms the shortcut: the reduced form
// evaluated at the fitted covariance, S_e replaced by Sigma_b,
//   Q(d, r) = (1/2) [ K Ddot_d ](r, r) / m_r.
// The two agree in every row of Q whose entries Sigma_b and S_e agree on: at
// every support of the GGM layer and for the parameters of every node of a DAG
// without immoralities. For the parameters of a collider with non-adjacent
// parents, and of its descendants, the shortcut is a first-order approximation
// in every response column of the parent block, parents included, and needs a
// correctly specified DAG; for the parameters of every other node it equals the
// reduced form.
//
// m_r is the average of the squared centered residuals Ec, which is the
// average the first stage fixes when the residuals of a response average zero.
// It differs from one by the finite-sample bias correction of the first
// stage, of order 1/n; with it the correction term sums to zero in the sample,
// and a quantity the first stage fixes exactly, as the innovation variance of
// a root node of a DAG with Gaussian intercept-only margins, gets a zero
// contribution.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Cross-sensitivity factor of the graph layer (C++ kernel)
//'
//' Forms the shortcut factor Q0(d, r) = (1/2) [ K Ddot_d ](r, r), which is the
//' reduced form factor (1/2) [ K Ddot_d K S_e ](r, r) evaluated at the fitted
//' covariance (S_e replaced by the fitted Sigma_b = K^-1), before the division
//' by the average squared residual m_r of each response.
//'
//' @param K the R x R precision matrix of the fitted layer.
//' @param Dlist list of P derivatives dSigma_b/dtheta_d, each R x R.
//' @return the P x R matrix Q0, whose entry (d, r) is the factor of parameter
//'   d and response r before division by m_r.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_two_stage_q(const arma::mat& K, const Rcpp::List& Dlist) {
  const arma::uword R = K.n_rows;
  const arma::uword P = Dlist.size();
  arma::mat Q(P, R, arma::fill::zeros);
  for (arma::uword d = 0; d < P; ++d) {
    arma::mat Dd = as<arma::mat>(Dlist[d]);
    if (Dd.n_rows != R || Dd.n_cols != R)
      stop("every derivative must be square with one row per response.");
    // only the diagonal of K * Dd is needed
    for (arma::uword r = 0; r < R; ++r)
      Q(d, r) = 0.5 * arma::dot(K.row(r), Dd.col(r));
  }
  return Q;
}

//' Corrected per-unit contributions of the graph layer (C++ kernel)
//'
//' Returns c_i = psi_i - Q (r_i o r_i - m), with m the observed average of the
//' squared residuals of each response and
//' Q(d, r) = (1/2) [ K Ddot_d ](r, r) / m_r.
//' Summed over units these vanish, and their outer product is the middle of the
//' corrected sandwich.
//'
//' @param Psi n x P matrix of Pearson contributions.
//' @param K the R x R precision matrix of the fitted layer.
//' @param Dlist list of P derivatives dSigma_b/dtheta_d.
//' @param Ec n x R matrix of centered Pearson residuals; m is the column mean
//'   of its squares, which must be positive.
//' @return list with \code{scores} (n x P), the factor \code{Q} (P x R) applied,
//'   already divided by m_r, and \code{centre}, the vector of the R averages m.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_two_stage_scores(const arma::mat& Psi, const arma::mat& K,
                                const Rcpp::List& Dlist, const arma::mat& Ec) {
  const arma::uword n = Ec.n_rows;
  const arma::uword R = Ec.n_cols;
  if (Psi.n_rows != n) stop("Psi must have one row per unit.");
  if (K.n_rows != R || K.n_cols != R)
    stop("K must be square with one row per response.");
  arma::mat Q = cpp_two_stage_q(K, Dlist);
  if (Q.n_rows != Psi.n_cols)
    stop("Dlist and Psi disagree on the number of parameters.");

  arma::mat Sq = Ec % Ec;                       // squared residuals
  arma::rowvec m = arma::mean(Sq, 0);
  if (!m.is_finite() || arma::any(m <= 0))
    stop("every response needs a positive average squared residual.");
  // influence of the dispersion, (r^2 - m_r) / m_r: column r is divided by m_r
  Q.each_row() /= m;
  Sq.each_row() -= m;
  arma::mat C = Psi - Sq * Q.t();
  return Rcpp::List::create(_["scores"] = C, _["Q"] = Q,
                            _["centre"] = m.t());
}
