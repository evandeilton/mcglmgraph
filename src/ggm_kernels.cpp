// ggm_kernels.cpp: Gaussian Graphical Model layer kernels (RcppArmadillo)
// Part of the mcglmgraph package. License: GPL-3.
// Author: Jose Evandeilton Lopes. MCGLM engine concepts: Wagner Hugo Bonat.
//
// Implements the precision-based parametrization of the between-response
// covariance: K_b = sum_d eta_d Z_d, Sigma_b = K_b^-1, with analytical
// derivatives dSigma_b/deta_d = -Sigma_b Z_d Sigma_b.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Build the between-response precision matrix (C++ kernel)
//'
//' Computes K_b = sum_d eta_d Z_d and symmetrizes the result.
//' @param eta numeric vector of length D (diagonal then edge parameters).
//' @param bases numeric cube R x R x D of symmetric base matrices Z_d.
//' @return The R x R symmetric precision matrix K_b = sum_d eta_d Z_d.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_build_precision(const arma::vec& eta, const arma::cube& bases) {
  const arma::uword R = bases.n_rows;
  const arma::uword D = bases.n_slices;
  if (eta.n_elem != D) stop("length(eta) must equal the number of bases.");
  arma::mat K(R, R, arma::fill::zeros);
  for (arma::uword d = 0; d < D; ++d) K += eta(d) * bases.slice(d);
  return 0.5 * (K + K.t());
}

//' Invert a symmetric positive-definite precision matrix (C++ kernel)
//'
//' Inverts the symmetric part of K by Cholesky-based inversion and stops with
//' an error when K is not positive definite.
//' @param K symmetric positive-definite matrix.
//' @return The covariance matrix Sigma_b = K^-1.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_precision_to_cov(const arma::mat& K) {
  arma::mat S;
  bool ok = arma::inv_sympd(S, 0.5 * (K + K.t()));
  if (!ok) stop("Precision matrix is not positive definite.");
  return S;
}

//' Analytical derivatives of Sigma_b with respect to eta (C++ kernel)
//'
//' Since Sigma_b = K_b^-1 and dK_b/deta_d = Z_d, the derivative is
//' dSigma_b/deta_d = -Sigma_b Z_d Sigma_b; each matrix is symmetrized.
//' @param Sigma_b the R x R covariance matrix.
//' @param bases numeric cube R x R x D of base matrices Z_d.
//' @return list of D symmetric matrices, each -Sigma_b Z_d Sigma_b.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_precision_derivatives(const arma::mat& Sigma_b,
                                     const arma::cube& bases) {
  const arma::uword D = bases.n_slices;
  Rcpp::List out(D);
  for (arma::uword d = 0; d < D; ++d) {
    arma::mat Dd = -Sigma_b * bases.slice(d) * Sigma_b;
    out[d] = 0.5 * (Dd + Dd.t());
  }
  return out;
}

//' Test symmetric positive-definiteness via the smallest eigenvalue (C++ kernel)
//'
//' The symmetric part of K is declared positive definite when its smallest
//' eigenvalue exceeds tol.
//' @param K square matrix; only its symmetric part is used.
//' @param tol numeric tolerance on the smallest eigenvalue (default 1e-10).
//' @return list with is_spd (logical) and min_eig (the smallest eigenvalue).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_check_spd(const arma::mat& K, double tol = 1e-10) {
  arma::vec ev = arma::eig_sym(0.5 * (K + K.t()));
  double mn = ev.min();
  return Rcpp::List::create(_["is_spd"] = (mn > tol), _["min_eig"] = mn);
}

//' Backtracking line search preserving positive-definiteness (C++ kernel)
//'
//' Largest step alpha among 1, 1/2, ..., 2^-max_halvings such that
//' K_b(eta + alpha*d) has smallest eigenvalue above tol. When no admissible
//' step is found, is_spd is false, eta_new is the point tried with the smallest
//' step 2^-max_halvings (it fails the test), halvings is max_halvings + 1 and
//' step_size is 2^-(max_halvings + 1), half the step that produced eta_new.
//' @param eta current parameter vector.
//' @param d Newton direction.
//' @param bases numeric cube R x R x D of base matrices Z_d.
//' @param max_halvings maximum number of step halvings (default 20).
//' @param tol tolerance on the smallest eigenvalue (default 1e-10).
//' @return list with eta_new (last point tried), step_size, halvings (number
//'   of halvings performed) and is_spd (whether eta_new passed the test).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_spd_line_search(const arma::vec& eta, const arma::vec& d,
                               const arma::cube& bases, int max_halvings = 20,
                               double tol = 1e-10) {
  double alpha = 1.0;
  int h = 0;
  bool ok = false;
  arma::vec eta_new = eta;
  for (h = 0; h <= max_halvings; ++h) {
    eta_new = eta + alpha * d;
    arma::mat K = cpp_build_precision(eta_new, bases);
    arma::vec ev = arma::eig_sym(K);
    if (ev.min() > tol) { ok = true; break; }
    alpha *= 0.5;
  }
  return Rcpp::List::create(_["eta_new"] = eta_new, _["step_size"] = alpha,
                            _["halvings"] = h, _["is_spd"] = ok);
}

//' Partial correlation matrix from a precision matrix (C++ kernel)
//'
//' Entry (i, j) is -K_ij / sqrt(K_ii K_jj), with unit diagonal. The diagonal
//' of K is assumed positive; this is not checked.
//' @param K symmetric positive-definite precision matrix.
//' @return R x R matrix of partial correlations.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_partial_cor(const arma::mat& K) {
  const arma::uword R = K.n_rows;
  arma::vec dd = arma::sqrt(K.diag());
  arma::mat P(R, R, arma::fill::zeros);
  for (arma::uword i = 0; i < R; ++i) {
    for (arma::uword j = 0; j < R; ++j) {
      P(i, j) = (i == j) ? 1.0 : -K(i, j) / (dd(i) * dd(j));
    }
  }
  return P;
}

//' Delta-method variance of the between-response partial correlations
//'
//' For the pair (i, j) the partial correlation is
//' rho_ij = -K_ij / sqrt(K_ii K_jj). Differentiating it with respect to a
//' covariance parameter theta_d gives
//'   g_d = -Kdot_d(i,j)/sqrt(K_ii K_jj)
//'         - (rho_ij/2) [ Kdot_d(i,i)/K_ii + Kdot_d(j,j)/K_jj ],
//' where Kdot_d = dK/dtheta_d. The delta method gives the variance
//' Var(rho_ij) = g' V g, with V the estimated covariance of the parameter
//' vector.
//'
//' The derivatives of the precision are obtained from the derivatives of the
//' covariance by Kdot_d = -K Sdot_d K, an identity that holds for any
//' parametrization. In the GGM layer it returns exactly the canonical base
//' Z_d; in the DAG layer it supplies the corresponding directions for
//' (gamma, delta) without changing the formula.
//'
//' K is assumed positive definite and is not checked. For a pair with
//' K_ii K_jj <= 0 the kernel returns rho = 0, a non-finite gradient and a
//' variance of zero. A quadratic form g' V g that is negative or not finite is
//' returned as zero.
//'
//' @param K the R x R precision matrix of the fitted layer.
//' @param Sdot list of D matrices dSigma_b/dtheta_d.
//' @param pairs integer matrix q x 2 with the response pairs (i, j), 1-based.
//' @param V estimated D x D covariance matrix of the parameter vector.
//' @return list with pcor (length q), variance (length q) and gradient (the
//'   q x D matrix of gradients g).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_pcor_delta(const arma::mat& K, const Rcpp::List& Sdot,
                          const arma::umat& pairs, const arma::mat& V) {
  const arma::uword D = Sdot.size();
  const arma::uword q = pairs.n_rows;
  if (V.n_rows != D || V.n_cols != D)
    stop("V must be square with one row per parameter.");

  // Kdot_d = -K Sdot_d K
  std::vector<arma::mat> Kdot(D);
  for (arma::uword d = 0; d < D; ++d) {
    arma::mat Sd = as<arma::mat>(Sdot[d]);
    arma::mat Kd = -K * Sd * K;
    Kdot[d] = 0.5 * (Kd + Kd.t());
  }

  arma::vec rho(q), var(q);
  arma::mat G(q, D, arma::fill::zeros);
  for (arma::uword e = 0; e < q; ++e) {
    arma::uword i = pairs(e, 0) - 1, j = pairs(e, 1) - 1;
    double kii = K(i, i), kjj = K(j, j);
    double den = std::sqrt(kii * kjj);
    double r = (den > 0.0) ? -K(i, j) / den : 0.0;
    rho(e) = r;
    for (arma::uword d = 0; d < D; ++d) {
      double g = -Kdot[d](i, j) / den
        - 0.5 * r * (Kdot[d](i, i) / kii + Kdot[d](j, j) / kjj);
      G(e, d) = g;
    }
    arma::rowvec ge = G.row(e);
    double v = arma::as_scalar(ge * V * ge.t());
    var(e) = (v > 0.0) ? v : 0.0;
  }
  return Rcpp::List::create(_["pcor"] = rho, _["variance"] = var,
                            _["gradient"] = G);
}
