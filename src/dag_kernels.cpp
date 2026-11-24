// dag_kernels.cpp: Directed Acyclic Graph layer kernels (RcppArmadillo)
// Part of the mcglmgraph package. License: GPL-3.
//
// Implements the modified Cholesky (LDL) parametrization of the
// between-response covariance: Sigma_b = T^-1 D (T^-1)^T, with the symmetric
// derivative in the structural coefficients:
//   dSigma_b/dgamma_ji = T^-1 E_ji Sigma_b + Sigma_b E_ij (T^-1)^T.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "equilibrated_solve.h"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

//' Build the unit lower-triangular structural matrix T (C++ kernel)
//'
//' Places -gamma_e at position (j, i) of the identity for each edge e = (i, j),
//' i < j, so that T is unit lower triangular with T[j, i] = -gamma_ji.
//'
//' @param gamma numeric vector of structural coefficients, one per edge, in the
//'   row order of edge_idx.
//' @param R integer number of responses.
//' @param edge_idx integer matrix |E| x 2 with edges (i, j), i < j, 1-based.
//' @return the R x R unit lower-triangular matrix T with T(j,i) = -gamma.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_build_T(const arma::vec& gamma, int R, const arma::umat& edge_idx) {
  arma::mat T(R, R, arma::fill::eye);
  const arma::uword E = edge_idx.n_rows;
  if (gamma.n_elem != E) stop("length(gamma) must equal nrow(edge_idx).");
  for (arma::uword e = 0; e < E; ++e) {
    arma::uword i = edge_idx(e, 0) - 1;
    arma::uword j = edge_idx(e, 1) - 1;
    T(j, i) = -gamma(e);
  }
  return T;
}

//' Covariance from the LDL parametrization (C++ kernel)
//'
//' Sigma_b = T^-1 D (T^-1)^T with D = diag(delta), symmetrized.
//'
//' @param gamma numeric vector of structural coefficients.
//' @param delta numeric vector of R positive innovation variances.
//' @param R integer number of responses.
//' @param edge_idx integer matrix |E| x 2 with edges (i, j), 1-based.
//' @return list with Sigma_b, T_inv and T.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_cholesky_to_cov(const arma::vec& gamma, const arma::vec& delta,
                               int R, const arma::umat& edge_idx) {
  arma::mat T = cpp_build_T(gamma, R, edge_idx);
  arma::mat Tinv = arma::inv(arma::trimatl(T));
  arma::mat D = arma::diagmat(delta);
  arma::mat Sigma = Tinv * D * Tinv.t();
  Sigma = 0.5 * (Sigma + Sigma.t());
  return Rcpp::List::create(_["Sigma_b"] = Sigma, _["T_inv"] = Tinv, _["T"] = T);
}

//' Precision from the LDL parametrization (C++ kernel)
//' @inheritParams cpp_cholesky_to_cov
//' @return the R x R precision matrix K_b = T^T D^-1 T.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_cholesky_to_precision(const arma::vec& gamma, const arma::vec& delta,
                                    int R, const arma::umat& edge_idx) {
  arma::mat T = cpp_build_T(gamma, R, edge_idx);
  arma::mat Dinv = arma::diagmat(1.0 / delta);
  arma::mat K = T.t() * Dinv * T;
  return 0.5 * (K + K.t());
}

//' Analytical LDL derivatives of Sigma_b (C++ kernel)
//'
//' The derivative in the innovation variance delta_k is
//' T^-1 E_kk (T^-1)^T, and the derivative in the structural coefficient
//' gamma_ji is T^-1 E_ji Sigma_b + Sigma_b E_ij (T^-1)^T, which is symmetric
//' because the second term is the transpose of the first. E_ab has a one in
//' position (a, b) and zeros elsewhere.
//' @inheritParams cpp_cholesky_to_cov
//' @return list with D_delta (R matrices) and D_gamma (|E| matrices).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_cholesky_derivatives(const arma::vec& gamma, const arma::vec& delta,
                                    int R, const arma::umat& edge_idx) {
  Rcpp::List cc = cpp_cholesky_to_cov(gamma, delta, R, edge_idx);
  arma::mat Sigma = cc["Sigma_b"];
  arma::mat Tinv = cc["T_inv"];
  const arma::uword E = edge_idx.n_rows;

  Rcpp::List D_delta(R), D_gamma(E);
  // delta derivatives: T^-1 E_kk (T^-1)^T
  for (int k = 0; k < R; ++k) {
    arma::vec ek(R, arma::fill::zeros); ek(k) = 1.0;
    arma::mat Dk = Tinv.col(k) * Tinv.col(k).t();
    D_delta[k] = 0.5 * (Dk + Dk.t());
  }
  // gamma derivatives: T^-1 E_ji Sigma + (that)^T
  for (arma::uword e = 0; e < E; ++e) {
    arma::uword i = edge_idx(e, 0) - 1;
    arma::uword j = edge_idx(e, 1) - 1;
    // T^-1 E_ji Sigma : E_ji has a single 1 at (j,i), so
    // T^-1 E_ji = (column j of T^-1) placed in column i; times Sigma:
    arma::mat term1(R, R, arma::fill::zeros);
    term1 = Tinv.col(j) * Sigma.row(i);   // (R x 1)(1 x R) = R x R
    arma::mat Dg = term1 + term1.t();
    D_gamma[e] = Dg;
  }
  return Rcpp::List::create(_["D_delta"] = D_delta, _["D_gamma"] = D_gamma);
}

// ---------------------------------------------------------------------------
// Closed-form layer: root, sensitivity and per-unit scores.
//
// Because |T| = 1, the Gaussian working criterion of the residual layer
// separates by node,
//   l = -(n/2) sum_j { log delta_j + q_j(gamma_j)/delta_j },
//   q_j = S_jj - 2 gamma_j' S_{P_j,j} + gamma_j' S_{P_jP_j} gamma_j,
// so the Pearson estimating equations of the directed layer are R independent
// least-squares problems and the root is available in closed form. The
// sensitivity inherits the same structure, one block per node, and the unit
// contribution of node j is the score of its own nodal regression.
// ---------------------------------------------------------------------------

namespace {

// Edge positions and parent indices of node j, in the order of edge_idx.
void node_parents(const arma::umat& edge_idx, arma::uword j,
                  std::vector<arma::uword>& parents,
                  std::vector<arma::uword>& positions) {
  parents.clear();
  positions.clear();
  for (arma::uword e = 0; e < edge_idx.n_rows; ++e) {
    if (edge_idx(e, 1) - 1 == j) {
      parents.push_back(edge_idx(e, 0) - 1);
      positions.push_back(e);
    }
  }
}

}  // namespace

//' Closed-form root of the directed layer (C++ kernel)
//'
//' Solves the Pearson estimating equations of the modified Cholesky layer node
//' by node: gamma_j = S[P_j, P_j]^-1 S[P_j, j] and
//' delta_j = S[j, j] - S[j, P_j] gamma_j. A node with no parent keeps
//' delta_j = S_jj and contributes no coefficient. Every delta_j is lifted to
//' floor_rel * S_jj when it falls below that floor.
//'
//' @param Shat R x R sample covariance of the Pearson residuals, in the
//'   positional (elimination) order of the nodes.
//' @param R integer number of responses.
//' @param edge_idx integer |E| x 2 matrix of edges (i, j), i < j, 1-based, in
//'   positional order.
//' @param floor_rel positivity floor of delta_j, relative to S_jj (to max |S_kk|
//'   when S_jj is not positive), so that it follows the scale of each response.
//' @param rcond_tol a parent block S[P_j, P_j] whose equilibrated form
//'   D^-1 S[P_j, P_j] D^-1, D = diag(S[P_j, P_j])^(1/2), has a reciprocal
//'   condition number not above this tolerance is near singular: its
//'   coefficients are not identified and are taken from a solve with a ridge
//'   of 1e-8 on the unit diagonal of the equilibrated block.
//' @return list with gamma (one per edge, in the order of edge_idx), delta
//'   (length R), the number of nodes whose normal equations needed a ridge,
//'   the number of innovation variances lifted to the positivity floor, and
//'   ill_conditioned, the positional (1-based) nodes whose parent block is
//'   near singular.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_dag_closed_root(const arma::mat& Shat, int R,
                               const arma::umat& edge_idx,
                               double floor_rel = 1e-10,
                               double rcond_tol = 3.666852862501036e-11) {
  if (Shat.n_rows != (arma::uword)R || Shat.n_cols != (arma::uword)R)
    stop("Shat must be square with one row per response.");
  const arma::uword E = edge_idx.n_rows;
  arma::vec gamma(E, arma::fill::zeros);
  arma::vec delta(R);
  int n_ridge = 0, n_floored = 0;
  std::vector<int> ill;

  const double scale = std::max(arma::abs(Shat.diag()).max(), 1e-300);

  std::vector<arma::uword> par, pos;
  for (arma::uword j = 0; j < (arma::uword)R; ++j) {
    node_parents(edge_idx, j, par, pos);
    delta(j) = Shat(j, j);
    // floor relative to the node's own variance, hence equivariant to scale
    const double floor_abs = floor_rel *
      ((std::isfinite(Shat(j, j)) && Shat(j, j) > 0.0) ? Shat(j, j) : scale);
    if (par.empty()) {
      if (delta(j) < floor_abs) { delta(j) = floor_abs; ++n_floored; }
      continue;
    }
    arma::uvec P(par);
    arma::uvec jj(1); jj(0) = j;
    arma::mat A = Shat.submat(P, P);
    A = 0.5 * (A + A.t());
    arma::vec b = Shat.submat(P, jj);

    // Solved on the equilibrated block, so the near-singularity test and the
    // fallback ridge do not depend on the scale of the responses.
    arma::vec g;
    bool ill_j = false, ridged = false;
    const bool ok = mcglmgraph::equilibrated_solve(A, b, rcond_tol, 1e-8, true,
                                                   g, ill_j, ridged);
    if (ill_j) ill.push_back((int)j + 1);
    if (ridged) ++n_ridge;
    if (!ok || !g.is_finite())
      stop("The normal equations of a node are not solvable.");

    for (std::size_t t = 0; t < pos.size(); ++t) gamma(pos[t]) = g(t);
    delta(j) = Shat(j, j) - arma::dot(b, g);
    if (!(delta(j) > floor_abs)) { delta(j) = floor_abs; ++n_floored; }
  }
  return Rcpp::List::create(_["gamma"] = gamma, _["delta"] = delta,
                            _["n_ridge"] = n_ridge, _["n_floored"] = n_floored,
                            _["ill_conditioned"] = Rcpp::IntegerVector(
                              ill.begin(), ill.end()));
}

//' Block-diagonal sensitivity of the directed layer (C++ kernel)
//'
//' Assembles J_n on the accumulated scale from its closed form: the block of
//' the coefficients of node j is
//' n Sigma[P_j, P_j] / delta_j, that of its innovation variance is
//' n / (2 delta_j^2), and every other block is zero, including the one
//' between coefficients of distinct nodes.
//'
//' @param Sigma R x R fitted covariance of the layer.
//' @param delta vector of R innovation variances.
//' @param n number of experimental units.
//' @param edge_idx integer |E| x 2 matrix of edges, 1-based.
//' @return the (|E| + R) x (|E| + R) sensitivity in the parameter order
//'   (gamma, delta).
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_dag_sensitivity(const arma::mat& Sigma, const arma::vec& delta,
                              double n, const arma::umat& edge_idx) {
  const arma::uword R = Sigma.n_rows;
  if (delta.n_elem != R) stop("delta must have one value per response.");
  if (delta.min() <= 0.0) stop("delta must be strictly positive.");
  const arma::uword E = edge_idx.n_rows;
  arma::mat J(E + R, E + R, arma::fill::zeros);

  std::vector<arma::uword> par, pos;
  for (arma::uword j = 0; j < R; ++j) {
    node_parents(edge_idx, j, par, pos);
    if (!par.empty()) {
      arma::uvec P(par), Q(pos);
      arma::mat blk = Sigma.submat(P, P);
      J.submat(Q, Q) = (n / delta(j)) * 0.5 * (blk + blk.t());
    }
    J(E + j, E + j) = n / (2.0 * delta(j) * delta(j));
  }
  return J;
}

//' Per-unit Pearson scores of the directed layer (C++ kernel)
//'
//' The unit contribution of the coefficients of node j is
//' eps_ij r[i, P_j] / delta_j, with eps_ij = r[i, j] - gamma_j' r[i, P_j], and
//' that of its innovation variance is (eps_ij^2 - delta_j) / (2 delta_j^2).
//' Summed over units these are the estimating functions whose root is
//' cpp_dag_closed_root(), and their outer product is the variability of the
//' sandwich.
//'
//' @param Ec n x R matrix of centered Pearson residuals.
//' @param gamma coefficients in the order of edge_idx.
//' @param delta vector of R innovation variances.
//' @param edge_idx integer |E| x 2 matrix of edges, 1-based.
//' @return list with scores (n x (|E| + R), in the order (gamma, delta)) and
//'   innovations (n x R).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_dag_scores(const arma::mat& Ec, const arma::vec& gamma,
                          const arma::vec& delta, const arma::umat& edge_idx) {
  const arma::uword R = Ec.n_cols;
  const arma::uword n = Ec.n_rows;
  if (delta.n_elem != R) stop("delta must have one value per response.");
  if (delta.min() <= 0.0) stop("delta must be strictly positive.");
  const arma::uword E = edge_idx.n_rows;
  if (gamma.n_elem != E) stop("gamma must have one value per edge.");

  arma::mat Psi(n, E + R, arma::fill::zeros);
  arma::mat eps(n, R);

  std::vector<arma::uword> par, pos;
  for (arma::uword j = 0; j < R; ++j) {
    node_parents(edge_idx, j, par, pos);
    arma::vec e_j = Ec.col(j);
    if (!par.empty()) {
      arma::uvec P(par), Q(pos);
      arma::mat Xp = Ec.cols(P);
      e_j -= Xp * gamma.elem(Q);
      Psi.cols(Q) = Xp.each_col() % (e_j / delta(j));
    }
    eps.col(j) = e_j;
    Psi.col(E + j) = (e_j % e_j - delta(j)) / (2.0 * delta(j) * delta(j));
  }
  return Rcpp::List::create(_["scores"] = Psi, _["innovations"] = eps);
}

//' Pearson estimating function of the directed layer in closed form (C++ kernel)
//'
//' Evaluates the estimating function on the accumulated scale directly from the
//' sample covariance, without forming the generic traces:
//'   psi of gamma_j = (n/delta_j) (S[P_j, j] - S[P_j, P_j] gamma_j),
//'   psi of delta_j = (n/2) (q_j - delta_j) / delta_j^2,
//' with q_j = S[j, j] - 2 gamma_j' S[P_j, j] + gamma_j' S[P_j, P_j] gamma_j. It
//' vanishes at the root of cpp_dag_closed_root() and is what certifies it.
//'
//' @param Shat R x R sample covariance, on the ordering of the support.
//' @param gamma coefficients in the order of edge_idx.
//' @param delta vector of R innovation variances.
//' @param n number of experimental units.
//' @param edge_idx integer |E| x 2 matrix of edges, 1-based.
//' @return the estimating function in the order (gamma, delta).
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_dag_pearson_closed(const arma::mat& Shat, const arma::vec& gamma,
                                 const arma::vec& delta, double n,
                                 const arma::umat& edge_idx) {
  const arma::uword R = Shat.n_rows;
  if (delta.n_elem != R) stop("delta must have one value per response.");
  if (delta.min() <= 0.0) stop("delta must be strictly positive.");
  const arma::uword E = edge_idx.n_rows;
  if (gamma.n_elem != E) stop("gamma must have one value per edge.");

  arma::vec psi(E + R, arma::fill::zeros);
  std::vector<arma::uword> par, pos;
  for (arma::uword j = 0; j < R; ++j) {
    node_parents(edge_idx, j, par, pos);
    double q = Shat(j, j);
    if (!par.empty()) {
      arma::uvec P(par), Q(pos), jj(1);
      jj(0) = j;
      arma::vec g = gamma.elem(Q);
      arma::mat A = Shat.submat(P, P);
      arma::vec b = Shat.submat(P, jj);
      arma::vec resid = b - A * g;
      psi.elem(Q) = (n / delta(j)) * resid;
      q += -2.0 * arma::dot(g, b) + arma::as_scalar(g.t() * A * g);
    }
    psi(E + j) = 0.5 * n * (q - delta(j)) / (delta(j) * delta(j));
  }
  return psi;
}

//' Inverse of the block-diagonal sensitivity of the directed layer (C++ kernel)
//'
//' The sensitivity of the directed layer is block diagonal with one block per
//' node, so its inverse is too: the block of node j is delta_j Sigma[P_j, P_j]^-1
//' / n and the entry of its innovation variance is 2 delta_j^2 / n. Inverting
//' block by block costs the sum of the cubes of the parent-set sizes instead of
//' the cube of the number of parameters.
//'
//' @inheritParams cpp_dag_sensitivity
//' @return list with the inverse and n_singular, the number of nodes whose
//'   parent block could not be inverted; a nonzero value means the caller must
//'   fall back on the general route.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_dag_sensitivity_inverse(const arma::mat& Sigma,
                                       const arma::vec& delta, double n,
                                       const arma::umat& edge_idx) {
  const arma::uword R = Sigma.n_rows;
  if (delta.n_elem != R) stop("delta must have one value per response.");
  if (delta.min() <= 0.0) stop("delta must be strictly positive.");
  const arma::uword E = edge_idx.n_rows;
  arma::mat Jinv(E + R, E + R, arma::fill::zeros);
  int n_singular = 0;

  std::vector<arma::uword> par, pos;
  for (arma::uword j = 0; j < R; ++j) {
    node_parents(edge_idx, j, par, pos);
    if (!par.empty()) {
      arma::uvec P(par), Q(pos);
      arma::mat blk = Sigma.submat(P, P);
      blk = 0.5 * (blk + blk.t());
      arma::mat inv_blk;
      if (arma::inv_sympd(inv_blk, blk) && inv_blk.is_finite())
        Jinv.submat(Q, Q) = (delta(j) / n) * inv_blk;
      else
        ++n_singular;
    }
    Jinv(E + j, E + j) = 2.0 * delta(j) * delta(j) / n;
  }
  return Rcpp::List::create(_["inverse"] = Jinv,
                            _["n_singular"] = n_singular);
}
