// selection_kernels.cpp: graph structure selection kernels (RcppArmadillo)
// Part of the mcglmgraph package. License: GPL-3.
//
// Second-order Pearson estimating functions for the between-response
// covariance, edge score/Wald statistics, soft-thresholding and the
// enumeration of unshielded triples for CPDAG discovery.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

//' Pearson estimating function and sensitivity for Sigma_b (C++ kernel)
//'
//' For n i.i.d. units with sample covariance Shat, computes the second-order
//' Pearson estimating function psi and the sensitivity J for the covariance
//' parameters, given the analytical derivatives Dlist of Sigma:
//' psi_d = -0.5 n [tr(S^-1 D_d) - tr(S^-1 D_d S^-1 Shat)],
//' J_ab = 0.5 n tr(S^-1 D_a S^-1 D_b),
//' with S = Sigma. Sigma must be symmetric positive definite.
//' @param Sigma current R x R symmetric positive definite covariance matrix.
//' @param Shat R x R sample covariance of the Pearson residuals.
//' @param n number of units.
//' @param Dlist list of P derivative matrices dSigma/dtheta_d.
//' @return list with psi (length P) and J (P x P).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_pearson_between(const arma::mat& Sigma, const arma::mat& Shat,
                               double n, const Rcpp::List& Dlist) {
  const int P = Dlist.size();
  arma::mat Sinv = arma::inv_sympd(0.5 * (Sigma + Sigma.t()));
  arma::vec psi(P, arma::fill::zeros);
  std::vector<arma::mat> SinvD(P);
  for (int d = 0; d < P; ++d) {
    arma::mat Dd = as<arma::mat>(Dlist[d]);
    arma::mat SD = Sinv * Dd;
    SinvD[d] = SD;
    double t1 = arma::trace(SD);
    double t2 = arma::trace(SD * Sinv * Shat);
    psi(d) = -0.5 * n * (t1 - t2);
  }
  arma::mat J(P, P, arma::fill::zeros);
  for (int a = 0; a < P; ++a) {
    for (int b = a; b < P; ++b) {
      double v = 0.5 * n * arma::trace(SinvD[a] * SinvD[b]);
      J(a, b) = v;
      J(b, a) = v;
    }
  }
  return Rcpp::List::create(_["psi"] = psi, _["J"] = J);
}

//' Soft-thresholding operator (C++ kernel)
//' @param x numeric vector.
//' @param thr non-negative threshold.
//' @return sign(x) * max(|x| - thr, 0).
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_soft_threshold(const arma::vec& x, double thr) {
  arma::vec a = arma::abs(x) - thr;
  a.transform([](double v) { return v > 0.0 ? v : 0.0; });
  return arma::sign(x) % a;
}

//' Per-edge Wald statistics (C++ kernel)
//' @param eta full parameter vector.
//' @param Vcov estimated covariance matrix of eta (model-based, J^-1).
//' @param edge_pos 0-based positions of edge parameters in eta.
//' @return vector of Wald statistics eta_e^2 / Vcov_ee; NA where the variance
//'   is not finite and positive, since the statistic is then undefined.
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_wald_edges(const arma::vec& eta, const arma::mat& Vcov,
                         const arma::uvec& edge_pos) {
  const arma::uword q = edge_pos.n_elem;
  arma::vec W(q);
  for (arma::uword e = 0; e < q; ++e) {
    arma::uword p = edge_pos(e);
    double v = Vcov(p, p);
    // an undefined statistic is NA, never zero, so it cannot rank as weakest
    W(e) = (std::isfinite(v) && v > 0.0) ? (eta(p) * eta(p) / v) : NA_REAL;
  }
  return W;
}

//' Enumerate unshielded triples of an undirected graph (C++ kernel)
//'
//' Returns triples (i, k, j) such that i-k and k-j are edges but i-j is not,
//' the candidate v-structures of CPDAG discovery.
//' @param adj integer adjacency matrix (R x R, 0/1, symmetric).
//' @return integer matrix with one row per triple, columns (i, k, j), 1-based.
//' @keywords internal
// [[Rcpp::export]]
arma::umat cpp_unshielded_triples(const arma::umat& adj) {
  const int R = adj.n_rows;
  std::vector<arma::uword> vi, vk, vj;
  for (int k = 0; k < R; ++k) {
    for (int i = 0; i < R; ++i) {
      if (i == k || adj(i, k) == 0) continue;
      for (int j = i + 1; j < R; ++j) {
        if (j == k || adj(j, k) == 0) continue;
        if (adj(i, j) == 0) {
          vi.push_back(i + 1); vk.push_back(k + 1); vj.push_back(j + 1);
        }
      }
    }
  }
  const arma::uword T = vi.size();
  arma::umat out(T, 3);
  for (arma::uword t = 0; t < T; ++t) {
    out(t, 0) = vi[t]; out(t, 1) = vk[t]; out(t, 2) = vj[t];
  }
  return out;
}

//' Score statistics for candidate (absent) edges (C++ kernel)
//'
//' For each candidate edge e=(i,j) of the linear-precision parametrization, the
//' Pearson score is psi_e = n (Sigma - Shat)[i,j] and the sensitivity is
//' J_ee = n (Sigma_ii Sigma_jj + Sigma_ij^2); the statistic is
//' S_e = psi_e^2 / J_ee, with S_e = 0 when J_ee is not positive. It ignores the
//' estimation of the other parameters, so it ranks candidates but is not
//' calibrated as chi-square(1); the nuisance-adjusted statistic is formed in R
//' by mc_edge_scores().
//' @param Sigma R x R fitted covariance matrix.
//' @param Shat R x R sample covariance of the Pearson residuals.
//' @param n number of units.
//' @param edges integer matrix q x 2 of candidate edges (i, j), 1-based.
//' @return numeric vector of score statistics.
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_score_edges(const arma::mat& Sigma, const arma::mat& Shat,
                          double n, const arma::umat& edges) {
  const arma::uword q = edges.n_rows;
  arma::vec S(q);
  for (arma::uword e = 0; e < q; ++e) {
    arma::uword i = edges(e, 0) - 1, j = edges(e, 1) - 1;
    double psi = n * (Sigma(i, j) - Shat(i, j));
    double Jee = n * (Sigma(i, i) * Sigma(j, j) + Sigma(i, j) * Sigma(i, j));
    S(e) = (Jee > 0.0) ? psi * psi / Jee : 0.0;
  }
  return S;
}
