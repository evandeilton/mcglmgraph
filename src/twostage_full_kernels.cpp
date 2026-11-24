// twostage_full_kernels.cpp: the full two-stage correction, for fits the
// shortcut does not reach (RcppArmadillo).
// Part of the mcglmgraph package. License: GPL-3.
//
// Stacking the first-stage and the graph estimating functions gives the
// influence function of the second stage,
//   IF_i = H^-1 ( psi_2i - B A^-1 psi_1i ),
// with A the sensitivity of the first stage and B = -E( d psi_2i / d alpha' )
// the cross sensitivity. The shortcut in twostage_kernels.cpp replaces
// B A^-1 psi_1i by a closed form that needs second-stage quantities only; that
// closed form rests on a single multiplicative dispersion per response, an
// identity within-response structure and constant prior weights, and on a
// first stage that moves the residual covariance as a diagonal rescaling.
// Outside those conditions the three pieces have to be formed, and that is what
// this file does.
//
// The first stage is the MCGLM of mcglm. With the R responses stacked response
// by response, the quasi-score and the Pearson functions are
//   psi_beta     = D' C^-1 W (y - mu),
//   psi_lambda_l = (y - mu)' C^-1 Cdot_l W C^-1 (y - mu) - tr(C^-1 Cdot_l W),
// with Cdot_l = dC / d lambda_l. When C is block diagonal over a partition of
// the rows (the rows of one unit, or of one cluster, across the R responses),
// both are sums over the blocks, and the block sums are the unit contributions
// psi_1u.
//
// The expected sensitivity is E(d psi_1 / d alpha') under a correctly specified
// first stage, E(y - mu) = 0 and E{(y - mu)(y - mu)'} = C:
//   S_bb  = -D' C^-1 W D,
//   S_ll' = -tr( C^-1 Cdot_l W C^-1 Cdot_l' ),
//   S_lb  = -tr( C^-1 Cdot_l W C^-1 Cdot_b ),   Cdot_b = dC / d beta_b,
//   S_bl  = 0,
// the last because the quasi-score is linear in y - mu and E(y - mu) = 0. mcglm
// reports the first two and sets S_lb to zero, since the covariance it builds at
// the end of the fit carries no derivative in beta. S_lb is not zero when the
// variance depends on the mean. Under the conditions of the closed form, with a
// correctly specified mean, the regression columns of B A^-1 are zero, and with
// S_lb set to zero they are not: the regression part of the first stage stops
// canceling from the propagation term.
//
// The Pearson residual of the second stage is r_ri = (y_ri - mu_ri) / s_ri, with
// s_ri^2 = C_(ri),(ri). Its derivative in the first-stage parameters is
//   d r_ri / d alpha = -{ s_ri^-1 d mu_ri / d alpha + r_ri d log s_ri / d alpha },
// with d log s_ri / d alpha = (1/2) [dC / d alpha]_(ri),(ri) / s_ri^2. Writing
// g_ri(alpha) for the expression in braces and e_i for the centered residuals,
//   B_n(d, .) = sum_i e_i' Sigma_b^-1 Sdot_d Sigma_b^-1 g_i,
// on the accumulated scale. Centering enters through e_i only, because the mean
// of the derivative multiplies residuals that sum to zero.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

namespace {

arma::vec sp_diag(const arma::sp_mat& M) {
  return arma::vec(arma::mat(M.diag()));
}

// columns of M multiplied by w: M W with W = diag(w)
arma::sp_mat scale_cols(arma::sp_mat M, const arma::vec& w) {
  for (arma::sp_mat::iterator it = M.begin(); it != M.end(); ++it)
    (*it) *= w(it.col());
  return M;
}

}  // namespace

//' First stage of an MCGLM, unit by unit (C++ kernel)
//'
//' Forms, for every unit of a partition of the rows, the sum over its rows of
//' the quasi-score and of the Pearson estimating functions of the first stage;
//' the sensitivity of the first stage from its definition, the block of the
//' covariance parameters on the regression parameters included; and the
//' derivative of the Pearson residuals in the first-stage parameters.
//'
//' @param invC the \eqn{nR \times nR} inverse of the first-stage covariance.
//' @param C the first-stage covariance.
//' @param Dblocks list of \eqn{R} dense \eqn{n \times p_r} derivatives of the
//'   means in the regression parameters of each response.
//' @param DC list of the derivatives of \eqn{C} in the covariance parameters.
//' @param DCbeta list, one entry per regression parameter, of the derivatives
//'   of \eqn{C} in that parameter, or \code{NULL} where the variance does not
//'   depend on the mean.
//' @param res raw residuals \eqn{y - \mu}, stacked response by response.
//' @param w prior weights of the first stage, one per stacked row.
//' @param unit zero-based unit of each observational row, length \eqn{n}.
//' @param n_units number of units.
//' @return A list with \code{psi} (\code{n_units} by \eqn{p_1}), the
//'   sensitivity \code{S} (\eqn{p_1 \times p_1}, on the sign convention of
//'   mcglm, \eqn{E(\partial\psi/\partial\alpha^\top)}), \code{G} (\eqn{nR} by
//'   \eqn{p_1}, the derivative of the Pearson residuals with the sign removed),
//'   \code{r} (the Pearson residuals), and the numbers \code{n_beta} of
//'   regression parameters and \code{n_cov} of covariance parameters, with
//'   \eqn{p_1} their sum.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_first_stage_units(const arma::sp_mat& invC, const arma::sp_mat& C,
                                 const Rcpp::List& Dblocks,
                                 const Rcpp::List& DC, const Rcpp::List& DCbeta,
                                 const arma::vec& res, const arma::vec& w,
                                 const arma::uvec& unit,
                                 const arma::uword n_units) {
  const arma::uword R = Dblocks.size();
  const arma::uword N = res.n_elem;
  if (R == 0 || N % R != 0)
    stop("the stacked residuals do not split into R responses.");
  const arma::uword n = N / R;
  if (unit.n_elem != n) stop("unit must have one entry per observational row.");
  if (invC.n_rows != N || invC.n_cols != N || C.n_rows != N || C.n_cols != N)
    stop("C and its inverse must match the stacked residuals.");
  if (w.n_elem != N) stop("w must have one weight per stacked row.");
  if (n_units == 0 || arma::any(unit >= n_units)) stop("unit index out of range.");

  // block-diagonal derivative of the stacked means
  std::vector<arma::mat> Dr(R);
  arma::uword pb = 0;
  for (arma::uword r = 0; r < R; ++r) {
    Dr[r] = as<arma::mat>(Dblocks[r]);
    if (Dr[r].n_rows != n) stop("every mean derivative must have one row per unit.");
    pb += Dr[r].n_cols;
  }
  const arma::uword pc = DC.size();
  if ((arma::uword) DCbeta.size() != pb)
    stop("DCbeta must have one entry per regression parameter.");
  arma::sp_mat D(N, pb);
  {
    arma::uword off = 0;
    for (arma::uword r = 0; r < R; ++r) {
      if (Dr[r].n_cols)
        D.submat(r * n, off, (r + 1) * n - 1, off + Dr[r].n_cols - 1) =
          arma::sp_mat(Dr[r]);
      off += Dr[r].n_cols;
    }
  }

  arma::uvec urow(N);
  for (arma::uword a = 0; a < N; ++a) urow(a) = unit(a % n);

  arma::mat psi(n_units, pb + pc, arma::fill::zeros);

  // quasi-score, D' C^-1 W res, split by unit
  arma::vec z = invC * (w % res);
  for (arma::sp_mat::const_iterator it = D.begin(); it != D.end(); ++it)
    psi(urow(it.row()), it.col()) += (*it) * z(it.row());

  // Pearson functions, zz' Cdot_l W zz - tr(C^-1 Cdot_l W) with zz = C^-1 res,
  // split by unit. diag(C^-1 Cdot_l) is read from the column sums of the
  // elementwise product, since C^-1 is symmetric.
  arma::vec zz = invC * res;
  arma::vec wzz = w % zz;
  std::vector<arma::sp_mat> Pc(pc);                // C^-1 Cdot_l
  std::vector<arma::sp_mat> PW(pc);                // C^-1 Cdot_l W
  std::vector<arma::sp_mat> dc(pc);
  for (arma::uword l = 0; l < pc; ++l) {
    dc[l] = as<arma::sp_mat>(DC[l]);
    if (dc[l].n_rows != N || dc[l].n_cols != N)
      stop("every covariance derivative must match the stacked residuals.");
    arma::vec q = zz % (dc[l] * wzz);
    arma::vec tr = arma::vec(arma::mat(arma::sum(invC % dc[l], 0)).t());
    for (arma::uword a = 0; a < N; ++a)
      psi(urow(a), pb + l) += q(a) - tr(a) * w(a);
    Pc[l] = invC * dc[l];
    PW[l] = scale_cols(Pc[l], w);
  }

  // derivatives of C in beta, as C^-1 Cdot_b
  std::vector<arma::sp_mat> Qb(pb);
  std::vector<bool> has_b(pb, false);
  std::vector<arma::vec> dgb(pb);
  for (arma::uword k = 0; k < pb; ++k) {
    if (Rf_isNull(DCbeta[k])) continue;
    arma::sp_mat Mk = as<arma::sp_mat>(DCbeta[k]);
    if (Mk.n_rows != N || Mk.n_cols != N)
      stop("every derivative in beta must match the stacked residuals.");
    Qb[k] = invC * Mk;
    dgb[k] = sp_diag(Mk);
    has_b[k] = true;
  }

  // sensitivity from its definition, on mcglm's sign convention
  arma::mat S(pb + pc, pb + pc, arma::fill::zeros);
  arma::sp_mat WD = D;                             // W D: rows scaled by w
  for (arma::sp_mat::iterator it = WD.begin(); it != WD.end(); ++it)
    (*it) *= w(it.row());
  S.submat(0, 0, pb - 1, pb - 1) = -arma::mat(D.t() * (invC * WD));
  for (arma::uword l = 0; l < pc; ++l) {
    for (arma::uword m = 0; m <= l; ++m) {
      // tr(X Y) = sum(X % Y'), for the products of the two derivatives
      double v = -arma::accu(PW[l] % Pc[m].t());
      S(pb + l, pb + m) = v;
      S(pb + m, pb + l) = v;
    }
    for (arma::uword k = 0; k < pb; ++k)
      if (has_b[k]) S(pb + l, k) = -arma::accu(PW[l] % Qb[k].t());
  }

  // derivative of the Pearson residuals, sign removed
  arma::vec s2 = sp_diag(C);
  if (arma::any(s2 <= 0)) stop("non-positive marginal variances in C.");
  arma::vec sdv = arma::sqrt(s2);
  arma::vec r = res / sdv;
  arma::mat G(N, pb + pc, arma::fill::zeros);
  for (arma::sp_mat::const_iterator it = D.begin(); it != D.end(); ++it)
    G(it.row(), it.col()) += (*it) / sdv(it.row());
  for (arma::uword k = 0; k < pb; ++k)
    if (has_b[k]) G.col(k) += r % (dgb[k] / (2.0 * s2));
  for (arma::uword l = 0; l < pc; ++l)
    G.col(pb + l) = r % (sp_diag(dc[l]) / (2.0 * s2));

  return Rcpp::List::create(_["psi"] = psi, _["S"] = S, _["G"] = G, _["r"] = r,
                            _["n_beta"] = (int) pb, _["n_cov"] = (int) pc);
}

//' Cross sensitivity between the graph layer and the first stage (C++ kernel)
//'
//' Forms \eqn{\mathbf B_n(d,\cdot)=\sum_i\mathbf e_i^\top\boldsymbol\Sigma_b^{-1}
//' \dot{\boldsymbol\Sigma}_d\boldsymbol\Sigma_b^{-1}\mathbf g_i}, the
//' accumulated cross sensitivity \eqn{-\partial\mathbf U_2/\partial
//' \boldsymbol\alpha^\top}, from the centered residuals and the derivative of the
//' Pearson residuals returned by \code{cpp_first_stage_units}.
//'
//' @param Ec \eqn{n \times R} centered Pearson residuals.
//' @param Sinv the inverse of the fitted \eqn{\boldsymbol\Sigma_b}.
//' @param Dlist list of the \eqn{P} derivatives of \eqn{\boldsymbol\Sigma_b}.
//' @param G \eqn{nR \times p_1} derivative of the Pearson residuals, stacked
//'   response by response.
//' @return The \eqn{P \times p_1} matrix \eqn{\mathbf B_n}.
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_two_stage_cross(const arma::mat& Ec, const arma::mat& Sinv,
                              const Rcpp::List& Dlist, const arma::mat& G) {
  const arma::uword n = Ec.n_rows;
  const arma::uword R = Ec.n_cols;
  const arma::uword P = Dlist.size();
  if (G.n_rows != n * R) stop("G must have one row per stacked observation.");
  if (Sinv.n_rows != R || Sinv.n_cols != R) stop("Sinv must be R by R.");
  arma::mat B(P, G.n_cols, arma::fill::zeros);
  for (arma::uword d = 0; d < P; ++d) {
    arma::mat Dd = as<arma::mat>(Dlist[d]);
    if (Dd.n_rows != R || Dd.n_cols != R)
      stop("every derivative must be square with one row per response.");
    arma::mat F = Ec * (Sinv * Dd * Sinv);      // row i is e_i' A_d
    B.row(d) = arma::vectorise(F).t() * G;      // column-major, response-major
  }
  return B;
}

//' Largest absolute difference between two sparse matrices (C++ kernel)
//'
//' @param A,B sparse matrices of the same dimensions.
//' @return \eqn{\max_{ab}|A_{ab}-B_{ab}|}, or \code{Inf} when the dimensions
//'   differ.
//' @keywords internal
// [[Rcpp::export]]
double cpp_sp_max_abs_diff(const arma::sp_mat& A, const arma::sp_mat& B) {
  if (A.n_rows != B.n_rows || A.n_cols != B.n_cols)
    return R_PosInf;
  arma::sp_mat Dm = A - B;
  double m = 0.0;
  for (arma::sp_mat::const_iterator it = Dm.begin(); it != Dm.end(); ++it)
    m = std::max(m, std::abs(*it));
  return m;
}

//' Totals of the first-stage estimating functions (C++ kernel)
//'
//' The quasi-score and the Pearson functions of an MCGLM summed over all rows,
//' \eqn{\mathbf D^\top\mathbf C^{-1}\mathbf W(\mathbf y-\boldsymbol\mu)} and
//' \eqn{(\mathbf y-\boldsymbol\mu)^\top\mathbf C^{-1}\dot{\mathbf C}_l\mathbf W
//' \mathbf C^{-1}(\mathbf y-\boldsymbol\mu)-\operatorname{tr}(\mathbf C^{-1}
//' \dot{\mathbf C}_l\mathbf W)}, without the sensitivity: what a numerical
//' derivative of the first stage needs at every step.
//'
//' @inheritParams cpp_first_stage_units
//' @return The vector of the \eqn{p_1} totals.
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_first_stage_totals(const arma::sp_mat& invC, const Rcpp::List& Dblocks,
                                 const Rcpp::List& DC, const arma::vec& res,
                                 const arma::vec& w) {
  const arma::uword R = Dblocks.size();
  const arma::uword N = res.n_elem;
  if (R == 0 || N % R != 0)
    stop("the stacked residuals do not split into R responses.");
  const arma::uword n = N / R;
  if (invC.n_rows != N || invC.n_cols != N || w.n_elem != N)
    stop("C^-1 and the weights must match the stacked residuals.");
  arma::vec z = invC * (w % res);
  std::vector<double> out;
  for (arma::uword r = 0; r < R; ++r) {
    arma::mat Dr = as<arma::mat>(Dblocks[r]);
    if (Dr.n_rows != n) stop("every mean derivative must have one row per unit.");
    arma::vec s = Dr.t() * z.subvec(r * n, (r + 1) * n - 1);
    out.insert(out.end(), s.begin(), s.end());
  }
  arma::vec zz = invC * res;
  arma::vec wzz = w % zz;
  for (arma::uword l = 0; l < (arma::uword) DC.size(); ++l) {
    arma::sp_mat Dl = as<arma::sp_mat>(DC[l]);
    if (Dl.n_rows != N || Dl.n_cols != N)
      stop("every covariance derivative must match the stacked residuals.");
    double q = arma::dot(zz, Dl * wzz);
    arma::rowvec tr = arma::rowvec(arma::mat(arma::sum(invC % Dl, 0)));
    out.push_back(q - arma::dot(tr, w));
  }
  return arma::vec(out);
}
