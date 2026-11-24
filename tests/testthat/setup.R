# =============================================================================
# SHARED FIXTURES AND HELPERS
# =============================================================================
# Numerical-differentiation helpers used to verify the analytical derivatives
# and the estimating system. Central differences with step 1e-5: the truncation
# error is of order eps^2 and the rounding error of order
# .Machine$double.eps / eps, both close to 1e-11 at this choice.

EPS_NUM <- 1e-5

## Numerical derivative of a matrix-valued function along coordinate k.
num_deriv_matrix <- function(f, par, k, eps = EPS_NUM) {
  pp <- par; pp[k] <- pp[k] + eps
  pm <- par; pm[k] <- pm[k] - eps
  (f(pp) - f(pm)) / (2 * eps)
}

## Numerical gradient of a scalar-valued function.
num_gradient <- function(f, par, eps = EPS_NUM) {
  vapply(seq_along(par), function(k) {
    pp <- par; pp[k] <- pp[k] + eps
    pm <- par; pm[k] <- pm[k] - eps
    (f(pp) - f(pm)) / (2 * eps)
  }, numeric(1))
}

## Numerical Jacobian of a vector-valued function.
num_jacobian <- function(f, par, eps = EPS_NUM) {
  cols <- lapply(seq_along(par), function(k) {
    pp <- par; pp[k] <- pp[k] + eps
    pm <- par; pm[k] <- pm[k] - eps
    (f(pp) - f(pm)) / (2 * eps)
  })
  do.call(cbind, cols)
}

## Scale-robust discrepancy: the denominator prevents an entry close to zero
## from inflating the relative error.
rel_gap <- function(analytic, numeric_value) {
  max(abs(analytic - numeric_value)) / max(1, max(abs(analytic)))
}

## Gaussian working log-criterion of the residual layer,
## up to the constant, written as an explicit function of Sigma_b.
working_ploglik <- function(Sigma_b, Shat, n) {
  -0.5 * n * (as.numeric(determinant(Sigma_b, logarithm = TRUE)$modulus) +
                sum(solve(Sigma_b) * Shat))
}

## ---- GGM configurations: dimension, support and admissible eta -------------
ggm_configs <- list(
  complete_R3 = list(
    R = 3L, edges = NULL,
    eta = c(2.0, 2.5, 3.0, -0.5, -0.3, -0.4)
  ),
  chain_R4 = list(
    R = 4L, edges = rbind(c(1, 2), c(2, 3), c(3, 4)),
    eta = c(2.0, 2.2, 2.4, 2.6, -0.5, -0.4, -0.3)
  ),
  sparse_R5 = list(
    R = 5L, edges = rbind(c(1, 2), c(1, 3), c(4, 5)),
    eta = c(1.8, 2.0, 2.2, 2.4, 2.6, -0.6, -0.35, -0.45)
  )
)

## ---- DAG configurations: support, gamma and delta --------------------------
dag_configs <- list(
  complete_R3 = list(
    R = 3L, edges = rbind(c(1, 2), c(1, 3), c(2, 3)),
    gamma = c(0.7, 0.3, 0.5), delta = c(1.0, 0.64, 0.49)
  ),
  chain_R4 = list(
    R = 4L, edges = rbind(c(1, 2), c(2, 3), c(3, 4)),
    gamma = c(0.6, -0.4, 0.5), delta = c(1.2, 0.8, 0.7, 0.9)
  ),
  hub_R5 = list(
    R = 5L, edges = rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5)),
    gamma = c(0.5, -0.3, 0.4, 0.2), delta = c(1.0, 0.9, 0.8, 0.7, 0.6)
  )
)

## Sigma_b of the GGM layer as a function of the parameter vector.
ggm_sigma <- function(eta, bases) {
  mc_precision_to_covariance(mc_build_precision(eta, bases))$Sigma_b
}

## Sigma_b of the DAG layer as a function of stacked (gamma, delta).
dag_sigma <- function(par, R, edges) {
  E <- nrow(edges)
  mc_cholesky_to_covariance(par[seq_len(E)], par[E + seq_len(R)],
                            R, edges)$Sigma_b
}

## Derivative list of the DAG layer in the same order as the (gamma, delta)
## parameter vector, which is the order the internal engine uses.
dag_derivatives <- function(par, R, edges) {
  E <- nrow(edges)
  der <- mc_cholesky_derivatives(par[seq_len(E)], par[E + seq_len(R)],
                                 R, edges)
  c(der$D_gamma, der$D_delta)
}

## Acceptance checks that run a Monte Carlo take minutes rather than seconds.
## They run when MCGLMGRAPH_LONG_TESTS is "true" and are skipped otherwise, so
## the routine suite stays fast; their recorded output lives with the checks of
## the release.
skip_long <- function() {
  testthat::skip_if_not(identical(Sys.getenv("MCGLMGRAPH_LONG_TESTS"), "true"),
                        "long check; set MCGLMGRAPH_LONG_TESTS=true to run it")
}
