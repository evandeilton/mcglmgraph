# =============================================================================
# NUMERICAL VERIFICATION OF THE ANALYTICAL DERIVATIVES AND ESTIMATING SYSTEM
# =============================================================================
# The closed-form derivatives are the central computational contribution of the
# package. Checking a single coordinate per layer would let through a sign
# error confined to a diagonal base or to one specific edge orientation, so
# these tests sweep EVERY coordinate, in more than one
# configuration of dimension and support, and also verify the two most used
# quantities of the engine: the estimating function and the sensitivity.

# ─────────────────────────────────────────────────────────────────────────────
# GGM LAYER DERIVATIVES
# ─────────────────────────────────────────────────────────────────────────────

test_that("GGM derivatives match numerical differentiation in every coordinate", {
  for (nm in names(ggm_configs)) {
    cfg <- ggm_configs[[nm]]
    bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
    eta <- cfg$eta
    expect_true(mc_check_spd(mc_build_precision(eta, bases))$is_spd,
                info = paste("configuration", nm, "must start from an SPD K_b"))

    Sigma <- ggm_sigma(eta, bases)
    D <- mc_precision_derivatives(Sigma, bases)
    expect_equal(length(D), length(eta))

    for (k in seq_along(eta)) {
      num <- num_deriv_matrix(function(p) ggm_sigma(p, bases), eta, k)
      expect_lt(rel_gap(D[[k]], num), 1e-6)
      ## the derivative is symmetric by construction
      expect_lt(max(abs(D[[k]] - t(D[[k]]))), 1e-12)
    }
  }
})

test_that("GGM diagonal derivatives are rank one and negative semidefinite", {
  ## dSigma/deta_kk = -Sigma E_kk Sigma is minus an outer product, hence
  ## negative semidefinite of rank one. It is the testable counterpart of the
  ## low-rank form at a diagonal
  ## coordinate.
  cfg <- ggm_configs$chain_R4
  bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
  Sigma <- ggm_sigma(cfg$eta, bases)
  D <- mc_precision_derivatives(Sigma, bases)
  for (k in seq_len(cfg$R)) {
    ev <- eigen(D[[k]], symmetric = TRUE, only.values = TRUE)$values
    expect_lt(max(ev), 1e-10)
    expect_equal(sum(abs(ev) > 1e-8), 1L)
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# DAG LAYER DERIVATIVES
# ─────────────────────────────────────────────────────────────────────────────

test_that("DAG derivatives match numerical differentiation in every coordinate", {
  for (nm in names(dag_configs)) {
    cfg <- dag_configs[[nm]]
    par <- c(cfg$gamma, cfg$delta)
    D <- dag_derivatives(par, cfg$R, cfg$edges)
    expect_equal(length(D), length(par))

    for (k in seq_along(par)) {
      num <- num_deriv_matrix(function(p) dag_sigma(p, cfg$R, cfg$edges),
                              par, k)
      expect_lt(rel_gap(D[[k]], num), 1e-6)
      expect_lt(max(abs(D[[k]] - t(D[[k]]))), 1e-12)
    }
  }
})

test_that("delta derivatives are rank one and positive semidefinite", {
  ## dSigma/ddelta_k = T^-1 E_kk T^-T = (column k of T^-1) squared.
  cfg <- dag_configs$chain_R4
  E <- nrow(cfg$edges)
  D <- dag_derivatives(c(cfg$gamma, cfg$delta), cfg$R, cfg$edges)
  for (k in seq_len(cfg$R)) {
    Dk <- D[[E + k]]
    ev <- eigen(Dk, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(ev), -1e-10)
    expect_equal(sum(abs(ev) > 1e-8), 1L)
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# ESTIMATING FUNCTION
# ─────────────────────────────────────────────────────────────────────────────
# psi is the gradient of the Gaussian working criterion with respect to the
# layer parameters. The check uses Shat DIFFERENT from Sigma, so that the
# gradient is not trivially zero.

test_that("psi is the gradient of the working criterion in the GGM layer", {
  set.seed(20260803)
  for (nm in names(ggm_configs)) {
    cfg <- ggm_configs[[nm]]
    bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
    n <- 250

    ## Shat perturbed around Sigma(eta), preserving symmetry and definiteness
    Sig0 <- ggm_sigma(cfg$eta, bases)
    P <- matrix(rnorm(cfg$R^2, sd = 0.03), cfg$R, cfg$R)
    Shat <- Sig0 + 0.5 * (P + t(P))
    Shat <- Shat + diag(0.05, cfg$R)
    expect_gt(min(eigen(Shat, symmetric = TRUE, only.values = TRUE)$values), 0)

    crit <- function(p) working_ploglik(ggm_sigma(p, bases), Shat, n)
    grad_num <- num_gradient(crit, cfg$eta)

    pj <- mcglmgraph:::cpp_pearson_between(
      Sig0, Shat, n, mc_precision_derivatives(Sig0, bases)
    )
    expect_lt(rel_gap(as.numeric(pj$psi), grad_num), 1e-6)
  }
})

test_that("psi is the gradient of the working criterion in the DAG layer", {
  set.seed(20260804)
  for (nm in names(dag_configs)) {
    cfg <- dag_configs[[nm]]
    par <- c(cfg$gamma, cfg$delta)
    n <- 300

    Sig0 <- dag_sigma(par, cfg$R, cfg$edges)
    P <- matrix(rnorm(cfg$R^2, sd = 0.03), cfg$R, cfg$R)
    Shat <- Sig0 + 0.5 * (P + t(P)) + diag(0.05, cfg$R)
    expect_gt(min(eigen(Shat, symmetric = TRUE, only.values = TRUE)$values), 0)

    crit <- function(p) working_ploglik(dag_sigma(p, cfg$R, cfg$edges),
                                        Shat, n)
    grad_num <- num_gradient(crit, par)

    pj <- mcglmgraph:::cpp_pearson_between(
      Sig0, Shat, n, dag_derivatives(par, cfg$R, cfg$edges)
    )
    expect_lt(rel_gap(as.numeric(pj$psi), grad_num), 1e-6)
  }
})

test_that("psi vanishes when Shat equals the fitted Sigma", {
  cfg <- ggm_configs$chain_R4
  bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
  Sig0 <- ggm_sigma(cfg$eta, bases)
  pj <- mcglmgraph:::cpp_pearson_between(
    Sig0, Sig0, 100, mc_precision_derivatives(Sig0, bases)
  )
  expect_lt(max(abs(as.numeric(pj$psi))), 1e-8)
})

# ─────────────────────────────────────────────────────────────────────────────
# SENSITIVITY
# ─────────────────────────────────────────────────────────────────────────────
# J is the expectation of minus the derivative of psi. The identity with the
# observed version holds at the point where Shat = Sigma(theta), which is where
# it is evaluated: there the term depending on Shat - Sigma vanishes.

test_that("J equals minus the Jacobian of psi at the root, GGM layer", {
  for (nm in names(ggm_configs)) {
    cfg <- ggm_configs[[nm]]
    bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
    n <- 200
    Shat <- ggm_sigma(cfg$eta, bases)   # evaluation at the root

    psi_of <- function(p) {
      Sig <- ggm_sigma(p, bases)
      as.numeric(mcglmgraph:::cpp_pearson_between(
        Sig, Shat, n, mc_precision_derivatives(Sig, bases))$psi)
    }
    Jnum <- -num_jacobian(psi_of, cfg$eta)
    Jana <- mcglmgraph:::cpp_pearson_between(
      Shat, Shat, n, mc_precision_derivatives(Shat, bases))$J

    expect_lt(rel_gap(Jana, Jnum), 1e-5)
    expect_lt(max(abs(Jana - t(Jana))), 1e-10)
    expect_gt(min(eigen(Jana, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
})

test_that("J equals minus the Jacobian of psi at the root, DAG layer", {
  for (nm in names(dag_configs)) {
    cfg <- dag_configs[[nm]]
    par <- c(cfg$gamma, cfg$delta)
    n <- 200
    Shat <- dag_sigma(par, cfg$R, cfg$edges)

    psi_of <- function(p) {
      Sig <- dag_sigma(p, cfg$R, cfg$edges)
      as.numeric(mcglmgraph:::cpp_pearson_between(
        Sig, Shat, n, dag_derivatives(p, cfg$R, cfg$edges))$psi)
    }
    Jnum <- -num_jacobian(psi_of, par)
    Jana <- mcglmgraph:::cpp_pearson_between(
      Shat, Shat, n, dag_derivatives(par, cfg$R, cfg$edges))$J

    expect_lt(rel_gap(Jana, Jnum), 1e-5)
    expect_lt(max(abs(Jana - t(Jana))), 1e-10)
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# SATURATED GGM-DAG EQUIVALENCE
# ─────────────────────────────────────────────────────────────────────────────

test_that("saturated GGM and DAG recover the same residual covariance", {
  skip_if_not_installed("MASS")
  set.seed(20260805)
  n <- 400
  K <- diag(4) * 2
  K[1, 2] <- K[2, 1] <- -0.6; K[2, 3] <- K[3, 2] <- -0.5
  K[3, 4] <- K[4, 3] <- -0.4
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                  x = rnorm(n))
  lp <- list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x)

  fit_ggm <- mcglm_graph(lp, between_graph = mc_between_graph("ggm"), data = d)
  fit_dag <- mcglm_graph(lp, between_graph = mc_between_graph("dag"), data = d)

  expect_true(fit_ggm$graph$converged)
  expect_true(fit_dag$graph$converged)
  ## Both sweep the whole cone: the saturated solution is the same and equals
  ## the sample covariance of the residuals.
  expect_lt(max(abs(fit_ggm$graph$Sigma_b - fit_dag$graph$Sigma_b)), 1e-6)
  expect_lt(max(abs(fit_ggm$graph$Sigma_b - fit_ggm$graph$Shat)), 1e-6)
  expect_lt(max(abs(fit_dag$graph$Sigma_b - fit_dag$graph$Shat)), 1e-6)
})
