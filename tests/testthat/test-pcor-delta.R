# =============================================================================
# DELTA-METHOD VARIANCE OF THE PARTIAL CORRELATIONS
# =============================================================================
# Fidelity to the theory is checked at three levels: the gradient against
# numerical differentiation; the identity Kdot = -K Sdot K against the closed
# form of the statement, where Kdot is the canonical base Z_d itself; and the
# variance against the quadratic form g' V g computed independently.

# ─────────────────────────────────────────────────────────────────────────────
# GRADIENT
# ─────────────────────────────────────────────────────────────────────────────

## Partial correlation of the pair (i, j) as a function of the parameters.
pcor_of_par_ggm <- function(par, bases, i, j) {
  K <- mc_build_precision(par, bases)
  -K[i, j] / sqrt(K[i, i] * K[j, j])
}

test_that("the partial-correlation gradient matches numerical differentiation, GGM", {
  for (nm in names(ggm_configs)) {
    cfg <- ggm_configs[[nm]]
    bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
    eta <- cfg$eta
    K <- mc_build_precision(eta, bases)
    Sigma <- mc_precision_to_covariance(K)$Sigma_b
    Sdot <- mc_precision_derivatives(Sigma, bases)

    pairs <- t(utils::combn(cfg$R, 2L))
    storage.mode(pairs) <- "integer"
    V <- diag(length(eta))          # irrelevant for the gradient
    out <- mcglmgraph:::cpp_pcor_delta(K, Sdot, pairs, V)

    for (e in seq_len(nrow(pairs))) {
      i <- pairs[e, 1L]; j <- pairs[e, 2L]
      gnum <- num_gradient(function(p) pcor_of_par_ggm(p, bases, i, j), eta)
      expect_lt(rel_gap(out$gradient[e, ], gnum), 1e-6)
    }
  }
})

test_that("the gradient also matches numerical differentiation in the DAG layer", {
  for (nm in names(dag_configs)) {
    cfg <- dag_configs[[nm]]
    par <- c(cfg$gamma, cfg$delta)
    E <- nrow(cfg$edges)
    pcor_of_par_dag <- function(p, i, j) {
      K <- mc_cholesky_to_precision(p[seq_len(E)], p[E + seq_len(cfg$R)],
                                    cfg$R, cfg$edges)
      -K[i, j] / sqrt(K[i, i] * K[j, j])
    }
    K <- mc_cholesky_to_precision(cfg$gamma, cfg$delta, cfg$R, cfg$edges)
    Sdot <- dag_derivatives(par, cfg$R, cfg$edges)

    pairs <- t(utils::combn(cfg$R, 2L))
    storage.mode(pairs) <- "integer"
    out <- mcglmgraph:::cpp_pcor_delta(K, Sdot, pairs, diag(length(par)))

    for (e in seq_len(nrow(pairs))) {
      i <- pairs[e, 1L]; j <- pairs[e, 2L]
      gnum <- num_gradient(function(p) pcor_of_par_dag(p, i, j), par)
      expect_lt(rel_gap(out$gradient[e, ], gnum), 1e-6)
    }
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# FIDELITY TO THE CLOSED FORM OF THE STATEMENT
# ─────────────────────────────────────────────────────────────────────────────

test_that("in the GGM layer, -K Sdot K recovers the canonical base exactly", {
  ## This is the step that licenses using the same expression in both layers:
  ## the gradient is written with Z_d, and in the GGM layer the
  ## identity returns Z_d.
  for (nm in names(ggm_configs)) {
    cfg <- ggm_configs[[nm]]
    bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
    K <- mc_build_precision(cfg$eta, bases)
    Sigma <- mc_precision_to_covariance(K)$Sigma_b
    Sdot <- mc_precision_derivatives(Sigma, bases)
    for (d in seq_along(cfg$eta)) {
      Kdot <- -K %*% Sdot[[d]] %*% K
      expect_lt(max(abs(Kdot - bases$bases[[d]])), 1e-8)
    }
  }
})

test_that("the gradient agrees with the closed form written in terms of Z_d", {
  cfg <- ggm_configs$sparse_R5
  bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
  K <- mc_build_precision(cfg$eta, bases)
  Sigma <- mc_precision_to_covariance(K)$Sigma_b
  Sdot <- mc_precision_derivatives(Sigma, bases)
  pairs <- t(utils::combn(cfg$R, 2L)); storage.mode(pairs) <- "integer"
  out <- mcglmgraph:::cpp_pcor_delta(K, Sdot, pairs, diag(length(cfg$eta)))

  ## Literal reimplementation of eq. grad-pcor in R, using the bases Z_d.
  for (e in seq_len(nrow(pairs))) {
    i <- pairs[e, 1L]; j <- pairs[e, 2L]
    rho <- -K[i, j] / sqrt(K[i, i] * K[j, j])
    g_lit <- vapply(seq_along(cfg$eta), function(d) {
      Z <- bases$bases[[d]]
      -Z[i, j] / sqrt(K[i, i] * K[j, j]) -
        (rho / 2) * (Z[i, i] / K[i, i] + Z[j, j] / K[j, j])
    }, numeric(1))
    expect_lt(rel_gap(out$gradient[e, ], g_lit), 1e-10)
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# VARIANCE
# ─────────────────────────────────────────────────────────────────────────────

test_that("the variance is the quadratic form g' V g", {
  cfg <- ggm_configs$chain_R4
  bases <- mc_precision_bases(cfg$R, edges = cfg$edges)
  K <- mc_build_precision(cfg$eta, bases)
  Sigma <- mc_precision_to_covariance(K)$Sigma_b
  Sdot <- mc_precision_derivatives(Sigma, bases)
  set.seed(20260806)
  A <- matrix(rnorm(length(cfg$eta)^2), length(cfg$eta))
  V <- crossprod(A) / 100                      # symmetric positive definite
  pairs <- t(utils::combn(cfg$R, 2L)); storage.mode(pairs) <- "integer"
  out <- mcglmgraph:::cpp_pcor_delta(K, Sdot, pairs, V)
  for (e in seq_len(nrow(pairs))) {
    g <- out$gradient[e, ]
    expect_equal(out$variance[e], as.numeric(t(g) %*% V %*% g),
                 tolerance = 1e-10)
  }
})

# ─────────────────────────────────────────────────────────────────────────────
# PUBLIC INTERFACE
# ─────────────────────────────────────────────────────────────────────────────

test_that("mc_partial_cor_delta returns the documented contract", {
  skip_if_not_installed("MASS")
  set.seed(20260807)
  n <- 300
  K <- diag(4) * 2
  K[1, 2] <- K[2, 1] <- -0.7; K[2, 3] <- K[3, 2] <- -0.5
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                  x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = d)

  res <- mc_partial_cor_delta(f)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), choose(4, 2))
  expect_true(all(c("response1", "response2", "present", "pcor", "se",
                    "statistic", "p_value", "lower", "upper") %in% names(res)))
  expect_identical(attr(res, "variability"), "model")
  expect_equal(dim(attr(res, "gradient")), c(choose(4, 2), f$graph$npar))

  ## the point estimates agree with the matrix from mc_partial_cor
  P <- mc_partial_cor(f)
  idx <- t(utils::combn(4, 2))
  expect_equal(res$pcor, P[idx], tolerance = 1e-12, ignore_attr = TRUE)

  ## the interval is symmetric around the estimate and contains it
  expect_true(all(res$lower <= res$pcor & res$pcor <= res$upper))
  expect_equal(res$upper - res$pcor, res$pcor - res$lower, tolerance = 1e-12)

  ## the empirical route returns positive, finite standard errors
  rob <- mc_partial_cor_delta(f, variability = "empirical")
  expect_identical(attr(rob, "variability"), "empirical")
  expect_true(all(is.finite(rob$se)))
  expect_true(all(rob$se >= 0))
})

test_that("pairs outside the support are structural zeros with zero standard error", {
  skip_if_not_installed("MASS")
  set.seed(20260808)
  n <- 300
  Y <- MASS::mvrnorm(n, rep(0, 4), diag(4))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                  x = rnorm(n))
  sparse <- mc_between_graph("ggm", edges = rbind(c(1, 2), c(3, 4)))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), between_graph = sparse,
                   data = d)
  res <- mc_partial_cor_delta(f)

  outside <- !res$present
  expect_true(any(outside))
  expect_true(all(abs(res$pcor[outside]) < 1e-12))
  expect_true(all(res$se[outside] < 1e-12))
  expect_true(all(is.na(res$statistic[outside])))
  ## on the edges of the support the standard error is strictly positive
  expect_true(all(res$se[res$present] > 0))
})

test_that("DAG precision support follows moralization rather than arrow adjacency", {
  spec <- mc_between_graph("dag", edges = rbind(c(1, 3), c(2, 3)),
                           ordering = 1:4)
  theta <- c(.6, -.5, 1, 1, 1, 1)
  built <- mc_build_sigma_between_graph(theta, 4, spec, inverse = TRUE)
  fit <- structure(
    list(graph = list(
      type = "dag", R = 4L, spec = spec, par = theta,
      K_b = built$Kb, Dlist = built$D_Sigmab,
      vcov = diag(length(theta)), npar = length(theta),
      response_names = paste0("y", 1:4)
    )),
    class = "mcglm_graph"
  )

  res <- mc_partial_cor_delta(fit)
  co_parents <- res$response1 == "y1" & res$response2 == "y2"
  outside_moral <- res$response1 == "y1" & res$response2 == "y4"

  ## 1 and 2 have no arrow between them, but are joined in the moral graph
  ## because both are parents of 3; their precision entry is structurally free.
  expect_true(res$present[co_parents])
  expect_gt(abs(res$pcor[co_parents]), 1e-8)
  expect_gt(res$se[co_parents], 0)

  ## Node 4 is isolated even after moralization, hence this remains an exact
  ## precision zero with zero gradient and standard error.
  expect_false(res$present[outside_moral])
  expect_equal(res$pcor[outside_moral], 0, tolerance = 1e-12)
  expect_equal(res$se[outside_moral], 0, tolerance = 1e-12)
})

test_that("mc_partial_cor_delta validates its arguments", {
  skip_if_not_installed("MASS")
  set.seed(20260809)
  Y <- MASS::mvrnorm(150, rep(0, 3), diag(3))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(150))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  expect_error(mc_partial_cor_delta(f, level = 1), "level")
  expect_error(mc_partial_cor_delta(f, level = c(0.9, 0.95)), "level")
  expect_error(mc_partial_cor_delta("not a fit"))
})

test_that("the clustering arguments reach the variability of the delta method", {
  skip_if_not_installed("MASS")
  set.seed(29); n <- 300
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  cl <- rep(seq_len(n / 2), each = 2L)

  plain <- mc_partial_cor_delta(f, variability = "empirical")
  clust <- mc_partial_cor_delta(f, variability = "empirical", cluster = cl)
  # the point estimates cannot move; only the standard errors may
  expect_equal(plain$pcor, clust$pcor, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(plain$se, clust$se)))

  # one cluster per row is the unclustered form
  solo <- mc_partial_cor_delta(f, variability = "empirical",
                               cluster = seq_len(n))
  expect_equal(plain$se, solo$se, tolerance = 1e-10)

  # the small-sample factor scales the standard errors by sqrt(G/(G-1))
  G <- n / 2
  ss <- mc_partial_cor_delta(f, variability = "empirical", cluster = cl,
                             small_sample = TRUE)
  expect_equal(ss$se, clust$se * sqrt(G / (G - 1)), tolerance = 1e-9)

  # the model-based route refuses them, as everywhere else
  expect_error(mc_partial_cor_delta(f, cluster = cl), "empirical variability")
})

test_that("the two-stage form leaves the partial correlations where they were", {
  skip_if_not_installed("MASS")
  set.seed(31); n <- 800
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.8
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  emp <- mc_partial_cor_delta(f, variability = "empirical")
  two <- mc_partial_cor_delta(f, variability = "two_stage")
  expect_equal(emp$pcor, two$pcor, tolerance = 1e-12)
  # the gradients annihilate the cross sensitivity, so the standard errors
  # agree to the order the proposition promises
  expect_equal(emp$se, two$se, tolerance = 0.01 * max(emp$se))
  expect_identical(attr(two, "variability"), "two_stage")
})

## The three forms on the partial correlations. The two empirical ones agree to
## first order under conditions (a) to (c); the model-based one agrees with
## them only when V = J. The quadratic forms are rebuilt here from the
## definitions of J and of the unit scores, without the package's kernels.
pcor_forms_theory <- function(f) {
  g <- f$graph
  Sinv <- solve(g$Sigma_b)
  P <- length(g$Dlist)
  J <- matrix(0, P, P)
  for (a in seq_len(P)) for (b in seq_len(P))
    J[a, b] <- 0.5 * g$n * sum(diag(Sinv %*% g$Dlist[[a]] %*% Sinv %*%
                                      g$Dlist[[b]]))
  Ec <- scale(g$residuals, center = TRUE, scale = FALSE)
  Psi <- sapply(seq_len(P), function(d) {
    A <- Sinv %*% g$Dlist[[d]] %*% Sinv
    -0.5 * (sum(diag(Sinv %*% g$Dlist[[d]])) - rowSums((Ec %*% A) * Ec))
  })
  G <- attr(mc_partial_cor_delta(f), "gradient")
  Jinv <- solve(J)
  list(model = sqrt(diag(G %*% Jinv %*% t(G))),
       empirical = sqrt(diag(G %*% Jinv %*% crossprod(Psi) %*% Jinv %*% t(G))))
}

test_that("with Gaussian residuals the three forms agree on the partial correlations", {
  skip_if_not_installed("MASS")
  set.seed(47); n <- 4000
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7; K[2, 3] <- K[3, 2] <- -0.5
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  th <- pcor_forms_theory(f)
  m <- mc_partial_cor_delta(f, "model")
  e <- mc_partial_cor_delta(f, "empirical")
  t2 <- mc_partial_cor_delta(f, "two_stage")
  expect_equal(m$se, unname(th$model), tolerance = 1e-8)
  expect_equal(e$se, unname(th$empirical), tolerance = 1e-8)
  # V = J holds for Gaussian residuals, so all three are the same to sampling
  # error, and the two empirical forms are the same exactly
  expect_true(all(abs(m$se / e$se - 1) < 0.06))
  expect_equal(t2$se, e$se, tolerance = 1e-8)
})

test_that("with gamma residuals the model-based form departs and the empirical ones do not", {
  skip_if_not_installed("MASS")
  # Gamma margins with a fixed shape under a Gaussian copula: a scale family,
  # so condition (c) holds, while the fourth moments are far from Gaussian and
  # V differs from J.
  set.seed(53); n <- 3000; x <- runif(n)
  Rl <- matrix(c(1, .6, .3, .6, 1, .5, .3, .5, 1), 3)
  U <- stats::pnorm(MASS::mvrnorm(n, rep(0, 3), Rl))
  shape <- 1.5
  mu <- cbind(exp(1 + .5 * x), exp(.5 - .3 * x), rep(exp(.8), n))
  Y <- sapply(1:3, function(r) stats::qgamma(U[, r], shape = shape,
                                             rate = shape / mu[, r]))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = x)
  f <- suppressMessages(mcglm_graph(
    list(y1 ~ x, y2 ~ x, y3 ~ x), data = d, link = rep("log", 3),
    variance = rep("tweedie", 3), power_fixed = rep(TRUE, 3),
    control_initial = list(regression = list(c(1, .5), c(.5, -.3), c(.8, 0)),
                           power = list(2, 2, 2),
                           tau = list(1 / shape, 1 / shape, 1 / shape),
                           rho = c(.5, .3, .4))))
  expect_true(mc_two_stage_available(f)$ok)
  th <- pcor_forms_theory(f)
  m <- mc_partial_cor_delta(f, "model")
  e <- mc_partial_cor_delta(f, "empirical")
  t2 <- mc_partial_cor_delta(f, "two_stage")
  expect_equal(m$se, unname(th$model), tolerance = 1e-8)
  expect_equal(e$se, unname(th$empirical), tolerance = 1e-8)
  # the model-based form is off by far more than sampling error ...
  expect_gt(max(abs(m$se / e$se - 1)), 0.08)
  # ... and the two empirical forms still coincide, as the gradients of the
  # partial correlations annihilate the cross sensitivity
  expect_equal(t2$se, e$se, tolerance = 1e-8)
})
