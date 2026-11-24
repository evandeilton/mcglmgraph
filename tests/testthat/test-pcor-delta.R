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

test_that("pairs outside the support are structural zeros without inference", {
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
  expect_true(all(res$pcor[outside] == 0))
  for (col in c("se", "statistic", "p_value", "lower", "upper"))
    expect_true(all(is.na(res[[col]][outside])), info = col)
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

  ## Node 4 is isolated even after moralization: a structural zero with no
  ## inference.
  expect_false(res$present[outside_moral])
  expect_identical(res$pcor[outside_moral], 0)
  expect_true(is.na(res$se[outside_moral]))
  expect_true(is.na(res$p_value[outside_moral]))
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


# ─────────────────────────────────────────────────────────────────────────────
# UNDEFINED PAIRS AND NEGATIVE VARIANCES IN THE KERNEL
# ─────────────────────────────────────────────────────────────────────────────

kernel_fixture <- function() {
  bases <- mc_precision_bases(3L)
  eta <- c(2, 2, 2, -.5, -.3, .2)
  K <- mc_build_precision(eta, bases)
  Sdot <- mc_precision_derivatives(mc_precision_to_covariance(K)$Sigma_b, bases)
  pairs <- t(utils::combn(3L, 2L)); storage.mode(pairs) <- "integer"
  list(K = K, Sdot = Sdot, pairs = pairs, V = diag(length(eta)), eta = eta)
}

test_that("an undefined pair gives NA and leaves the other pairs unchanged", {
  fx <- kernel_fixture()
  ref <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, fx$V)
  for (bad in c(0, -1, Inf, NaN)) {
    K <- fx$K
    K[1, 1] <- bad                   # pairs (1,2) and (1,3) are undefined
    out <- mcglmgraph:::cpp_pcor_delta(K, fx$Sdot, fx$pairs, fx$V)
    und <- c(1L, 2L)
    expect_true(all(is.na(out$pcor[und])))
    expect_true(all(is.na(out$variance[und])))
    expect_true(all(is.na(out$gradient[und, ])))
    ## the defined pair equals the same pair computed alone on that K
    alone <- mcglmgraph:::cpp_pcor_delta(K, fx$Sdot, fx$pairs[3L, , drop = FALSE],
                                         fx$V)
    expect_identical(out$pcor[3L], as.numeric(alone$pcor))
    expect_identical(out$variance[3L], as.numeric(alone$variance))
    expect_identical(out$gradient[3L, ], alone$gradient[1L, ])
  }
  ## a clean K is untouched
  expect_false(anyNA(ref$pcor)); expect_false(anyNA(ref$variance))
  expect_false(anyNA(ref$gradient))
})

test_that("rounding-size negative variance is 0 and a clear one is NA", {
  fx <- kernel_fixture()
  ref <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, fx$V)
  g <- ref$gradient[1L, ]
  ## rank-one V = a a' with a orthogonal to g: the quadratic form is zero up to
  ## rounding, whatever its sign
  a <- c(g[2L], -g[1L], rep(0, length(g) - 2L))
  o1 <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, tcrossprod(a))
  expect_false(anyNA(o1$variance))
  expect_true(all(o1$variance >= 0))
  expect_lt(o1$variance[1L], 1e-14)
  ## V = -I is clearly not positive semidefinite
  o2 <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs,
                                    -diag(length(g)))
  expect_true(all(is.na(o2$variance)))
  expect_equal(o2$pcor, ref$pcor)
  ## a covariance that is positive semidefinite gives g' V g >= 0
  expect_true(all(ref$variance >= 0))
  ## non-finite V gives NA, not an error
  Vn <- fx$V; Vn[1, 1] <- NA_real_
  o3 <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, Vn)
  expect_true(any(is.na(o3$variance)))
})

# ─────────────────────────────────────────────────────────────────────────────
# STRUCTURAL ZEROS IN THE PUBLIC FUNCTION, AGAINST THE DEFINITION
# ─────────────────────────────────────────────────────────────────────────────

# Reference delta-method table from the definition rho = -K_ij / sqrt(K_ii K_jj)
# with a numerical gradient, independent of the package's own gradient.
pcor_reference <- function(par, precision_of, V, pairs, level = 0.95) {
  t(vapply(seq_len(nrow(pairs)), function(e) {
    i <- pairs[e, 1L]; j <- pairs[e, 2L]
    f <- function(p) {
      K <- precision_of(p); -K[i, j] / sqrt(K[i, i] * K[j, j])
    }
    g <- numDeriv::grad(f, par)
    se <- sqrt(drop(crossprod(g, V %*% g)))
    rho <- f(par)
    c(pcor = rho, se = se, statistic = rho / se,
      p_value = 2 * stats::pnorm(-abs(rho / se)),
      lower = rho - stats::qnorm(1 - (1 - level) / 2) * se,
      upper = rho + stats::qnorm(1 - (1 - level) / 2) * se)
  }, numeric(6)))
}

test_that("structural-zero rows of a sparse GGM carry 0 and NA inference", {
  skip_if_not_installed("MASS"); skip_if_not_installed("numDeriv")
  set.seed(20261124)
  n <- 300
  Y <- MASS::mvrnorm(n, rep(0, 4), diag(4))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                  x = rnorm(n))
  edges <- rbind(c(1, 2), c(3, 4))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x),
                   between_graph = mc_between_graph("ggm", edges = edges),
                   data = d)
  res <- mc_partial_cor_delta(f)
  bases <- mc_precision_bases(4L, edges = edges)
  par <- f$graph$par
  expect_equal(mc_build_precision(par, bases), unname(as.matrix(f$graph$K_b)),
               ignore_attr = TRUE)
  V <- as.matrix(mc_vcov_graph(f, type = "model"))
  pairs <- t(utils::combn(4L, 2L))
  ref <- pcor_reference(par, function(p) mc_build_precision(p, bases), V, pairs)

  abs_row <- !res$present
  expect_equal(sum(abs_row), 4L)
  expect_true(all(res$pcor[abs_row] == 0))
  for (col in c("se", "statistic", "p_value", "lower", "upper"))
    expect_true(all(is.na(res[[col]][abs_row])), info = col)
  for (col in colnames(ref))
    expect_equal(res[[col]][!abs_row], ref[!abs_row, col],
                 tolerance = 1e-5, info = col)

  ## gradient rows: NA for structural zeros, the numerical gradient otherwise
  G <- attr(res, "gradient")
  expect_true(all(is.na(G[abs_row, ])))
  for (e in which(!abs_row)) {
    i <- pairs[e, 1L]; j <- pairs[e, 2L]
    g_num <- numDeriv::grad(function(p) {
      K <- mc_build_precision(p, bases); -K[i, j] / sqrt(K[i, i] * K[j, j])
    }, par)
    expect_equal(G[e, ], g_num, tolerance = 1e-7)
  }

  ## the Fisher interval keeps the same structure
  fi <- mc_partial_cor_delta(f, interval = "fisher")
  expect_true(all(is.na(fi$lower[abs_row])) && all(is.na(fi$upper[abs_row])))
  expect_false(anyNA(fi$lower[!abs_row]))
})

test_that("structural-zero rows of a DAG carry 0 and NA inference", {
  skip_if_not_installed("numDeriv")
  spec <- mc_between_graph("dag", edges = rbind(c(1, 3), c(2, 3)),
                           ordering = 1:4)
  theta <- c(.6, -.5, 1, 1, 1, 1)
  built <- mc_build_sigma_between_graph(theta, 4, spec, inverse = TRUE)
  fit <- structure(
    list(graph = list(type = "dag", R = 4L, spec = spec, par = theta,
                      K_b = built$Kb, Dlist = built$D_Sigmab,
                      vcov = diag(length(theta)), npar = length(theta),
                      response_names = paste0("y", 1:4))),
    class = "mcglm_graph")
  res <- mc_partial_cor_delta(fit)
  pairs <- t(utils::combn(4L, 2L))
  precision_of <- function(p)
    mc_build_sigma_between_graph(p, 4, spec, inverse = TRUE)$Kb
  ref <- pcor_reference(theta, precision_of, as.matrix(mc_vcov_graph(fit)),
                        pairs)
  zero <- !res$present
  expect_equal(sum(zero), 3L)            # pairs (1,4), (2,4), (3,4)
  ## the married pair (1,2) is in the moral graph and keeps its inference
  expect_true(res$present[1L])
  expect_false(anyNA(res[1L, c("pcor", "se", "statistic", "p_value")]))
  G <- attr(res, "gradient")
  expect_true(all(is.na(G[zero, ])))
  expect_false(anyNA(G[!zero, ]))
  expect_true(all(res$pcor[zero] == 0))
  for (col in c("se", "statistic", "p_value", "lower", "upper"))
    expect_true(all(is.na(res[[col]][zero])), info = col)
  for (col in colnames(ref))
    expect_equal(res[[col]][!zero], ref[!zero, col], tolerance = 1e-5,
                 info = col)
})

test_that("the rounding tolerance is relative to the scale of the quadratic form", {
  fx <- kernel_fixture()
  ref <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, fx$V)
  g <- ref$gradient[1L, ]
  a <- c(g[2L], -g[1L], rep(0, length(g) - 2L))      # orthogonal to g
  P <- tcrossprod(a)
  ## V = P - eps I has g' V g = -eps |g|^2; eps sets the ratio to the scale
  ratio_of <- function(eps) {
    V <- P - eps * diag(length(g))
    -drop(g %*% V %*% g) / drop(abs(g) %*% abs(V) %*% abs(g))
  }
  unit <- drop(abs(g) %*% abs(P) %*% abs(g)) / sum(g^2)
  eps_small <- 1e-12 * unit; eps_large <- 1e-8 * unit   # around the 1e-10 limit
  expect_gt(ratio_of(eps_small), 1e-13)
  expect_lt(ratio_of(eps_small), 1e-11)
  expect_gt(ratio_of(eps_large), 1e-9)
  for (s in c(1e-200, 1, 1e200)) {
    v_small <- mcglmgraph:::cpp_pcor_delta(
      fx$K, fx$Sdot, fx$pairs, s * (P - eps_small * diag(length(g))))$variance
    v_large <- mcglmgraph:::cpp_pcor_delta(
      fx$K, fx$Sdot, fx$pairs, s * (P - eps_large * diag(length(g))))$variance
    expect_identical(v_small[1L], 0, info = s)
    expect_true(is.na(v_large[1L]), info = s)
    ## a negative definite V is never rounding, however small
    v_neg <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs,
                                         -s * fx$V)$variance
    expect_true(all(is.na(v_neg)), info = s)
    ## and a positive definite one keeps its variance, however small
    v_pos <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs,
                                         s * fx$V)$variance
    expect_equal(v_pos, s * ref$variance, tolerance = 1e-12, info = s)
    expect_true(all(v_pos > 0), info = s)
  }
})

test_that("non-finite inputs of the kernel give NA and never a number", {
  fx <- kernel_fixture()
  ref <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, fx$V)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    ## in V: the estimate and the gradient do not depend on V
    Vb <- fx$V; Vb[2, 3] <- Vb[3, 2] <- bad
    o <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, Vb)
    expect_true(all(is.na(o$variance)), info = bad)
    expect_identical(o$pcor, ref$pcor)
    expect_identical(o$gradient, ref$gradient)
    ## in one derivative: the estimate stays, gradient column and variance go
    Sb <- fx$Sdot; Sb[[2L]][1, 1] <- bad
    o <- mcglmgraph:::cpp_pcor_delta(fx$K, Sb, fx$pairs, fx$V)
    expect_identical(o$pcor, ref$pcor)
    expect_true(all(is.na(o$gradient[, 2L])), info = bad)
    expect_true(all(is.na(o$variance)), info = bad)
  }
  ## an off-diagonal NaN of K reaches every output of that pair
  K <- fx$K; K[1, 2] <- K[2, 1] <- NaN
  o <- mcglmgraph:::cpp_pcor_delta(K, fx$Sdot, fx$pairs[1L, , drop = FALSE],
                                   fx$V)
  expect_true(is.na(o$pcor)); expect_true(is.na(o$variance))
})

test_that("the gradient of defined pairs matches numerical differentiation", {
  skip_if_not_installed("numDeriv")
  fx <- kernel_fixture()
  bases <- mc_precision_bases(3L)
  out <- mcglmgraph:::cpp_pcor_delta(fx$K, fx$Sdot, fx$pairs, fx$V)
  for (e in seq_len(nrow(fx$pairs))) {
    i <- fx$pairs[e, 1L]; j <- fx$pairs[e, 2L]
    f <- function(p) {
      K <- mc_build_precision(p, bases); -K[i, j] / sqrt(K[i, i] * K[j, j])
    }
    expect_equal(out$pcor[e], f(fx$eta), tolerance = 1e-14)
    expect_equal(out$gradient[e, ], numDeriv::grad(f, fx$eta), tolerance = 1e-8)
  }
})

test_that("mc_partial_cor_delta carries an undefined pair and a bad covariance as NA", {
  spec <- mc_between_graph("dag", edges = rbind(c(1, 3), c(2, 3)),
                           ordering = 1:4)
  theta <- c(.6, -.5, 1, 1, 1, 1)
  built <- mc_build_sigma_between_graph(theta, 4, spec, inverse = TRUE)
  make_fit <- function(K, V)
    structure(
      list(graph = list(type = "dag", R = 4L, spec = spec, par = theta,
                        K_b = K, Dlist = built$D_Sigmab, vcov = V,
                        npar = length(theta),
                        response_names = paste0("y", 1:4))),
      class = "mcglm_graph")
  inference <- c("se", "statistic", "p_value", "lower", "upper")
  ok <- mc_partial_cor_delta(make_fit(built$Kb, diag(length(theta))))

  ## K_11 = 0: the present pairs (1,2) and (1,3) are undefined, (2,3) is not
  K <- built$Kb; K[1, 1] <- 0
  for (ci in c("wald", "fisher")) {
    res <- expect_silent(
      mc_partial_cor_delta(make_fit(K, diag(length(theta))), interval = ci))
    und <- 1:2
    expect_true(all(res$present[und]))
    expect_true(all(is.na(unlist(res[und, c("pcor", inference)]))), info = ci)
    expect_true(all(is.na(attr(res, "gradient")[und, ])))
    expect_false(anyNA(res[4L, c("pcor", inference)]))        # pair (2,3)
    ## structural zeros stay 0, whatever the diagonal
    expect_identical(res$pcor[!res$present], c(0, 0, 0))
  }

  ## a covariance that is not positive semidefinite: estimate kept, no inference
  for (ci in c("wald", "fisher")) {
    res <- expect_silent(
      mc_partial_cor_delta(make_fit(built$Kb, -diag(length(theta))),
                           interval = ci))
    expect_identical(res$pcor, ok$pcor)
    expect_true(all(is.na(unlist(res[inference]))), info = ci)
  }
})
