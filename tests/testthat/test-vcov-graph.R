# Exported accessor for the covariance of the graph-layer parameters.
#
# Every check below rebuilds the target from its definition,
# in plain R and without reusing the package's own kernels, and compares:
#   sensitivity      J_ab = (n/2) tr(S^-1 D_a S^-1 D_b)        (eq. sens-graph)
#   unit score       psi_di = -(1/2)[tr(S^-1 D_d)
#                                    - r_i' S^-1 D_d S^-1 r_i] (eq. pearson-entre)
#   model-based      Var = J_n^-1
#   sandwich         Var = J_n^-1 (n V) J_n^-1, V = n^-1 sum_i psi_i psi_i'

## Strip every attribute but the dimensions: the accessor deliberately carries
## provenance attributes, and they must not take part in a numeric comparison.
bare <- function(M) {
  M <- as.matrix(M)
  attributes(M) <- list(dim = dim(M))
  M
}

make_fit <- function(n = 400, seed = 11) {
  set.seed(seed)
  K <- diag(3) * 2
  K[1, 2] <- K[2, 1] <- -0.7
  K[2, 3] <- K[3, 2] <- -0.5
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

## Sensitivity from its definition, one entry at a time.
theory_J <- function(Sigma, Dlist, n) {
  Sinv <- solve(Sigma)
  P <- length(Dlist)
  J <- matrix(0, P, P)
  for (a in seq_len(P)) for (b in seq_len(P))
    J[a, b] <- 0.5 * n * sum(diag(Sinv %*% Dlist[[a]] %*% Sinv %*% Dlist[[b]]))
  J
}

## Per-unit Pearson scores from their definition, one unit at a time.
theory_psi <- function(E, Sigma, Dlist) {
  Ec <- scale(E, center = TRUE, scale = FALSE)
  Sinv <- solve(Sigma)
  n <- nrow(Ec); P <- length(Dlist)
  Psi <- matrix(0, n, P)
  for (d in seq_len(P)) {
    A <- Sinv %*% Dlist[[d]] %*% Sinv
    tr <- sum(diag(Sinv %*% Dlist[[d]]))
    for (i in seq_len(n))
      Psi[i, d] <- -0.5 * (tr - drop(crossprod(Ec[i, ], A %*% Ec[i, ])))
  }
  Psi
}

test_that("the model-based form is the inverse of the sensitivity of the theory", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  J_th <- theory_J(g$Sigma_b, g$Dlist, g$n)
  expect_equal(bare(g$J), J_th, tolerance = 1e-10)

  V <- mc_vcov_graph(f, type = "model")
  expect_equal(bare(V), solve(J_th), tolerance = 1e-8)
  expect_identical(attr(V, "type"), "model")
  expect_identical(attr(V, "variability"), "model")
  expect_identical(dimnames(V), list(names(g$par), names(g$par)))
  expect_equal(bare(V), bare(t(V)), tolerance = 1e-12)
})

test_that("the sandwich is J^-1 (n V) J^-1 with V built unit by unit", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  Psi <- theory_psi(g$residuals, g$Sigma_b, g$Dlist)
  n <- g$n
  Vemp <- crossprod(Psi) / n
  Jinv <- solve(theory_J(g$Sigma_b, g$Dlist, n))
  target <- Jinv %*% (n * Vemp) %*% Jinv

  V <- mc_vcov_graph(f, type = "sandwich")
  expect_equal(bare(V), bare(0.5 * (target + t(target))), tolerance = 1e-8)
  expect_identical(attr(V, "variability"), "empirical")
  expect_true(all(diag(V) > 0))
})

test_that("the accumulated scale is the one the estimates live on", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  # J is accumulated (proportional to n), so J^-1 shrinks like 1/n and the
  # standard errors are comparable with the coefficients themselves.
  f2 <- make_fit(n = 1600, seed = 11)
  se1 <- sqrt(diag(mc_vcov_graph(f, "model")))
  se2 <- sqrt(diag(mc_vcov_graph(f2, "model")))
  expect_lt(max(se2 / se1), 0.75)          # halving twice: expect about 0.5
  expect_gt(min(se2 / se1), 0.25)
})

test_that("Gaussian residuals put the sandwich close to the model-based form", {
  skip_if_not_installed("MASS")
  # Under Gaussian fourth moments V = J, so the two forms agree up to Monte
  # Carlo error. This is the Bartlett-type identity the working model assumes.
  f <- make_fit(n = 4000, seed = 3)
  Vm <- mc_vcov_graph(f, "model")
  Vs <- mc_vcov_graph(f, "sandwich")
  ratio <- sqrt(diag(Vs)) / sqrt(diag(Vm))
  expect_true(all(ratio > 0.8 & ratio < 1.25))
})

test_that("one cluster per row reproduces the unclustered sandwich", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  V0 <- mc_vcov_graph(f, "sandwich")
  V1 <- mc_vcov_graph(f, "sandwich", cluster = seq_len(f$graph$n))
  expect_equal(bare(V0), bare(V1), tolerance = 1e-12)
  expect_identical(attr(V1, "n_clusters"), f$graph$n)
})

test_that("clustering aggregates the scores before squaring them", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  cl <- rep(seq_len(g$n / 2), each = 2L)
  Psi <- theory_psi(g$residuals, g$Sigma_b, g$Dlist)
  Pc <- rowsum(Psi, cl, reorder = FALSE)
  Vemp <- crossprod(Pc) / g$n
  Jinv <- solve(theory_J(g$Sigma_b, g$Dlist, g$n))
  target <- Jinv %*% (g$n * Vemp) %*% Jinv
  V <- mc_vcov_graph(f, "sandwich", cluster = cl)
  expect_equal(bare(V), bare(0.5 * (target + t(target))), tolerance = 1e-8)

  # small_sample multiplies the clustered variability by G/(G-1)
  G <- g$n / 2
  Vs <- mc_vcov_graph(f, "sandwich", cluster = cl, small_sample = TRUE)
  expect_equal(bare(Vs), bare(V) * G / (G - 1), tolerance = 1e-10)
  expect_true(attr(Vs, "small_sample"))
})

test_that("the joint Wald test is a quadratic form in this very matrix", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  for (robust in c(FALSE, TRUE)) {
    V <- mc_vcov_graph(f, type = if (robust) "sandwich" else "model")
    ep <- attr(V, "edge_pos")
    theta <- g$par[ep]
    W <- drop(t(theta) %*% solve(V[ep, ep, drop = FALSE]) %*% theta)
    tst <- mc_test_independence(f, robust = robust)
    expect_equal(unname(tst$statistic), W, tolerance = 1e-9)
    expect_equal(unname(tst$parameter), length(ep))
    expect_equal(unname(tst$p.value),
                 pchisq(W, length(ep), lower.tail = FALSE), tolerance = 1e-12)
  }
})

test_that("the DAG layer is served by the same accessor", {
  skip_if_not_installed("MASS")
  set.seed(5); n <- 400
  x <- rnorm(n)
  y1 <- 1 + 0.5 * x + rnorm(n)
  y2 <- 0.6 * y1 + rnorm(n)
  y3 <- 0.4 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag"), data = d)
  g <- f$graph
  V <- mc_vcov_graph(f, "model")
  expect_equal(bare(V), solve(theory_J(g$Sigma_b, g$Dlist, g$n)),
               tolerance = 1e-8)
  expect_equal(nrow(V), length(g$par))
})

test_that("the contracts of the arguments are enforced", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200)
  expect_error(mc_vcov_graph(f, type = "robust"), "arg")
  expect_error(mc_vcov_graph(f, "model", cluster = rep(1:2, 100)),
               "empirical variability only")
  expect_error(mc_vcov_graph(f, "model", small_sample = TRUE),
               "empirical variability only")
  expect_error(mc_vcov_graph(f, "sandwich", small_sample = NA), "TRUE or FALSE")
  expect_error(mc_vcov_graph(f, "sandwich", cluster = 1:3),
               "one identifier per observational row")
  expect_error(mc_vcov_graph(f, "sandwich", cluster = rep(1L, f$graph$n)),
               "at least two independent clusters")
  expect_error(mc_vcov_graph(structure(list(), class = "lm")), "inherits")
})

test_that("an uncertified fit is refused rather than silently served", {
  fake <- structure(list(graph = list(
    inference_ok = FALSE, status = "score_not_converged",
    inversion_status = "ridge_regularized", rcond_J = 1e-18
  )), class = "mcglm_graph")
  expect_error(mc_vcov_graph(fake), "not certified")
  expect_error(mc_vcov_graph(fake, "sandwich"), "not certified")
})

test_that("the Wald test refuses the arguments its variance cannot use", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200)
  expect_error(mc_test_independence(f, robust = FALSE, cluster = rep(1:2, 100)),
               "only for robust = TRUE")
  expect_error(mc_test_independence(f, robust = FALSE, small_sample = TRUE),
               "only for robust = TRUE")
  # and it accepts them when the sandwich is in use
  cl <- rep(seq_len(f$graph$n / 2), each = 2L)
  tst <- mc_test_independence(f, robust = TRUE, cluster = cl,
                              small_sample = TRUE)
  expect_s3_class(tst, "htest")
  expect_true(is.finite(tst$statistic))
})

test_that("the small-sample flag says whether G/(G-1) was applied", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  n <- f$graph$n
  cl <- rep(seq_len(n / 2), each = 2L)
  G <- n / 2
  for (type in c("sandwich", "two_stage")) {
    plain <- mc_vcov_graph(f, type, cluster = cl)
    corr <- mc_vcov_graph(f, type, cluster = cl, small_sample = TRUE)
    # applied exactly once: the ratio is G/(G-1), not its square
    expect_equal(bare(corr), bare(plain) * G / (G - 1), tolerance = 1e-12,
                 info = type)
    expect_true(attr(corr, "small_sample"), info = type)
    expect_false(attr(plain, "small_sample"), info = type)
    # without clusters there is nothing to apply it to, and the flag says so
    solo <- mc_vcov_graph(f, type, small_sample = TRUE)
    expect_equal(bare(solo), bare(mc_vcov_graph(f, type)), tolerance = 0,
                 info = type)
    expect_false(attr(solo, "small_sample"), info = type)
  }
  expect_false(attr(mc_vcov_graph(f, "model"), "small_sample"))
})
