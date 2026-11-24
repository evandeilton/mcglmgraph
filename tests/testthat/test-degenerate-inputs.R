# Degenerate inputs of the graph layer: an empty DAG support, collinear parents,
# missing responses and an edge whose estimated variance is not positive.

degenerate_data <- function(n = 200, seed = 5) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y1 <- stats::rnorm(n) + 0.3 * x
  y2 <- 0.5 * y1 + stats::rnorm(n)
  y3 <- 0.6 * y1 + 0.4 * y2 + stats::rnorm(n)
  data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
}
degenerate_lp <- list(y1 ~ x, y2 ~ x, y3 ~ x)
fit_degenerate <- function(d, spec, ...)
  suppressMessages(mcglm_graph(degenerate_lp, between_graph = spec, data = d,
                               ...))

test_that("a DAG layer with no edge fits each node by its own variance", {
  d <- degenerate_data()
  for (solver in c("closed", "newton")) {
    f <- fit_degenerate(d, mc_between_graph("dag", matrix(0L, 0, 2)),
                        control = list(dag_solver = solver))
    g <- f$graph
    expect_true(g$converged, info = solver)
    expect_true(g$inference_ok, info = solver)
    expect_identical(names(g$par), paste0("delta_", 1:3), info = solver)
    expect_equal(unname(g$par), unname(diag(g$Shat)), tolerance = 1e-8,
                 info = solver)
    expect_equal(unname(g$Sigma_b), diag(unname(diag(g$Shat))),
                 tolerance = 1e-8, info = solver)
    expect_length(g$edge_pos, 0L)
  }
  expect_output(print(f), "edges: 0")
  expect_output(print(summary(f)), "Graph layer")
  for (ty in c("model", "sandwich", "two_stage")) {
    V <- mc_vcov_graph(f, ty)
    expect_identical(dim(V), c(3L, 3L), info = ty)
    expect_true(all(is.finite(V)), info = ty)
  }
  cf <- mc_dag_coefficients(f)
  expect_identical(nrow(cf$edges), 0L)
  expect_identical(nrow(cf$nodes), 3L)
  expect_true(all(is.finite(mc_compute_criterion(f, c("aic", "bic", "qic",
                                                       "gbic")))))
  expect_error(mc_test_independence(f), "no edge parameters")
})

test_that("collinear parents are reported, not certified", {
  d <- degenerate_data()
  spec <- mc_between_graph("dag", rbind(c(1, 3), c(2, 3)))
  set.seed(8)
  cases <- list(identical = d$y1, near = d$y1 + 1e-9 * stats::rnorm(nrow(d)))
  for (nm in names(cases)) {
    dd <- d
    dd$y2 <- cases[[nm]]
    out <- capture.output(
      expect_warning(f <- fit_degenerate(dd, spec), "did not converge"),
      type = "output")
    # no raw console warning from the linear algebra
    expect_false(any(grepl("warning", out, ignore.case = TRUE)), info = nm)
    g <- f$graph
    expect_match(g$status, "^parent_block_ill_conditioned: y3", info = nm)
    expect_false(g$converged, info = nm)
    expect_false(g$inference_ok, info = nm)
    expect_identical(g$linear_solve_status, "ridge_regularized", info = nm)
    expect_error(mc_vcov_graph(f, "sandwich"), "not certified", info = nm)
    root <- mc_dag_closed_root(g$Shat, 3, rbind(c(1, 3), c(2, 3)))
    expect_identical(root$status$ill_conditioned, 3L, info = nm)
    expect_gt(root$status$n_ridge, 0L)
  }
  # a well-conditioned parent block is certified
  f_ok <- fit_degenerate(d, spec)
  expect_identical(f_ok$graph$status, "closed_form")
  expect_true(f_ok$graph$inference_ok)
  root <- mc_dag_closed_root(f_ok$graph$Shat, 3, rbind(c(1, 3), c(2, 3)))
  expect_identical(root$status$ill_conditioned, integer(0))
})

test_that("the GES local score refuses a near-singular parent block", {
  d <- degenerate_data()
  S <- stats::cov(cbind(d$y1, d$y1 + 1e-9 * stats::rnorm(nrow(d)), d$y3))
  ls <- getFromNamespace("cpp_ges_local_score", "mcglmgraph")
  expect_identical(ls(S, nrow(d), 3L, c(1L, 2L), log(nrow(d))), -Inf)
  expect_true(is.finite(ls(S, nrow(d), 3L, 1L, log(nrow(d)))))
})

test_that("missing responses and covariates are refused with a clear error", {
  d <- degenerate_data()
  d1 <- d; d1$y2[c(3, 7)] <- NA
  expect_error(fit_degenerate(d1, mc_between_graph("ggm")),
               "response y2 has 2 rows with a missing or non-finite value")
  d2 <- d; d2$y3[5] <- Inf
  expect_error(fit_degenerate(d2, mc_between_graph("dag")),
               "response y3 has 1 row with a missing or non-finite value")
  d3 <- d; d3$x[1] <- NA
  expect_error(fit_degenerate(d3, mc_between_graph("ggm")),
               "covariates of y1 have 1 row with a missing value")
  off <- rep(list(rep(0, nrow(d))), 3)
  off[[2]][4] <- NA
  expect_error(fit_degenerate(d, mc_between_graph("ggm"), offset = off),
               "offset of y2 has 1 missing or non-finite value")
})

test_that("an edge with non-positive variance has an undefined Wald statistic", {
  wald <- getFromNamespace("cpp_wald_edges", "mcglmgraph")
  eta <- c(2, 2, 2, -0.5, 0.3, -0.2)
  V <- diag(c(1, 1, 1, 0.04, 0, -0.01))
  w <- wald(eta, V, 3:5)
  expect_equal(w[1], 0.25 / 0.04)
  expect_true(is.na(w[2]) && is.na(w[3]))
  V[4, 4] <- NaN
  expect_true(is.na(wald(eta, V, 3L)))

  # Propagated through the fit: the statistic and p-value are NA, and the
  # backward search does not take the edge as the weakest one.
  skip_if_not_installed("MASS")
  set.seed(1)
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
  Y <- MASS::mvrnorm(300, rep(0, 3), solve(K))
  dd <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = stats::rnorm(300))
  f <- fit_degenerate(dd, mc_between_graph("ggm"))
  s0 <- mc_edge_scores(f, "wald")
  weakest <- which.min(s0$stat)
  nm <- paste0("eta_", s0$i[weakest], s0$j[weakest])
  p <- match(nm, names(f$graph$par))
  f$graph$vcov[p, ] <- 0
  f$graph$vcov[, p] <- 0
  s <- mc_edge_scores(f, "wald")
  expect_true(is.na(s$stat[weakest]) && is.na(s$p_value[weakest]))
  expect_true(all(is.finite(s$stat[-weakest])))
  # The stepwise search refits its start, so the undefined statistic is
  # injected in the kernel: the weakest edge gets NA, and with k_confirm = 1
  # the search must take the next weakest one and say why.
  real_wald <- wald
  local_mocked_bindings(
    cpp_wald_edges = function(eta, Vcov, edge_pos) {
      w <- real_wald(eta, Vcov, edge_pos)
      w[which.min(w)] <- NA_real_
      w
    }, .package = "mcglmgraph")
  expect_warning(
    sel <- mc_graph_stepwise(f, direction = "backward", criterion = "bic",
                             k_confirm = 1L, max_steps = 1L),
    "undefined")
  ranked <- s0[order(s0$stat), ]
  expect_identical(nrow(sel$trace), 1L)
  expect_identical(c(sel$trace$edge_i, sel$trace$edge_j),
                   c(ranked$i[2], ranked$j[2]))
})

test_that("the closed DAG root is equivariant to a rescaling of the responses", {
  set.seed(21)
  n <- 400
  Z <- matrix(stats::rnorm(n * 4), n, 4)
  Y <- Z %*% chol(matrix(c(1, .4, .3, .2, .4, 1, .5, .3, .3, .5, 1, .4,
                           .2, .3, .4, 1), 4))
  S <- stats::cov(Y)
  edges <- rbind(c(1, 3), c(2, 3), c(1, 4), c(3, 4))
  r0 <- mc_dag_closed_root(S, 4, edges)
  expect_identical(r0$status$ill_conditioned, integer(0))
  scales <- list(c(1e-6, 1, 1e3, 1e6), c(1e6, 1e-6, 1, 1e-3),
                 c(1, 1, 1, 1e6))
  for (a in scales) {
    Sa <- S * outer(a, a)
    ra <- mc_dag_closed_root(Sa, 4, edges)
    expect_identical(ra$status$ill_conditioned, integer(0))
    expect_identical(ra$status$n_ridge, 0L)
    pa <- ra$edges[, "parent"]; ch <- ra$edges[, "child"]
    expect_equal(ra$gamma, r0$gamma * a[ch] / a[pa], tolerance = 1e-9)
    expect_equal(ra$delta, r0$delta * a^2, tolerance = 1e-9)
  }
})

test_that("GES is invariant to a rescaling of the residual covariance", {
  skip_if_not_installed("MASS")
  set.seed(22)
  n <- 500
  x <- stats::rnorm(n)
  z1 <- stats::rnorm(n); z2 <- stats::rnorm(n)
  z3 <- 0.8 * z1 + 0.7 * z2 + stats::rnorm(n); z4 <- 0.6 * z3 + stats::rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, y4 = z4, x = x)
  f <- suppressMessages(mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x),
                                    data = d))
  g0 <- mc_ges(f)
  a <- c(1e-6, 1, 1e3, 1e6)
  f$graph$Shat <- f$graph$Shat * outer(a, a)
  ga <- mc_ges(f)
  expect_identical(ga$cpdag, g0$cpdag)
  ls <- getFromNamespace("cpp_ges_local_score", "mcglmgraph")
  S <- f$graph$Shat
  expect_true(is.finite(ls(S, n, 3L, c(1L, 2L), log(n))))
})

test_that("collinear parents are flagged after one parent is rescaled", {
  d <- degenerate_data()
  set.seed(23)
  ls <- getFromNamespace("cpp_ges_local_score", "mcglmgraph")
  a <- c(1e6, 1, 1)
  for (y2 in list(d$y1, d$y1 + 1e-9 * stats::rnorm(nrow(d)))) {
    S <- stats::cov(cbind(d$y1, y2, d$y3)) * outer(a, a)
    root <- mc_dag_closed_root(S, 3, rbind(c(1, 3), c(2, 3)))
    expect_identical(root$status$ill_conditioned, 3L)
    expect_identical(ls(S, nrow(d), 3L, c(1L, 2L), log(nrow(d))), -Inf)
  }
})
