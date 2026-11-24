# Closed-form layer of the directed parametrization: root, sensitivity, per-unit
# scores and the scale-free coefficients.
#
# Targets rebuilt from their definitions:
#   root         gamma_j = S_{P_jP_j}^{-1} S_{P_j,j},  delta_j = S_jj - S_{j,P_j} gamma_j
#   sensitivity  J_{gamma_j gamma_j} = n Sigma_{P_jP_j}/delta_j,  J_{delta_j delta_j} = n/(2 delta_j^2)
#   unit score   psi_{gamma_j,i} = eps_ij r_{i,P_j}/delta_j,  eps_ij = r_ij - gamma_j' r_{i,P_j}
#   scale        gamma_jk -> gamma_jk a_j/a_k,  delta_j -> a_j^2 delta_j
#   invariants   gamma*_jk = gamma_jk sqrt(Sigma_kk/Sigma_jj),  R2 = 1 - delta_j/Sigma_jj

## A DAG whose parents are mutually adjacent (no immorality) and one with a
## v-structure, so that the two regimes of the fitted variances are separated.
perfect_edges <- rbind(c(1, 2), c(1, 3), c(2, 3))      # 1 -> 2 -> 3, 1 -> 3
vstruct_edges <- rbind(c(1, 3), c(2, 3))               # 1 -> 3 <- 2, 1 and 2 not adjacent

## Least-squares root written out in plain R, node by node.
ref_root <- function(S, edges) {
  R <- nrow(S)
  g <- numeric(nrow(edges)); d <- diag(S)
  for (j in seq_len(R)) {
    ix <- which(edges[, 2] == j)
    if (!length(ix)) next
    P <- edges[ix, 1]
    g[ix] <- solve(S[P, P, drop = FALSE], S[P, j])
    d[j] <- S[j, j] - sum(S[j, P] * g[ix])
  }
  list(gamma = g, delta = unname(d))
}

test_that("the root recovers the coefficients that generated the covariance", {
  gamma <- c(0.7, 0.3, 0.5); delta <- c(1, 0.64, 0.49)
  S <- mc_cholesky_to_covariance(gamma, delta, R = 3, edges = perfect_edges)$Sigma_b
  fit <- mc_dag_closed_root(S, R = 3, edges = perfect_edges)
  expect_equal(fit$gamma, gamma, tolerance = 1e-12)
  expect_equal(fit$delta, delta, tolerance = 1e-12)
  expect_equal(fit$Sigma_b, S, tolerance = 1e-12)
  expect_equal(fit$status$n_ridge, 0L)
  expect_equal(fit$status$n_floored, 0L)
})

test_that("the root is the node-by-node least-squares solution", {
  set.seed(9)
  X <- matrix(rnorm(500 * 4), 500, 4)
  X[, 2] <- 0.6 * X[, 1] + X[, 2]
  X[, 3] <- 0.4 * X[, 1] + 0.3 * X[, 2] + X[, 3]
  X[, 4] <- 0.5 * X[, 3] + X[, 4]
  S <- crossprod(scale(X, scale = FALSE)) / 500
  ed <- rbind(c(1, 2), c(1, 3), c(2, 3), c(3, 4))
  got <- mc_dag_closed_root(S, R = 4, edges = ed)
  want <- ref_root(S, ed)
  expect_equal(got$gamma, want$gamma, tolerance = 1e-12)
  expect_equal(got$delta, want$delta, tolerance = 1e-12)
  # and each block is the coefficient vector of the corresponding regression
  fit3 <- lm(X[, 3] ~ X[, 1] + X[, 2])
  expect_equal(unname(got$gamma[c(2, 3)]), unname(coef(fit3)[-1]),
               tolerance = 1e-10)
})

test_that("the root annihilates the Pearson estimating function", {
  set.seed(11)
  X <- matrix(rnorm(400 * 3), 400, 3)
  X[, 2] <- 0.6 * X[, 1] + X[, 2]
  X[, 3] <- 0.4 * X[, 2] + X[, 3]
  S <- crossprod(scale(X, scale = FALSE)) / 400
  for (ed in list(perfect_edges, vstruct_edges)) {
    fit <- mc_dag_closed_root(S, R = 3, edges = ed)
    der <- mc_cholesky_derivatives(fit$gamma, fit$delta, R = 3, edges = ed)
    pj <- mcglmgraph:::cpp_pearson_between(fit$Sigma_b, S, 400,
                                           c(der$D_gamma, der$D_delta))
    expect_lt(max(abs(pj$psi)) / 400, 1e-12)
  }
})

test_that("Newton and the closed form solve the same system", {
  skip_if_not_installed("MASS")
  set.seed(4); n <- 500; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + 0.3 * y1 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  args <- list(list(y1 ~ x, y2 ~ x, y3 ~ x),
               between_graph = mc_between_graph("dag"), data = d)
  fn <- do.call(mcglm_graph, c(args, list(control = list(dag_solver = "newton"))))
  fc <- do.call(mcglm_graph, args)          # the closed route is the default
  expect_identical(fn$graph$solver, "newton")
  expect_identical(fc$graph$solver, "closed")
  expect_identical(fc$graph$status, "closed_form")
  expect_true(fc$graph$converged)
  expect_equal(unname(fn$graph$par), unname(fc$graph$par), tolerance = 1e-8)
  expect_lt(fc$graph$score_norm_unit, fc$graph$score_tol)
  expect_equal(fn$graph$Sigma_b, fc$graph$Sigma_b, tolerance = 1e-8)
})

test_that("the closed-form sensitivity equals the generic assembly", {
  gamma <- c(0.7, 0.3, 0.5); delta <- c(1, 0.64, 0.49); n <- 250
  for (ed in list(perfect_edges, vstruct_edges)) {
    g <- if (identical(ed, vstruct_edges)) gamma[1:2] else gamma
    S <- mc_cholesky_to_covariance(g, delta, R = 3, edges = ed)$Sigma_b
    der <- mc_cholesky_derivatives(g, delta, R = 3, edges = ed)
    Dl <- c(der$D_gamma, der$D_delta)
    J_gen <- mcglmgraph:::cpp_pearson_between(S, S, n, Dl)$J
    J_cf <- mc_dag_sensitivity(S, delta, n = n, R = 3, edges = ed)
    expect_equal(J_cf, J_gen, tolerance = 1e-10)
  }
})

test_that("the sensitivity has the block structure the theory gives it", {
  # two nodes with parents, disjoint parent sets except for node 3
  ed <- rbind(c(1, 2), c(1, 3), c(2, 3), c(3, 4))
  gamma <- c(0.5, 0.3, 0.4, 0.6); delta <- c(1, 0.8, 0.7, 0.9); n <- 300
  S <- mc_cholesky_to_covariance(gamma, delta, R = 4, edges = ed)$Sigma_b
  J <- mc_dag_sensitivity(S, delta, n = n, R = 4, edges = ed)
  E <- nrow(ed)
  # innovation variances: diagonal, with the closed value
  expect_equal(diag(J)[E + seq_len(4)], n / (2 * delta^2), tolerance = 1e-12)
  expect_equal(J[E + seq_len(4), E + seq_len(4)],
               diag(n / (2 * delta^2)), tolerance = 1e-12)
  # no coupling between coefficients and innovation variances
  expect_equal(max(abs(J[seq_len(E), E + seq_len(4)])), 0)
  # node 3 has parents 1 and 2, at edge positions 2 and 3
  expect_equal(J[c(2, 3), c(2, 3)], (n / delta[3]) * S[c(1, 2), c(1, 2)],
               tolerance = 1e-12, ignore_attr = TRUE)
  # coefficients of distinct nodes do not interact: edge 1 belongs to node 2
  expect_equal(max(abs(J[1, c(2, 3, 4)])), 0)
  expect_equal(max(abs(J[4, c(1, 2, 3)])), 0)
})

test_that("the per-unit scores agree with the generic Pearson contributions", {
  set.seed(13); n <- 300
  X <- matrix(rnorm(n * 3), n, 3)
  X[, 2] <- 0.6 * X[, 1] + X[, 2]
  X[, 3] <- 0.4 * X[, 2] + X[, 3]
  Ec <- scale(X, center = TRUE, scale = FALSE)
  S <- crossprod(Ec) / n
  ed <- perfect_edges
  fit <- mc_dag_closed_root(S, R = 3, edges = ed)
  sc <- mcglmgraph:::cpp_dag_scores(Ec, fit$gamma, fit$delta, ed)

  # (a) against the definition eps_ij r_{i,P_j} / delta_j
  eps3 <- Ec[, 3] - Ec[, c(1, 2)] %*% fit$gamma[c(2, 3)]
  expect_equal(sc$innovations[, 3], as.numeric(eps3), tolerance = 1e-12)
  expect_equal(sc$scores[, 2], Ec[, 1] * as.numeric(eps3) / fit$delta[3],
               tolerance = 1e-12)
  expect_equal(sc$scores[, nrow(ed) + 3],
               (as.numeric(eps3)^2 - fit$delta[3]) / (2 * fit$delta[3]^2),
               tolerance = 1e-12)
  # (b) against the generic trace form used by the sandwich
  der <- mc_cholesky_derivatives(fit$gamma, fit$delta, R = 3, edges = ed)
  Dl <- c(der$D_gamma, der$D_delta)
  Sinv <- solve(fit$Sigma_b)
  gen <- vapply(Dl, function(Dd) {
    A <- Sinv %*% Dd %*% Sinv
    -0.5 * (sum(diag(Sinv %*% Dd)) - rowSums((Ec %*% A) * Ec))
  }, numeric(n))
  expect_equal(unname(sc$scores), unname(gen), tolerance = 1e-9)
  # (c) the scores sum to zero at the root
  expect_lt(max(abs(colSums(sc$scores))) / n, 1e-12)
})

test_that("a DAG without immorality reproduces the sample covariance", {
  set.seed(17); n <- 400
  X <- matrix(rnorm(n * 3), n, 3)
  X[, 3] <- 0.5 * X[, 1] + 0.5 * X[, 2] + X[, 3]
  Ec <- scale(X, center = TRUE, scale = FALSE)
  S <- crossprod(Ec) / n
  perfect <- mc_dag_closed_root(S, R = 3, edges = perfect_edges)
  expect_equal(perfect$Sigma_b, S, tolerance = 1e-10)   # saturated, chordal
  vst <- mc_dag_closed_root(S, R = 3, edges = vstruct_edges)
  # the two nonadjacent parents keep their model-implied covariance, and the
  # variance of the common child differs from its sample value
  expect_equal(vst$Sigma_b[1, 2], 0, tolerance = 1e-12)
  expect_gt(abs(vst$Sigma_b[3, 3] - S[3, 3]), 1e-6)
  expect_equal(diag(vst$Sigma_b)[1:2], diag(S)[1:2], tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("the node sandwich is the heteroscedasticity-consistent covariance", {
  set.seed(19); n <- 500
  X <- matrix(rnorm(n * 3), n, 3)
  X[, 3] <- 0.5 * X[, 1] + 0.5 * X[, 2] + abs(X[, 1]) * X[, 3]  # heteroscedastic
  Ec <- scale(X, center = TRUE, scale = FALSE)
  S <- crossprod(Ec) / n

  node_sandwich <- function(ed) {
    fit <- mc_dag_closed_root(S, R = 3, edges = ed)
    sc <- mcglmgraph:::cpp_dag_scores(Ec, fit$gamma, fit$delta, ed)
    J <- mc_dag_sensitivity(fit$Sigma_b, fit$delta, n = n, R = 3, edges = ed)
    Jinv <- solve(J)
    V <- Jinv %*% crossprod(sc$scores) %*% Jinv
    ix <- which(ed[, 2] == 3)
    P <- ed[ix, 1]
    u <- Ec[, 3] - Ec[, P, drop = FALSE] %*% fit$gamma[ix]
    meat <- crossprod(Ec[, P, drop = FALSE] * as.numeric(u)) / n
    bread <- solve(S[P, P, drop = FALSE])
    list(pkg = V[ix, ix, drop = FALSE], hc0 = bread %*% meat %*% bread / n)
  }

  # no immorality: the two coincide exactly, because Sigma_{P_jP_j} = S_{P_jP_j}
  perf <- node_sandwich(perfect_edges)
  expect_equal(unname(perf$pkg), unname(perf$hc0), tolerance = 1e-9)
  # with a v-structure the bread is the model-implied block and they differ
  vst <- node_sandwich(vstruct_edges)
  expect_gt(max(abs(vst$pkg - vst$hc0)) / max(abs(vst$hc0)), 1e-3)
})

test_that("the fit is equivariant under a diagonal rescaling", {
  # Proposition on scale: A Sigma A has factors A T A^-1 and A D A
  set.seed(23)
  ed <- rbind(c(1, 2), c(1, 3), c(2, 3), c(2, 4), c(3, 4))
  gamma <- c(0.5, -0.3, 0.4, 0.2, 0.6); delta <- c(1, 0.8, 0.7, 0.9)
  S <- mc_cholesky_to_covariance(gamma, delta, R = 4, edges = ed)$Sigma_b
  a <- c(0.4, 3, 1.3, 0.7)
  Sa <- diag(a) %*% S %*% diag(a)

  f0 <- mc_dag_closed_root(S, R = 4, edges = ed)
  fa <- mc_dag_closed_root(Sa, R = 4, edges = ed)
  # gamma_jk -> gamma_jk a_j / a_k  (j the child, k the parent)
  expect_equal(fa$gamma, f0$gamma * a[ed[, 2]] / a[ed[, 1]], tolerance = 1e-10)
  expect_equal(fa$delta, f0$delta * a^2, tolerance = 1e-10)
  expect_equal(fa$Sigma_b, diag(a) %*% f0$Sigma_b %*% diag(a), tolerance = 1e-10)

  # the standardized quantities are invariant
  gstar <- function(f, S) f$gamma * sqrt(diag(S)[ed[, 1]] / diag(S)[ed[, 2]])
  r2 <- function(f, S) 1 - f$delta / diag(S)
  expect_equal(gstar(fa, fa$Sigma_b), gstar(f0, f0$Sigma_b), tolerance = 1e-10)
  expect_equal(r2(fa, fa$Sigma_b), r2(f0, f0$Sigma_b), tolerance = 1e-10)
  # the support is untouched
  expect_equal(which(abs(fa$gamma) > 1e-12), which(abs(f0$gamma) > 1e-12))
})

test_that("the Wald statistic of an edge is invariant under rescaling", {
  set.seed(29); n <- 400
  ed <- rbind(c(1, 2), c(1, 3), c(2, 3))
  X <- matrix(rnorm(n * 3), n, 3)
  X[, 2] <- 0.6 * X[, 1] + X[, 2]
  X[, 3] <- 0.4 * X[, 2] + 0.2 * X[, 1] + X[, 3]
  Ec <- scale(X, center = TRUE, scale = FALSE)
  a <- c(0.5, 2, 1.4)
  z_of <- function(M) {
    S <- crossprod(M) / n
    fit <- mc_dag_closed_root(S, R = 3, edges = ed)
    J <- mc_dag_sensitivity(fit$Sigma_b, fit$delta, n = n, R = 3, edges = ed)
    sc <- mcglmgraph:::cpp_dag_scores(M, fit$gamma, fit$delta, ed)
    Jinv <- solve(J)
    Vs <- Jinv %*% crossprod(sc$scores) %*% Jinv
    list(model = fit$gamma / sqrt(diag(Jinv)[seq_len(nrow(ed))]),
         sandwich = fit$gamma / sqrt(diag(Vs)[seq_len(nrow(ed))]))
  }
  z0 <- z_of(Ec)
  za <- z_of(Ec %*% diag(a))
  expect_equal(za$model, z0$model, tolerance = 1e-9)
  expect_equal(za$sandwich, z0$sandwich, tolerance = 1e-9)
})

test_that("mc_dag_coefficients reports the invariants and their inference", {
  skip_if_not_installed("MASS")
  set.seed(4); n <- 400; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + 0.3 * y1 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag"), data = d,
                   control = list(dag_solver = "closed"))
  co <- mc_dag_coefficients(f)
  expect_s3_class(co, "mcg_dag_coef")
  g <- f$graph
  E <- nrow(co$edges)
  expect_equal(co$edges$gamma, unname(g$par[seq_len(E)]), tolerance = 1e-12)
  Sd <- diag(g$Sigma_b)
  expect_equal(co$edges$gamma_std,
               co$edges$gamma * sqrt(Sd[co$edges$from] / Sd[co$edges$to]),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(co$nodes$r_squared, 1 - co$nodes$delta / Sd[co$nodes$position],
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(co$edges$se,
               sqrt(diag(mc_vcov_graph(f, "model")))[seq_len(E)],
               tolerance = 1e-12, ignore_attr = TRUE)
  # the complete DAG on three responses has no immorality, so R^2 of a node
  # without parents is zero and the diagonal is the sample one
  expect_equal(co$nodes$r_squared[1], 0, tolerance = 1e-12)
  expect_true(all(co$nodes$r_squared >= 0 & co$nodes$r_squared < 1))
  expect_equal(nrow(mc_dag_coefficients(f, "empirical")$edges), E)
})

test_that("the closed-form API enforces its contracts", {
  S <- diag(3)
  expect_error(mc_dag_closed_root(S[1:2, 1:2], R = 3), "3 by 3")
  bad <- S; bad[1, 2] <- 0.5
  expect_error(mc_dag_closed_root(bad, R = 3), "symmetric")
  expect_error(mc_dag_sensitivity(S, c(1, 1), n = 10, R = 3), "exactly 3")
  expect_error(mc_dag_sensitivity(S, c(1, 1, -1), n = 10, R = 3), "positive")
  expect_error(mc_dag_sensitivity(S, c(1, 1, 1), n = 0, R = 3), "positive")
  skip_if_not_installed("MASS")
  set.seed(1); n <- 150
  d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = rnorm(n))
  fg <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  expect_error(mc_dag_coefficients(fg), "applies to a DAG layer")
  expect_error(mcglm_graph(list(y1 ~ x, y2 ~ x), data = d,
                           control = list(dag_solver = "bfgs")), "arg")
})

test_that("the blockwise inverse of the sensitivity matches the general solve", {
  # supplied in lexicographic order, which is the order the package reads
  ed <- rbind(c(1, 2), c(1, 3), c(1, 5), c(2, 3), c(2, 4), c(3, 4), c(4, 5))
  gamma <- c(0.5, -0.3, 0.1, 0.4, 0.2, 0.6, -0.25)
  delta <- c(1, 0.8, 0.7, 0.9, 0.6); n <- 350
  S <- mc_cholesky_to_covariance(gamma, delta, R = 5, edges = ed)$Sigma_b
  J <- mc_dag_sensitivity(S, delta, n = n, R = 5, edges = ed)
  bi <- mcglmgraph:::cpp_dag_sensitivity_inverse(S, delta, n, ed)
  expect_equal(bi$n_singular, 0L)
  expect_equal(bi$inverse, solve(J), tolerance = 1e-9, ignore_attr = TRUE)
  # the closed entries: delta_j Sigma_{P_jP_j}^-1 / n and 2 delta_j^2 / n
  E <- nrow(ed)
  expect_equal(diag(bi$inverse)[E + seq_len(5)], 2 * delta^2 / n,
               tolerance = 1e-12)
  ix <- which(ed[, 2] == 4)                     # node 4 has parents 2 and 3
  expect_equal(bi$inverse[ix, ix],
               (delta[4] / n) * solve(S[c(2, 3), c(2, 3)]),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the two routes agree on the parameters and on their variance", {
  set.seed(31); R <- 8; n <- 400
  L <- diag(R)
  for (j in 2:R) for (i in 1:(j - 1)) L[j, i] <- 0.25 / (j - i)
  X <- matrix(rnorm(n * R), n, R) %*% t(solve(L))
  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
  spec <- mc_between_graph("dag")
  a <- mcglmgraph:::.fit_between(S, n, spec, dag_solver = "newton")
  b <- mcglmgraph:::.fit_between(S, n, spec, dag_solver = "closed")
  expect_equal(a$par, b$par, tolerance = 1e-10)
  expect_equal(a$J, b$J, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(a$vcov, b$vcov, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(a$Sigma_b, b$Sigma_b, tolerance = 1e-10)
  expect_true(b$inference_ok)
  expect_identical(attr(b$vcov, "inversion_status"), "exact_blockwise")
  expect_identical(b$solver, "closed")
  expect_identical(a$solver, "newton")
})

test_that("a non-identity ordering is honored by the closed route", {
  set.seed(37); R <- 4; n <- 300
  X <- matrix(rnorm(n * R), n, R)
  X[, 2] <- 0.5 * X[, 1] + X[, 2]
  X[, 4] <- 0.4 * X[, 2] + X[, 4]
  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
  ordering <- c(3L, 1L, 2L, 4L)
  spec <- mc_between_graph("dag", ordering = ordering)
  a <- mcglmgraph:::.fit_between(S, n, spec, dag_solver = "newton")
  b <- mcglmgraph:::.fit_between(S, n, spec, dag_solver = "closed")
  expect_equal(a$par, b$par, tolerance = 1e-9)
  expect_equal(a$Sigma_b, b$Sigma_b, tolerance = 1e-10)
  # the fitted covariance is expressed in response order under both routes
  expect_equal(diag(b$Sigma_b), diag(S), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("edges are read in lexicographic order, whatever order they arrive in", {
  shuffled <- rbind(c(2, 3), c(1, 2), c(1, 3))
  sorted <- rbind(c(1, 2), c(1, 3), c(2, 3))
  expect_identical(mcglmgraph:::.validate_edges(shuffled, 3),
                   matrix(as.integer(sorted), ncol = 2L))
  # a coefficient vector therefore follows the sorted order, and the two calls
  # below describe the same model only because gamma is permuted to match
  g_sorted <- c(0.7, 0.3, 0.5)
  g_shuffled <- g_sorted                         # already in sorted order
  A <- mc_build_T_matrix(g_sorted, R = 3, edges = sorted)
  B <- mc_build_T_matrix(g_shuffled, R = 3, edges = shuffled)
  expect_identical(A, B)
  # the closed root reports the order it used, so the caller can align
  S <- mc_cholesky_to_covariance(g_sorted, c(1, 0.8, 0.6), R = 3,
                                 edges = shuffled)$Sigma_b
  fit <- mc_dag_closed_root(S, R = 3, edges = shuffled)
  expect_identical(unname(fit$edges), matrix(as.integer(sorted), ncol = 2L))
  expect_equal(fit$gamma, g_sorted, tolerance = 1e-12)
})

## ---- arcs given directly ---------------------------------------------------

test_that("arbitrary acyclic arcs are fitted through a derived ordering", {
  set.seed(41); R <- 5; n <- 300
  X <- matrix(rnorm(n * R), n, R)
  X[, 2] <- 0.5 * X[, 3] + X[, 2]
  X[, 1] <- 0.4 * X[, 2] + X[, 1]
  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
  # arcs that run against the natural order: 3 -> 2 -> 1
  arcs <- rbind(c(3L, 2L), c(2L, 1L), c(4L, 3L), c(5L, 4L))
  got <- mc_dag_closed_root(S, R = R, edges = arcs, directed = TRUE)
  rank_of <- match(seq_len(R), got$ordering)
  expect_true(all(rank_of[arcs[, 1]] < rank_of[arcs[, 2]]))  # parents first
  expect_identical(unname(got$edges), unname(arcs[order(rank_of[arcs[, 1]],
                                                        rank_of[arcs[, 2]]), ]))

  # the same fit through an explicit ordering and unordered pairs
  same <- mc_dag_closed_root(S, R = R, edges = arcs, directed = TRUE,
                             ordering = got$ordering)
  expect_equal(same$gamma, got$gamma, tolerance = 1e-12)
  expect_equal(same$Sigma_b, got$Sigma_b, tolerance = 1e-12)

  # delta is indexed by response, alongside the diagonal of Sigma_b
  expect_equal(length(got$delta), R)
  expect_true(all(got$delta > 0))
  expect_lt(got$delta[1], S[1, 1])       # node 1 has a parent, so it explains
  expect_equal(got$delta[got$ordering[1]], S[got$ordering[1],
                                             got$ordering[1]],
               tolerance = 1e-12)        # the first node has none
})

test_that("Markov equivalent members give the same fitted covariance", {
  # 1 -> 2 -> 3, 2 -> 1 and 2 -> 3, 3 -> 2 -> 1: same skeleton, no collider
  set.seed(43); n <- 400
  X <- matrix(rnorm(n * 3), n, 3)
  X[, 2] <- 0.6 * X[, 1] + X[, 2]
  X[, 3] <- 0.5 * X[, 2] + X[, 3]
  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n
  members <- list(rbind(c(1L, 2L), c(2L, 3L)),
                  rbind(c(2L, 1L), c(2L, 3L)),
                  rbind(c(2L, 1L), c(3L, 2L)))
  crit <- function(Sg) -0.5 * n * (determinant(Sg)$modulus[1] +
                                     sum(diag(solve(Sg, S))))
  base <- mc_dag_closed_root(S, 3, members[[1]], directed = TRUE)
  for (m in members[-1]) {
    alt <- mc_dag_closed_root(S, 3, m, directed = TRUE)
    expect_equal(alt$Sigma_b, base$Sigma_b, tolerance = 1e-10)
    expect_equal(crit(alt$Sigma_b), crit(base$Sigma_b), tolerance = 1e-10)
    expect_equal(length(alt$gamma), length(base$gamma))
  }
  # the collider 1 -> 2 <- 3 is not in the class and gives a different fit
  coll <- mc_dag_closed_root(S, 3, rbind(c(1L, 2L), c(3L, 2L)), directed = TRUE)
  expect_gt(max(abs(coll$Sigma_b - base$Sigma_b)), 1e-6)
})

test_that("a cyclic or incompatible request is refused", {
  S <- diag(3)
  cyc <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L))
  expect_error(mc_dag_closed_root(S, 3, cyc, directed = TRUE),
               "directed cycle")
  arcs <- rbind(c(1L, 2L), c(2L, 3L))
  expect_error(mc_dag_closed_root(S, 3, arcs, directed = TRUE,
                                  ordering = c(3L, 2L, 1L)),
               "not compatible with the arcs")
  expect_error(mc_dag_closed_root(S, 3, rbind(c(1L, 2L), c(2L, 1L)),
                                  directed = TRUE), "appears twice")
  expect_error(mc_dag_closed_root(S, 3, rbind(c(1L, 1L)), directed = TRUE),
               "self-loops")
  expect_error(mc_dag_closed_root(S, 3, arcs, directed = NA), "TRUE or FALSE")
  expect_error(mc_dag_closed_root(S, 3, ordering = c(1, 1, 2)), "permutation")
})

test_that("the directed layer takes the closed route unless told otherwise", {
  skip_if_not_installed("MASS")
  set.seed(51); n <- 300; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.5 * y1 + rnorm(n); y3 <- 0.3 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag"), data = d)
  expect_identical(f$graph$solver, "closed")
  expect_identical(f$graph$status, "closed_form")
  expect_true(f$graph$converged)
  expect_true(f$graph$inference_ok)
  expect_identical(attr(f$graph$vcov, "inversion_status"), "exact_blockwise")
  # the undirected layer is untouched by the default
  g <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  expect_identical(g$graph$solver, "newton")
  expect_identical(g$graph$status, "converged")
  # and a refit of the graph block inherits the route from the control
  h <- mcglmgraph:::.refit_between(f, mc_between_graph("dag",
                                                       edges = rbind(c(1, 2))))
  expect_identical(h$graph$solver, "closed")
})
