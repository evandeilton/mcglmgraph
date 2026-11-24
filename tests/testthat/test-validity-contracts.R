# Contracts introduced when the public API was tightened: conflict-aware
# orientation with validity gates, the disabled preconditioned proximal step,
# the KKT certificate of the L1 route, the nuisance-adjusted edge score and the
# diagnostics that make a stopping rule auditable.

chain_skeleton <- function() {
  R <- 4L
  adj <- matrix(0L, R, R)
  for (e in list(c(1, 2), c(2, 3), c(3, 4))) {
    adj[e[1], e[2]] <- 1L; adj[e[2], e[1]] <- 1L
  }
  sep <- vector("list", R * R)
  put <- function(i, j, S) {
    sep[[(i - 1L) * R + j]] <<- S; sep[[(j - 1L) * R + i]] <<- S
  }
  # neither center separates its pair, so both triples claim a collider and the
  # two claims require opposite directions on the edge 2-3
  put(1, 3, integer(0)); put(2, 4, integer(0)); put(1, 4, integer(0))
  list(adjacency = adj, sepset = sep)
}

test_that("conflicting colliders orient nothing and are reported", {
  res <- mcglmgraph:::.orient_vstructures(chain_skeleton())
  expect_equal(nrow(res$conflicts), 1L)
  expect_equal(sort(c(res$conflicts$from, res$conflicts$to)), c(2L, 3L))
  expect_true(all(res$vstructures$is_collider))
  expect_false(any(res$vstructures$applied))
  expect_equal(res$amat, t(res$amat))          # everything left undirected
})

test_that("validity gates recognize cycles and consistent extensions", {
  cyc <- matrix(0L, 3, 3); cyc[1, 2] <- 1L; cyc[2, 3] <- 1L; cyc[3, 1] <- 1L
  expect_true(mcglmgraph:::.has_directed_cycle(cyc))
  dag <- matrix(0L, 3, 3); dag[1, 2] <- 1L; dag[2, 3] <- 1L
  expect_false(mcglmgraph:::.has_directed_cycle(dag))
  chain <- matrix(c(0L,1L,0L, 1L,0L,1L, 0L,1L,0L), 3, 3, byrow = TRUE)
  expect_true(mcglmgraph:::.has_consistent_extension(chain))
  coll <- matrix(c(0L,1L,0L, 0L,0L,0L, 0L,1L,0L), 3, 3, byrow = TRUE)
  expect_true(mcglmgraph:::.has_consistent_extension(coll))
})

test_that("the pipeline earns the CPDAG label only through its gates", {
  set.seed(11); n <- 1500
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  p <- mc_graph_pipeline(f)
  expect_true(all(p$gates))
  expect_true(p$valid_cpdag)
  expect_identical(p$label, "CPDAG")
  expect_setequal(names(p$gates),
                  c("no_orientation_conflict", "unregularized_separation",
                    "complete_separation",
                    "orientation_within_skeleton", "skeleton_preserved",
                    "acyclic", "vstructures_preserved",
                    "consistent_extension"))
})

test_that("the preconditioned proximal step is refused", {
  skip_if_not_installed("MASS")
  set.seed(2); n <- 200
  Y <- MASS::mvrnorm(n, rep(0, 3), diag(3))
  d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  expect_error(mc_graph_penalized(f, lambda = 1, precondition = TRUE),
               "disabled")
})

test_that("L1 convergence carries a KKT certificate", {
  skip_if_not_installed("MASS")
  set.seed(9); n <- 400
  K <- diag(4) * 2
  K[1,2] <- K[2,1] <- -0.6; K[2,3] <- K[3,2] <- -0.6; K[3,4] <- K[4,3] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
  d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], y4 = Y[,4], x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = d)
  pen <- mc_graph_penalized(f, lambda = 2)
  expect_true(is.finite(pen$kkt))
  if (isTRUE(pen$converged)) expect_lte(pen$kkt_unit, pen$kkt_tol)
  # a KKT tolerance that cannot be met: the iteration goes on to max_iter, as
  # the L1-penalized algorithm prescribes, with the parameter
  # gate passed and no certificate
  strict <- mc_graph_penalized(f, lambda = 2, kkt_tol = 0)
  expect_false(isTRUE(strict$converged))
  expect_identical(strict$status, "kkt_not_certified")
  expect_identical(strict$iterations, 300L)
})

test_that("absent-edge scores are adjusted for the nuisance parameters", {
  skip_if_not_installed("MASS")
  set.seed(9); n <- 400
  K <- diag(4) * 2
  K[1,2] <- K[2,1] <- -0.6; K[2,3] <- K[3,2] <- -0.6; K[3,4] <- K[4,3] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
  d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], y4 = Y[,4], x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = d)
  sel <- mc_graph_stepwise(f, criterion = "gbic")
  es <- mc_edge_scores(sel$fit)
  absent <- es[!es$present, , drop = FALSE]
  skip_if(nrow(absent) == 0L, "the selected graph is complete")
  # Var(u_e) = J_ee - J_{e,lam} J_{lam,lam}^{-1} J_{lam,e} <= J_ee, so the
  # efficient statistic cannot fall below the unadjusted one
  expect_true(all(absent$stat >= absent$stat_unadjusted - 1e-8))
  expect_true(all(is.finite(absent$variance)))
  rob <- mc_edge_scores(sel$fit, variability = "empirical")
  expect_identical(attr(rob, "variability"), "empirical")
  expect_true(all(is.finite(rob$stat[!rob$present])))
})

test_that("the fit reports the diagnostics of its own stopping rule", {
  skip_if_not_installed("MASS")
  set.seed(4); n <- 300
  Y <- MASS::mvrnorm(n, rep(0, 3), diag(3))
  d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  expect_true(is.finite(f$graph$score_norm))
  expect_true(is.finite(f$graph$rcond_J))
  expect_true(nzchar(f$graph$stopping_rule))
  if (isTRUE(f$graph$converged)) expect_lt(f$graph$score_norm_unit, 1e-4)
  cr <- mc_compute_criterion(f, c("bic", "gbic"))
  expect_identical(attr(cr, "criterion_status"), "godambe")
  expect_true(is.finite(attr(cr, "effective_dimension")))
})

test_that("the residual degrees of freedom of the separation test are settable", {
  S <- diag(3); S[1, 2] <- S[2, 1] <- 0.4; S[2, 3] <- S[3, 2] <- 0.4
  a <- mcglmgraph:::.ci_pcor(S, 1, 3, 2, n = 400)
  b <- mcglmgraph:::.ci_pcor(S, 1, 3, 2, n = 400, df_residual = 60)
  expect_equal(a$df, 400 - 1 - 3)
  expect_equal(b$df, 60 - 1 - 2)
  expect_equal(a$pcor, b$pcor)
  expect_identical(a$inversion_status, "exact")
})

## ---- one edge contract, whichever door the caller came through -------------
## The fitting path reaches the kernels through .dag_setup() and .check_edges(),
## and a silent coercion there would let a user ask for one edge and be fitted
## another. These tests pin the strict contract on the path that mcglm_graph()
## actually takes, not only on the low-level LDL wrappers.

test_that("malformed edges are refused identically by both layers", {
  bad <- list(
    "non-integer index"  = rbind(c(1, 2.7)),
    "non-finite index"   = rbind(c(1, NA_real_)),
    "out of range"       = rbind(c(1, 9)),
    "self-loop"          = rbind(c(2, 2)),
    "duplicated pair"    = rbind(c(1, 2), c(1, 2)),
    "reversed duplicate" = rbind(c(1, 2), c(2, 1))
  )
  for (nm in names(bad)) {
    expect_error(mc_precision_bases(3, bad[[nm]]), info = paste("ggm:", nm))
    expect_error(mc_build_T_matrix(rep(0.1, nrow(bad[[nm]])), 3, bad[[nm]]),
                 info = paste("dag:", nm))
  }
})

test_that("the fitting path refuses what the low-level API refuses", {
  skip_if_not_installed("MASS")
  set.seed(11)
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
  Y <- MASS::mvrnorm(60, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(60))
  lp <- list(y1 ~ x, y2 ~ x, y3 ~ x)
  fit_with <- function(type, edges)
    mcglm_graph(lp, between_graph = mc_between_graph(type, edges), data = d,
                link = rep("identity", 3), variance = rep("constant", 3),
                power_fixed = rep(TRUE, 3))
  for (type in c("ggm", "dag")) {
    # A non-integer index is refused, never truncated to a different edge.
    expect_error(fit_with(type, rbind(c(1, 2.7))))
    expect_error(fit_with(type, rbind(c(1, 2), c(2, 1))))
    expect_error(fit_with(type, rbind(c(1, NA_real_))))
  }
})

test_that("R is validated rather than coerced", {
  expect_error(mc_precision_bases(2.7))
  expect_error(mc_precision_bases(c(3, 4)))
  expect_error(mc_build_T_matrix(numeric(0), 2.7,
                                 matrix(integer(0), ncol = 2)))
})

test_that("ordering is refused on the unordered GGM layer", {
  expect_error(mc_between_graph("ggm", ordering = c(3, 2, 1)))
  expect_error(mc_between_graph("dag", ordering = c(1, 1, 2)))
  expect_identical(mc_between_graph("dag", ordering = c(3, 1, 2))$ordering,
                   c(3L, 1L, 2L))
})
