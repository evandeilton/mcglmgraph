test_that("Fisher residual degrees of freedom and failures are explicit", {
  S <- diag(3)
  S[1, 2] <- S[2, 1] <- 0.3
  raw <- mcglmgraph:::.ci_pcor(S, 1, 2, integer(), n = 100)
  residual <- mcglmgraph:::.ci_pcor(
    S, 1, 2, integer(), n = 100, df_residual = 80
  )
  expect_equal(raw$df, 97)
  expect_equal(residual$df, 78)
  expect_error(
    mcglmgraph:::.ci_pcor(
      S, 1, 2, S = 3L, n = 100, df_residual = 3
    ),
    "Insufficient Fisher degrees of freedom"
  )

  singular <- matrix(1, 3, 3)
  reg <- mcglmgraph:::.ci_pcor(singular, 1, 2, integer(), n = 100)
  expect_identical(reg$inversion_status, "ridge_regularized")
  expect_gt(reg$inversion_ridge, 0)
  skel <- mc_pc_skeleton(singular, n = 100, max_ord = 0)
  expect_gt(skel$diagnostics$n_regularized, 0L)
  expect_identical(skel$diagnostics$df_contract, "n - |S| - 3")
})

test_that("df_residual is propagated through v-structure discovery", {
  skip_if_not_installed("MASS")
  set.seed(120); n <- 250
  d <- data.frame(
    y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = rnorm(n)
  )
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  expect_no_error(mc_find_vstructures(f, df_residual = 230))
  p <- mc_graph_pipeline(f, df_residual = 230)
  expect_equal(p$df_residual, 230)
  expect_identical(
    p$separation_diagnostics$df_contract,
    "df_residual - |S| - 2"
  )
})

test_that("Meek propagation leaves simultaneous opposite requests unresolved", {
  A <- matrix(0L, 4, 4)
  A[1, 2] <- A[2, 1] <- 1L
  A[3, 1] <- 1L
  A[4, 2] <- 1L
  out <- mc_orient_meek(A)
  expect_equal(out[1, 2], 1L)
  expect_equal(out[2, 1], 1L)
  conflicts <- attr(out, "meek_conflicts", exact = TRUE)
  expect_equal(nrow(conflicts), 1L)
  expect_equal(unname(as.integer(conflicts[1, ])), c(1L, 2L))
})

test_that("score orientation reports its own skeleton and preserves the PC one", {
  S <- matrix(0.8, 3, 3); diag(S) <- 1
  fake <- structure(list(graph = list(
    Shat = S, n = 10000, response_names = c("y1", "y2", "y3")
  )), class = "mcglm_graph")
  score_cpdag <- matrix(0L, 3, 3)
  score_cpdag[1, 2] <- 1L
  local_mocked_bindings(
    .bnlearn_available = function() TRUE,
    .orient_bnlearn = function(skel_adj, fit, method) score_cpdag,
    .package = "mcglmgraph"
  )
  out <- mc_graph_pipeline(fake, orient_method = "hc", max_ord = 0)
  expect_true(all(out$pc_skeleton[upper.tri(out$pc_skeleton)] == 1L))
  expect_equal(out$skeleton, 1L * ((out$cpdag == 1L) |
                                    (t(out$cpdag) == 1L)))
  expect_equal(unname(out$edges), matrix(c(1L, 2L), ncol = 2L))
  expect_false(out$gates[["skeleton_preserved"]])
  expect_false(out$valid_cpdag)
})

test_that("matrix inversion never regularizes without a status", {
  exact <- mcglmgraph:::.safe_solve(diag(2))
  expect_identical(attr(exact, "inversion_status"), "exact")
  expect_equal(attr(exact, "inversion_ridge"), 0)

  regularized <- mcglmgraph:::.safe_solve(matrix(0, 2, 2))
  expect_identical(
    attr(regularized, "inversion_status"), "ridge_regularized"
  )
  expect_gt(attr(regularized, "inversion_ridge"), 0)
})

test_that("penalized attachment replaces every fit diagnostic", {
  fake <- structure(list(graph = list(
    R = 2L, n = 100, Shat = diag(2), score_norm = 999,
    score_norm_unit = 999, rcond_J = 0, stopping_rule = "stale",
    rcond_tol = .Machine$double.eps^(2 / 3)
  )), class = "mcglm_graph")
  spec <- mc_between_graph("ggm", matrix(integer(), ncol = 2L))
  pen <- list(
    eta = c(1, 1, 0), edges = matrix(integer(), ncol = 2L),
    iterations = 2L, converged = TRUE, status = "converged",
    kkt_unit = 0, kkt_tol = 1e-4
  )
  out <- mcglmgraph:::.attach_penalized_fit(fake, spec, pen)
  expect_equal(out$graph$score_norm, 0)
  expect_equal(out$graph$score_norm_unit, 0)
  expect_true(is.finite(out$graph$rcond_J))
  expect_false(identical(out$graph$stopping_rule, "stale"))
  expect_false(out$graph$inference_ok)
  expect_identical(out$graph$linear_solve_status,
                   "not_applicable_penalized")
})

test_that("inferential APIs refuse explicitly uncertified fits", {
  fake <- structure(list(graph = list(
    type = "ggm", inference_ok = FALSE, status = "score_not_converged",
    inversion_status = "ridge_regularized", rcond_J = 0
  )), class = "mcglm_graph")
  expect_error(mc_test_independence(fake), "not certified")
  expect_error(mc_partial_cor_delta(fake), "not certified")
  expect_error(mc_edge_scores(fake), "not certified")
})

test_that("the identity structure is recognized before Matrix is attached", {
  # The default within-response structure is a ddiMatrix. Matrix is in Imports,
  # so in a session that has not attached it the S4 method for dim() is not
  # dispatched and dim() answers NULL on such an object, which would send
  # .is_identity_structure() down the FALSE branch before it could reach the
  # ddiMatrix case, so mc_graph_stability(refit = "full") must work on the
  # first call of a session. The check runs in a
  # fresh process, since any earlier call in this one would already have loaded
  # the namespace.
  skip_if_not_installed("callr")
  out <- callr::r(function() {
    set.seed(11)
    n <- 60; x <- stats::rnorm(n)
    d <- data.frame(y1 = stats::rnorm(n), y2 = stats::rnorm(n),
                    y3 = stats::rnorm(n), x = x)
    f <- mcglmgraph::mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
    M <- f$matrix_pred[[1L]][[1L]]
    list(attached   = "package:Matrix" %in% search(),
         dim_null   = is.null(dim(M)),
         recognized = mcglmgraph:::.identity_matrix_pred(f))
  })
  expect_false(out$attached)          # the premise of the test
  expect_true(out$recognized)         # the contract
})
