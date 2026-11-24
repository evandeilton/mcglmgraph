test_that("LDL covariance is positive definite for arbitrary gamma", {
  edges <- rbind(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4))
  S <- mc_cholesky_to_covariance(c(0.7, 0.3, 0.5, 0.2, 0.1, 0.4),
                                 c(1, 0.8, 0.6, 0.5), R = 4, edges = edges)$Sigma_b
  expect_true(min(eigen(S, only.values = TRUE)$values) > 0)
})

test_that("DAG analytical derivatives match numerical (symmetric structural derivative)", {
  edges <- rbind(c(1,2), c(1,3), c(2,3))
  gamma <- c(0.7, 0.3, 0.5); delta <- c(1, 0.64, 0.49); R <- 3
  der <- mc_cholesky_derivatives(gamma, delta, R, edges)
  eps <- 1e-6
  # gamma derivative
  gp <- gamma; gp[3] <- gp[3] + eps; gm <- gamma; gm[3] <- gm[3] - eps
  Sp <- mc_cholesky_to_covariance(gp, delta, R, edges)$Sigma_b
  Sm <- mc_cholesky_to_covariance(gm, delta, R, edges)$Sigma_b
  expect_lt(max(abs(der$D_gamma[[3]] - (Sp - Sm) / (2 * eps))), 1e-5)
  # delta derivative
  dp <- delta; dp[2] <- dp[2] + eps; dm <- delta; dm[2] <- dm[2] - eps
  Sp <- mc_cholesky_to_covariance(gamma, dp, R, edges)$Sigma_b
  Sm <- mc_cholesky_to_covariance(gamma, dm, R, edges)$Sigma_b
  expect_lt(max(abs(der$D_delta[[2]] - (Sp - Sm) / (2 * eps))), 1e-5)
})

test_that("precision is the inverse of the LDL covariance", {
  edges <- rbind(c(1,2), c(1,3), c(2,3))
  S <- mc_cholesky_to_covariance(c(.7,.3,.5), c(1,.64,.49), 3, edges)$Sigma_b
  K <- mc_cholesky_to_precision(c(.7,.3,.5), c(1,.64,.49), 3, edges)
  expect_lt(max(abs(K %*% S - diag(3))), 1e-10)
})

test_that("DAG ordering permutes the structural system consistently", {
  theta <- c(.4, .2, .3, 1, .8, .6)
  natural <- mc_build_sigma_between_graph(
    theta, 3, mc_between_graph("dag", ordering = 1:3), inverse = TRUE)
  reversed <- mc_build_sigma_between_graph(
    theta, 3, mc_between_graph("dag", ordering = 3:1), inverse = TRUE)
  expect_equal(reversed$Sigmab, natural$Sigmab[3:1, 3:1], tolerance = 1e-12)
  expect_equal(reversed$Kb, natural$Kb[3:1, 3:1], tolerance = 1e-12)
})

test_that("invalid DAG ordering is rejected", {
  expect_error(
    mc_build_sigma_between_graph(rep(1, 6), 3,
                                 mc_between_graph("dag", ordering = c(1, 1, 3))),
    "permutation"
  )
})

test_that("LDL wrappers reject duplicated or malformed edge coordinates", {
  duplicated <- rbind(c(1, 2), c(2, 1))
  expect_error(mc_build_T_matrix(c(.2, .7), 2, duplicated), "duplicated")
  expect_error(
    mc_cholesky_derivatives(c(.2, .7), c(1, 1), 2, duplicated),
    "duplicated"
  )

  expect_error(mc_build_T_matrix(.2, 2, matrix(1:3, nrow = 1)),
               "exactly two columns")
  expect_error(mc_build_T_matrix(.2, 2, matrix(c(1, Inf), nrow = 1)),
               "finite integer indices")
  expect_error(mc_build_T_matrix(.2, 2, matrix(c(1, 1.5), nrow = 1)),
               "finite integer indices")
})

test_that("LDL wrappers validate dimensions, finiteness and positive delta", {
  edge <- matrix(c(1, 2), nrow = 1)

  expect_error(mc_build_T_matrix(numeric(0), 2, edge), "gamma")
  expect_error(mc_build_T_matrix(Inf, 2, edge), "gamma")
  expect_error(mc_cholesky_to_covariance(c(.2, .3), c(1, 1), 2, edge),
               "gamma")

  expect_error(mc_build_D_matrix(c(1, 2), 3), "exactly 3")
  expect_error(mc_build_D_matrix(c(1, NA_real_), 2), "finite positive")
  expect_error(mc_build_D_matrix(c(1, 0), 2), "strictly positive")
  expect_error(mc_cholesky_to_covariance(.2, c(1, Inf), 2, edge),
               "finite positive")
  expect_error(mc_cholesky_to_precision(.2, c(1, -1), 2, edge),
               "strictly positive")
  expect_error(mc_cholesky_derivatives(.2, c(1, 0), 2, edge),
               "strictly positive")

  expect_error(mc_build_D_matrix(1, 0), "positive integer")
  expect_error(mc_build_T_matrix(.2, 2.5, edge), "positive integer")
  expect_error(mc_build_T_matrix(.2, c(2, 2), edge), "positive integer")

  expect_equal(mc_build_D_matrix(2, 1), matrix(2, 1, 1))
  expect_equal(mc_build_T_matrix(numeric(0), 1), matrix(1, 1, 1))
})

test_that("DAG initialization regularizes a singular sample covariance", {
  par <- mcglmgraph:::.init_par(
    mc_between_graph("dag", ordering = 1:3), 3L, matrix(1, 3L, 3L))
  info <- attr(par, "init_info", exact = TRUE)
  expect_true(all(is.finite(par)))
  expect_true(all(tail(par, 3L) > 0))
  expect_true(is.list(info))
  expect_gt(info$spd_adjustment, 0)
  expect_lte(info$min_eigen_before, info$eigen_floor)
})
