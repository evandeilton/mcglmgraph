# =============================================================================
# ALIAS mcglmgraph() / mcglm_graph()
# =============================================================================
# The two names denote the same object. This test exists to prevent a future
# edit from turning them into separate implementations, which is the only way
# an alias bound by assignment can become a problem.

test_that("both names point to the same function", {
  expect_identical(mcglmgraph, mcglm_graph)
  expect_identical(formals(mcglmgraph), formals(mcglm_graph))
  expect_identical(body(mcglmgraph), body(mcglm_graph))
})

test_that("both names produce the same fit and the same class", {
  skip_if_not_installed("MASS")
  set.seed(20260810)
  n <- 200
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  lp <- list(y1 ~ x, y2 ~ x, y3 ~ x)

  a <- mcglm_graph(lp, data = d)
  b <- mcglmgraph(lp, data = d)

  ## The class stays "mcglm_graph" under both names: it is the identifier
  ## recorded inside the archived fits and the one the S3 methods are
  ## registered on.
  expect_s3_class(a, "mcglm_graph")
  expect_s3_class(b, "mcglm_graph")
  expect_identical(class(a), class(b))

  expect_equal(a$graph$par, b$graph$par)
  expect_equal(a$graph$Sigma_b, b$graph$Sigma_b)
  expect_equal(a$graph$K_b, b$graph$K_b)
  expect_equal(a$graph$vcov, b$graph$vcov)

  ## the reader functions accept an object coming from either name
  expect_equal(mc_partial_cor(a), mc_partial_cor(b))
  expect_equal(mc_compute_criterion(a, "gbic"), mc_compute_criterion(b, "gbic"))
})
