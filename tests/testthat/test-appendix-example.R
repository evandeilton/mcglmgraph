# Locks a complete worked example of the workflow. The example must keep
# running literally on the public API:
# mcglm_graph(between_graph = ...) -> mc_graph_stepwise() -> sel$fit -> readers.

test_that("the worked example runs on the public API", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("mcglm")

  set.seed(2026)
  n <- 200
  z <- MASS::mvrnorm(n, mu = rep(0, 3),
                     Sigma = matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3, 3))
  dados <- data.frame(x = rnorm(n), y1 = NA, y2 = NA, y3 = NA)
  dados$y1 <- 1 + 0.8 * dados$x + z[, 1]
  dados$y2 <- rpois(n, exp(0.5 + 0.3 * dados$x + 0.4 * z[, 2]))
  dados$y3 <- rpois(n, exp(0.7 - 0.2 * dados$x + 0.4 * z[, 3]))

  Z0 <- mcglm::mc_id(dados)

  fit <- mcglm_graph(
    linear_pred   = c(y1 ~ x, y2 ~ x, y3 ~ x),
    matrix_pred   = list(Z0, Z0, Z0),
    link          = c("identity", "log", "log"),
    variance      = c("constant", "tweedie", "tweedie"),
    between_graph = mc_between_graph("ggm"),
    data          = dados
  )
  expect_s3_class(fit, "mcglm_graph")
  expect_true(fit$graph$converged)

  sel <- mc_graph_stepwise(fit, criterion = "gbic")
  expect_s3_class(sel, "mcg_selection")
  expect_s3_class(sel$fit, "mcglm_graph")

  P <- mc_partial_cor(sel$fit)
  expect_equal(dim(P), c(3L, 3L))
  expect_equal(diag(P), rep(1, 3), ignore_attr = TRUE)

  skip_if(sel$n_edges == 0L, "empty support has no edge parameter to test")
  tt <- mc_test_independence(sel$fit, robust = TRUE)
  expect_s3_class(tt, "htest")
  expect_true(is.finite(tt$p.value))
})

test_that("the misspecified forms of the appendix example fail, as documented", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("mcglm")

  set.seed(7)
  n <- 120
  Y <- MASS::mvrnorm(n, rep(0, 3), diag(3))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))

  # graph_layer is not an argument of mcglm_graph(); it reaches mcglm() via ...
  expect_error(
    mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), graph_layer = "ggm", data = d)
  )

  # the readers require the refitted mcglm_graph object, not the selection
  fit <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  sel <- mc_graph_stepwise(fit, criterion = "gbic")
  expect_error(mc_partial_cor(sel))
  expect_error(mc_test_independence(sel))
})
