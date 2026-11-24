# Residual degrees of freedom of the first stage.
#
# Targets, in plain R:
#   exact case   nu_e = n - rank(X) when every response is Gaussian with
#                identity link and constant variance and shares the column space
#   reduction    with the intercept alone, nu_e - |S| - 2 = n - |S| - 3
#   size         after projecting out h covariates, the partial correlation of
#                the residuals behaves as one of n - h observations, so the
#                Fisher test calibrated with nu_e is closer to nominal than the
#                one calibrated with n

test_that("the exact case is n minus the rank of the common design", {
  skip_if_not_installed("MASS")
  set.seed(1); n <- 60; x <- rnorm(n); g <- gl(3, 20)
  d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x, g = g)
  f <- mcglm_graph(list(y1 ~ x + g, y2 ~ x + g, y3 ~ x + g), data = d)
  rd <- mc_residual_df(f)
  expect_s3_class(rd, "mcg_residual_df")
  expect_true(rd$exact)
  expect_identical(rd$df, n - qr(stats::model.matrix(~ x + g, d))$rank)
  expect_identical(rd$rank, 4L)
  expect_length(rd$why, 0L)
  expect_output(print(rd), "56")
})

test_that("outside the exact case nothing is invented", {
  skip_if_not_installed("MASS")
  set.seed(2); n <- 80; x <- rnorm(n)
  d <- data.frame(y1 = rpois(n, exp(1 + .3 * x)), y2 = rpois(n, 3),
                  y3 = rpois(n, 4), x = x)
  f <- suppressMessages(mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d,
                                    link = rep("log", 3),
                                    variance = rep("tweedie", 3),
                                    power_fixed = rep(TRUE, 3)))
  rd <- mc_residual_df(f)
  expect_false(rd$exact)
  expect_true(is.na(rd$df))
  expect_match(paste(rd$why, collapse = " "), "constant variance")
  expect_match(paste(rd$why, collapse = " "), "identity link")
  expect_error(mc_graph_pipeline(f, df_residual = "auto"), "not exact")

  # Gaussian responses with different designs do not share a column space
  d2 <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x,
                   z = rnorm(n))
  f2 <- mcglm_graph(list(y1 ~ x, y2 ~ z, y3 ~ x), data = d2)
  rd2 <- mc_residual_df(f2)
  expect_false(rd2$exact)
  expect_match(paste(rd2$why, collapse = " "), "same column space")
})

test_that("with the intercept alone the contract is the default one", {
  skip_if_not_installed("MASS")
  set.seed(3); n <- 200
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.6 * z1 + 0.6 * z2 + rnorm(n)
  z4 <- 0.5 * z3 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, y4 = z4)
  f <- mcglm_graph(list(y1 ~ 1, y2 ~ 1, y3 ~ 1, y4 ~ 1), data = d)
  expect_identical(mc_residual_df(f)$df, n - 1)
  a <- mc_graph_pipeline(f)
  b <- mc_graph_pipeline(f, df_residual = "auto")
  expect_identical(a$cpdag, b$cpdag)
  expect_identical(a$separation_diagnostics$n_tests,
                   b$separation_diagnostics$n_tests)
  va <- mc_find_vstructures(f)
  vb <- mc_find_vstructures(f, df_residual = "auto")
  expect_identical(va$p_value, vb$p_value)
})

test_that("the functions that receive a covariance take the object, not auto", {
  skip_if_not_installed("MASS")
  set.seed(4); n <- 50; x <- rnorm(n)
  d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  S <- f$graph$Shat
  expect_error(mc_pc_skeleton(S, n, df_residual = "auto"),
               "needs the fitted model")
  a <- mc_pc_skeleton(S, n, df_residual = mc_residual_df(f))
  b <- mc_pc_skeleton(S, n, df_residual = n - 2)
  expect_identical(a$diagnostics$min_df, b$diagnostics$min_df)
  expect_identical(a$adjacency, b$adjacency)
  expect_error(mc_pc_skeleton(S, n, df_residual = list(1)), "must be NULL")
})

test_that("calibrated with nu_e the Fisher test is closer to nominal", {
  skip_on_cran()
  # three independent Gaussian responses, n = 30 units, eight covariates in
  # the mean: the partial correlation of responses 1 and 2 given 3 is zero
  set.seed(5); n <- 30; p <- 8; B <- 3000
  X <- cbind(1, matrix(rnorm(n * p), n, p))
  h <- qr(X)$rank
  rej <- matrix(FALSE, B, 2)
  for (b in seq_len(B)) {
    Y <- matrix(rnorm(n * 3), n, 3)
    E <- stats::lm.fit(X, Y)$residuals
    S <- crossprod(E) / n
    t_n <- mcglmgraph:::.ci_pcor(S, 1, 2, 3, n)
    t_e <- mcglmgraph:::.ci_pcor(S, 1, 2, 3, n, df_residual = n - h)
    rej[b, ] <- c(t_n$p_value, t_e$p_value) < 0.05
  }
  size <- colMeans(rej)
  se <- sqrt(0.05 * 0.95 / B)
  # the approximation by n rejects too often; nu_e is within Monte Carlo error
  expect_gt(size[1], 0.05 + 3 * se)
  expect_lt(abs(size[2] - 0.05), 3 * se)
  expect_lt(abs(size[2] - 0.05), abs(size[1] - 0.05))
})
