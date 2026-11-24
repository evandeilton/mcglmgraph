# Scale-free quantities of the directed layer and their intervals.
#
# Targets, in plain R:
#   gamma*_jk = gamma_jk (Sigma_kk / Sigma_jj)^(1/2), R2_j = 1 - delta_j / Sigma_jj,
#   with Sigma = T^-1 D T^-T built from (gamma, delta)
#   gradients  by numDeriv on those two functions of the parameters
#   variance   g' V g, the delta method
#   coverage   nominal in a Gaussian chain, by Monte Carlo

dag_fit <- function(n = 400, seed = 4, edges = NULL) {
  set.seed(seed); x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + 0.3 * y1 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
              between_graph = mc_between_graph("dag", edges = edges), data = d)
}

## gamma* and R2 as functions of the parameter vector, through the layer's own
## parametrization Sigma = T^-1 D T^-T.
scale_free_of <- function(fit) {
  g <- fit$graph
  R <- g$R
  dag <- mcglmgraph:::.dag_setup(g$spec, R)
  eo <- dag$edges_original
  E <- nrow(eo)
  function(par) {
    S <- mc_cholesky_to_covariance(par[seq_len(E)], par[E + seq_len(R)], R,
                                   dag$edges_pos)$Sigma_b
    S <- S[dag$inverse, dag$inverse]          # back to response order
    gs <- par[seq_len(E)] * sqrt(S[cbind(eo[, 1], eo[, 1])] /
                                   S[cbind(eo[, 2], eo[, 2])])
    r2 <- 1 - par[E + seq_len(R)] /
      S[cbind(dag$ordering, dag$ordering)]
    c(gs, r2)
  }
}

test_that("the gradients of gamma* and R2 are those of the definitions", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  for (f in list(dag_fit(), dag_fit(edges = rbind(c(1, 2), c(2, 3))))) {
    co <- mc_dag_coefficients(f)
    E <- nrow(co$edges)
    fn <- scale_free_of(f)
    par <- as.numeric(f$graph$par)
    expect_equal(fn(par)[seq_len(E)], co$edges$gamma_std, tolerance = 1e-12)
    expect_equal(fn(par)[-seq_len(E)], co$nodes$r_squared, tolerance = 1e-12)
    Jn <- numDeriv::jacobian(fn, par)
    expect_lt(max(abs(Jn[seq_len(E), ] - attr(co, "gradient_std"))), 1e-7)
    expect_lt(max(abs(Jn[-seq_len(E), ] - attr(co, "gradient_r2"))), 1e-7)
  }
})

test_that("the standard errors are the delta method on the chosen covariance", {
  skip_if_not_installed("MASS")
  f <- dag_fit()
  for (v in c("model", "empirical", "two_stage")) {
    co <- mc_dag_coefficients(f, variability = v)
    V <- as.matrix(mc_vcov_graph(f, switch(v, model = "model",
                                           empirical = "sandwich",
                                           two_stage = "two_stage")))
    Gs <- attr(co, "gradient_std"); Gr <- attr(co, "gradient_r2")
    expect_equal(co$edges$se_std, sqrt(diag(Gs %*% V %*% t(Gs))),
                 tolerance = 1e-10, info = v)
    expect_equal(co$nodes$se_r_squared, sqrt(diag(Gr %*% V %*% t(Gr))),
                 tolerance = 1e-10, info = v)
    q <- stats::qnorm(0.975)
    expect_equal(co$edges$upper_std - co$edges$gamma_std, q * co$edges$se_std,
                 tolerance = 1e-12)
  }
})

test_that("the fixed-factor interval is the one of gamma carried through the factor", {
  skip_if_not_installed("MASS")
  f <- dag_fit()
  co <- mc_dag_coefficients(f, interval = "fixed_factor")
  fac <- co$edges$gamma_std / co$edges$gamma
  expect_equal(co$edges$lower_std, co$edges$lower * fac, tolerance = 1e-12)
  expect_equal(co$edges$upper_std, co$edges$upper * fac, tolerance = 1e-12)
  # the delta interval is wider, since it also carries the variances in fac
  cd <- mc_dag_coefficients(f)
  expect_true(all(cd$edges$se_std > 0))
  expect_identical(cd$interval, "delta")
})

test_that("the delta intervals cover at the nominal rate in a Gaussian chain", {
  skip_on_cran()
  skip_long()
  skip_if_not_installed("MASS")
  # population: y1 -> y2 -> y3 with unit innovations
  g12 <- 0.6; g23 <- 0.5
  S <- mc_cholesky_to_covariance(c(g12, 0, g23), c(1, 1, 1), R = 3)$Sigma_b
  truth_std <- c(g12 * sqrt(S[1, 1] / S[2, 2]), 0, g23 * sqrt(S[2, 2] / S[3, 3]))
  truth_r2 <- 1 - c(1, 1, 1) / diag(S)
  B <- 500
  hit <- t(vapply(seq_len(B), function(b) {
    set.seed(5000 + b); n <- 300; x <- rnorm(n)
    y1 <- rnorm(n); y2 <- g12 * y1 + rnorm(n); y3 <- g23 * y2 + rnorm(n)
    d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
    f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                     between_graph = mc_between_graph("dag"), data = d)
    co <- mc_dag_coefficients(f, variability = "empirical")
    c(co$edges$lower_std <= truth_std & truth_std <= co$edges$upper_std,
      co$nodes$lower_r_squared <= truth_r2 & truth_r2 <= co$nodes$upper_r_squared)
  }, logical(6)))
  cover <- colMeans(hit)
  se <- sqrt(0.95 * 0.05 / B)
  # the first response has no parent, so its R2 is zero by construction, with
  # a zero standard error and a degenerate interval: it is not a quantity to
  # cover, and it is left out
  expect_true(all(abs(cover[-4] - 0.95) < 3.5 * se))
})

test_that("a response without parents has R2 zero and nothing to estimate", {
  skip_if_not_installed("MASS")
  co <- mc_dag_coefficients(dag_fit())
  root <- co$nodes$n_parents == 0L
  expect_true(any(root))
  expect_equal(co$nodes$r_squared[root], rep(0, sum(root)), tolerance = 1e-12)
  expect_equal(co$nodes$se_r_squared[root], rep(0, sum(root)), tolerance = 1e-12)
})
