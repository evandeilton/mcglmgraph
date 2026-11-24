# Methods inherited from mcglm.
#
# An mcglm_graph object is an mcglm object with a graph block attached, so the
# mcglm methods have to answer on it as they do on the first-stage fit. They
# reach the coefficients through coef() with the mcglm types, which the method
# of this class has to pass on.

inherit_fit <- function(n = 150, seed = 3) {
  set.seed(seed)
  x <- rnorm(n)
  y1 <- 1 + 0.5 * x + rnorm(n)
  y2 <- 0.4 * y1 + 0.3 * x + rnorm(n)
  y3 <- 0.2 * y2 + rnorm(n)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
              data = data.frame(y1, y2, y3, x))
}

test_that("coef() answers the mcglm types as coef.mcglm does", {
  f <- inherit_fit()
  m <- f
  class(m) <- setdiff(class(f), "mcglm_graph")
  for (tp in c("beta", "tau", "power", "correlation"))
    expect_identical(coef(f, type = tp), coef(m, type = tp))
  expect_identical(coef(f, std.error = TRUE, response = 2, type = "beta"),
                   coef(m, std.error = TRUE, response = 2, type = "beta"))
  expect_identical(coef(f, type = c("beta", "tau")),
                   coef(m, type = c("beta", "tau")))
  expect_identical(coef(f), coef(m))
})

test_that("coef() keeps the types of the graph class", {
  f <- inherit_fit()
  expect_identical(coef(f, type = "graph"), f$graph$par)
  expect_identical(coef(f, type = "regression"), f$Regression)
  expect_identical(coef(f, type = "dispersion"), f$Covariance)
  expect_identical(coef(f, type = "gr"), f$graph$par)
  expect_error(coef(f, type = "nonsense"), "type must be one of")
  expect_error(coef(f, type = c("graph", "beta")), "type must be one of")
  expect_error(coef(f, type = 1), "type must be one of")
})

test_that("summary() prints the first stage and the graph layer", {
  f <- inherit_fit()
  m <- f
  class(m) <- setdiff(class(f), "mcglm_graph")
  out <- expect_no_warning(capture.output(summary(f)))
  ref <- capture.output(summary(m))
  # every line of the mcglm summary is there, followed by the graph block
  expect_identical(out[seq_along(ref)], ref)
  expect_true(any(grepl("Regression:", out, fixed = TRUE)))
  expect_true(any(grepl("=== Graph layer (GGM) ===", out, fixed = TRUE)))
})

test_that("anova() is inherited", {
  f <- inherit_fit()
  m <- f
  class(m) <- setdiff(class(f), "mcglm_graph")
  expect_identical(capture.output(a <- anova(f)), capture.output(b <- anova(m)))
  expect_identical(a, b)
})

test_that("a failing first-stage summary is reported, not hidden", {
  f <- inherit_fit()
  f$beta_names <- NULL               # what summary.mcglm reads first
  expect_warning(capture.output(summary(f)), "first-stage summary")
})
