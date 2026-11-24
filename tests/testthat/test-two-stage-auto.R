# Route of method = "auto" for the two-stage correction. The shortcut evaluates
# the reduced form at the fitted covariance: exact at every node of a DAG
# without immoralities, first-order only at a collider with non-adjacent
# parents and its descendants. On such a DAG "auto" takes the full form.

auto_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
  x <- stats::rnorm(n)
  y3 <- 0.6 * z[, 1] + 0.5 * z[, 2] + stats::rnorm(n)
  y4 <- 0.5 * y3 + stats::rnorm(n)
  data.frame(y1 = z[, 1] + 0.3 * x, y2 = z[, 2], y3 = y3, y4 = y4, x = x)
}
auto_fit <- function(d, spec)
  suppressMessages(mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x),
                               between_graph = spec, data = d))
bare_v <- function(M) {
  M <- as.matrix(M)
  attributes(M) <- list(dim = dim(M))
  M
}

test_that("on a DAG with an immorality auto takes the full form", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- auto_data()
  f <- auto_fit(d, mc_between_graph("dag", rbind(c(1, 3), c(2, 3), c(3, 4))))
  imm <- getFromNamespace(".dag_immoralities", "mcglmgraph")
  expect_identical(imm(f), 3L)
  Va <- mc_vcov_graph(f, "two_stage")
  Vf <- mc_vcov_graph(f, "two_stage", method = "full")
  Vs <- expect_silent(mc_vcov_graph(f, "two_stage", method = "shortcut"))
  expect_identical(attr(Va, "two_stage_method"), "full")
  expect_identical(bare_v(Va), bare_v(Vf))
  expect_gt(max(abs(bare_v(Vs) - bare_v(Vf))) / max(abs(bare_v(Vf))), 1e-3)
  # the gap is the collider, not the first-stage sensitivity: the shortcut
  # differs from the full form with the expected sensitivity as well
  Ve <- mc_vcov_graph(f, "two_stage", method = "full", sensitivity = "expected")
  expect_gt(max(abs(bare_v(Vs) - bare_v(Ve))) / max(abs(bare_v(Ve))), 1e-3)
  expect_identical(mc_two_stage_available(f)$method, "full")
})

test_that("auto keeps the shortcut with a warning when the full form is unavailable", {
  skip_if_not_installed("MASS")
  d <- auto_data()
  f <- auto_fit(d, mc_between_graph("dag", rbind(c(1, 3), c(2, 3), c(3, 4))))
  # A fit without its first-stage components cannot rebuild the first stage.
  f$list_X <- NULL
  expect_false(mc_two_stage_available(f)$full_ok)
  expect_warning(Va <- mc_vcov_graph(f, "two_stage"), "first-order")
  expect_identical(attr(Va, "two_stage_method"), "shortcut")
  expect_silent(mc_vcov_graph(f, "two_stage", method = "shortcut"))
  expect_identical(mc_two_stage_available(f)$method, "shortcut")
})

test_that("on a DAG without immorality auto is the shortcut, which is exact", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- auto_data()
  imm <- getFromNamespace(".dag_immoralities", "mcglmgraph")
  supports <- list(chain = rbind(c(1, 2), c(2, 3), c(3, 4)),
                   triangle = rbind(c(1, 3), c(2, 3), c(1, 2)))
  for (nm in names(supports)) {
    f <- auto_fit(d, mc_between_graph("dag", supports[[nm]]))
    expect_identical(imm(f), integer(0), info = nm)
    Va <- expect_silent(mc_vcov_graph(f, "two_stage"))
    Vs <- mc_vcov_graph(f, "two_stage", method = "shortcut")
    expect_identical(attr(Va, "two_stage_method"), "shortcut", info = nm)
    expect_identical(bare_v(Va), bare_v(Vs), info = nm)
    # With the observed sensitivity of the equations the first stage solves,
    # the full form is the reduced form at the data; without immoralities the
    # shortcut equals it, up to the convergence tolerance of the first stage.
    Vo <- mc_vcov_graph(f, "two_stage", method = "full")
    expect_lt(max(abs(bare_v(Vs) - bare_v(Vo))) / max(abs(bare_v(Vo))), 1e-5,
              label = nm)
  }
})

test_that("the GGM layer keeps its route", {
  skip_if_not_installed("MASS")
  d <- auto_data()
  f <- auto_fit(d, mc_between_graph("ggm", rbind(c(1, 3), c(2, 3), c(3, 4))))
  imm <- getFromNamespace(".dag_immoralities", "mcglmgraph")
  expect_identical(imm(f), integer(0))
  Va <- expect_silent(mc_vcov_graph(f, "two_stage"))
  expect_identical(attr(Va, "two_stage_method"), "shortcut")
  expect_identical(bare_v(Va),
                   bare_v(mc_vcov_graph(f, "two_stage", method = "shortcut")))
  expect_identical(mc_two_stage_available(f)$method, "shortcut")
})

# ---- the variance function decides the route ---------------------------------

# The route stated from how each fit is built: the full form when a response
# has a variance function other than "constant" or the DAG has an immorality,
# the shortcut otherwise.
route_expected <- function(variance, immoral)
  if (any(variance != "constant") || immoral) "full" else "shortcut"

test_that("auto reads the variance functions and the immorality of the fit", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  imm <- mc_between_graph("dag", rbind(c(1, 3), c(2, 3)))
  chain <- mc_between_graph("dag", rbind(c(1, 2), c(2, 3)))
  cases <- list(
    list(kind = "gaussian", graph = mc_between_graph("ggm"),
         variance = rep("constant", 3), immoral = FALSE),
    list(kind = "gaussian", graph = chain,
         variance = rep("constant", 3), immoral = FALSE),
    list(kind = "gaussian", graph = imm,
         variance = rep("constant", 3), immoral = TRUE),
    list(kind = "tweedie", graph = mc_between_graph("ggm"),
         variance = rep("tweedie", 3), immoral = FALSE),
    list(kind = "mixed", graph = mc_between_graph("ggm"),
         variance = c("constant", "tweedie", "tweedie"), immoral = FALSE),
    list(kind = "binomialP", graph = mc_between_graph("ggm"),
         variance = rep("binomialP", 3), immoral = FALSE))
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    f <- route_fit(cs$kind, between_graph = cs$graph)
    want <- route_expected(cs$variance, cs$immoral)
    expect_identical(unlist(f$variance), cs$variance, info = k)
    Va <- mc_vcov_graph(f, "two_stage")
    Vm <- mc_vcov_graph(f, "two_stage", method = want)
    expect_identical(attr(Va, "two_stage_method"), want, info = k)
    expect_identical(mc_two_stage_available(f)$method, want, info = k)
    expect_identical(bare_v(Va), bare_v(Vm), info = k)
  }
})

test_that("a forced route is honored whatever the variance function", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- route_fit("tweedie")
  Vs <- expect_silent(mc_vcov_graph(f, "two_stage", method = "shortcut"))
  Vf <- mc_vcov_graph(f, "two_stage", method = "full")
  expect_identical(attr(Vs, "two_stage_method"), "shortcut")
  expect_identical(attr(Vf, "two_stage_method"), "full")
  expect_gt(max(abs(bare_v(Vs) - bare_v(Vf))), 0)
  g <- route_fit("gaussian")
  expect_identical(
    attr(mc_vcov_graph(g, "two_stage", method = "full"), "two_stage_method"),
    "full")
})

test_that("a mean-dependent fit without its first-stage components takes the shortcut", {
  skip_if_not_installed("MASS")
  f <- route_fit("tweedie")
  f$list_X <- NULL
  av <- mc_two_stage_available(f)
  expect_false(av$full_ok)
  expect_true(av$shortcut_ok)
  expect_identical(av$method, "shortcut")
  Va <- expect_silent(mc_vcov_graph(f, "two_stage"))
  expect_identical(attr(Va, "two_stage_method"), "shortcut")
  expect_identical(bare_v(Va), bare_v(mc_vcov_graph(f, "two_stage",
                                                    method = "shortcut")))
  expect_error(mc_vcov_graph(f, "two_stage", method = "full"),
               "not available|does not apply")
})

test_that("variance labels that cannot be read count as depending on the mean", {
  skip_if_not_installed("MASS")
  dep <- getFromNamespace(".two_stage_mean_dependent", "mcglmgraph")
  f <- route_fit("gaussian")
  expect_false(dep(f))
  g <- f; g$variance <- list("constant", "constant")   # wrong length
  expect_true(dep(g))
  g <- f; g$variance <- NULL
  expect_true(dep(g))
  g <- f; g$variance <- list(1, 2, 3)
  expect_true(dep(g))
  g <- f; g$variance <- list("constant", NA_character_, "constant")
  expect_true(dep(g))
})

test_that("the name of the variance function decides, not its power", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- route_fit("tweedie0")
  expect_identical(mc_two_stage_available(f)$method, "full")
  Vf <- mc_vcov_graph(f, "two_stage")
  Vs <- mc_vcov_graph(f, "two_stage", method = "shortcut")
  expect_identical(attr(Vf, "two_stage_method"), "full")
  # with the power fixed at zero the variance is constant: the routes agree
  expect_lt(max(abs(bare_v(Vf) - bare_v(Vs))) / max(abs(bare_v(Vf))), 1e-6)
})

test_that("the full form is the derivative of the second stage for binomial margins", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  f <- route_fit("binomialP")
  g <- f$graph
  V <- mc_vcov_graph(f, "two_stage")
  expect_identical(attr(V, "two_stage_method"), "full")
  B <- attr(V, "B")
  X <- cbind(1, f$data$x)
  Y <- as.matrix(f$data[, c("y1", "y2", "y3")])
  Sinv <- solve(g$Sigma_b)
  # residuals as a function of alpha = (beta, rho, tau), in plain R: logit
  # link and v = tau {mu (1 - mu)}^p / N, with the power fixed at its start
  pw <- unlist(f$list_initial$power)
  resid_at <- function(a) {
    b <- matrix(a[1:6], 2)
    tau <- a[10:12]
    sapply(1:3, function(r) {
      mu <- stats::plogis(drop(X %*% b[, r]))
      (Y[, r] - mu) / sqrt(tau[r] * (mu * (1 - mu))^pw[r] / 10)
    })
  }
  U2 <- function(a) {
    Ec <- scale(resid_at(a), center = TRUE, scale = FALSE)
    vapply(g$Dlist, function(Dd) {
      A <- Sinv %*% Dd %*% Sinv
      -0.5 * sum(sum(diag(Sinv %*% Dd)) - rowSums((Ec %*% A) * Ec))
    }, numeric(1))
  }
  a0 <- c(f$Regression, f$Covariance)
  expect_equal(unname(resid_at(a0)), unname(as.matrix(g$residuals)),
               tolerance = 1e-8)
  Bnum <- -numDeriv::jacobian(U2, a0)
  expect_lt(max(abs(Bnum - B)) / max(abs(B)), 1e-6)
})
