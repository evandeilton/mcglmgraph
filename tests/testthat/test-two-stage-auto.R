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
