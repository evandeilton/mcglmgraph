# Coverage of the theory:
# Godambe-adjusted BIC, selection-frequency aggregation and the robust
# sandwich Wald test.

make_chain <- function(n = 400, seed = 7) {
  set.seed(seed)
  K <- diag(4) * 2
  K[1, 2] <- K[2, 1] <- -0.6; K[2, 3] <- K[3, 2] <- -0.6; K[3, 4] <- K[4, 3] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K)); x <- rnorm(n)
  data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4], x = x)
}

test_that("GBIC is computed with the effective dimension and has the alias clbic", {
  skip_if_not_installed("MASS")
  d <- make_chain()
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
  cr <- mc_compute_criterion(f, c("aic", "bic", "qic", "gbic", "clbic"))
  expect_true(all(c("aic", "bic", "qic", "gbic", "clbic") %in% names(cr)))
  expect_true(is.finite(cr["gbic"]))
  expect_equal(unname(cr["gbic"]), unname(cr["clbic"]))
  base <- -2 * mc_graph_ploglik(f)
  pen_qic <- cr["qic"] - base       # = 2 * effdim
  pen_gb <- cr["gbic"] - base       # = log(n) * effdim
  expect_gt(unname(pen_gb), unname(pen_qic))           # log(400) > 2
  # effective dimension recovered consistently from both penalties
  expect_equal(unname(pen_gb / log(nrow(d))), unname(pen_qic / 2), tolerance = 1e-8)
})

test_that("stepwise and lambda selection accept GBIC and its alias clbic", {
  skip_if_not_installed("MASS")
  d <- make_chain()
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
  sel <- mc_graph_stepwise(f, criterion = "gbic")
  expect_s3_class(sel, "mcg_selection")
  expect_equal(sel$criterion, "gbic")
  sel_alias <- mc_graph_stepwise(f, criterion = "clbic")
  expect_equal(sel_alias$criterion, "clbic")
})

test_that("selection-frequency aggregation returns a valid frequency matrix", {
  skip_if_not_installed("MASS")
  d <- make_chain()
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
  st <- mc_graph_stability(f, B = 25, pi_thr = 0.5)
  expect_s3_class(st, "mcg_stability")
  expect_equal(dim(st$selection_prob), c(4L, 4L))
  expect_true(all(st$selection_prob >= 0 & st$selection_prob <= 1))
  # adjacent (true) edges selected more often than the mediated pair (1,3)
  expect_gte(st$selection_prob[1, 2], st$selection_prob[1, 3])
  expect_gte(st$selection_prob[3, 4], st$selection_prob[1, 4])
})

test_that("selection frequencies use valid subsamples as their denominator", {
  set.seed(91)
  fake <- structure(list(graph = list(
    type = "ggm", residuals = matrix(rnorm(60), ncol = 2L),
    R = 2L, response_names = c("y1", "y2")
  )), class = "mcglm_graph")
  calls <- 0L
  local_mocked_bindings(
    mc_graph_stepwise = function(fit, ...) {
      calls <<- calls + 1L
      if (calls == 1L) stop("falha deliberada")
      fit$graph$converged <- TRUE
      edges <- if (calls == 2L) matrix(c(1L, 2L), ncol = 2L) else
        matrix(integer(), ncol = 2L)
      structure(list(edges = edges, fit = fit), class = "mcg_selection")
    },
    .package = "mcglmgraph"
  )
  st <- mc_graph_stability(fake, B = 3L, pi_thr = 0.4)
  expect_equal(st$B_valid, 2L)
  expect_equal(st$B_failed, 1L)
  expect_equal(st$selection_prob[1, 2], 0.5)
  expect_equal(st$selection_prob_conservative[1, 2], 1 / 3)
  expect_equal(nrow(st$failures), 1L)
})

test_that("penalized fitting preserves positive definiteness and decreases its objective", {
  skip_if_not_installed("MASS")
  d <- make_chain()
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
  pen <- mc_graph_penalized(f, lambda = 1, max_iter = 100)
  expect_true(mc_check_spd(pen$K_b)$is_spd)
  expect_true(all(diff(pen$objective_path) <=
                    1e-8 * (1 + abs(head(pen$objective_path, -1)))))
  expect_true(is.finite(pen$objective))
  expect_true(pen$status %in% c("converged", "maximum_iterations",
                                "objective_line_search_failed",
                                "kkt_not_certified"))
})

test_that("robust Wald test runs and is calibrated under the Gaussian model", {
  skip_if_not_installed("MASS")
  d <- make_chain()
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
  t_model <- mc_test_independence(f)
  t_robust <- mc_test_independence(f, robust = TRUE)
  expect_s3_class(t_robust, "htest")
  expect_true(is.finite(t_robust$statistic))
  # under correct Gaussian specification V ~ J, so both statistics are close
  expect_lt(abs(log(t_robust$statistic / t_model$statistic)), 0.5)
})
