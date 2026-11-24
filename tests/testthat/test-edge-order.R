# The edge parameters of a fit are stored in canonical order (i < j, rows in
# lexicographic order). Every consumer of the support must read it in that same
# order, so nothing may depend on the order in which the edges were supplied.

edge_order_data <- function(n = 300, seed = 3) {
  set.seed(seed)
  z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
  x <- stats::rnorm(n)
  y3 <- 0.6 * z[, 1] + 0.5 * z[, 2] + stats::rnorm(n)
  y4 <- 0.5 * y3 + stats::rnorm(n)
  data.frame(y1 = z[, 1] + 0.3 * x, y2 = z[, 2], y3 = y3, y4 = y4, x = x)
}
edge_order_lp <- list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x)
fit_with_spec <- function(d, spec, ...)
  suppressMessages(mcglm_graph(edge_order_lp, between_graph = spec, data = d,
                               ...))

test_that("the specification stores its edges in canonical order", {
  sp <- mc_between_graph("ggm", rbind(c(4, 3), c(2, 1), c(3, 2)))
  expect_identical(sp$edges, rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L)))
  sp <- mc_between_graph("dag", rbind(c(3, 1), c(2, 1)), ordering = c(3, 2, 1))
  expect_identical(sp$edges, rbind(c(1L, 2L), c(1L, 3L)))
  expect_null(mc_between_graph("ggm")$edges)
  expect_error(mc_between_graph("ggm", rbind(c(1, 2), c(2, 1))), "duplicated")
  expect_error(mc_between_graph("dag", rbind(c(1, 4)), ordering = 1:3),
               "out of range")
})

test_that("GGM edge scores, Wald tests and criteria ignore the input order", {
  skip_if_not_installed("MASS")
  d <- edge_order_data()
  inputs <- list(sorted = rbind(c(1, 2), c(2, 3), c(3, 4)),
                 permuted = rbind(c(3, 4), c(1, 2), c(2, 3)),
                 reversed = rbind(c(4, 3), c(2, 1), c(3, 2)))
  fits <- lapply(inputs, function(e) fit_with_spec(d, mc_between_graph("ggm", e)))
  # A specification built by hand, bypassing mc_between_graph(), is
  # canonicalized by the fit as well.
  raw <- structure(list(type = "ggm", edges = inputs$permuted, ordering = NULL),
                   class = "mcg_graph_spec")
  fits$raw <- fit_with_spec(d, raw)
  ref <- fits$sorted
  s_ref <- mc_edge_scores(ref)
  # The Wald statistic of (3,4) is about 74.7 here; reading another edge's
  # parameter gives a different value, so equality is a real check.
  expect_gt(s_ref$stat[s_ref$i == 3 & s_ref$j == 4], 70)
  for (nm in names(fits)[-1L]) {
    f <- fits[[nm]]
    expect_identical(f$graph$spec$edges, ref$graph$spec$edges, info = nm)
    expect_equal(f$graph$par, ref$graph$par, tolerance = 1e-12, info = nm)
    for (w in c("both", "wald", "score"))
      expect_equal(mc_edge_scores(f, w), mc_edge_scores(ref, w),
                   tolerance = 1e-10, info = paste(nm, w))
    expect_equal(mc_edge_scores(f, "score", variability = "empirical"),
                 mc_edge_scores(ref, "score", variability = "empirical"),
                 tolerance = 1e-10, info = nm)
    expect_equal(mc_vcov_graph(f, "sandwich"), mc_vcov_graph(ref, "sandwich"),
                 tolerance = 1e-12, info = nm)
    expect_equal(mc_wald_graph(f, "eta_34")$statistic,
                 mc_wald_graph(ref, "eta_34")$statistic, tolerance = 1e-12,
                 info = nm)
    expect_equal(mc_test_independence(f)$statistic,
                 mc_test_independence(ref)$statistic, tolerance = 1e-12,
                 info = nm)
    for (cr in c("aic", "bic", "qic", "gbic"))
      expect_equal(mc_compute_criterion(f, cr), mc_compute_criterion(ref, cr),
                   tolerance = 1e-12, info = paste(nm, cr))
  }
})

test_that("the Wald statistic of an edge reads that edge's own parameter", {
  skip_if_not_installed("MASS")
  d <- edge_order_data()
  f <- fit_with_spec(d, mc_between_graph("ggm", rbind(c(3, 4), c(1, 2), c(2, 3))))
  s <- mc_edge_scores(f, "wald")
  V <- f$graph$vcov
  for (k in seq_len(nrow(s))) {
    p <- match(paste0("eta_", s$i[k], s$j[k]), names(f$graph$par))
    expect_equal(s$stat[k], unname(f$graph$par[p]^2 / V[p, p]),
                 tolerance = 1e-12)
  }
  pos <- getFromNamespace(".edge_param_pos", "mcglmgraph")
  expect_identical(names(f$graph$par)[pos(f, 4, 3)], "eta_34")
  expect_error(pos(f, 1, 4), "not in the support")
})

test_that("forward and bidirectional stepwise rank candidates at the fitted point", {
  skip_if_not_installed("MASS")
  d <- edge_order_data()
  f <- fit_with_spec(d, mc_between_graph("ggm"))
  refit <- getFromNamespace(".refit_between", "mcglmgraph")
  for (dir in c("forward", "both")) {
    sel <- mc_graph_stepwise(f, direction = dir, criterion = "bic",
                             k_confirm = 1L, trace = TRUE)
    tr <- sel$trace
    expect_true(nrow(tr) >= 2L, info = dir)
    # the stored support is canonical
    ed <- sel$fit$graph$spec$edges
    expect_identical(ed, ed[order(ed[, 1L], ed[, 2L]), , drop = FALSE],
                     info = dir)
    # Replay: with k_confirm = 1 every accepted move is the best-ranked one
    # under the scores of a fresh fit of the current support, given sorted.
    # "forward" starts from the empty graph, "both" from the complete one
    cur <- if (dir == "forward") matrix(integer(0), 0, 2) else
      t(utils::combn(4L, 2L))
    for (k in seq_len(nrow(tr))) {
      cf <- refit(f, mc_between_graph("ggm", cur))
      e <- c(tr$edge_i[k], tr$edge_j[k])
      if (tr$action[k] == "add") {
        sc <- mc_edge_scores(cf, "score")
        top <- sc[which.max(sc$stat), ]
        cand_add <- c(top$i, top$j)
        # with direction = "both" the move may also be a removal; an add
        # must be the top-ranked absent edge
        expect_identical(as.integer(e), as.integer(cand_add),
                         info = paste(dir, "step", k))
        cur <- rbind(cur, e)
      } else {
        wd <- mc_edge_scores(cf, "wald")
        low <- wd[which.min(wd$stat), ]
        expect_identical(as.integer(e), as.integer(c(low$i, low$j)),
                         info = paste(dir, "step", k))
        cur <- cur[!(cur[, 1L] == e[1L] & cur[, 2L] == e[2L]), , drop = FALSE]
      }
      cur <- cur[order(cur[, 1L], cur[, 2L]), , drop = FALSE]
      expect_equal(tr$criterion[k],
                   as.numeric(mc_compute_criterion(
                     refit(f, mc_between_graph("ggm", cur)), "bic")),
                   tolerance = 1e-10, info = paste(dir, "step", k))
    }
    # the final fit scores exactly like a fresh fit of its sorted support
    fresh <- refit(f, mc_between_graph("ggm", cur))
    expect_equal(mc_edge_scores(sel$fit), mc_edge_scores(fresh),
                 tolerance = 1e-10, info = dir)
  }
})

test_that("a DAG support given in another row order fits the same layer", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  d <- edge_order_data()
  ord <- c(2L, 1L, 3L, 4L)
  inputs <- list(sorted = rbind(c(1, 3), c(2, 3), c(3, 4)),
                 permuted = rbind(c(3, 4), c(2, 3), c(1, 3)),
                 reversed = rbind(c(4, 3), c(3, 1), c(3, 2)))
  fits <- lapply(inputs, function(e)
    fit_with_spec(d, mc_between_graph("dag", e, ordering = ord)))
  ref <- fits$sorted
  for (nm in names(fits)[-1L]) {
    f <- fits[[nm]]
    expect_identical(f$graph$spec$edges, ref$graph$spec$edges, info = nm)
    expect_identical(names(f$graph$par), names(ref$graph$par), info = nm)
    expect_equal(f$graph$par, ref$graph$par, tolerance = 1e-12, info = nm)
    expect_equal(mc_vcov_graph(f, "sandwich"), mc_vcov_graph(ref, "sandwich"),
                 tolerance = 1e-12, info = nm)
    expect_equal(mc_vcov_graph(f, "two_stage"), mc_vcov_graph(ref, "two_stage"),
                 tolerance = 1e-10, info = nm)
    expect_equal(mc_dag_coefficients(f)$edges, mc_dag_coefficients(ref)$edges,
                 tolerance = 1e-12, info = nm)
    expect_equal(mc_compute_criterion(f, "gbic"),
                 mc_compute_criterion(ref, "gbic"), tolerance = 1e-12,
                 info = nm)
  }
})
