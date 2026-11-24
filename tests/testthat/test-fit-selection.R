make_chain_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  K <- diag(4) * 2
  K[1,2] <- K[2,1] <- -0.6; K[2,3] <- K[3,2] <- -0.6; K[3,4] <- K[4,3] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
  x <- rnorm(n)
  data.frame(y1 = Y[,1] + 0.5 * x, y2 = Y[,2], y3 = Y[,3],
             y4 = Y[,4] - 0.3 * x, x = x)
}

test_that("complete GGM reproduces the sample covariance and converges", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data()
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  expect_true(fit$graph$converged)
  expect_identical(fit$graph$status, "converged")
  expect_true(is.list(fit$graph$meta))
  expect_true(length(fit$graph$iter) == 1L &&
                is.finite(fit$graph$iter) && fit$graph$iter >= 1L)
  expect_lt(max(abs(fit$graph$Sigma_b - fit$graph$Shat)), 1e-6)
  sparse <- mcglmgraph:::.refit_between(
    fit, mc_between_graph("ggm", rbind(c(1, 2), c(2, 3), c(3, 4))))
  expect_true(sparse$graph$converged)
  expect_identical(sparse$graph$status, "converged")
  expect_equal(nrow(sparse$graph$meta$bases$edges), 3L)
  expect_equal(length(sparse$graph$edge_pos), 3L)
})

test_that("the fit is a genuine mcglm object with working methods", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data()
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  expect_s3_class(fit, "mcglm")            # inherits mcglm
  expect_false(is.null(coef(fit)))         # inherited mcglm coef works
  expect_false(is.null(vcov(fit)))         # inherited mcglm vcov works
  expect_equal(ncol(as.matrix(residuals(fit))), 4L)
  expect_false(is.null(mcglm::gof(fit)))   # inherited mcglm gof works
})

test_that("regression coefficients are invariant to the graph layer", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data()
  fg <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x),
                    between_graph = mc_between_graph("ggm"), data = dat)
  fd <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x),
                    between_graph = mc_between_graph("dag"), data = dat)
  expect_equal(coef(fg, "regression"), coef(fd, "regression"), tolerance = 1e-6)
})

test_that("stepwise selection recovers the chain structure", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data()
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  sel <- mc_graph_stepwise(fit, criterion = "bic")
  truth <- rbind(c(1,2), c(2,3), c(3,4))
  got <- sel$edges[order(sel$edges[,1], sel$edges[,2]), , drop = FALSE]
  expect_equal(unname(got), unname(truth))
})

test_that("GGM and DAG yield the same saturated covariance", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data()
  fg <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x),
                    between_graph = mc_between_graph("ggm"), data = dat)
  fd <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x),
                    between_graph = mc_between_graph("dag"), data = dat)
  expect_equal(fg$graph$Sigma_b, fd$graph$Sigma_b, tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("partial correlations recover second-order conditional separation", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data(n = 600)
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  P <- mc_partial_cor(fit)
  expect_gt(abs(P[1,2]), abs(P[1,3]))   # adjacent stronger than mediated
  expect_gt(abs(P[3,4]), abs(P[1,4]))
})

test_that("lambda selection accepts only converged path candidates", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data(n = 400)
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  sel <- mc_select_lambda(fit, n_lambda = 5L, criterion = "gbic")
  accepted <- is.finite(sel$lambda_path$criterion)
  expect_true(any(accepted))
  expect_true(all(sel$lambda_path$penalized_converged[accepted]))
  expect_true(all(sel$lambda_path$refit_converged[accepted]))
})

test_that("stepwise reports and rejects non-converged candidate refits", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data(n = 300)
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  sel <- mc_graph_stepwise(fit, criterion = "gbic", trace = FALSE)
  expect_true(isTRUE(sel$fit$graph$converged))
  expect_s3_class(sel$candidate_failures, "data.frame")
  expect_named(sel$candidate_failures,
               c("step", "action", "edge_i", "edge_j", "status"))
})

test_that("cluster sandwich aggregates Pearson contributions by cluster", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data(n = 300)
  fit <- mcglm_graph(list(y1~x, y2~x, y3~x, y4~x), data = dat)
  V_iid <- mcglmgraph:::.robust_vcov(fit)
  V_unique <- mcglmgraph:::.robust_vcov(
    fit, cluster = seq_len(nrow(dat)), small_sample = FALSE
  )
  expect_equal(V_unique, V_iid, tolerance = 1e-12)
  cl <- rep(seq_len(60), each = 5)
  V_cluster <- mcglmgraph:::.robust_vcov(
    fit, cluster = cl, small_sample = TRUE
  )
  expect_true(all(is.finite(V_cluster)))
  expect_equal(V_cluster, t(V_cluster), tolerance = 1e-12)
  expect_error(
    mc_test_independence(fit, cluster = cl),
    "robust = TRUE"
  )
  expect_s3_class(
    mc_test_independence(
      fit, robust = TRUE, cluster = cl, small_sample = TRUE
    ),
    "htest"
  )
})

test_that("stepwise refuses a directed fit instead of converting it", {
  skip_if_not_installed("MASS")
  dat <- make_chain_data()
  lp <- list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x)
  dag <- mcglm_graph(lp, between_graph = mc_between_graph("dag"), data = dat)
  expect_error(mc_graph_stepwise(dag, criterion = "bic"),
               "searches over supports of the GGM layer")
  # the undirected fit of the same data is served as before
  ggm <- mcglm_graph(lp, data = dat)
  expect_s3_class(mc_graph_stepwise(ggm, criterion = "bic", trace = FALSE),
                  "mcg_selection")
})
