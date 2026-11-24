test_that("sparse DAG inference uses the observed curvature under misspecification", {
  skip_if_not_installed("numDeriv")
  S <- matrix(c(1, 0.5, 0.7,
                0.5, 1, 0.7,
                0.7, 0.7, 2), 3, 3)
  edges <- rbind(c(1L, 3L), c(2L, 3L))
  spec <- mc_between_graph("dag", edges = edges)
  n <- 500
  bfit <- mcglmgraph:::.fit_between(S, n, spec)
  fit <- list(graph = c(list(type = "dag", Shat = S, n = n, R = 3L), bfit))
  H <- mcglmgraph:::.graph_inference_bread(fit)
  builder <- mcglmgraph:::.make_builder(spec, 3L)
  criterion <- function(par) {
    Sigma <- builder(par)$Sigma
    ld <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
    -n / 2 * (ld + sum(solve(Sigma) * S))
  }
  H_numeric <- -numDeriv::hessian(criterion, as.numeric(bfit$par))
  expect_equal(unname(H), unname(H_numeric), tolerance = 1e-4)
  expect_equal(H[1:2, 1:2], n * S[1:2, 1:2] / bfit$par[5],
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(unname(bfit$J[1, 2]), 0, tolerance = 1e-12)
  expect_gt(abs(H[1, 2] - bfit$J[1, 2]), 100)
})

test_that("DAG sandwich and Godambe criteria use the observed bread", {
  skip_if_not_installed("MASS")
  set.seed(6101)
  S <- matrix(c(1, 0.5, 0.7,
                0.5, 1, 0.7,
                0.7, 0.7, 2), 3, 3)
  Y <- MASS::mvrnorm(400, rep(0, 3), S)
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3])
  edges <- rbind(c(1L, 3L), c(2L, 3L))
  f <- mcglm_graph(list(y1 ~ 1, y2 ~ 1, y3 ~ 1),
                   between_graph = mc_between_graph("dag", edges = edges),
                   data = d)
  g <- f$graph
  expect_true(g$converged)
  H <- mcglmgraph:::.graph_inference_bread(f)
  Psi <- mcglmgraph:::.pearson_scores(g$residuals, g$Sigma_b, g$Dlist)
  Hinv <- solve(H)
  V_manual <- Hinv %*% crossprod(Psi) %*% Hinv
  V <- mc_vcov_graph(f, type = "sandwich")
  expect_equal(as.numeric(V), as.numeric(V_manual), tolerance = 1e-8)
  Ec <- scale(g$residuals, center = TRUE, scale = FALSE)
  corrected <- mcglmgraph:::cpp_two_stage_scores(
    Psi, as.matrix(g$K_b), g$Dlist, Ec)$scores
  V_two_stage <- mc_vcov_graph(f, type = "two_stage", method = "shortcut")
  expect_equal(as.numeric(V_two_stage),
               as.numeric(Hinv %*% crossprod(corrected) %*% Hinv),
               tolerance = 1e-8)

  cl <- rep(seq_len(g$n / 2L), each = 2L)
  Pc <- rowsum(Psi, cl, reorder = FALSE)
  Vunit <- crossprod(Pc) / g$n
  penalty <- sum(diag(solve(H / g$n, Vunit)))
  crit <- mc_compute_criterion(f, c("qic", "gbic", "bic"), cluster = cl)
  expect_equal(attr(crit, "effective_dimension"), penalty, tolerance = 1e-8)
  ll <- mc_graph_ploglik(f)
  expect_equal(unname(crit["qic"] + 2 * ll), 2 * penalty, tolerance = 1e-8)
  expect_equal(unname(crit["gbic"] + 2 * ll),
               log(length(unique(cl))) * penalty, tolerance = 1e-8)
  expect_equal(unname(crit["bic"] + 2 * ll),
               log(length(unique(cl))) * g$npar, tolerance = 1e-8)
})

test_that("clustered Godambe criteria pass through both GGM selectors", {
  set.seed(6102)
  n <- 160
  d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n))
  f <- mcglm_graph(list(y1 ~ 1, y2 ~ 1, y3 ~ 1), data = d)
  cl <- rep(seq_len(n / 2L), each = 2L)
  step <- mc_graph_stepwise(f, criterion = "gbic", cluster = cl,
                            max_steps = 1L)
  expect_s3_class(step, "mcg_selection")
  path <- mc_select_lambda(f, criterion = "gbic", cluster = cl,
                           n_lambda = 1L, include_zero = TRUE)
  full <- mcglmgraph:::.refit_between(f, mc_between_graph("ggm"))
  expect_equal(tail(path$lambda_path$criterion, 1),
               as.numeric(mc_compute_criterion(full, "gbic", cluster = cl)),
               tolerance = 1e-8)
})
