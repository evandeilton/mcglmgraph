# Resampling that refits both stages.
#
# The bootstrap is the yardstick for the analytical correction: refitting
# propagates the first stage by construction, so its standard errors are what
# mc_vcov_graph(type = "two_stage") approximates to first order and what the
# naive sandwich, which conditions on the residuals, does not.

make_fit <- function(n = 250, seed = 11, k12 = -0.8, k23 = -0.7) {
  set.seed(seed)
  K <- diag(3) * 2
  K[1, 2] <- K[2, 1] <- k12
  K[2, 3] <- K[3, 2] <- k23
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

test_that("refitting both stages reproduces the analytical correction", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 300, seed = 23)
  set.seed(5)
  bt <- mc_graph_bootstrap(f, B = 150)
  q <- length(f$graph$par)
  se_boot <- bt$se[seq_len(q)]
  se_two <- sqrt(diag(mc_vcov_graph(f, "two_stage")))
  se_naive <- sqrt(diag(mc_vcov_graph(f, "sandwich")))

  # the corrected form tracks the bootstrap and the naive one does not
  r_two <- se_boot / se_two
  expect_true(all(r_two > 0.7 & r_two < 1.4))
  expect_lt(mean(abs(log(r_two))), mean(abs(log(se_boot / se_naive))))

  # where the two analytical forms disagree, the bootstrap sides with the
  # corrected one. On a diagonal coefficient the ratio of the two is
  # sqrt(2 (K_kk - 1) / K_kk), below one when K_kk < 2 and above it when
  # K_kk > 2, so the direction is read off the fit rather than assumed.
  K <- as.matrix(f$graph$K_b)
  pred <- sqrt(2 * (diag(K) - 1) / diag(K))
  # the formula is asymptotic and both sides here carry the sampling noise of
  # an empirical fourth moment at this sample size; the tight check on it lives
  # in the two-stage tests, at n = 20000
  obs <- (se_two / se_naive)[1:3]
  expect_equal(unname(obs), unname(pred), tolerance = 0.15)
  moves <- which(abs(log(pred)) > 0.2)            # where it matters
  expect_gt(length(moves), 0L)
  for (k in moves)
    expect_lt(abs(log(se_boot[k] / se_two[k])),
              abs(log(se_boot[k] / se_naive[k])))
})

test_that("resampling the residuals alone reproduces the naive sandwich", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 300, seed = 23)
  set.seed(7)
  br <- mc_graph_bootstrap(f, B = 150, stage = "residual")
  q <- length(f$graph$par)
  r <- br$se[seq_len(q)] / sqrt(diag(mc_vcov_graph(f, "sandwich")))
  expect_true(all(r > 0.75 & r < 1.3))
  expect_identical(br$stage, "residual")
  # and it is visibly different from the two-stage bootstrap on the diagonal
  set.seed(7)
  bt <- mc_graph_bootstrap(f, B = 150)
  expect_lt(bt$se[1], br$se[1])
})

test_that("the resamples do not depend on the worker count", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 150)
  set.seed(99); a <- mc_graph_bootstrap(f, B = 20, workers = 1L)
  set.seed(99); b <- mc_graph_bootstrap(f, B = 20, workers = 2L)
  expect_equal(a$replicates, b$replicates, tolerance = 1e-12)
  expect_equal(a$se, b$se, tolerance = 1e-12)
})

test_that("the returned object holds what it claims", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 150)
  set.seed(3)
  bt <- mc_graph_bootstrap(f, B = 30, level = 0.9)
  expect_s3_class(bt, "mcg_bootstrap")
  nm <- c(names(f$graph$par), "pcor_y1_y2", "pcor_y1_y3", "pcor_y2_y3")
  expect_identical(names(bt$estimate), nm)
  expect_identical(colnames(bt$replicates), nm)
  expect_equal(nrow(bt$replicates), bt$B_valid)
  expect_equal(bt$B_valid + bt$B_failed, 30L)
  expect_equal(unname(bt$se), unname(apply(bt$replicates, 2, sd)),
               tolerance = 1e-12)
  expect_equal(unname(bt$bias),
               unname(colMeans(bt$replicates) - bt$estimate), tolerance = 1e-12)
  expect_equal(unname(bt$vcov), unname(cov(bt$replicates)), tolerance = 1e-12)
  # the percentile interval is the empirical quantile pair
  expect_equal(unname(bt$ci[1, ]),
               unname(quantile(bt$replicates[, 1], c(0.05, 0.95), names = FALSE)),
               tolerance = 1e-12)
  # the point estimate is the one of the original fit
  expect_equal(unname(bt$estimate[seq_along(f$graph$par)]),
               unname(f$graph$par), tolerance = 1e-12)
  txt <- capture.output(print(bt))
  expect_true(any(grepl("both stages refitted", txt)))
})

test_that("the support of the graph is held fixed across resamples", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200)
  sub <- .refit_between(f, mc_between_graph("ggm", rbind(c(1, 2), c(2, 3))))
  set.seed(4)
  bt <- mc_graph_bootstrap(sub, B = 20)
  # five parameters, three diagonal and two edges, in every resample
  expect_equal(ncol(bt$replicates), length(sub$graph$par) + 3L)
  expect_identical(names(bt$estimate)[seq_along(sub$graph$par)],
                   names(sub$graph$par))
})

test_that("whole clusters are resampled when the unit of dependence is one", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200)
  cl <- rep(seq_len(100), each = 2L)
  set.seed(8)
  bt <- mc_graph_bootstrap(f, B = 30, cluster = cl)
  expect_true(bt$clustered)
  expect_identical(bt$n_clusters, 100L)
  expect_equal(nrow(bt$replicates), bt$B_valid)
  expect_error(mc_graph_bootstrap(f, B = 5, cluster = 1:3),
               "one identifier per observational row")
  expect_error(mc_graph_bootstrap(f, B = 5, cluster = rep(1L, 200)),
               "at least two independent clusters")
})

test_that("a custom statistic is honored", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- make_fit(n = 150)
  set.seed(6)
  bt <- mc_graph_bootstrap(f, B = 20, statistic = function(z)
    c(logdet = as.numeric(determinant(z$graph$Sigma_b,
                                      logarithm = TRUE)$modulus)))
  expect_identical(colnames(bt$replicates), "logdet")
  expect_equal(unname(bt$estimate),
               as.numeric(determinant(f$graph$Sigma_b,
                                      logarithm = TRUE)$modulus),
               tolerance = 1e-12)
  expect_error(mc_graph_bootstrap(f, B = 5, statistic = function(z) 1),
               "named numeric vector")
  expect_error(mc_graph_bootstrap(f, B = 5, statistic = "mean"),
               "must be a function")
})

test_that("a fit whose rows are not exchangeable is resampled by clusters only", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 120)
  # a within-response structure that is not the identity describes dependence
  # among the very rows the resampling would permute
  g <- f
  M <- diag(120); M[1, 2] <- M[2, 1] <- 0.5
  g$matrix_pred[[1]] <- list(M)
  expect_error(mc_graph_bootstrap(g, B = 5), "needs cluster")
  # a cluster that separates rows the structure joins would cut it
  split_12 <- c(1L, 2L, rep(3:61, each = 2L))
  expect_error(mc_graph_bootstrap(g, B = 5, cluster = split_12),
               "connects rows of different clusters")
  expect_error(mc_cpdag_stability(g, B = 5, refit = "full",
                                  cluster = split_12),
               "connects rows of different clusters")
  # the residual route stays available, since it does not touch the rows
  expect_s3_class(mc_graph_bootstrap(g, B = 5, stage = "residual"),
                  "mcg_bootstrap")
})

test_that("a cluster drawn twice gives two blocks that do not touch", {
  # the structure of a resample, rebuilt in plain R: each drawn block keeps its
  # own submatrix and nothing joins two copies of the same cluster
  Z <- matrix(0, 5, 5)
  Z[1:3, 1:3] <- 0.3; Z[4:5, 4:5] <- 0.6; diag(Z) <- 1
  blocks <- list(1:3, 4:5, 1:3)
  got <- as.matrix(mcglmgraph:::.resample_structure(Z, blocks))
  want <- matrix(0, 8, 8)
  want[1:3, 1:3] <- Z[1:3, 1:3]
  want[4:5, 4:5] <- Z[4:5, 4:5]
  want[6:8, 6:8] <- Z[1:3, 1:3]
  expect_equal(unname(got), want, tolerance = 0)
  # the submatrix of the rows drawn would join the two copies instead
  rows <- unlist(blocks)
  expect_gt(max(abs(Z[rows, rows][1:3, 6:8])), 0)
  expect_equal(max(abs(got[1:3, 6:8])), 0)
  # on the identity the rebuilt structure is the identity of the resample
  I5 <- diag(5)
  expect_equal(unname(as.matrix(mcglmgraph:::.resample_structure(I5, blocks))),
               diag(8), tolerance = 0)
})

test_that("the unit draws follow a fixed sequence of calls to the generator", {
  n <- 40
  by_cluster <- split(seq_len(n), rep(1:10, each = 4))
  set.seed(3)
  old <- lapply(1:6, function(b)
    unlist(by_cluster[sample.int(10, 10, replace = TRUE)], use.names = FALSE))
  set.seed(3)
  new <- mcglmgraph:::.draw_units(6, n, by_cluster)
  expect_identical(lapply(new, `[[`, "rows"), old)
  expect_identical(lapply(new, function(z) unlist(z$blocks, use.names = FALSE)),
                   old)
  set.seed(4)
  old <- lapply(1:6, function(b) sample.int(n, n, replace = TRUE))
  set.seed(4)
  new <- mcglmgraph:::.draw_units(6, n)
  expect_identical(lapply(new, `[[`, "rows"), old)
  expect_true(all(vapply(new, function(z) is.null(z$blocks), logical(1))))
})

test_that("random-effect blocks nested in clusters are resampled whole", {
  skip_if_not_installed("MASS")
  set.seed(99); G <- 20; m <- 4
  cl <- rep(seq_len(G), each = m); blk <- factor(rep(seq_len(G * m / 2), each = 2))
  n <- G * m; x <- stats::runif(n)
  b <- stats::rnorm(G * m / 2, sd = 0.4)[blk]
  Rl <- matrix(c(1, .5, .5, 1), 2)
  U <- stats::pnorm(MASS::mvrnorm(n, rep(0, 2), Rl))
  mu <- exp(cbind(1 + .5 * x + b, .5 - .3 * x + b))
  d <- data.frame(y1 = stats::qpois(U[, 1], mu[, 1]),
                  y2 = stats::qpois(U[, 2], mu[, 2]), x = x, cl = cl, blk = blk)
  Z <- c(mcglm::mc_id(d), mcglm::mc_mixed(~ 0 + blk, data = d))
  f <- suppressMessages(mcglm_graph(
    list(y1 ~ x, y2 ~ x), data = d, matrix_pred = rep(list(Z), 2),
    link = rep("log", 2), variance = rep("poisson_tweedie", 2),
    power_fixed = rep(TRUE, 2),
    control_initial = list(regression = list(c(1, .5), c(.5, -.3)),
                           power = list(1, 1),
                           tau = list(c(.2, .1), c(.2, .1)), rho = .4)))
  expect_error(mc_graph_bootstrap(f, B = 3), "needs cluster")
  set.seed(5)
  bt <- mc_graph_bootstrap(f, B = 3, cluster = d$cl)
  expect_true(bt$clustered)
  expect_identical(bt$n_clusters, as.integer(G))
  expect_gte(bt$B_valid, 1L)
  # the selection frequencies draw half of the clusters, whole
  set.seed(6)
  st <- mc_graph_stability(f, B = 3, refit = "full", cluster = d$cl)
  expect_true(st$clustered)
  expect_equal(st$B_valid + st$B_failed, 3L)
  expect_error(mc_graph_stability(f, B = 3, refit = "full"), "needs cluster")
})

test_that("the contracts of the arguments are enforced", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 120)
  expect_error(mc_graph_bootstrap(f, B = 1), "at least two")
  expect_error(mc_graph_bootstrap(f, B = 5, level = 1), "strictly between")
  expect_error(mc_graph_bootstrap(f, B = 5, stage = "first"), "arg")
  expect_error(mc_graph_bootstrap(structure(list(), class = "lm")), "inherits")
})

test_that("the DAG layer is bootstrapped through the same route", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(15); n <- 250; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag"), data = d)
  set.seed(2)
  bt <- mc_graph_bootstrap(f, B = 30)
  expect_equal(ncol(bt$replicates), length(f$graph$par) + 3L)
  expect_true(all(bt$se > 0))
  expect_equal(bt$B_failed, 0L)
})

test_that("a resample carries its own first stage", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200, seed = 31)
  ns <- asNamespace("mcglmgraph")
  set.seed(8)
  rows <- sample.int(200, 200, replace = TRUE)
  fb <- ns$.refit_on_rows(f, rows, "both", f$graph$spec, ns$.graph_control(f))
  ref <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = f$data[rows, ])

  # both stages are those of the resample, as a fit on its rows gives them
  expect_s3_class(fb, c("mcglm_graph", "mcglm"), exact = TRUE)
  expect_false(isTRUE(all.equal(fb$Regression, f$Regression)))
  expect_equal(fb$Regression, ref$Regression, tolerance = 1e-6)
  # the two first stages start from different values, so they agree to the
  # tolerance of the first-stage algorithm
  expect_equal(fb$Covariance, ref$Covariance, tolerance = 1e-4)
  expect_equal(fb$graph$par, ref$graph$par, tolerance = 1e-6)
  expect_identical(names(fb$graph), names(ref$graph))

  # so a statistic of the first stage has a spread
  set.seed(9)
  bt <- mc_graph_bootstrap(f, B = 20, statistic = function(z)
    c(slope = unname(z$Regression[2L]), eta = unname(z$graph$par[1L])))
  expect_gt(bt$se[["slope"]], 0)
  expect_gt(bt$se[["eta"]], 0)
  # of the size of the model-based standard error of that slope
  se_model <- sqrt(diag(as.matrix(stats::vcov(f))))[2L]
  expect_lt(abs(log(bt$se[["slope"]] / se_model)), log(2))
})

test_that("every route builds the same graph block", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200, seed = 32)
  ns <- asNamespace("mcglmgraph")
  rf <- ns$.refit_between(f, f$graph$spec)
  expect_identical(names(rf$graph), names(f$graph))
  expect_equal(rf$graph, f$graph)
  rows <- seq_len(200)
  rr <- ns$.refit_on_rows(f, rows, "residual", f$graph$spec,
                          ns$.graph_control(f))
  expect_identical(names(rr$graph), names(f$graph))
  expect_equal(rr$graph, f$graph)
  # a sparser support refreshes the diagnostics along with the parameters
  sp <- ns$.refit_on_rows(f, seq(1L, 200L, by = 2L), "residual",
                          mc_between_graph("ggm", rbind(c(1, 2))),
                          ns$.graph_control(f))
  expect_identical(sp$graph$n, 100L)
  expect_identical(sp$graph$npar, 4L)
  expect_identical(length(sp$graph$convergence_gates), 3L)
  expect_identical(sp$graph$meta$E, 1L)
})
