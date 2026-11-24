# Extended BIC, warm starts and the unpenalized end of the L1 path, and
# stability selection with error control.
#
# Targets, in plain R:
#   EBIC      -2 l + log(n) npar + 4 gamma |E| log R (Foygel and Drton, 2010)
#   path      from the diagonal model, the first edge to enter the L1 path is the
#             one with the largest |S_ij|, since the KKT condition of an edge at
#             zero is n |Sigma_ij - S_ij| <= lambda and Sigma is diagonal there
#   bound     E(V) <= q^2 / ((2 pi - 1) |E_max|) (Meinshausen and Buehlmann,
#             2010, Theorem 1)

sel_fit <- function(n = 300, seed = 1, R = 4) {
  set.seed(seed)
  K <- diag(R) * 2
  for (k in seq_len(R - 1L)) K[k, k + 1L] <- K[k + 1L, k] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, R), solve(K))
  d <- as.data.frame(Y); names(d) <- paste0("y", seq_len(R))
  d$x <- rnorm(n)
  mcglm_graph(lapply(paste0("y", seq_len(R)),
                     function(v) stats::reformulate("x", response = v)),
              data = d)
}

## ---- EBIC -------------------------------------------------------------------

test_that("EBIC is BIC plus 4 gamma |E| log R", {
  skip_if_not_installed("MASS")
  f <- sel_fit()
  g <- f$graph
  ll <- -0.5 * g$n * (g$R * log(2 * pi) +
                        as.numeric(determinant(g$Sigma_b)$modulus) +
                        sum(solve(g$Sigma_b) * g$Shat))
  for (gam in c(0, 0.25, 0.5, 1)) {
    want <- -2 * ll + log(g$n) * g$npar + 4 * gam * length(g$edge_pos) * log(g$R)
    expect_equal(as.numeric(mc_compute_criterion(f, "ebic", ebic_gamma = gam)),
                 want, tolerance = 1e-10)
  }
  expect_identical(unname(mc_compute_criterion(f, "ebic", ebic_gamma = 0)),
                   unname(mc_compute_criterion(f, "bic")))
  # without `which` the output is the one it was before EBIC existed
  expect_identical(names(mc_compute_criterion(f)),
                   c("aic", "bic", "qic", "gbic", "clbic"))
  expect_error(mc_compute_criterion(f, "ebic", ebic_gamma = -1), "ebic_gamma")
})

test_that("the selectors accept EBIC and reduce to BIC at gamma zero", {
  skip_if_not_installed("MASS")
  f <- sel_fit(R = 5)
  a <- mc_graph_stepwise(f, criterion = "bic", trace = FALSE)
  b <- mc_graph_stepwise(f, criterion = "ebic", ebic_gamma = 0, trace = FALSE)
  expect_identical(a$edges, b$edges)
  big <- mc_graph_stepwise(f, criterion = "ebic", ebic_gamma = 5, trace = FALSE)
  expect_lte(big$n_edges, a$n_edges)
  l0 <- mc_select_lambda(f, criterion = "bic", n_lambda = 12)
  l1 <- mc_select_lambda(f, criterion = "ebic", ebic_gamma = 0, n_lambda = 12)
  expect_identical(l0$edges, l1$edges)
  expect_equal(l0$lambda_path$criterion, l1$lambda_path$criterion,
               tolerance = 1e-12)
})

## ---- warm starts and lambda = 0 ---------------------------------------------

test_that("warm starts reach the same certified supports as independent fits", {
  skip_if_not_installed("MASS")
  f <- sel_fit(R = 5, seed = 3)
  cold <- mc_select_lambda(f, n_lambda = 15)
  warm <- mc_select_lambda(f, n_lambda = 15, warm_start = TRUE)
  expect_true(warm$warm_start)
  both <- cold$lambda_path$penalized_converged & warm$lambda_path$penalized_converged
  expect_true(any(both))
  expect_identical(cold$lambda_path$n_edges[both], warm$lambda_path$n_edges[both])
  expect_equal(cold$lambda_path$criterion[both], warm$lambda_path$criterion[both],
               tolerance = 1e-8)
  expect_identical(cold$edges, warm$edges)
  expect_warning(mc_select_lambda(f, n_lambda = 5, warm_start = TRUE,
                                  workers = 2L), "sequence")
})

test_that("the start of a penalized fit is a starting value and nothing more", {
  skip_if_not_installed("MASS")
  f <- sel_fit(R = 4, seed = 5)
  S <- f$graph$Shat
  lam <- 0.6 * f$graph$n * max(abs(S[upper.tri(S)]))
  a <- mc_graph_penalized(f, lambda = lam, max_iter = 5000L)
  b <- mc_graph_penalized(f, lambda = lam, start = a$eta)
  expect_true(a$converged && b$converged)
  expect_identical(a$edges, b$edges)
  expect_lt(b$iterations, a$iterations)
  expect_error(mc_graph_penalized(f, lam, start = 1:3), "finite values")
  bad <- c(1, 1, 1, 1, 5, 0, 0, 0, 0, 0)
  expect_error(mc_graph_penalized(f, lam, start = bad), "positive definite")
})

test_that("the unpenalized end of the path is the complete fit", {
  skip_if_not_installed("MASS")
  f <- sel_fit(R = 4, seed = 7)
  s <- mc_select_lambda(f, n_lambda = 10, include_zero = TRUE)
  last <- s$lambda_path[nrow(s$lambda_path), ]
  expect_identical(last$lambda, 0)
  expect_identical(last$penalized_status, "unpenalized")
  expect_identical(last$n_edges, choose(4, 2))
  full <- mcglmgraph:::.refit_between(f, mc_between_graph("ggm"))
  expect_equal(last$criterion, as.numeric(mc_compute_criterion(full, "bic")),
               tolerance = 1e-12)
  # without it the path is the geometric grid alone
  s0 <- mc_select_lambda(f, n_lambda = 10)
  expect_identical(s0$lambda_path, s$lambda_path[-nrow(s$lambda_path), ])
})

## ---- stability selection with error control ---------------------------------

test_that("the bound is the one of Meinshausen and Buehlmann", {
  expect_equal(mcglmgraph:::.stability_bound(3, 0.75, 10), 9 / (0.5 * 10))
  expect_equal(mcglmgraph:::.stability_bound(2, 0.9, 45), 4 / (0.8 * 45))
  expect_true(is.na(mcglmgraph:::.stability_bound(2, 0.5, 45)))
})

test_that("the first edge to enter the path is the largest sample covariance", {
  skip_if_not_installed("MASS")
  f <- sel_fit(R = 5, seed = 9)
  S <- f$graph$Shat
  ut <- which(upper.tri(S), arr.ind = TRUE)
  first <- ut[which.max(abs(S[ut])), ]
  got <- mcglmgraph:::.first_q_edges(f, 1L)
  expect_identical(as.integer(got[1, ]), as.integer(sort(first)))
  three <- mcglmgraph:::.first_q_edges(f, 3L)
  expect_equal(nrow(three), 3L)
  expect_identical(as.integer(three[1, ]), as.integer(got[1, ]))
})

test_that("q turns the frequencies into stability selection", {
  skip_if_not_installed("MASS")
  f <- sel_fit(R = 5, seed = 11)
  set.seed(2)
  st <- mc_graph_stability(f, B = 10, q = 3, pi_thr = 0.8)
  expect_identical(st$q, 3L)
  expect_equal(st$error_bound, 9 / (0.6 * 10))
  expect_equal(st$B_valid, 10L)
  # every subsample selected exactly q edges
  expect_equal(sum(st$selection_prob[upper.tri(st$selection_prob)]), 3)
  # complementary pairs: 2B subsamples, frequencies on a grid of 1/(2B)
  set.seed(2)
  cp <- mc_graph_stability(f, B = 6, q = 3, pi_thr = 0.8,
                           sampling = "complementary")
  expect_identical(cp$n_subsamples, 12L)
  pr <- cp$selection_prob[upper.tri(cp$selection_prob)]
  expect_equal(pr * 12, round(pr * 12), tolerance = 1e-12)
  expect_error(mc_graph_stability(f, B = 5, q = 3, method = "stepwise"),
               "applies to method")
  expect_error(mc_graph_stability(f, B = 5, q = 11), "between 1 and")
  expect_output(print(st), "Bound on the expected number of false edges")
})

test_that("under independence the false selections stay below the bound", {
  skip_on_cran()
  skip_long()
  skip_if_not_installed("MASS")
  # every edge is null, so every selected edge is false
  q <- 3L; pi_thr <- 0.75; R <- 5
  bound <- q^2 / ((2 * pi_thr - 1) * choose(R, 2))
  V <- vapply(seq_len(60), function(r) {
    set.seed(700 + r); n <- 200
    d <- as.data.frame(matrix(rnorm(n * R), n, R))
    names(d) <- paste0("y", seq_len(R)); d$x <- rnorm(n)
    f <- mcglm_graph(lapply(paste0("y", seq_len(R)),
                            function(v) stats::reformulate("x", response = v)),
                     data = d)
    st <- mc_graph_stability(f, B = 20, q = q, pi_thr = pi_thr)
    nrow(st$stable_edges)
  }, numeric(1))
  expect_lte(mean(V), bound)
})

## ---- L1-penalized algorithm: the start and the stopping rule of the iteration

kkt_case <- function() {
  set.seed(1); R <- 5; n <- 150
  K <- diag(R) * 2
  for (k in 1:(R - 1)) K[k, k + 1] <- K[k + 1, k] <- -0.7
  Y <- MASS::mvrnorm(n, rep(0, R), solve(K))
  d <- as.data.frame(Y); names(d) <- paste0("y", 1:R); d$x <- rnorm(n)
  f <- mcglm_graph(lapply(paste0("y", 1:R),
                          function(v) stats::reformulate("x", response = v)),
                   data = d)
  S <- f$graph$Shat
  list(f = f, lambda = 0.2 * n * max(abs(S[upper.tri(S)])),
       inverse_start = c(diag(solve(S)), rep(0, choose(R, 2))))
}

test_that("the iteration starts from the diagonal fit of the empty graph", {
  skip_if_not_installed("MASS")
  kc <- kkt_case()
  f <- kc$f; S <- f$graph$Shat
  # the start is 1 / S_kk, which is the root of the diagonal fit
  d0 <- mcglmgraph:::.refit_between(f, mc_between_graph("ggm", matrix(0L, 0, 2)))
  expect_equal(unname(d0$graph$par), unname(1 / diag(S)), tolerance = 1e-8)
  # the objective at the start is the one of K = diag(1 / S_kk):
  # (n / 2) {tr(K S) - log|K|} = (n / 2) {R + sum log S_kk}
  p0 <- mc_graph_penalized(f, kc$lambda, max_iter = 1L)
  expect_equal(p0$objective_path[1],
               0.5 * f$graph$n * (nrow(S) + sum(log(diag(S)))),
               tolerance = 1e-12)
  # a convex problem: the certified solutions agree whatever the start
  a <- mc_graph_penalized(f, kc$lambda)
  b <- mc_graph_penalized(f, kc$lambda, start = kc$inverse_start)
  expect_true(a$converged && b$converged)
  expect_identical(a$edges, b$edges)
  expect_lt(max(abs(a$eta - b$eta)), 1e-4)
})

test_that("a small step does not stop the iteration while KKT fails", {
  skip_if_not_installed("MASS")
  kc <- kkt_case()
  # from the inverse of S as start, the parameter gate passes at the fifth
  # iteration with a KKT residual above the tolerance, so a rule based on the
  # step alone would stop there; with max_iter = 5 the status is
  # kkt_not_certified
  early <- mc_graph_penalized(kc$f, kc$lambda, tol = 1e-2,
                              start = kc$inverse_start, max_iter = 5L)
  expect_identical(early$status, "kkt_not_certified")
  expect_gt(early$kkt_unit, early$kkt_tol)
  # the rule of the L1-penalized algorithm goes on and certifies
  full <- mc_graph_penalized(kc$f, kc$lambda, tol = 1e-2,
                             start = kc$inverse_start)
  expect_identical(full$status, "converged")
  expect_lte(full$kkt_unit, full$kkt_tol)
  expect_gt(full$iterations, 5L)
})

test_that("a response without residual variance is refused on the penalty path", {
  skip_if_not_installed("MASS")
  f <- sel_fit(n = 120)
  # degenerate residual second moment: constant Pearson residuals for y2
  f$graph$Shat[2, ] <- f$graph$Shat[, 2] <- 0
  expect_error(mc_graph_penalized(f, lambda = 1), "y2")
  expect_error(mc_select_lambda(f, n_lambda = 5), "y2")
  expect_error(mc_graph_penalized(f, lambda = 1), "residual variance")
  # subsample path of stability selection goes through the same guard
  expect_error(mcglmgraph:::.first_q_edges(f, q = 2, n_lambda = 5),
               "subsample path.*y2")
  # without response names the position is reported; Inf is refused too
  dimnames(f$graph$Shat) <- NULL
  expect_error(mc_graph_penalized(f, lambda = 1), "response\\(s\\) 2 ")
  f$graph$Shat[2, 2] <- Inf
  expect_error(mc_graph_penalized(f, lambda = 1), "non-finite")
})

# ---- step acceptance of the penalized iteration -----------------------------

test_that("the penalized objective never increases and the path is certified", {
  skip_if_not_installed("MASS")
  set.seed(41); n <- 400; R <- 8
  K <- diag(R) * 2
  for (k in seq_len(R - 1L)) K[k, k + 1L] <- K[k + 1L, k] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, R), solve(K))
  colnames(Y) <- paste0("y", seq_len(R))
  f <- mcglm_graph(lapply(seq_len(R), function(k)
    stats::as.formula(paste0("y", k, " ~ x"))), data = data.frame(Y, x = rnorm(n)))

  # a penalty in the middle of the path, where a step of fixed scale overshoots
  pen <- mc_graph_penalized(f, lambda = 5)
  expect_true(pen$converged)
  expect_identical(pen$status, "converged")
  expect_lt(pen$iterations, 300L)
  expect_lte(pen$kkt_unit, 1e-4)
  # to the rounding of the objective itself
  expect_true(all(diff(pen$objective_path) <=
                    1e-12 * (1 + abs(pen$objective_path[1L]))))
  # a tolerance a hundred times below the default is reached as well
  tight <- mc_graph_penalized(f, lambda = 5, tol = 1e-8, max_iter = 5000L)
  expect_true(tight$converged)

  # the support is the set of coefficients that are not exactly zero
  edge_coef <- pen$eta[-seq_len(R)]
  expect_identical(pen$n_selected, sum(edge_coef != 0))
  expect_identical(nrow(pen$edges), pen$n_selected)

  # every point of the default grid carries the certificate
  sel <- mc_select_lambda(f, n_lambda = 15)
  expect_true(all(sel$lambda_path$penalized_converged))
  expect_identical(unname(sel$edges), cbind(seq_len(R - 1L), 2:R))
})

test_that("the path passes its tolerances to the penalized fit", {
  skip_if_not_installed("MASS")
  f <- sel_fit()
  # one iteration certifies the empty graph at the largest penalty and no more
  one <- mc_select_lambda(f, n_lambda = 8, max_iter = 1L)
  expect_true(one$lambda_path$penalized_converged[1L])
  expect_false(all(one$lambda_path$penalized_converged))
  expect_true(all(one$lambda_path$penalized_status[
    !one$lambda_path$penalized_converged] == "maximum_iterations"))
  # no relative change is below a tolerance of zero, so nothing is certified
  expect_error(mc_select_lambda(f, n_lambda = 5, tol = 0), "No L1 candidate")
  expect_error(mc_select_lambda(f, n_lambda = 5, max_iter = 0), "max_iter")
  expect_error(mc_select_lambda(f, n_lambda = 5, tol = -1), "tol")
  expect_error(mc_select_lambda(f, n_lambda = 5, kkt_tol = NA), "kkt_tol")
})
