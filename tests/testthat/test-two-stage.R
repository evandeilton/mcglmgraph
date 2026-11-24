# Propagation of the first-stage uncertainty into the graph layer.
#
# Targets rebuilt from their definitions:
#   Q          Q[d, r] = -(1/2) [ Z_d Sigma_b ](r, r)
#   shortcut   IF_i = H^-1 { psi_i - Q (r_i o r_i - m) / m }, m the average
#              squared residual of each response, which the first stage fixes
#   geometry   H^-1 Q[, r] is the derivative of eta along the rescaling
#              K -> S^-1 K S^-1 with S = I + eps e_r e_r', so it carries
#              -K_rr in the row of eta_rr, -K_rj/2 in the rows of the edges
#              incident to r, and zero elsewhere
#   invariance g_ij' H^-1 Q = 0 for every pair, so the delta-method variance of
#              a partial correlation is unchanged by the correction
#   size       with Gaussian residuals the variance of eta_kk goes from
#              2 K_kk^2 / n to 4 K_kk (K_kk - 1) / n

bare <- function(M) {
  M <- as.matrix(M)
  attributes(M) <- list(dim = dim(M))
  M
}

make_fit <- function(n = 600, seed = 11, k12 = -0.8, k23 = 0) {
  set.seed(seed)
  K <- diag(3) * 2
  K[1, 2] <- K[2, 1] <- k12
  K[2, 3] <- K[3, 2] <- k23
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

## The canonical directions Z_d of the precision, from the derivatives of the
## covariance: Z_d = -K Ddot_d K.
theory_Z <- function(K, Dlist) lapply(Dlist, function(Dd) -K %*% Dd %*% K)

test_that("Q is the diagonal of Z_d Sigma_b, halved and signed", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  Q <- mcglmgraph:::cpp_two_stage_q(as.matrix(g$K_b), g$Dlist)
  Z <- theory_Z(as.matrix(g$K_b), g$Dlist)
  want <- t(vapply(Z, function(Zd) -0.5 * diag(Zd %*% g$Sigma_b), numeric(3)))
  expect_equal(bare(Q), bare(want), tolerance = 1e-10)

  # in the GGM layer the canonical directions are the bases themselves, so the
  # rows have the closed form the proposition states
  b <- mc_precision_bases(3)
  for (d in seq_along(b$bases))
    expect_equal(Z[[d]], b$bases[[d]], tolerance = 1e-9, ignore_attr = TRUE)
  S <- g$Sigma_b
  for (k in 1:3) {                       # diagonal parameter eta_kk
    row <- rep(0, 3); row[k] <- -0.5 * S[k, k]
    expect_equal(unname(Q[k, ]), row, tolerance = 1e-10)
  }
  ed <- b$edges                          # edge parameters
  for (e in seq_len(nrow(ed))) {
    i <- ed[e, 1]; j <- ed[e, 2]
    row <- rep(0, 3); row[i] <- row[j] <- -0.5 * S[i, j]
    expect_equal(unname(Q[3 + e, ]), row, tolerance = 1e-10)
  }
})

test_that("H^-1 Q is the direction of a diagonal rescaling of the precision", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  Q <- mcglmgraph:::cpp_two_stage_q(as.matrix(g$K_b), g$Dlist)
  Hinv <- g$n * solve(g$J)               # J is accumulated: J = n H
  HQ <- Hinv %*% Q
  K <- as.matrix(g$K_b)
  b <- mc_precision_bases(3)
  ed <- b$edges
  for (r in 1:3) {
    want <- numeric(nrow(HQ))
    want[r] <- -K[r, r]                                  # row of eta_rr
    for (e in seq_len(nrow(ed))) {
      i <- ed[e, 1]; j <- ed[e, 2]
      if (i == r) want[3 + e] <- -0.5 * K[r, j]
      if (j == r) want[3 + e] <- -0.5 * K[r, i]
    }
    expect_equal(unname(HQ[, r]), want, tolerance = 1e-7)
  }

  # the same vector, read as a derivative: eta(K_eps) with K_eps the rescaled
  # precision, differentiated numerically at eps = 0, equals 2 H^-1 Q[, r]
  eps <- 1e-5
  for (r in 1:3) {
    eta_of <- function(e) {
      Se <- diag(3); Se[r, r] <- 1 + e
      Ke <- solve(Se) %*% K %*% solve(Se)
      c(diag(Ke), Ke[cbind(ed[, 1], ed[, 2])])
    }
    num <- (eta_of(eps) - eta_of(-eps)) / (2 * eps)
    expect_equal(unname(num), unname(2 * HQ[, r]), tolerance = 1e-5)
  }
})

test_that("the correction leaves the partial correlations untouched", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  Q <- mcglmgraph:::cpp_two_stage_q(as.matrix(g$K_b), g$Dlist)
  Hinv <- g$n * solve(g$J)
  del <- mc_partial_cor_delta(f)
  G <- attr(del, "gradient")             # q x D gradients of the pcors
  # g_ij' H^-1 Q = 0 for every pair: the proposition that makes the partial
  # correlations free of the first stage to first order
  expect_lt(max(abs(G %*% Hinv %*% Q)), 1e-6 * max(1, max(abs(G))))

  # so the delta-method standard errors agree under the two sandwiches
  a <- mc_partial_cor_delta(f, variability = "empirical")
  Vt <- mc_vcov_graph(f, "two_stage")
  se_two <- sqrt(pmax(diag(G %*% Vt %*% t(G)), 0))
  expect_equal(a$se, se_two, tolerance = 5e-3 * max(a$se))
})

test_that("an absent edge is not touched by the correction", {
  skip_if_not_installed("MASS")
  # K_23 = 0 by construction, so the correction to eta_23, which is
  # proportional to K_23, vanishes
  f <- make_fit(n = 2000, seed = 5, k12 = -0.8, k23 = 0)
  Vs <- mc_vcov_graph(f, "sandwich")
  Vt <- mc_vcov_graph(f, "two_stage")
  nm <- names(f$graph$par)
  r <- sqrt(diag(Vt)) / sqrt(diag(Vs))
  expect_equal(unname(r[nm == "eta_23"]), 1, tolerance = 0.02)
  expect_equal(unname(r[nm == "eta_13"]), 1, tolerance = 0.02)
  expect_gt(unname(r[nm == "eta_12"]), 1)        # the present edge does move
})

test_that("the corrected contributions sum to zero", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  Ec <- scale(g$residuals, center = TRUE, scale = FALSE)
  Psi <- mcglmgraph:::.pearson_scores(g$residuals, g$Sigma_b, g$Dlist)
  ts <- mcglmgraph:::cpp_two_stage_scores(Psi, as.matrix(g$K_b), g$Dlist, Ec)
  expect_lt(max(abs(colSums(ts$scores))) / g$n, 1e-12)
  # and they are the naive ones minus Q (r o r - m)
  m <- colMeans(Ec^2)
  want <- Psi - sweep(Ec^2, 2L, m) %*% t(ts$Q)
  expect_equal(bare(ts$scores), bare(want), tolerance = 1e-10)
  expect_equal(as.numeric(ts$centre), unname(m), tolerance = 1e-12)
  # the center is the first-stage bias correction, 1 - k_r/n, not exactly one
  expect_equal(as.numeric(ts$centre), rep(1 - 2 / g$n, 3), tolerance = 1e-8)
})

test_that("the sandwich is recovered when Q is set aside", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  g <- f$graph
  Psi <- mcglmgraph:::.pearson_scores(g$residuals, g$Sigma_b, g$Dlist)
  Jinv <- solve(g$J)
  naive <- Jinv %*% crossprod(Psi) %*% Jinv
  expect_equal(bare(mc_vcov_graph(f, "sandwich")),
               bare(0.5 * (naive + t(naive))), tolerance = 1e-9)
})

test_that("the size of the correction matches the closed form", {
  skip_if_not_installed("MASS")
  # With Gaussian residuals, and with the diagonal of Sigma_b already fixed at
  # one by the first stage, the variance of eta_kk goes from 2 K_kk^2 / n to
  # 4 K_kk (K_kk - 1) / n, so the ratio of standard errors is
  # sqrt(2 (K_kk - 1) / K_kk). Response 3 is independent of the others, so its
  # K_kk is one and the correction wipes its variance out.
  f <- make_fit(n = 20000, seed = 17, k12 = -0.8)
  g <- f$graph
  K <- as.matrix(g$K_b)
  Vs <- mc_vcov_graph(f, "sandwich")
  Vt <- mc_vcov_graph(f, "two_stage")
  ratio <- sqrt(diag(Vt)) / sqrt(diag(Vs))
  want <- unname(sqrt(2 * (diag(K) - 1) / diag(K)))
  for (k in 1:3)
    expect_equal(unname(ratio[k]), want[k], tolerance = 0.08,
                 info = paste("response", k))
  expect_lt(unname(ratio[3]), 0.1)          # K_33 = 1: nothing left
})

test_that("a fit the shortcut does not cover is refused, with a reason", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 300)

  # a Poisson-Tweedie response whose power is not fixed at one: the variance
  # mu + tau mu^p is additive, so the first stage does not fix the average of
  # the squared residuals and the shortcut has nothing to stand on
  g1 <- f
  g1$variance <- as.list(c("poisson_tweedie", "constant", "constant"))
  g1$power_fixed <- list(TRUE, TRUE, TRUE)
  g1$list_initial$power <- list(2, 0, 0)
  av <- mc_two_stage_available(g1)
  expect_false(av$ok)
  expect_match(paste(av$why, collapse = " "), "Poisson-Tweedie")
  expect_error(mc_vcov_graph(g1, "two_stage", method = "shortcut"),
               "does not apply to this fit")

  # the same family with the power fixed at one does qualify
  g2 <- g1
  g2$list_initial$power <- list(1, 0, 0)
  expect_true(mc_two_stage_available(g2)$ok)

  # more than one covariance component means more than one dispersion parameter
  g3 <- f
  g3$matrix_pred[[2]] <- c(g3$matrix_pred[[2]], g3$matrix_pred[[2]])
  av3 <- mc_two_stage_available(g3)
  expect_false(av3$ok)
  expect_match(paste(av3$why, collapse = " "), "more than one dispersion")

  # an unknown variance function is refused rather than assumed
  g4 <- f
  g4$variance <- as.list(c("something_else", "constant", "constant"))
  expect_false(mc_two_stage_available(g4)$ok)

  # and the conditional forms remain available throughout
  expect_true(is.matrix(mc_vcov_graph(g1, "sandwich")))
  expect_true(is.matrix(mc_vcov_graph(g1, "model")))
})

test_that("the shortcut needs constant prior weights", {
  skip_if_not_installed("MASS")
  set.seed(7); n <- 200
  K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3])
  fit_w <- function(w) suppressMessages(mcglm_graph(
    list(y1 ~ 1, y2 ~ 1, y3 ~ 1), data = d, weights = w,
    control_algorithm = list(tol = 1e-10, max_iter = 100)))

  # weights that vary: the first stage fixes a weighted average, so the
  # shortcut is refused and auto takes the full form
  set.seed(8)
  fu <- fit_w(rep(list(stats::runif(n, 0.5, 2)), 3))
  av <- mc_two_stage_available(fu)
  expect_false(av$shortcut_ok)
  expect_true(av$full_ok)
  expect_match(paste(av$why_shortcut, collapse = " "), "prior weights")
  expect_identical(av$method, "full")
  expect_identical(attr(mc_vcov_graph(fu, "two_stage"), "two_stage_method"),
                   "full")
  expect_error(mc_vcov_graph(fu, "two_stage", method = "shortcut"),
               "does not apply to this fit")

  # a constant weight w0 keeps the shortcut. The bias correction of mcglm does
  # not scale with w0, so the first stage fixes mean(r^2) at 1 - p/(n w0^2)
  # instead of 1 - p/n; the shortcut reads that average and stays exact
  fc <- fit_w(rep(list(rep(2, n)), 3))
  f1 <- fit_w(NULL)
  expect_true(mc_two_stage_available(fc)$shortcut_ok)
  expect_identical(mc_two_stage_available(fc)$method, "shortcut")
  expect_equal(unname(colMeans(fc$graph$residuals^2)), rep(1 - 1 / (4 * n), 3),
               tolerance = 1e-8)
  expect_equal(unname(colMeans(f1$graph$residuals^2)), rep(1 - 1 / n, 3),
               tolerance = 1e-8)
  Vs <- mc_vcov_graph(fc, "two_stage", method = "shortcut")
  Vo <- mc_vcov_graph(fc, "two_stage", method = "full")
  expect_lt(max(abs(bare(Vs) - bare(Vo))) / max(abs(bare(Vo))), 1e-6)
})

test_that("counts fitted with a multiplicative variance qualify", {
  skip_if_not_installed("MASS")
  set.seed(37); n <- 400; x <- rnorm(n)
  d <- data.frame(y1 = rpois(n, exp(1 + 0.3 * x)),
                  y2 = rpois(n, exp(1 - 0.2 * x)),
                  y3 = rpois(n, 6), x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d,
                   link = rep("log", 3), variance = rep("tweedie", 3),
                   power_fixed = rep(TRUE, 3))
  expect_true(mc_two_stage_available(f)$ok)
  V <- mc_vcov_graph(f, "two_stage")
  expect_true(is.matrix(V))
  expect_true(all(diag(V) >= 0))
})

test_that("the two-stage form honors the clustering arguments", {
  skip_if_not_installed("MASS")
  f <- make_fit()
  n <- f$graph$n
  V0 <- mc_vcov_graph(f, "two_stage")
  V1 <- mc_vcov_graph(f, "two_stage", cluster = seq_len(n))
  expect_equal(bare(V0), bare(V1), tolerance = 1e-12)
  cl <- rep(seq_len(n / 2), each = 2L)
  Vc <- mc_vcov_graph(f, "two_stage", cluster = cl)
  Vs <- mc_vcov_graph(f, "two_stage", cluster = cl, small_sample = TRUE)
  G <- n / 2
  expect_equal(bare(Vs), bare(Vc) * G / (G - 1), tolerance = 1e-10)
  expect_identical(attr(Vc, "n_clusters"), as.integer(G))
  expect_identical(attr(Vc, "type"), "two_stage")
})

test_that("the DAG layer is served by the same correction", {
  skip_if_not_installed("MASS")
  set.seed(41); n <- 500; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag"), data = d)
  expect_true(mc_two_stage_available(f)$ok)
  Vt <- mc_vcov_graph(f, "two_stage")
  expect_equal(dim(Vt), rep(length(f$graph$par), 2L))
  expect_true(all(diag(Vt) >= 0))
  # Q is built from the derivatives of this layer, so it is not the GGM one
  Q <- attr(Vt, "Q")
  expect_equal(dim(Q), c(length(f$graph$par), 3L))
})

test_that("a single dispersion parameter is necessary and not sufficient", {
  skip_if_not_installed("MASS")
  f <- make_fit(n = 200)
  expect_true(mc_two_stage_available(f)$ok)

  # One covariance component that is not a multiple of the identity keeps the
  # parameter count at one and still breaks the factorization over units that
  # the identity behind the shortcut rests on.
  g <- f
  M <- diag(f$graph$n); M[1, 2] <- M[2, 1] <- 0.5
  g$matrix_pred[[2]] <- list(M)
  av <- mc_two_stage_available(g)
  expect_false(av$ok)
  expect_match(paste(av$why, collapse = " "), "not a multiple of the identity")
  expect_error(mc_vcov_graph(g, "two_stage", method = "shortcut"),
               "does not apply to this fit")
  # the full form needs the units the structure is block diagonal over
  expect_error(mc_vcov_graph(g, "two_stage"), "pass cluster")
  expect_false(mc_two_stage_available(g)$full_ok)
  expect_true(mc_two_stage_available(g, cluster = c(1L, 1L,
                                                    3:f$graph$n))$full_ok)

  # the conditional forms stay available, since they assume nothing of the kind
  expect_true(is.matrix(mc_vcov_graph(g, "sandwich")))
})
