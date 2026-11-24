# The full two-stage correction.
#
# Targets rebuilt from the theory, in plain R, without the package's kernels:
#   units        psi_beta,i = D_i' C_i^-1 (y_i - mu_i),
#                psi_l,i    = e_i' C_i^-1 Cdot_l,i C_i^-1 e_i - tr(C_i^-1 Cdot_l,i),
#                with C_i the R x R block of unit i of the first-stage covariance
#   sensitivity  S = d E_alpha0{ psi(alpha) } / d alpha' at alpha0, by numDeriv
#   cross        B_n = -d/d alpha' sum_i psi_2i(alpha), residuals recomputed at
#                alpha, by numDeriv
#   influence    IF_u = H^-1 { psi_2u - B A^-1 psi_1u }, the two-stage
#                influence function
#   shortcut     on the fits the shortcut covers, the two forms agree
#   rescaling    with common Gaussian designs B has the rescaling form exactly;
#                with (c) false the partial correlations stop being free of it

bare <- function(M) {
  M <- as.matrix(M)
  attributes(M) <- list(dim = dim(M))
  M
}

## Block-diagonal matrix from a list of dense blocks.
bdiag_dense <- function(blocks) {
  nr <- vapply(blocks, nrow, integer(1)); nc <- vapply(blocks, ncol, integer(1))
  out <- matrix(0, sum(nr), sum(nc))
  ro <- c(0L, cumsum(nr)); co <- c(0L, cumsum(nc))
  for (k in seq_along(blocks))
    out[ro[k] + seq_len(nr[k]), co[k] + seq_len(nc[k])] <- blocks[[k]]
  out
}

## Gamma margins with a fixed shape under a Gaussian copula: a scale family,
## so condition (c) holds unless `cvary`, which lets the copula correlation of
## the first two responses move with x.
gamma_fit <- function(n, seed, cvary = FALSE, shape = 2, power = 2) {
  set.seed(seed); x <- stats::runif(n)
  if (!cvary) {
    Rl <- matrix(c(1, .5, .2, .5, 1, .4, .2, .4, 1), 3)
    U <- stats::pnorm(MASS::mvrnorm(n, rep(0, 3), Rl))
  } else {
    U <- t(vapply(seq_len(n), function(i) {
      a <- 0.1 + 0.7 * x[i]
      Rl <- matrix(c(1, a, .2, a, 1, .4, .2, .4, 1), 3)
      stats::pnorm(MASS::mvrnorm(1, rep(0, 3), Rl))
    }, numeric(3)))
  }
  mu <- cbind(exp(1 + .5 * x), exp(.5 - .3 * x), rep(exp(.8), n))
  Y <- sapply(1:3, function(r) stats::qgamma(U[, r], shape = shape,
                                             rate = shape / mu[, r]))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = x)
  suppressMessages(mcglm_graph(
    list(y1 ~ x, y2 ~ x, y3 ~ x), data = d, link = rep("log", 3),
    variance = rep("tweedie", 3), power_fixed = rep(TRUE, 3),
    control_initial = list(regression = list(c(1, .5), c(.5, -.3), c(.8, 0)),
                           power = list(power, power, power),
                           tau = list(1 / shape, 1 / shape, 1 / shape),
                           rho = c(.5, .2, .4))))
}

gauss_fit <- function(n, seed, intercept_only = FALSE) {
  set.seed(seed); x <- stats::runif(n)
  Rl <- matrix(c(1, .5, .2, .5, 1, .4, .2, .4, 1), 3)
  Z <- MASS::mvrnorm(n, rep(0, 3), Rl)
  d <- data.frame(y1 = 1 + x + Z[, 1], y2 = 2 - x + 1.5 * Z[, 2],
                  y3 = .5 * Z[, 3], x = x)
  lp <- if (intercept_only) list(y1 ~ 1, y2 ~ 1, y3 ~ 1) else
    list(y1 ~ x, y2 ~ x, y3 ~ x)
  mcglm_graph(lp, data = d)
}

## Poisson-Tweedie counts through a Poisson-gamma mixture: y | Z ~ Poisson(mu Z),
## Z ~ Gamma(k, k) with k = mu^(2 - p) / tau, so Var y = mu + tau mu^p. The
## power is estimated in the first stage, which is the case the shortcut
## refuses.
pt_fit <- function(n, seed) {
  set.seed(seed); x <- stats::runif(n)
  Rl <- matrix(c(1, .5, .2, .5, 1, .4, .2, .4, 1), 3)
  U <- stats::pnorm(MASS::mvrnorm(n, rep(0, 3), Rl))
  mu <- cbind(exp(1 + .5 * x), exp(.5 - .3 * x), rep(exp(.8), n))
  Y <- sapply(1:3, function(r) {
    k <- mu[, r]^0.5 / 0.5
    stats::rpois(n, mu[, r] * stats::qgamma(U[, r], shape = k, rate = k))
  })
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = x)
  suppressMessages(mcglm_graph(
    list(y1 ~ x, y2 ~ x, y3 ~ x), data = d, link = rep("log", 3),
    variance = rep("poisson_tweedie", 3), power_fixed = rep(FALSE, 3),
    control_initial = list(regression = list(c(1, .5), c(.5, -.3), c(.8, 0)),
                           power = list(1.5, 1.5, 1.5), tau = list(.5, .5, .5),
                           rho = c(.4, .2, .3)),
    control_algorithm = list(tuning = 0.5, max_iter = 300, tol = 1e-6)))
}

## A structure like that of the Hunting application: the identity plus a
## random-effect block per pair of rows, the blocks nested in clusters.
blocked_fit <- function(seed, G = 40, m = 6) {
  set.seed(seed)
  Rl <- matrix(c(1, .5, .2, .5, 1, .4, .2, .4, 1), 3)
  cl <- rep(seq_len(G), each = m); blk <- rep(seq_len(G * m / 2), each = 2)
  n <- G * m; x <- stats::runif(n)
  b <- MASS::mvrnorm(G * m / 2, rep(0, 3), 0.15 * Rl)[blk, ]
  U <- stats::pnorm(MASS::mvrnorm(n, rep(0, 3), Rl))
  mu <- exp(cbind(1 + .5 * x, .5 - .3 * x, rep(.8, n)) + b)
  Y <- sapply(1:3, function(r) stats::qpois(U[, r], mu[, r]))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = x, cl = cl,
                  blk = factor(blk))
  Z <- c(mcglm::mc_id(d), mcglm::mc_mixed(~ 0 + blk, data = d))
  suppressMessages(mcglm_graph(
    list(y1 ~ x, y2 ~ x, y3 ~ x), data = d, matrix_pred = rep(list(Z), 3),
    link = rep("log", 3), variance = rep("poisson_tweedie", 3),
    power_fixed = rep(TRUE, 3),
    control_initial = list(regression = list(c(1, .5), c(.5, -.3), c(.8, 0)),
                           power = list(1, 1, 1),
                           tau = list(c(.2, .1), c(.2, .1), c(.2, .1)),
                           rho = c(.4, .2, .3))))
}

## The first stage rebuilt with mcglm's own builder at a given parameter
## vector, for the numerical derivatives below.
first_stage_at <- function(fit, alpha) {
  R <- length(fit$linear_pred)
  inf <- fit$Information
  pb <- sum(unlist(inf$n_betas))
  li <- mcglm:::mc_updateBeta(fit$list_initial, alpha[seq_len(pb)],
                              information = inf, n_resp = R)
  li <- mcglm:::mc_updateCov(li, covariance = alpha[-seq_len(pb)],
                             list_power_fixed = fit$power_fixed,
                             information = inf, n_resp = R)
  ml <- Map(mcglm:::mc_link_function, beta = li$regression,
            offset = fit$offset, X = fit$list_X, link = fit$link)
  Cf <- mcglm:::mc_build_C(list_mu = ml, list_Ntrial = fit$Ntrial,
                           rho = li$rho, list_tau = li$tau,
                           list_power = li$power, list_Z = fit$matrix_pred,
                           list_sparse = fit$sparse,
                           list_variance = fit$variance,
                           list_covariance = fit$covariance,
                           list_power_fixed = fit$power_fixed,
                           compute_C = TRUE)
  list(mu = unlist(lapply(ml, function(m) m$mu)),
       D = lapply(ml, function(m) as.matrix(m$D)),
       C = as.matrix(Cf$C), invC = as.matrix(Cf$inv_C),
       DC = lapply(Cf$D_C, as.matrix))
}

## Totals of the first-stage estimating functions at alpha, in plain R with
## unit weights; with `correct`, the Pearson functions carry the bias
## correction tr{ (D'C^-1 D)^-1 D'C^-1 Cdot_l C^-1 D } that mcglm adds.
totals_dense <- function(fit, a, y, correct) {
  s <- first_stage_at(fit, a)
  D <- bdiag_dense(s$D)
  z <- s$invC %*% (y - s$mu)
  U <- c(drop(t(D) %*% z),
         vapply(s$DC, function(Cd) drop(t(z) %*% Cd %*% z) -
                  sum(diag(s$invC %*% Cd)), numeric(1)))
  U <- unname(U)
  if (!correct) return(U)
  F <- s$invC %*% D
  iJ <- solve(t(D) %*% F)
  cr <- unname(vapply(s$DC, function(Cd)
    sum(diag(iJ %*% t(F) %*% Cd %*% F)), numeric(1)))
  U + c(rep(0, ncol(D)), cr)
}

test_that("the first stage is rebuilt as mcglm fitted it", {
  skip_if_not_installed("MASS")
  f <- gamma_fit(400, 3)
  g <- f$graph
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(g$n), g$n)
  # the Pearson residuals the graph layer was fitted on
  expect_equal(fs$r, as.numeric(g$residuals), tolerance = 1e-12)
  # the two blocks mcglm reports, against the inverse it stores
  JIS <- as.matrix(f$joint_inv_sensitivity)
  pb <- fs$n_beta
  expect_equal(bare(solve(fs$S[1:pb, 1:pb])), bare(JIS[1:pb, 1:pb]),
               tolerance = 1e-8)
  expect_equal(bare(solve(fs$S[-(1:pb), -(1:pb)])),
               bare(JIS[-(1:pb), -(1:pb)]), tolerance = 1e-8)
  # the block mcglm leaves at zero is not zero here
  expect_equal(max(abs(JIS[-(1:pb), 1:pb])), 0)
  expect_gt(max(abs(fs$S[-(1:pb), 1:pb])), 1)
})

test_that("the unit contributions are the estimating functions of each unit", {
  skip_if_not_installed("MASS")
  f <- gamma_fit(60, 5)
  g <- f$graph
  n <- g$n; R <- 3
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(n), n)
  alpha <- c(f$Regression, f$Covariance)
  st <- first_stage_at(f, alpha)
  res <- as.numeric(f$residuals)
  pb <- length(f$Regression)
  Dfull <- bdiag_dense(st$D)
  want <- matrix(0, n, ncol(fs$psi))
  for (i in seq_len(n)) {
    rows <- i + (seq_len(R) - 1L) * n          # unit i across the responses
    Ci_inv <- solve(st$C[rows, rows])
    ei <- res[rows]
    want[i, seq_len(pb)] <- t(Dfull[rows, , drop = FALSE]) %*% Ci_inv %*% ei
    for (l in seq_along(st$DC)) {
      Cd <- st$DC[[l]][rows, rows]
      want[i, pb + l] <- drop(t(ei) %*% Ci_inv %*% Cd %*% Ci_inv %*% ei) -
        sum(diag(Ci_inv %*% Cd))
    }
  }
  expect_equal(bare(fs$psi), bare(want), tolerance = 1e-9)
})

test_that("the sensitivity is the derivative of the expected estimating function", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  f <- gamma_fit(40, 9, power = 1.5)
  g <- f$graph
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(g$n), g$n)
  a0 <- c(f$Regression, f$Covariance)
  s0 <- first_stage_at(f, a0)
  pb <- length(f$Regression)
  # E_{alpha0} psi(alpha): the quasi-score has mean D' C^-1 (mu0 - mu) and the
  # Pearson functions tr{ W_l (C0 + (mu - mu0)(mu - mu0)' - C) }
  expected_psi <- function(a) {
    s <- first_stage_at(f, a)
    dm <- s0$mu - s$mu
    Dfull <- bdiag_dense(s$D)
    qb <- drop(t(Dfull) %*% s$invC %*% dm)
    M <- s0$C + tcrossprod(dm) - s$C
    ql <- vapply(s$DC, function(Cd)
      sum(diag(s$invC %*% Cd %*% s$invC %*% M)), numeric(1))
    c(qb, ql)
  }
  Snum <- numDeriv::jacobian(expected_psi, a0)
  expect_lt(max(abs(Snum - fs$S)) / max(abs(fs$S)), 1e-6)
  # the cross block in particular
  expect_lt(max(abs(Snum[-(1:pb), 1:pb] - fs$S[-(1:pb), 1:pb])) /
              max(abs(fs$S[-(1:pb), 1:pb])), 1e-6)
})

test_that("the observed sensitivity is the Jacobian of the first-stage totals", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  f <- gamma_fit(60, 5, power = 1.5)
  a0 <- c(f$Regression, f$Covariance)
  y <- as.numeric(f$fitted) + as.numeric(f$residuals)
  expect_true(f$con$correct)
  # the totals the first stage solves, with dense matrices: the Pearson
  # functions carry the bias correction tr{(D'C^-1 D)^-1 D'C^-1 Cdot_l C^-1 D}
  Snum <- numDeriv::jacobian(function(a) totals_dense(f, a, y, TRUE), a0)
  Sobs <- mcglmgraph:::.first_stage_observed(f)
  expect_lt(max(abs(Sobs - Snum)) / max(abs(Snum)), 1e-6)
  # it is the matrix the full form uses by default
  V <- mc_vcov_graph(f, "two_stage", method = "full")
  expect_identical(attr(V, "sensitivity"), "observed")
  expect_equal(bare(attr(V, "A")), bare(-Sobs), tolerance = 1e-12)
  Ve <- mc_vcov_graph(f, "two_stage", method = "full", sensitivity = "expected")
  expect_identical(attr(Ve, "sensitivity"), "expected")
})

test_that("B is the Jacobian of the averaged second-stage function", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  f <- gamma_fit(300, 13, power = 1.5)
  g <- f$graph
  V <- mc_vcov_graph(f, "two_stage", method = "full")
  B <- attr(V, "B")
  X <- cbind(1, f$data$x)
  Y <- as.matrix(f$data[, c("y1", "y2", "y3")])
  Sinv <- solve(g$Sigma_b)
  # the residuals as a function of alpha = (beta, rho, tau), in plain R: log
  # link and v = tau mu^p, with the power fixed; the correlations of the first
  # stage do not enter
  resid_at <- function(a) {
    b <- matrix(a[1:6], 2)
    tau <- a[10:12]
    sapply(1:3, function(r) {
      mu <- exp(drop(X %*% b[, r]))
      (Y[, r] - mu) / sqrt(tau[r] * mu^1.5)
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
               tolerance = 1e-10)
  Bnum <- -numDeriv::jacobian(U2, a0)
  expect_lt(max(abs(Bnum - B)) / max(abs(B)), 1e-6)
  # the block of the correlations is zero: they do not enter the scale
  expect_equal(max(abs(B[, 7:9])), 0)
})

test_that("on the fits the shortcut covers, the full form agrees with it", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  for (f in list(gauss_fit(20000, 101), gamma_fit(20000, 101))) {
    Vf <- mc_vcov_graph(f, "two_stage", method = "full")
    Vs <- mc_vcov_graph(f, "two_stage", method = "shortcut")
    ratio <- sqrt(diag(Vf)) / sqrt(diag(Vs))
    expect_true(all(abs(ratio - 1) < 0.01))
    expect_identical(attr(Vf, "two_stage_method"), "full")
    expect_identical(attr(Vs, "two_stage_method"), "shortcut")
  }
})

test_that("without the cross block the regression part does not cancel", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  # The same pieces, but with the block-diagonal sensitivity mcglm reports: the
  # diagonal coefficients drift far from the shortcut, which is what the cross
  # block is there to prevent.
  f <- gamma_fit(20000, 11)
  g <- f$graph
  Vs <- mc_vcov_graph(f, "two_stage", method = "shortcut")
  Vf <- mc_vcov_graph(f, "two_stage", method = "full")
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(g$n), g$n)
  psi1 <- sweep(fs$psi, 2L, colMeans(fs$psi))
  A0 <- -solve(as.matrix(f$joint_inv_sensitivity))
  M0 <- t(solve(t(A0), t(attr(Vf, "B"))))
  C0 <- mcglmgraph:::.pearson_scores(g$residuals, g$Sigma_b, g$Dlist) -
    psi1 %*% t(M0)
  V0 <- solve(g$J) %*% crossprod(C0) %*% solve(g$J)
  r0 <- sqrt(diag(V0)) / sqrt(diag(Vs))
  rf <- sqrt(diag(Vf)) / sqrt(diag(Vs))
  expect_gt(max(abs(r0[1:3] - 1)), 0.2)
  expect_lt(max(abs(rf - 1)), 0.01)
})

test_that("with common Gaussian designs B has the rescaling form exactly", {
  skip_if_not_installed("MASS")
  f <- gauss_fit(500, 23)
  g <- f$graph
  V <- mc_vcov_graph(f, "two_stage", method = "full")
  B <- attr(V, "B")
  Q <- mcglmgraph:::cpp_two_stage_q(as.matrix(g$K_b), g$Dlist)
  tau <- f$Covariance[4:6]
  # no mean term: the residuals are orthogonal to the common design ...
  expect_lt(max(abs(B[, 1:6])), 1e-9 * max(abs(B)))
  # ... and the scale term is n Q[, r] / tau_r, the diagonal rescaling
  for (r in 1:3)
    expect_equal(unname(B[, 6 + 3 + r]), unname(g$n * Q[, r] / tau[r]),
                 tolerance = 1e-9)
  # so the gradients of the partial correlations annihilate H^-1 B
  G <- attr(mc_partial_cor_delta(f), "gradient")
  HB <- solve(g$J) %*% B
  expect_lt(max(abs(G %*% HB)) / (max(abs(G)) * max(abs(HB))), 1e-10)
})

test_that("with the residual correlation moving with x, the partial correlations are no longer free of B", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  orth <- function(f) {
    g <- f$graph
    B <- attr(mc_vcov_graph(f, "two_stage", method = "full"), "B")
    HB <- solve(g$J) %*% B
    G <- attr(mc_partial_cor_delta(f), "gradient")
    max(abs(G %*% HB)) / (max(abs(G)) * max(abs(HB)))
  }
  held <- orth(gamma_fit(20000, 7))
  broken <- orth(gamma_fit(20000, 7, cvary = TRUE))
  # with (c) the residual is sampling error, of order n^-1/2; without it the
  # population value is not zero
  expect_lt(held, 5e-3)
  expect_gt(broken, 1.5e-2)
})

test_that("clusters of an identity structure sum the unit contributions", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- gamma_fit(400, 17)
  n <- f$graph$n
  V0 <- mc_vcov_graph(f, "two_stage", method = "full")
  V1 <- mc_vcov_graph(f, "two_stage", method = "full", cluster = seq_len(n))
  expect_equal(bare(V0), bare(V1), tolerance = 1e-10)
  cl <- rep(seq_len(n / 4), each = 4L)
  Vc <- mc_vcov_graph(f, "two_stage", method = "full", cluster = cl)
  Vss <- mc_vcov_graph(f, "two_stage", method = "full", cluster = cl,
                       small_sample = TRUE)
  G <- n / 4
  expect_equal(bare(Vss), bare(Vc) * G / (G - 1), tolerance = 1e-12)
  expect_true(attr(Vss, "small_sample"))
  # rows of one cluster are summed before squaring: rebuild from the row scores
  g <- f$graph
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(n), n)
  psi1 <- sweep(fs$psi, 2L, colMeans(fs$psi))
  M <- t(solve(t(attr(V0, "A")), t(attr(V0, "B"))))
  Cr <- mcglmgraph:::.pearson_scores(g$residuals, g$Sigma_b, g$Dlist) -
    psi1 %*% t(M)
  Cc <- rowsum(Cr, cl)
  want <- solve(g$J) %*% crossprod(Cc) %*% solve(g$J)
  expect_equal(bare(Vc), bare(0.5 * (want + t(want))), tolerance = 1e-10)
})

test_that("auto takes the full form for a variance that depends on the mean and the shortcut for constant ones", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- gamma_fit(300, 19)
  Va <- mc_vcov_graph(f, "two_stage")
  Vf <- mc_vcov_graph(f, "two_stage", method = "full")
  expect_identical(attr(Va, "two_stage_method"), "full")
  expect_equal(bare(Va), bare(Vf), tolerance = 0)
  av <- mc_two_stage_available(f)
  expect_true(av$shortcut_ok && av$full_ok)
  expect_identical(av$method, "full")
  expect_identical(av$ok, av$shortcut_ok)

  # constant variance functions: the shortcut, bit for bit
  g <- gauss_fit(300, 19)
  Vg <- mc_vcov_graph(g, "two_stage")
  Vgs <- mc_vcov_graph(g, "two_stage", method = "shortcut")
  expect_identical(attr(Vg, "two_stage_method"), "shortcut")
  expect_equal(bare(Vg), bare(Vgs), tolerance = 0)
  expect_identical(mc_two_stage_available(g)$method, "shortcut")

  # a power estimated in the first stage takes the shortcut out of play
  p <- pt_fit(300, 11)
  av <- mc_two_stage_available(p)
  expect_false(av$shortcut_ok)
  expect_true(av$full_ok)
  expect_identical(av$method, "full")
  Vp <- mc_vcov_graph(p, "two_stage")
  expect_identical(attr(Vp, "two_stage_method"), "full")
  expect_true(all(diag(Vp) > 0))
  expect_error(mc_vcov_graph(p, "two_stage", method = "shortcut"),
               "does not apply to this fit")
  expect_error(mc_vcov_graph(p, "sandwich", method = "full"),
               "has no meaning")
})

test_that("a graph layer refitted on other residuals is refused", {
  skip_if_not_installed("MASS")
  f <- gamma_fit(200, 31)
  h <- f
  h$graph$residuals <- h$graph$residuals[c(2:200, 1), ]
  expect_error(mc_vcov_graph(h, "two_stage", method = "full"),
               "not fitted on the Pearson residuals")
})

## ---- estimated power, and structures that are not the identity -------------

test_that("with the power estimated, B is still the Jacobian of the second stage", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  f <- pt_fit(300, 11)
  g <- f$graph
  expect_false(mc_two_stage_available(f)$shortcut_ok)
  V <- mc_vcov_graph(f, "two_stage")
  expect_identical(attr(V, "two_stage_method"), "full")
  B <- attr(V, "B")
  X <- cbind(1, f$data$x)
  Y <- as.matrix(f$data[, c("y1", "y2", "y3")])
  Sinv <- solve(g$Sigma_b)
  # alpha = (beta, rho, (p_r, tau_r) per response); the Poisson-Tweedie scale
  # is sqrt(mu + tau mu^p)
  resid_at <- function(a) {
    b <- matrix(a[1:6], 2)
    sapply(1:3, function(r) {
      p <- a[10 + 2 * (r - 1)]; tau <- a[11 + 2 * (r - 1)]
      mu <- exp(drop(X %*% b[, r]))
      (Y[, r] - mu) / sqrt(mu + tau * mu^p)
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
               tolerance = 1e-10)
  Bnum <- -numDeriv::jacobian(U2, a0)
  expect_lt(max(abs(Bnum - B)) / max(abs(B)), 1e-6)
  # the columns of the powers are not zero: the power moves the scale
  expect_gt(max(abs(B[, c(10, 12, 14)])), 0)
  # the blocks of the sensitivity mcglm reports, powers included
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(g$n), g$n)
  JIS <- as.matrix(f$joint_inv_sensitivity)
  expect_equal(bare(solve(fs$S[1:6, 1:6])), bare(JIS[1:6, 1:6]),
               tolerance = 1e-8)
  expect_equal(bare(solve(fs$S[-(1:6), -(1:6)])), bare(JIS[-(1:6), -(1:6)]),
               tolerance = 1e-8)
})

test_that("with the power estimated, the sensitivity is still the derivative of the expected function", {
  skip_on_cran()
  skip_long()
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  f <- pt_fit(300, 11)
  g <- f$graph
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(g$n), g$n)
  a0 <- c(f$Regression, f$Covariance)
  s0 <- first_stage_at(f, a0)
  expected_psi <- function(a) {
    s <- first_stage_at(f, a)
    dm <- s0$mu - s$mu
    qb <- drop(t(bdiag_dense(s$D)) %*% s$invC %*% dm)
    M <- s0$C + tcrossprod(dm) - s$C
    ql <- vapply(s$DC, function(Cd)
      sum(diag(s$invC %*% Cd %*% s$invC %*% M)), numeric(1))
    c(qb, ql)
  }
  Snum <- numDeriv::jacobian(expected_psi, a0, method.args = list(r = 2))
  expect_lt(max(abs(Snum - fs$S)) / max(abs(fs$S)), 1e-6)
})

test_that("blocks nested in clusters are the units of the full form", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- blocked_fit(1, G = 20, m = 4)
  cl <- f$data$cl
  av <- mc_two_stage_available(f, cluster = cl)
  expect_false(av$shortcut_ok)
  expect_true(av$full_ok)
  expect_false(mc_two_stage_available(f)$full_ok)
  # a grouping that splits a block is refused
  bad <- rep(seq_len(40), each = 2L)[c(2:80, 1)]
  expect_false(mc_two_stage_available(f, cluster = bad)$full_ok)
  expect_error(mc_vcov_graph(f, "two_stage", cluster = bad),
               "connects rows of different clusters")
  V <- mc_vcov_graph(f, "two_stage", cluster = cl)
  expect_identical(attr(V, "two_stage_method"), "full")
  expect_true(all(diag(V) > 0))
  # the blocks themselves are units as well: pairs of rows
  Vb <- mc_vcov_graph(f, "two_stage", cluster = as.integer(f$data$blk))
  expect_true(all(diag(Vb) > 0))
})

## ---- the equations the first stage solves ---------------------------------

## Gaussian margins with intercept only, on a chain DAG 1 -> 2 -> 3 or a GGM,
## the first stage solved to a tight tolerance with or without mcglm's bias
## correction.
chain_fit <- function(n, seed, correct = TRUE, type = "dag") {
  set.seed(seed)
  S0 <- matrix(c(1, .5, .25, .5, 1, .5, .25, .5, 1), 3)
  Y <- MASS::mvrnorm(n, c(0, 2, 0), S0)
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3])
  bg <- if (type == "dag")
    mc_between_graph("dag", edges = rbind(c(1L, 2L), c(2L, 3L)),
                     ordering = 1:3) else mc_between_graph("ggm")
  suppressMessages(mcglm_graph(
    list(y1 ~ 1, y2 ~ 1, y3 ~ 1), between_graph = bg, data = d,
    control_algorithm = list(correct = correct, tol = 1e-10, max_iter = 100)))
}

test_that("the observed sensitivity differentiates the equations as solved", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  for (cr in c(TRUE, FALSE)) {
    f <- chain_fit(100, 1, correct = cr)
    expect_identical(f$con$correct, cr)
    a0 <- c(f$Regression, f$Covariance)
    y <- as.numeric(f$fitted) + as.numeric(f$residuals)
    w <- rep(1, length(y))
    # the package totals are the plain-R ones, with or without the correction
    for (k in c(TRUE, FALSE))
      expect_equal(mcglmgraph:::.first_stage_totals(f, a0, y, w, k),
                   totals_dense(f, a0, y, k), tolerance = 1e-10)
    # the equations solved vanish at the estimate; the other set does not
    U <- totals_dense(f, a0, y, cr)
    expect_lt(max(abs(U)), 1e-6)
    expect_gt(max(abs(totals_dense(f, a0, y, !cr))), 0.5)
    # the observed sensitivity is their Jacobian
    Snum <- numDeriv::jacobian(function(a) totals_dense(f, a, y, cr), a0)
    Sobs <- mcglmgraph:::.first_stage_observed(f)
    expect_lt(max(abs(Sobs - Snum)) / max(abs(Snum)), 1e-7)
    V <- mc_vcov_graph(f, "two_stage", method = "full")
    expect_equal(bare(attr(V, "A")), bare(-Sobs), tolerance = 1e-12)
  }
})

test_that("the unit contributions sum to minus the bias correction", {
  skip_if_not_installed("MASS")
  # the correction is one term for the sample; centering the contributions
  # adds it back, so the sandwich does not see it
  f <- chain_fit(100, 1)
  g <- f$graph
  a0 <- c(f$Regression, f$Covariance)
  y <- as.numeric(f$fitted) + as.numeric(f$residuals)
  fs <- mcglmgraph:::.first_stage_system(f, seq_len(g$n), g$n)
  corr <- totals_dense(f, a0, y, TRUE) - totals_dense(f, a0, y, FALSE)
  expect_equal(unname(colSums(fs$psi)), -corr, tolerance = 1e-6)
  expect_gt(max(abs(corr)), 0.5)
})

test_that("a root innovation variance fixed by the first stage has no variance", {
  skip_if_not_installed("MASS")
  # With Gaussian intercept-only margins the root variance is the average
  # squared Pearson residual, which the corrected first stage fixes at
  # (n - 1) / n in every sample, and without the correction at one.
  for (cr in c(TRUE, FALSE)) {
    f <- chain_fit(100, 1, correct = cr)
    d1 <- unname(f$graph$par["delta_1"])
    expect_equal(d1, if (cr) 99 / 100 else 1, tolerance = 1e-10)
    Vs <- mc_vcov_graph(f, "sandwich")
    Vc <- mc_vcov_graph(f, "two_stage", method = "shortcut")
    Vf <- mc_vcov_graph(f, "two_stage", method = "full")
    sw <- sqrt(Vs["delta_1", "delta_1"])
    expect_lt(sqrt(max(Vc["delta_1", "delta_1"], 0)), 1e-8 * sw)
    expect_lt(sqrt(max(Vf["delta_1", "delta_1"], 0)), 1e-8 * sw)
  }
})

test_that("the shortcut factor is the dispersion derivative times its influence", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("numDeriv")
  # The first stage fixes mean(r_r^2) at m_r, and r_r scales as tau_r^-1/2,
  # so tau_r-hat = mean(e_r^2) / m_r and unit i moves it, on the accumulated
  # scale, by tau_r (r_ri^2 - m_r) / m_r. With B(d, r) = -dU_2d / dtau_r the
  # propagation term is sum_r B(d, r) tau_r (r_ri^2 - m_r) / (n m_r), so the
  # factor applied to r_ri^2 - m_r is Q(d, r) = B(d, r) tau_r / (n m_r).
  f <- chain_fit(100, 3)
  g <- f$graph
  Y <- as.matrix(f$data[, c("y1", "y2", "y3")])
  e <- sweep(Y, 2L, colMeans(Y))
  tau <- f$Covariance[4:6]
  m <- colMeans(as.matrix(g$residuals)^2)
  expect_equal(unname(m), rep(99 / 100, 3), tolerance = 1e-10)
  Sinv <- solve(g$Sigma_b)
  U2 <- function(t) {
    Ec <- scale(sweep(e, 2L, sqrt(t), "/"), center = TRUE, scale = FALSE)
    vapply(g$Dlist, function(Dd) {
      A <- Sinv %*% Dd %*% Sinv
      -0.5 * sum(sum(diag(Sinv %*% Dd)) - rowSums((Ec %*% A) * Ec))
    }, numeric(1))
  }
  B <- -numDeriv::jacobian(U2, tau)
  want <- sweep(B, 2L, tau / (g$n * m), "*")
  Vc <- mc_vcov_graph(f, "two_stage", method = "shortcut")
  expect_equal(bare(attr(Vc, "Q")), bare(want), tolerance = 1e-7)
  # the factor before the division is the one of cpp_two_stage_q
  Q0 <- mcglmgraph:::cpp_two_stage_q(as.matrix(g$K_b), g$Dlist)
  expect_equal(bare(attr(Vc, "Q")), bare(sweep(Q0, 2L, m, "/")),
               tolerance = 1e-12)
})

test_that("without immorality the shortcut and the full form agree", {
  skip_if_not_installed("MASS")
  # Gaussian intercept-only margins: the full form, with the observed
  # sensitivity of the equations as solved, and the shortcut give the same
  # standard errors on a chain DAG and on a GGM
  for (type in c("dag", "ggm")) {
    f <- chain_fit(100, 5, type = type)
    sc <- sqrt(pmax(diag(mc_vcov_graph(f, "two_stage", method = "shortcut")), 0))
    fo <- sqrt(pmax(diag(mc_vcov_graph(f, "two_stage", method = "full")), 0))
    expect_lt(max(abs(sc - fo)) / max(sc), 1e-6)
  }
})

## Monte Carlo acceptance checks: the full form against the sampling standard
## deviation of the estimates over replicated data sets, and against the
## bootstrap that refits both stages. They take minutes.
## The ratio is of the median standard error to the Monte Carlo deviation: in a
## few replicates the first stage converges to a weakly identified root, where
## the observed sensitivity is near singular and the propagated variance is far
## too large, and a root mean square would be driven by those alone.
mc_ratio <- function(fits, cluster_of = NULL) {
  ok <- Filter(function(z) !inherits(z, "try-error"), fits)
  P <- do.call(rbind, lapply(ok, function(f) f$graph$par))
  se <- do.call(rbind, lapply(ok, function(f) sqrt(diag(mc_vcov_graph(
    f, "two_stage", cluster = if (is.null(cluster_of)) NULL else cluster_of(f))))))
  list(n = length(ok),
       ratio = apply(se, 2, stats::median) / apply(P, 2, stats::sd))
}

## The first seed from `from` whose Poisson-Tweedie fit converges.
first_pt_fit <- function(n, from) {
  for (s in from + 0:50) {
    f <- tryCatch(pt_fit(n, s), error = function(e) NULL)
    if (!is.null(f)) return(f)
  }
  stop("no converging seed")
}

test_that("with the power estimated, the full form tracks the sampling spread", {
  skip_on_cran()
  skip_long()
  skip_if_not_installed("MASS")
  cores <- max(1L, min(8L, parallel::detectCores() - 1L))
  fits <- parallel::mclapply(1:400, function(s) try(pt_fit(1000, 20000 + s),
                                                    silent = TRUE),
                             mc.cores = cores)
  mr <- mc_ratio(fits)
  expect_gte(mr$n, 100L)
  expect_true(all(mr$ratio > 0.8 & mr$ratio < 1.25))
  # and against the bootstrap of both stages on one data set
  f <- first_pt_fit(1000, 20001)
  set.seed(3)
  bt <- mc_graph_bootstrap(f, B = 200, workers = cores)
  P <- length(f$graph$par)
  r <- sqrt(diag(mc_vcov_graph(f, "two_stage"))) / bt$se[seq_len(P)]
  expect_true(all(r > 0.7 & r < 1.4))
})

test_that("with blocks nested in clusters, the full form tracks the sampling spread", {
  skip_on_cran()
  skip_long()
  skip_if_not_installed("MASS")
  cores <- max(1L, min(8L, parallel::detectCores() - 1L))
  fits <- parallel::mclapply(1:400, function(s) try(blocked_fit(1000 + s),
                                                    silent = TRUE),
                             mc.cores = cores)
  mr <- mc_ratio(fits, cluster_of = function(f) f$data$cl)
  expect_gte(mr$n, 300L)
  expect_true(all(mr$ratio > 0.85 & mr$ratio < 1.15))
  # the bootstrap by clusters of E2 is the other yardstick; with 40 clusters
  # it runs somewhat wide of the sampling spread, so the band is wider
  f <- blocked_fit(99)
  set.seed(8)
  bt <- mc_graph_bootstrap(f, B = 200, cluster = f$data$cl, workers = cores)
  P <- length(f$graph$par)
  r <- sqrt(diag(mc_vcov_graph(f, "two_stage", cluster = f$data$cl))) /
    bt$se[seq_len(P)]
  expect_true(all(r > 0.7 & r < 1.3))
})
