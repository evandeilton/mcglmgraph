# Wald tests of linear hypotheses on the graph layer.
#
# Targets rebuilt from the theory, in plain R:
#   statistic  W = n (L th - c)' [L J^-1 V J^-1 L']^-1 (L th - c), with J and V
#              per unit, J_ab = (1/2) tr(S^-1 D_a S^-1 D_b)
#              and V = n^-1 sum_i psi_i psi_i'
#   invariance W(L, c) = W(A L, A c) for any invertible A
#   rank       redundant rows add no degree of freedom
#   size       under H0 with a fixed support, W is chi-squared(rank L)

bare <- function(M) {
  M <- as.matrix(M)
  attributes(M) <- list(dim = dim(M))
  M
}

wald_fit <- function(n = 400, seed = 11, R = 3) {
  set.seed(seed)
  K <- diag(R) * 2
  K[1, 2] <- K[2, 1] <- -0.7
  if (R >= 4) K[3, 4] <- K[4, 3] <- -0.6
  Y <- MASS::mvrnorm(n, rep(0, R), solve(K))
  d <- as.data.frame(Y); names(d) <- paste0("y", seq_len(R))
  d$x <- rnorm(n)
  mcglm_graph(lapply(paste0("y", seq_len(R)),
                     function(v) stats::reformulate("x", response = v)),
              data = d)
}

## J and V per unit, from their definitions.
theory_JV <- function(f) {
  g <- f$graph
  Sinv <- solve(g$Sigma_b)
  P <- length(g$Dlist)
  J <- matrix(0, P, P)
  for (a in seq_len(P)) for (b in seq_len(P))
    J[a, b] <- 0.5 * sum(diag(Sinv %*% g$Dlist[[a]] %*% Sinv %*% g$Dlist[[b]]))
  Ec <- scale(g$residuals, center = TRUE, scale = FALSE)
  Psi <- sapply(seq_len(P), function(d) {
    A <- Sinv %*% g$Dlist[[d]] %*% Sinv
    -0.5 * (sum(diag(Sinv %*% g$Dlist[[d]])) - rowSums((Ec %*% A) * Ec))
  })
  list(J = J, V = crossprod(Psi) / g$n)
}

test_that("the statistic is eq. wald on the per-unit scale", {
  skip_if_not_installed("MASS")
  f <- wald_fit()
  g <- f$graph
  th <- theory_JV(f)
  L <- rbind(c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 1, -1))
  cc <- c(-0.5, 0.1)
  d <- L %*% g$par - cc
  Jinv <- solve(th$J)
  W_th <- g$n * drop(t(d) %*% solve(L %*% Jinv %*% th$V %*% Jinv %*% t(L)) %*% d)
  w <- mc_wald_graph(f, L, cc, type = "sandwich")
  expect_equal(unname(w$statistic), W_th, tolerance = 1e-8)
  expect_identical(unname(w$parameter), 2L)
  expect_equal(w$p.value, stats::pchisq(W_th, 2, lower.tail = FALSE),
               tolerance = 1e-8)
  # model-based: V replaced by J
  W_m <- g$n * drop(t(d) %*% solve(L %*% Jinv %*% t(L)) %*% d)
  expect_equal(unname(mc_wald_graph(f, L, cc, type = "model")$statistic), W_m,
               tolerance = 1e-8)
})

test_that("selecting every edge is the joint test of independence", {
  skip_if_not_installed("MASS")
  f <- wald_fit()
  ep <- f$graph$edge_pos
  L <- matrix(0, length(ep), length(f$graph$par))
  L[cbind(seq_along(ep), ep)] <- 1
  for (type in c("model", "sandwich", "two_stage")) {
    a <- mc_wald_graph(f, L, type = type)
    b <- mc_test_independence(f, type = type)
    expect_identical(unname(a$statistic), unname(b$statistic), info = type)
    expect_identical(unname(a$parameter), unname(b$parameter), info = type)
  }
  # robust = TRUE and type = "sandwich" name the same covariance
  expect_identical(mc_test_independence(f, robust = TRUE)$statistic,
                   mc_test_independence(f, type = "sandwich")$statistic)
  expect_identical(mc_test_independence(f)$statistic,
                   mc_test_independence(f, type = "model")$statistic)
  expect_error(mc_test_independence(f, robust = TRUE, type = "model"),
               "disagree")
  # the names select the same parameters
  nm <- names(f$graph$par)[ep]
  expect_identical(mc_wald_graph(f, nm)$statistic, mc_wald_graph(f, L)$statistic)
})

test_that("the statistic is invariant to recombining the rows", {
  skip_if_not_installed("MASS")
  f <- wald_fit(R = 4)
  D <- length(f$graph$par)
  set.seed(3)
  L <- matrix(stats::rnorm(3 * D), 3, D)
  cc <- stats::rnorm(3, sd = 0.1)
  A <- matrix(stats::rnorm(9), 3, 3)
  w1 <- mc_wald_graph(f, L, cc)
  w2 <- mc_wald_graph(f, A %*% L, drop(A %*% cc))
  expect_equal(unname(w1$statistic), unname(w2$statistic), tolerance = 1e-9)
})

test_that("redundant rows add no degree of freedom", {
  skip_if_not_installed("MASS")
  f <- wald_fit()
  L <- rbind(c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 0, 1, 0))
  L3 <- rbind(L, L[1, ] + L[2, ])
  w2 <- mc_wald_graph(f, L)
  w3 <- mc_wald_graph(f, L3)
  expect_identical(unname(w3$parameter), 2L)
  expect_equal(unname(w3$statistic), unname(w2$statistic), tolerance = 1e-9)
  # a right-hand side the redundant row contradicts is refused
  expect_error(mc_wald_graph(f, L3, c(0, 0, 1)), "not consistent")
})

test_that("the arguments are checked", {
  skip_if_not_installed("MASS")
  f <- wald_fit(n = 200)
  expect_error(mc_wald_graph(f, c("eta_12", "eta_99")), "does not have")
  expect_error(mc_wald_graph(f, rep(1, 4)), "one column per graph parameter")
  expect_error(mc_wald_graph(f, "eta_12", rhs = c(0, 1)), "one per row")
  expect_error(mc_wald_graph(f, matrix(0, 0, 6)), "no rows")
})

test_that("the directed layer is tested through the same statistic", {
  skip_if_not_installed("MASS")
  set.seed(4); n <- 400; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag"), data = d)
  co <- mc_dag_coefficients(f, variability = "empirical")
  # one arc at a time, the statistic is the square of the reported z
  for (e in seq_len(nrow(co$edges))) {
    L <- numeric(length(f$graph$par)); L[e] <- 1
    w <- mc_wald_graph(f, L, type = "sandwich")
    expect_equal(unname(w$statistic), co$edges$statistic[e]^2,
                 tolerance = 1e-10)
  }
})

test_that("under H0 with a fixed support the size is nominal", {
  skip_on_cran()
  skip_long()
  skip_if_not_installed("MASS")
  # Gaussian residuals, n = 500, eta_13 = eta_23 = 0 in the population; the
  # two-degree-of-freedom test of both is referred to chi-squared(2)
  B <- 400
  rej <- vapply(seq_len(B), function(b) {
    f <- wald_fit(n = 500, seed = 1000 + b)
    w <- mc_wald_graph(f, c("eta_13", "eta_23"), type = "sandwich")
    w$p.value < 0.05
  }, logical(1))
  se <- sqrt(0.05 * 0.95 / B)
  expect_lt(abs(mean(rej) - 0.05), 3 * se)
})

test_that("the adequacy of a support is the Wald test of its omitted edges", {
  skip_if_not_installed("MASS")
  f <- wald_fit(n = 500, R = 4)          # edges 1-2 and 3-4 in the population
  chain <- mcglmgraph:::.refit_between(
    f, mc_between_graph("ggm", edges = rbind(c(1, 2), c(3, 4))))
  a <- mc_test_adequacy(chain)
  expect_identical(unname(a$parameter), 4L)
  expect_true(a$descriptive_after_selection)
  expect_equal(nrow(a$omitted), 4L)
  # the same number, formed by hand on the saturated fit
  sat <- mcglmgraph:::.refit_between(f, mc_between_graph("ggm"))
  nm <- c("eta_13", "eta_14", "eta_23", "eta_24")
  expect_identical(unname(a$statistic),
                   unname(mc_wald_graph(sat, nm)$statistic))
  # the true support is compatible with the data; one missing edge is not
  expect_gt(a$p.value, 0.001)
  lean <- mcglmgraph:::.refit_between(
    f, mc_between_graph("ggm", edges = rbind(c(3, 4))))
  expect_lt(mc_test_adequacy(lean)$p.value, 1e-6)
  # nothing to test on the saturated graph itself
  expect_error(mc_test_adequacy(sat), "nothing to test")
})

test_that("the adequacy of a directed support is tested along its ordering", {
  skip_if_not_installed("MASS")
  set.seed(8); n <- 500; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.5 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph(
                     "dag", edges = rbind(c(1, 2), c(2, 3))), data = d)
  a <- mc_test_adequacy(f)
  expect_identical(unname(a$parameter), 1L)          # the arc 1 -> 3
  sat <- mcglmgraph:::.refit_between(f, mc_between_graph("dag"))
  L <- numeric(length(sat$graph$par)); L[2] <- 1     # gamma of the pair (1, 3)
  expect_identical(unname(a$statistic), unname(mc_wald_graph(sat, L)$statistic))
  expect_gt(a$p.value, 0.001)
})
