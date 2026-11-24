# Greedy equivalence search over Markov equivalence classes.
#
# A search over classes is well posed only if the score cannot tell two members
# of a class apart. That is checked first, since everything else rests on it.

collider_fit <- function(n = 800, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.8 * z1 + 0.8 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

dag_fit <- function(B, n = 3000, seed = 2) {
  set.seed(seed)
  R <- nrow(B)
  IB <- solve(diag(R) - B)
  Y <- MASS::mvrnorm(n, rep(0, R), IB %*% t(IB))
  d <- as.data.frame(Y); names(d) <- paste0("y", seq_len(R))
  d$x <- rnorm(n)
  mcglm_graph(lapply(seq_len(R), function(i)
    stats::as.formula(paste0("y", i, " ~ x"))), data = d)
}

# Total penalized score of a DAG, summing the local contributions.
dag_score <- function(S, n, dag, penalty) {
  sum(vapply(seq_len(nrow(dag)), function(j)
    mcglmgraph:::cpp_ges_local_score(
      S, n, j, which(dag[, j] == 1L & dag[j, ] == 0L), penalty),
    numeric(1)))
}

test_that("the score cannot tell two members of a class apart", {
  skip_if_not_installed("MASS")
  f <- collider_fit()
  S <- f$graph$Shat; n <- f$graph$n
  set.seed(4)
  for (k in 1:12) {
    R <- 4L
    A <- matrix(0L, R, R)
    for (a in seq_len(R - 1L)) for (b in (a + 1L):R)
      if (runif(1) < 0.6) A[a, b] <- 1L
    cp <- mc_true_cpdag(t(A))
    mem <- mc_class_members(cp)
    if (length(mem) < 2L) next
    S4 <- stats::cov(MASS::mvrnorm(500, rep(0, R), diag(R)))
    sc <- vapply(mem, function(D) dag_score(S4, 500, D, log(500)), numeric(1))
    expect_equal(sc, rep(sc[1L], length(sc)), tolerance = 1e-8)
  }
})

test_that("the local score is the penalized residual variance of a regression", {
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 400)
  S <- f$graph$Shat; n <- f$graph$n
  for (j in 1:3) for (pa in list(integer(0), setdiff(1:3, j), j %% 3 + 1)) {
    got <- mcglmgraph:::cpp_ges_local_score(S, n, j, as.integer(pa), log(n))
    resid <- if (!length(pa)) S[j, j] else
      S[j, j] - drop(S[j, pa] %*% solve(S[pa, pa, drop = FALSE], S[pa, j]))
    want <- -0.5 * n * log(resid) - 0.5 * log(n) * (length(pa) + 1)
    expect_equal(got, want, tolerance = 1e-10)
  }
})

test_that("a DAG consistent with the partially directed graph is produced", {
  set.seed(6)
  for (k in 1:20) {
    R <- sample(3:5, 1)
    A <- matrix(0L, R, R)
    for (a in seq_len(R - 1L)) for (b in (a + 1L):R)
      if (runif(1) < 0.6) A[a, b] <- 1L
    cp <- mc_true_cpdag(t(A))
    res <- mcglmgraph:::cpp_pdag_to_dag(cp)
    expect_true(res$ok)
    dag <- matrix(as.integer(res$dag), R, R)
    expect_true(all((dag == 1L) + (t(dag) == 1L) <= 1L))   # fully directed
    expect_identical(unname((dag + t(dag)) > 0), unname((cp + t(cp)) > 0))
    # every arc the class compelled is kept
    comp <- which(cp == 1L & t(cp) == 0L, arr.ind = TRUE)
    for (e in seq_len(nrow(comp)))
      expect_equal(dag[comp[e, 1], comp[e, 2]], 1L)
    # and it belongs to the class it came from
    expect_identical(unname(mcglmgraph:::.dag_to_cpdag(dag)), unname(cp[, ]))
  }
})

test_that("the collider is recovered and the output is a valid class", {
  skip_if_not_installed("MASS")
  f <- collider_fit()
  out <- mc_ges(f)
  expect_s3_class(out, "mcg_ges")
  T <- diag(3); T[3, 1] <- T[3, 2] <- -0.8
  want <- mc_true_cpdag(T)
  expect_identical(unname(out$cpdag), unname(want[, ]))
  # the output is a fixed point of Meek and admits a consistent extension
  fixed <- mc_orient_meek(out$cpdag)
  attributes(fixed) <- list(dim = dim(fixed))
  expect_identical(fixed, unname(out$cpdag))
  expect_true(mcglmgraph:::.has_consistent_extension(out$cpdag))
  expect_true(is.finite(out$score))
  expect_gt(out$steps[["forward"]], 0L)
})

test_that("the class of a chain is recovered, undirected", {
  skip_if_not_installed("MASS")
  B <- matrix(0, 3, 3); B[2, 1] <- 0.8; B[3, 2] <- 0.8
  f <- dag_fit(B)
  out <- mc_ges(f)
  # a chain has no collider, so its class leaves both edges undirected
  expect_identical(unname(out$cpdag), unname(mc_true_cpdag(-B)[, ]))
  expect_true(all((out$cpdag == 1L) == (t(out$cpdag) == 1L)))
})

test_that("the search raises the score at every step it takes", {
  skip_if_not_installed("MASS")
  f <- collider_fit()
  out <- mc_ges(f)
  expect_true(all(out$trace$gain > 0))
  empty <- matrix(0L, 3, 3)
  s0 <- mcglmgraph:::.class_score(f$graph$Shat, f$graph$n, empty, log(f$graph$n))
  expect_gt(out$score, s0)
  # starting from the answer leaves it alone
  again <- mc_ges(f, start = out$cpdag)
  expect_identical(unname(again$cpdag), unname(out$cpdag))
})

test_that("a lighter penalty keeps more edges", {
  skip_if_not_installed("MASS")
  B <- matrix(0, 4, 4)
  B[2, 1] <- 0.7; B[3, 2] <- 0.7; B[4, 1] <- 0.25
  f <- dag_fit(B, n = 400, seed = 8)
  e <- function(o) sum(o$cpdag == 1L | t(o$cpdag) == 1L) / 2
  heavy <- mc_ges(f, penalty = 2 * log(f$graph$n))
  light <- mc_ges(f, penalty = 0)
  expect_gte(e(light), e(heavy))
  expect_equal(e(light), 6)                  # no penalty: the complete graph
})

test_that("the search agrees with the constraint-based route on clear cases", {
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 2000, seed = 5)
  expect_identical(unname(mc_ges(f)$cpdag), unname(mc_graph_pipeline(f)$cpdag))
})

test_that("the contracts of the arguments are enforced", {
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 200)
  expect_error(mc_ges(f, penalty = -1), "non-negative")
  expect_error(mc_ges(f, max_set = -1), "non-negative")
  expect_error(mc_ges(f, max_steps = 0), "positive integer")
  expect_error(mc_ges(f, start = matrix(0L, 2, 2)), "square 0/1 adjacency")
  expect_error(mc_ges(f, start = matrix(2L, 3, 3)), "square 0/1 adjacency")
  expect_error(mc_ges(structure(list(), class = "lm")), "inherits")
})
