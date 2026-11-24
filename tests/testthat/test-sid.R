# Structural intervention distance.
#
# The measure is tied to its meaning here rather than to an implementation: in a
# linear Gaussian model generated from the true graph, the population regression
# of j on i and on the parents of i in the estimate returns the total causal
# effect of i on j exactly for the pairs the distance does not count. That is
# what a valid adjustment set is, and it is checked pair by pair.

random_dag <- function(R, p = 0.5) {
  A <- matrix(0L, R, R)
  for (a in seq_len(R - 1L)) for (b in (a + 1L):R)
    if (runif(1) < p) A[a, b] <- 1L
  o <- sample(R)
  A[o, o]
}

# Covariance and total effects of a linear model on the DAG A.
linear_model <- function(A, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  R <- nrow(A)
  B <- matrix(0, R, R)
  B[A == 1L] <- runif(sum(A), 0.4, 0.9) * sample(c(-1, 1), sum(A), TRUE)
  B <- t(B)                       # B[v, p] is the coefficient of p in v
  IB <- solve(diag(R) - B)
  list(Sigma = IB %*% diag(R) %*% t(IB), total = IB)   # total[j, i]
}

# Coefficient of i in the population regression of j on (i, S).
reg_coef <- function(Sigma, j, i, S) {
  idx <- c(i, S)
  drop(solve(Sigma[idx, idx, drop = FALSE], Sigma[idx, j])[1L])
}

test_that("a graph is at distance zero from itself", {
  set.seed(5)
  for (k in 1:200) {
    A <- random_dag(sample(3:6, 1), runif(1, 0.2, 0.8))
    expect_equal(as.numeric(mc_sid(A, A)), 0)
  }
})

test_that("it counts exactly the pairs the adjustment gets wrong", {
  # The criterion is graphical, the check is numerical: in a linear Gaussian
  # model on the true graph, the population regression of j on i and on the
  # parents of i in the estimate returns the total causal effect of i on j
  # exactly for the pairs the distance does not count. The two must agree pair
  # by pair, in both directions.
  set.seed(7)
  tested <- 0L
  for (k in 1:80) {
    R <- sample(3:5, 1)
    G <- random_dag(R, runif(1, 0.3, 0.7))
    H <- random_dag(R, runif(1, 0.3, 0.7))
    lm_ <- linear_model(G)
    bad <- attr(mc_sid(G, H), "pairs")
    for (i in seq_len(R)) {
      S <- which(H[, i] == 1L)
      for (j in setdiff(seq_len(R), i)) {
        want <- lm_$total[j, i]                 # the true total effect
        got <- reg_coef(lm_$Sigma, j, i, S)
        expect_equal(bad[i, j] == 1L, abs(got - want) >= 1e-8,
                     info = sprintf("R=%d i=%d j=%d", R, i, j))
        tested <- tested + 1L
      }
    }
  }
  expect_gt(tested, 800L)
})

test_that("a descendant off every causal path does not invalidate the set", {
  # The back-door criterion refuses any adjustment set containing a descendant
  # of the treatment. The generalized criterion refuses only descendants of the
  # nodes on proper causal paths to the response, and this graph separates the
  # two: 1 -> 4 makes 4 a descendant of 1, but 4 lies on no path from 1 to 3.
  G <- matrix(0L, 4, 4)
  G[1, 2] <- G[2, 3] <- G[1, 4] <- 1L
  H <- matrix(0L, 4, 4)
  H[4, 1] <- 1L                                  # the estimate adjusts 1 for 4
  bad <- attr(mc_sid(G, H), "pairs")
  expect_equal(unname(bad[1, 3]), 0L)            # still valid for the effect on 3
  lm_ <- linear_model(G, seed = 3)
  expect_equal(reg_coef(lm_$Sigma, 3, 1, 4), lm_$total[3, 1], tolerance = 1e-9)
})

test_that("hand-worked cases come out as expected", {
  # truth 1 -> 2 -> 3
  G <- matrix(0L, 3, 3); G[1, 2] <- G[2, 3] <- 1L
  # the empty graph adjusts for nothing: the marginal association of 1 with 2
  # and of 1 with 3 and of 2 with 3 is the total effect, so only the pairs whose
  # association is confounded are counted
  d_empty <- mc_sid(G, matrix(0L, 3, 3))
  expect_equal(as.numeric(d_empty), 3)
  # reversing the first arc: 2 -> 1, 2 -> 3
  H <- matrix(0L, 3, 3); H[2, 1] <- H[2, 3] <- 1L
  expect_equal(as.numeric(mc_sid(G, H)), 3)
  # the distance is not symmetric: the empty graph sends every pair to the
  # unadjusted association, which is wrong under the chain, while under an empty
  # truth the chain's parent sets adjust for variables that do not matter
  empty <- matrix(0L, 3, 3)
  expect_equal(as.numeric(mc_sid(empty, G)), 0)
  expect_false(identical(as.numeric(mc_sid(G, empty)),
                         as.numeric(mc_sid(empty, G))))
  # it never exceeds the number of ordered pairs
  expect_lte(as.numeric(d_empty), 3 * 2)
  expect_gte(as.numeric(d_empty), 0)
})

test_that("a partially directed estimate gives an interval containing its members", {
  G <- matrix(0L, 3, 3); G[1, 2] <- G[2, 3] <- 1L
  cp <- mc_true_cpdag(t(G))          # the class of the chain
  out <- mc_sid(G, cp)
  expect_s3_class(out, "mcg_sid")
  expect_equal(out$n_members, 3L)
  # the truth is one of the members, so the lower end is zero
  expect_equal(out$lower, 0)
  expect_gt(out$upper, 0)
  expect_equal(length(out$values), 3L)
  expect_equal(min(out$values), out$lower)
  expect_equal(max(out$values), out$upper)
  # every member is a DAG whose own distance is the value reported
  mem <- mc_class_members(cp)
  for (m in seq_along(mem))
    expect_equal(as.numeric(mc_sid(G, mem[[m]])), out$values[m])
  txt <- capture.output(print(out))
  expect_true(any(grepl("Interval", txt)))
})

test_that("a CPDAG whose class contains the truth reaches zero", {
  set.seed(11)
  for (k in 1:25) {
    R <- sample(3:5, 1)
    G <- random_dag(R, runif(1, 0.3, 0.7))
    cp <- mc_true_cpdag(t(G))
    out <- mc_sid(G, cp)
    if (inherits(out, "mcg_sid")) expect_equal(out$lower, 0) else
      expect_equal(as.numeric(out), 0)
  }
})

test_that("an mcg_cpdag object is accepted directly", {
  skip_if_not_installed("MASS")
  set.seed(13); n <- 1500; x <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.8 * z1 + 0.8 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  cp <- mc_graph_pipeline(f)
  G <- matrix(0L, 3, 3); G[1, 3] <- G[2, 3] <- 1L
  out <- mc_sid(G, cp)
  # the collider is recovered, so the estimate is the true DAG itself
  expect_equal(as.numeric(out), 0)
})

test_that("the contracts of the arguments are enforced", {
  G <- matrix(0L, 3, 3); G[1, 2] <- 1L
  cyc <- matrix(0L, 3, 3); cyc[1, 2] <- cyc[2, 3] <- cyc[3, 1] <- 1L
  expect_error(mc_sid(cyc, G), "acyclic")
  und <- matrix(0L, 3, 3); und[1, 2] <- und[2, 1] <- 1L
  expect_error(mc_sid(und, G), "undirected edge")
  expect_error(mc_sid(G, matrix(0L, 4, 4)), "different number of responses")
  expect_error(mc_sid(matrix(2L, 3, 3), G), "zeros and ones")
})
