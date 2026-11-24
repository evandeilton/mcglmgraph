# The unit-diagonal constraint of the directed layer, and the two validity
# gates that moved to C++.

test_that("the recursion produces a covariance with unit diagonal", {
  cfgs <- list(
    list(g = c(0.6, 0.6), R = 3L, ed = rbind(c(1, 2), c(2, 3))),
    list(g = c(0.5, -0.4, 0.3), R = 3L, ed = rbind(c(1, 2), c(1, 3), c(2, 3))),
    list(g = c(0.4, 0.3, -0.2, 0.25), R = 4L,
         ed = rbind(c(1, 2), c(1, 3), c(2, 4), c(3, 4)))
  )
  for (cf in cfgs) {
    u <- mc_dag_unit_diagonal(cf$g, R = cf$R, edges = cf$ed)
    expect_true(u$admissible)
    expect_equal(diag(u$Sigma_b), rep(1, cf$R), tolerance = 1e-12)
    # the closed form of the constraint, parent set by parent set
    for (j in seq_len(cf$R)) {
      P <- cf$ed[cf$ed[, 2] == j, 1]
      g_j <- cf$g[cf$ed[, 2] == j]
      want <- if (!length(P)) 1 else
        1 - drop(t(g_j) %*% u$Sigma_b[P, P, drop = FALSE] %*% g_j)
      expect_equal(u$delta[j], want, tolerance = 1e-10)
    }
    # and the proportion of variance the parents account for
    expect_equal(u$r_squared, 1 - u$delta, tolerance = 1e-12)
  }
})

test_that("the admissible set contains the origin and is a proper subset", {
  ed <- rbind(c(1, 3), c(2, 3))
  expect_true(mc_dag_unit_diagonal(c(0, 0), R = 3, edges = ed)$admissible)
  # two strong parents leave the child with negative innovation variance
  bad <- mc_dag_unit_diagonal(c(0.9, 0.9), R = 3, edges = ed)
  expect_false(bad$admissible)
  expect_lt(bad$min_delta, 0)
  expect_true(all(is.na(bad$to_unit_innovations)))
  # the boundary is where the child keeps nothing of its own
  edge <- mc_dag_unit_diagonal(c(sqrt(0.5), sqrt(0.5)), R = 3, edges = ed)
  expect_equal(edge$min_delta, 0, tolerance = 1e-12)
  expect_false(edge$admissible)
})

test_that("every orbit has one representative with unit innovations", {
  g <- c(0.5, -0.4, 0.3)
  ed <- rbind(c(1, 2), c(1, 3), c(2, 3))
  u <- mc_dag_unit_diagonal(g, R = 3, edges = ed)
  a <- u$to_unit_innovations
  A <- diag(a)
  Sa <- A %*% u$Sigma_b %*% A
  # fitting the rescaled covariance under the same support gives unit deltas
  fit <- mc_dag_closed_root(Sa, R = 3, edges = ed)
  expect_equal(fit$delta, rep(1, 3), tolerance = 1e-10)
  # and the coefficients moved as the scale proposition says
  expect_equal(fit$gamma, g * a[ed[, 2]] / a[ed[, 1]], tolerance = 1e-10)
})

test_that("the constrained covariance is recovered by the closed-form fit", {
  g <- c(0.45, 0.3, -0.25)
  ed <- rbind(c(1, 2), c(1, 3), c(2, 3))
  u <- mc_dag_unit_diagonal(g, R = 3, edges = ed)
  fit <- mc_dag_closed_root(u$Sigma_b, R = 3, edges = ed)
  expect_equal(fit$gamma, g, tolerance = 1e-10)
  expect_equal(fit$delta, u$delta, tolerance = 1e-10)
})

test_that("the contracts of the arguments are enforced", {
  expect_error(mc_dag_unit_diagonal(c(0.5), R = 3,
                                    edges = rbind(c(1, 2), c(2, 3))),
               "exactly 2 finite")
  expect_error(mc_dag_unit_diagonal(c(0.5, NA), R = 3,
                                    edges = rbind(c(1, 2), c(2, 3))),
               "finite")
  expect_error(mc_dag_unit_diagonal(numeric(0), R = 0), "positive integer")
})

## ---- the gates, against plain-R versions ----------------------------------

ref_cycle <- function(A) {
  R <- nrow(A); removed <- rep(FALSE, R)
  repeat {
    live <- which(!removed)
    if (!length(live)) break
    dir <- (A == 1L) & (t(A) == 0L)
    indeg <- colSums(dir[live, live, drop = FALSE])
    src <- live[indeg == 0]
    if (!length(src)) break
    removed[src] <- TRUE
  }
  any(!removed)
}

ref_extension <- function(A) {
  R <- nrow(A); alive <- rep(TRUE, R)
  for (step in seq_len(R)) {
    found <- FALSE
    for (x in which(alive)) {
      out <- which(alive & A[x, ] == 1L & A[, x] == 0L)
      if (length(out)) next
      und <- which(alive & A[x, ] == 1L & A[, x] == 1L)
      nb <- which(alive & (A[x, ] == 1L | A[, x] == 1L))
      okx <- TRUE
      for (y in und) {
        others <- setdiff(nb, y)
        if (length(others) &&
            any(A[y, others] == 0L & A[others, y] == 0L)) { okx <- FALSE; break }
      }
      if (okx) { alive[x] <- FALSE; found <- TRUE; break }
    }
    if (!found) break
  }
  !any(alive)
}

test_that("the gates in C++ agree with the plain-R versions", {
  set.seed(17)
  for (k in 1:300) {
    R <- sample(3:6, 1)
    A <- matrix(0L, R, R)
    for (a in seq_len(R - 1L)) for (b in (a + 1L):R) {
      u <- runif(1)
      if (u < 0.25) { A[a, b] <- 1L }
      else if (u < 0.5) { A[b, a] <- 1L }
      else if (u < 0.75) { A[a, b] <- A[b, a] <- 1L }
    }
    o <- sample(R); A <- A[o, o]
    expect_identical(mcglmgraph:::.has_directed_cycle(A), ref_cycle(A))
    expect_identical(mcglmgraph:::.has_consistent_extension(A), ref_extension(A))
  }
})

test_that("the gates answer the cases they exist for", {
  # a directed cycle
  cyc <- matrix(0L, 3, 3); cyc[1, 2] <- cyc[2, 3] <- cyc[3, 1] <- 1L
  expect_true(mcglmgraph:::.has_directed_cycle(cyc))
  expect_false(mcglmgraph:::.has_consistent_extension(cyc))
  # a CPDAG always extends
  cp <- matrix(0L, 3, 3); cp[1, 3] <- cp[2, 3] <- 1L
  expect_false(mcglmgraph:::.has_directed_cycle(cp))
  expect_true(mcglmgraph:::.has_consistent_extension(cp))
  # a partially directed graph with no extension: 1 -> 2 -> 3 -> 4 and 1 - 4.
  # Orienting 4 -> 1 closes a cycle; orienting 1 -> 4 creates the collider
  # 1 -> 4 <- 3 with 1 and 3 nonadjacent, which the graph did not have. Neither
  # direction is available, so no DAG extends it.
  bad <- matrix(0L, 4, 4)
  bad[1, 2] <- bad[2, 3] <- bad[3, 4] <- 1L
  bad[1, 4] <- bad[4, 1] <- 1L
  expect_false(mcglmgraph:::.has_consistent_extension(bad))
  expect_false(mcglmgraph:::.has_directed_cycle(bad))   # the arcs alone are fine
})
