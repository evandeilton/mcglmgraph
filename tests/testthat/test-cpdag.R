test_that("v-structure detection runs and returns a data frame", {
  skip_if_not_installed("MASS")
  set.seed(3); n <- 600
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  vs <- mc_find_vstructures(f)
  expect_s3_class(vs, "data.frame")
  expect_true(all(c("i", "k", "j", "is_collider") %in% names(vs)))
})

test_that("the CPDAG pipeline returns a square adjacency matrix", {
  skip_if_not_installed("MASS")
  set.seed(4); n <- 500
  z1 <- rnorm(n); z2 <- 0.8 * z1 + rnorm(n); z3 <- 0.8 * z1 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  pipe <- mc_graph_pipeline(f)
  expect_s3_class(pipe, "mcg_cpdag")
  expect_true(is.matrix(pipe$cpdag))
  expect_equal(dim(pipe$cpdag), c(3L, 3L))
})

test_that("the native PC route recovers a collider", {
  # y1 -> y3 <- y2 : the parents y1, y2 are marginally independent, so the
  # skeleton must keep y1-y3, y2-y3 and drop y1-y2, and y3 is a collider.
  set.seed(11); n <- 1500
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  p <- mc_graph_pipeline(f)
  A <- p$cpdag
  expect_equal(A[1, 2], 0); expect_equal(A[2, 1], 0)   # no y1-y2 edge
  expect_equal(A[1, 3], 1); expect_equal(A[3, 1], 0)   # y1 -> y3
  expect_equal(A[2, 3], 1); expect_equal(A[3, 2], 0)   # y2 -> y3
  expect_true(any(p$vstructures$is_collider))
})

test_that("the native PC route does NOT invent a collider in a chain", {
  # y1 -> y2 -> y3 : y2 is a mediator, not a collider. The unshielded triple
  # (1, 2, 3) must be left unoriented (Markov equivalence class of the chain).
  set.seed(12); n <- 1500
  z1 <- rnorm(n); z2 <- 0.8 * z1 + rnorm(n); z3 <- 0.8 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  p <- mc_graph_pipeline(f)
  A <- p$cpdag
  expect_equal(A[1, 3], 0); expect_equal(A[3, 1], 0)   # y1, y3 nonadjacent
  expect_false(any(p$vstructures$is_collider))         # no false collider
  expect_true(A[1, 2] == 1 && A[2, 1] == 1)            # y1 - y2 undirected
  expect_true(A[2, 3] == 1 && A[3, 2] == 1)            # y2 - y3 undirected
})

test_that("Meek rules orient a collider configuration", {
  amat <- matrix(0L, 3, 3)
  amat[1, 2] <- 1L                       # 1 -> 2
  amat[3, 2] <- 1L                       # 3 -> 2  (collider at 2)
  out <- mc_orient_meek(amat)
  expect_equal(out[1, 2], 1L); expect_equal(out[2, 1], 0L)
  expect_equal(out[3, 2], 1L); expect_equal(out[2, 3], 0L)
})

test_that("PC-stable skeleton is invariant to response relabeling", {
  K <- diag(4)
  K[1, 2] <- K[2, 1] <- -.35
  K[2, 3] <- K[3, 2] <- -.30
  K[3, 4] <- K[4, 3] <- -.25
  Sigma <- solve(K)
  A <- mc_pc_skeleton(Sigma, n = 10000, alpha = .001)$adjacency
  p <- c(4, 2, 1, 3)
  Ap <- mc_pc_skeleton(Sigma[p, p], n = 10000, alpha = .001)$adjacency
  expect_equal(Ap[order(p), order(p)], A)
})

test_that("mc_test_independence performs the Wald test", {
  skip_if_not_installed("MASS")
  set.seed(5); n <- 300
  K <- diag(3) * 2; K[1,2] <- K[2,1] <- -0.7
  Y <- MASS::mvrnorm(n, rep(0, 3), solve(K))
  d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], x = rnorm(n))
  f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
  a <- mc_test_independence(f)
  expect_s3_class(a, "htest")
  expect_true(a$p.value < 0.05)          # rejects independence (edge 1-2 present)
})
