# Comparing an estimated structure with a known one.
#
# The counts are checked against hand-worked cases and against their own
# definitions, recomputed here, rather than against the implementation.

bare <- function(M) { M <- as.matrix(M); attributes(M) <- list(dim = dim(M)); M }

test_that("the skeleton counts are the confusion table of the pairs", {
  # three responses, pairs in order (1,2), (1,3), (2,3)
  truth <- rbind(c(1, 2), c(2, 3))            # two edges
  est <- rbind(c(1, 2), c(1, 3))              # one hit, one false, one missed
  m <- mc_structure_metrics(est, truth, R = 3)
  expect_equal(unname(m[c("tp", "fp", "fn", "tn")]), c(1, 1, 1, 0))
  expect_equal(unname(m[["tpr"]]), 1 / 2)
  expect_equal(unname(m[["fpr"]]), 1 / 1)
  expect_equal(unname(m[["precision"]]), 1 / 2)
  expect_equal(unname(m[["f1"]]), 2 * 1 / (2 * 1 + 1 + 1))
  expect_equal(unname(m[["shd"]]), 2)
})

test_that("a structure compared with itself is recovered exactly", {
  set.seed(1)
  for (R in 3:6) {
    A <- matrix(0L, R, R)
    ut <- which(upper.tri(A), arr.ind = TRUE)
    keep <- ut[runif(nrow(ut)) < 0.6, , drop = FALSE]
    if (!nrow(keep)) next
    A[keep] <- 1L; A[keep[, c(2, 1), drop = FALSE]] <- 1L
    m <- mc_structure_metrics(A, A)
    expect_equal(unname(m[["f1"]]), 1)
    expect_equal(unname(m[["shd"]]), 0)
    expect_equal(unname(m[["fp"]]), 0)
    expect_equal(unname(m[["fn"]]), 0)
    cm <- mc_cpdag_metrics(A, A)
    expect_equal(unname(cm[["cpdag_shd"]]), 0)
    expect_equal(unname(cm[["cpdag_skeleton_f1"]]), 1)
  }
})

test_that("an empty denominator gives NA rather than a number", {
  empty <- matrix(0L, 3, 3)
  full <- matrix(1L, 3, 3); diag(full) <- 0L
  # nothing selected: precision undefined, recall zero
  m <- mc_structure_metrics(empty, full)
  expect_true(is.na(m[["precision"]]))
  expect_equal(unname(m[["tpr"]]), 0)
  # nothing to find: recall undefined, false positive rate one
  m2 <- mc_structure_metrics(full, empty)
  expect_true(is.na(m2[["tpr"]]))
  expect_equal(unname(m2[["fpr"]]), 1)
  expect_equal(unname(m2[["precision"]]), 0)
})

test_that("the CPDAG measures separate the skeleton from the orientation", {
  # truth: 1 -> 3 <- 2, a collider, so both arcs are compelled
  truth <- matrix(0L, 3, 3)
  truth[1, 3] <- truth[2, 3] <- 1L
  expect_equal(unname(mc_cpdag_metrics(truth, truth)[["cpdag_shd"]]), 0)

  # an estimate that finds the skeleton but leaves it undirected
  est <- matrix(0L, 3, 3)
  est[1, 3] <- est[3, 1] <- est[2, 3] <- est[3, 2] <- 1L
  m <- mc_cpdag_metrics(est, truth)
  expect_equal(unname(m[["cpdag_skeleton_f1"]]), 1)      # skeleton is right
  expect_equal(unname(m[["cpdag_shd"]]), 2)              # both pairs disagree
  expect_equal(unname(m[["cpdag_orient_recall"]]), 0)    # nothing oriented
  expect_true(is.na(m[["cpdag_vstructure_precision"]]))  # none declared
  expect_equal(unname(m[["cpdag_vstructure_recall"]]), 0)
  expect_true(is.na(m[["cpdag_overorient"]]))            # truth directs all
})

test_that("over-orientation counts directions the truth does not have", {
  # truth: a chain with no collider, so the CPDAG is fully undirected
  truth <- matrix(0L, 3, 3)
  truth[1, 2] <- truth[2, 1] <- truth[2, 3] <- truth[3, 2] <- 1L
  est <- matrix(0L, 3, 3)
  est[1, 2] <- 1L                       # one of the two edges invented
  est[2, 3] <- est[3, 2] <- 1L
  m <- mc_cpdag_metrics(est, truth)
  expect_equal(unname(m[["cpdag_overorient"]]), 1 / 2)
  expect_true(is.na(m[["cpdag_orient_recall"]]))         # nothing to recover
  expect_equal(unname(m[["cpdag_skeleton_f1"]]), 1)
  expect_equal(unname(m[["cpdag_shd"]]), 1)
})

test_that("v-structure precision and recall count the triples", {
  truth <- matrix(0L, 4, 4)
  truth[1, 3] <- truth[2, 3] <- 1L      # 1 -> 3 <- 2
  truth[3, 4] <- 1L                     # compelled by propagation
  est <- matrix(0L, 4, 4)
  est[1, 3] <- est[2, 3] <- 1L
  est[4, 3] <- 1L                       # a second, wrong collider at 3
  m <- mc_cpdag_metrics(est, truth)
  # the truth has one v-structure; the estimate declares three at node 3
  expect_equal(unname(m[["cpdag_vstructure_recall"]]), 1)
  expect_equal(unname(m[["cpdag_vstructure_precision"]]), 1 / 3)
})

test_that("the CPDAG of a known DAG has the arcs the class compels", {
  # a collider is compelled
  T1 <- diag(3); T1[3, 1] <- -0.5; T1[3, 2] <- -0.4
  cp1 <- mc_true_cpdag(T1)
  expect_equal(unname(cp1[1, 3]), 1L); expect_equal(unname(cp1[3, 1]), 0L)
  expect_equal(unname(cp1[2, 3]), 1L); expect_equal(unname(cp1[3, 2]), 0L)
  expect_equal(unname(cp1[1, 2]), 0L)

  # a chain has no collider, so nothing is compelled
  T2 <- diag(3); T2[2, 1] <- -0.5; T2[3, 2] <- -0.4
  cp2 <- mc_true_cpdag(T2)
  expect_true(all((cp2 == 1L) == (t(cp2) == 1L)))

  # a complete DAG has no immorality either: every parent set is adjacent
  T3 <- diag(4)
  for (j in 2:4) for (i in 1:(j - 1)) T3[j, i] <- -0.3
  cp3 <- mc_true_cpdag(T3)
  expect_true(all((cp3 == 1L) == (t(cp3) == 1L)))
  expect_equal(sum(cp3), 4 * 3)                 # every pair joined

  # tol decides what counts as an arc
  T4 <- diag(3); T4[2, 1] <- -1e-12
  expect_equal(sum(mc_true_cpdag(T4)), 2)
  expect_equal(sum(mc_true_cpdag(T4, tol = 1e-8)), 0)
})

test_that("the CPDAG of a DAG is what the discovery recovers in large samples", {
  skip_if_not_installed("MASS")
  set.seed(19); n <- 3000
  x <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n)
  z3 <- 0.8 * z1 + 0.8 * z2 + rnorm(n)          # 1 -> 3 <- 2
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
  cp <- mc_graph_pipeline(f)
  T <- diag(3); T[3, 1] <- -0.8; T[3, 2] <- -0.8
  want <- mc_true_cpdag(T)
  expect_identical(bare(cp$cpdag), bare(want))
  m <- mc_cpdag_metrics(cp$cpdag, want)
  expect_equal(unname(m[["cpdag_shd"]]), 0)
  expect_equal(unname(m[["cpdag_vstructure_recall"]]), 1)
})

test_that("the inputs may be adjacencies, indicators or edge lists alike", {
  A <- matrix(0L, 3, 3)
  A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1L
  ed <- rbind(c(1, 2), c(2, 3))
  ind <- c(1L, 0L, 1L)                   # pairs (1,2), (1,3), (2,3)
  a <- mc_structure_metrics(A, A)
  b <- mc_structure_metrics(ed, ed, R = 3)
  c3 <- mc_structure_metrics(ind, ind)
  expect_equal(unname(a), unname(b))
  expect_equal(unname(a), unname(c3))
  # a directed adjacency is reduced to its skeleton
  D <- matrix(0L, 3, 3); D[1, 2] <- 1L; D[2, 3] <- 1L
  expect_equal(unname(mc_structure_metrics(D, A)[["shd"]]), 0)
})

test_that("the contracts of the arguments are enforced", {
  A <- matrix(0L, 3, 3)
  expect_error(mc_cpdag_metrics(matrix(2L, 3, 3), A), "zeros and ones")
  expect_error(mc_cpdag_metrics(matrix(0L, 2, 3), A), "square")
  expect_error(mc_cpdag_metrics(matrix(0L, 4, 4), A),
               "different number of responses")
  expect_error(mc_structure_metrics(c(1L, 0L), c(1L, 0L, 1L)),
               "different number of pairs")
  expect_error(mc_true_cpdag(matrix(NA_real_, 2, 2)), "finite numeric")
  expect_error(mc_true_cpdag(diag(3), tol = -1), "non-negative")
})
