# Orientation of unshielded triples: Meek's rules against their definitions,
# and the conservative and majority rules against the mechanism they exist to
# defend against.
#
# Adjacency convention: amat[a, b] == 1 && amat[b, a] == 0 means a -> b.

is_dir <- function(A, a, b) A[a, b] == 1L && A[b, a] == 0L
is_undir <- function(A, a, b) A[a, b] == 1L && A[b, a] == 1L

## ---- Meek's three rules, one canonical instance each -----------------------

test_that("R1 orients the edge that would otherwise create a collider", {
  # 1 -> 2 - 3 with 1 and 3 nonadjacent: orienting 2 -> 3 avoids 1 -> 2 <- 3
  A <- matrix(0L, 3, 3)
  A[1, 2] <- 1L
  A[2, 3] <- A[3, 2] <- 1L
  out <- mc_orient_meek(A)
  expect_true(is_dir(out, 2, 3))
  expect_equal(attr(out, "meek_blocked"), 0L)
  # and it is withheld when that very triple is ambiguous
  held <- mc_orient_meek(A, ambiguous = rbind(c(1L, 2L, 3L)))
  expect_true(is_undir(held, 2, 3))
  expect_equal(attr(held, "meek_blocked"), 1L)
  # the triple is unordered in its endpoints
  held2 <- mc_orient_meek(A, ambiguous = rbind(c(3L, 2L, 1L)))
  expect_true(is_undir(held2, 2, 3))
})

test_that("R2 orients the edge that would otherwise close a cycle", {
  # 1 -> 3 -> 2 with 1 - 2 undirected: orienting 1 -> 2 avoids the cycle
  A <- matrix(0L, 3, 3)
  A[1, 3] <- 1L
  A[3, 2] <- 1L
  A[1, 2] <- A[2, 1] <- 1L
  out <- mc_orient_meek(A)
  expect_true(is_dir(out, 1, 2))
  # R2 rests on no unshielded triple, so no ambiguity can withhold it
  held <- mc_orient_meek(A, ambiguous = rbind(c(1L, 3L, 2L)))
  expect_true(is_dir(held, 1, 2))
  expect_equal(attr(held, "meek_blocked"), 0L)
})

test_that("R3 orients through two undirected legs into a common head", {
  # 1 - 3, 1 - 4 undirected, 3 -> 2, 4 -> 2, 3 and 4 nonadjacent, 1 - 2
  A <- matrix(0L, 4, 4)
  A[1, 2] <- A[2, 1] <- 1L
  A[1, 3] <- A[3, 1] <- 1L
  A[1, 4] <- A[4, 1] <- 1L
  A[3, 2] <- 1L
  A[4, 2] <- 1L
  out <- mc_orient_meek(A)
  expect_true(is_dir(out, 1, 2))
  # withheld when the unshielded triple at its center is ambiguous
  held <- mc_orient_meek(A, ambiguous = rbind(c(3L, 1L, 4L)))
  expect_true(is_undir(held, 1, 2))
  expect_gte(attr(held, "meek_blocked"), 1L)
})

test_that("the rules reach a fixed point and report the sweeps", {
  A <- matrix(0L, 4, 4)
  A[1, 2] <- 1L
  A[2, 3] <- A[3, 2] <- 1L
  A[3, 4] <- A[4, 3] <- 1L
  out <- mc_orient_meek(A)
  expect_true(is_dir(out, 2, 3))          # R1 on 1 -> 2 - 3
  expect_true(is_dir(out, 3, 4))          # R1 again on 2 -> 3 - 4
  expect_gte(attr(out, "meek_iterations"), 2L)
  expect_identical(mc_orient_meek(out), out[, , drop = FALSE],
                   ignore_attr = TRUE)    # idempotent at the fixed point
})

test_that("mc_orient_meek enforces its contracts", {
  A <- diag(0L, 3)
  expect_error(mc_orient_meek(matrix(2L, 2, 2)), "0/1")
  expect_error(mc_orient_meek(matrix(1L, 2, 3)), "square")
  expect_error(mc_orient_meek(A, ambiguous = rbind(c(1L, 2L))), "three columns")
  expect_error(mc_orient_meek(A, ambiguous = rbind(c(1L, 2L, 9L))),
               "node indices in 1:3")
})

## ---- the three rules for the triples ---------------------------------------

chain_sigma <- function(a = 0.4) {
  S <- matrix(c(1, a, a * a, a, 1, a, a * a, a, 1), 3, 3)
  dimnames(S) <- list(c("z1", "z2", "z3"), c("z1", "z2", "z3"))
  S
}

two_sepset_sigma <- function() {
  S <- matrix(c(1, 0.3, 0.6, 0.5,
                0.3, 1, 0.5, 0.6,
                0.6, 0.5, 1, 0.2,
                0.5, 0.6, 0.2, 1), 4, 4)
  dimnames(S) <- list(c("i", "j", "l", "m"), c("i", "j", "l", "m"))
  S
}

decide <- function(S, n, rule) {
  sk <- mc_pc_skeleton(S, n)
  sep <- if (rule == "standard") NULL else mc_separating_sets(S, n, sk$adjacency)
  mcglmgraph:::.orient_vstructures(sk, rule = rule, sep = sep)
}

test_that("the order-zero mechanism fools the standard rule and not the others", {
  # Chain 1 -> 2 -> 3. The marginal correlation of the extremes is a^2, the
  # order-zero test accepts it, the recorded separating set is empty, and the
  # center is therefore declared a collider that does not exist.
  S <- chain_sigma(0.4); n <- 100
  std <- decide(S, n, "standard")$vstructures
  expect_equal(nrow(std), 1L)
  expect_true(std$is_collider)          # the false collider
  expect_true(std$applied)

  for (rule in c("conservative", "majority")) {
    r <- decide(S, n, rule)$vstructures
    expect_false(r$is_collider)
    expect_true(r$ambiguous)
    expect_equal(r$n_sets, 2L)          # the empty set and the mediator
    expect_equal(r$prop_center, 0.5)
  }
})

test_that("a true collider is declared by all three rules", {
  S <- matrix(c(1, 0, 0.7, 0, 1, 0.7, 0.7, 0.7, 1.98), 3, 3)
  dimnames(S) <- list(c("a", "b", "c"), c("a", "b", "c"))
  n <- 600
  for (rule in c("standard", "conservative", "majority")) {
    r <- decide(S, n, rule)
    v <- r$vstructures
    expect_true(v$is_collider, info = rule)
    expect_true(v$applied, info = rule)
    expect_true(is_dir(r$amat, 1, 3) && is_dir(r$amat, 2, 3), info = rule)
  }
})

test_that("the conservative and majority rules do not depend on the labeling", {
  # i and j are separated by {l} and by {m} alike. Which one the skeleton
  # records is decided by the order of enumeration, so the standard rule
  # declares a collider at whichever of l, m comes second, and relabeling the
  # responses moves it. The other two rules see both sets and orient neither.
  S <- two_sepset_sigma(); n <- 200
  perm <- c(1L, 2L, 4L, 3L)                   # swap l and m

  canon <- function(rule, P) {
    Sp <- S[P, P]
    r <- decide(Sp, n, rule)
    am <- mc_orient_meek(r$amat,
                         ambiguous = if (nrow(r$ambiguous)) r$ambiguous else NULL)
    dimnames(am) <- list(colnames(Sp), colnames(Sp))
    nm <- sort(colnames(Sp))
    am[nm, nm]
  }

  std_a <- canon("standard", seq_len(4L))
  std_b <- canon("standard", perm)
  expect_false(identical(unname(std_a), unname(std_b)))   # order dependent

  for (rule in c("conservative", "majority")) {
    a <- canon(rule, seq_len(4L))
    b <- canon(rule, perm)
    expect_identical(unname(a), unname(b), info = rule)
    # and nothing is oriented: both triples are ambiguous
    expect_true(all((a == 1L) == (t(a) == 1L)), info = rule)
  }
})

test_that("an ambiguous triple withholds the propagation it would settle", {
  S <- chain_sigma(0.4); n <- 100
  r <- decide(S, n, "conservative")
  expect_equal(nrow(r$ambiguous), 1L)
  am <- mc_orient_meek(r$amat, ambiguous = r$ambiguous)
  expect_true(all((am == 1L) == (t(am) == 1L)))   # fully undirected output
})

test_that("the rule needs the collection when it is not the standard one", {
  S <- chain_sigma(); n <- 100
  sk <- mc_pc_skeleton(S, n)
  expect_error(mcglmgraph:::.orient_vstructures(sk, rule = "conservative"),
               "needs the separating sets")
  expect_error(mcglmgraph:::.orient_vstructures(sk, rule = "majority"),
               "needs the separating sets")
})

## ---- the pipeline ----------------------------------------------------------

make_chain_fit <- function(n = 400, seed = 21) {
  set.seed(seed)
  x <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- 0.4 * z1 + sqrt(1 - 0.16) * rnorm(n)
  z3 <- 0.4 * z2 + sqrt(1 - 0.16) * rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

test_that("the pipeline carries the rule and the ambiguity through", {
  f <- make_chain_fit()
  std <- mc_graph_pipeline(f)
  expect_identical(std$orient_rule, "standard")
  expect_equal(std$n_ambiguous, 0L)
  expect_null(std$sepsets)

  cons <- mc_graph_pipeline(f, orient_rule = "conservative")
  expect_identical(cons$orient_rule, "conservative")
  expect_s3_class(cons$sepsets, "mcg_sepsets")
  expect_equal(nrow(cons$vstructures), nrow(std$vstructures))
  expect_true(all(c("ambiguous", "n_sets", "prop_center") %in%
                    names(cons$vstructures)))
  # the skeleton is untouched by the rule: only the orientation can differ
  expect_identical(cons$pc_skeleton, std$pc_skeleton)
  expect_identical(cons$sepset, std$sepset)
  # and the validity gates keep their seven names
  expect_equal(length(std$gates), 7L)
  expect_identical(names(cons$gates), names(std$gates))
})

test_that("the default pipeline is unchanged by the new machinery", {
  f <- make_chain_fit(seed = 33)
  a <- mc_graph_pipeline(f)
  b <- mc_graph_pipeline(f, orient_rule = "standard")
  expect_identical(a$cpdag, b$cpdag)
  expect_identical(a$vstructures$is_collider, b$vstructures$is_collider)
  expect_identical(a$valid_cpdag, b$valid_cpdag)
})

test_that("mc_find_vstructures exposes the rule with the same columns", {
  f <- make_chain_fit()
  v <- mc_find_vstructures(f, orient_rule = "conservative")
  expect_true(all(c("i", "k", "j", "is_collider", "applied", "ambiguous",
                    "n_sets", "prop_center", "pcor", "p_value") %in% names(v)))
  vs <- mc_find_vstructures(f)
  expect_true(all(is.na(vs$n_sets)))          # no collection under "standard"
  expect_false(any(vs$ambiguous))
})
