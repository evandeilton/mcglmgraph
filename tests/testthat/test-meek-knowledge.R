# Meek's four rules and background knowledge.
#
# The rules are checked against their definition rather than against a second
# implementation. Given a DAG D, the maximally oriented partially directed graph
# under background knowledge K is the union of the DAGs Markov equivalent to D
# that satisfy K: an edge is directed exactly when every such DAG agrees on its
# direction. That union is computed here by enumerating all orientations of the
# skeleton, in plain R, and is what mc_orient_meek() must return.

## ---- the definition, by enumeration ----------------------------------------

orc_has_cycle <- function(A) {
  R <- nrow(A); indeg <- colSums(A); live <- rep(TRUE, R)
  repeat {
    src <- which(live & indeg == 0)
    if (!length(src)) break
    s <- src[1]; live[s] <- FALSE
    indeg[A[s, ] == 1] <- indeg[A[s, ] == 1] - 1
    indeg[s] <- -1
  }
  any(live)
}

## The immoralities of a DAG: i -> k <- j with i and j nonadjacent.
orc_vstructures <- function(A) {
  R <- nrow(A); adj <- (A + t(A)) > 0; out <- character(0)
  for (k in seq_len(R)) {
    pa <- which(A[, k] == 1)
    if (length(pa) < 2) next
    cb <- utils::combn(pa, 2)
    for (t in seq_len(ncol(cb)))
      if (!adj[cb[1, t], cb[2, t]])
        out <- c(out, paste0(paste(sort(cb[, t]), collapse = ":"), "|", k))
  }
  sort(unique(out))
}

orc_skeleton <- function(A) {
  S <- (A + t(A)) > 0
  which(S & upper.tri(S), arr.ind = TRUE)
}

orc_all_dags <- function(ed, R) {
  m <- nrow(ed); out <- list()
  for (mask in seq_len(2^m) - 1L) {
    A <- matrix(0L, R, R)
    for (e in seq_len(m)) {
      i <- ed[e, 1]; j <- ed[e, 2]
      if (bitwAnd(mask, bitwShiftL(1L, e - 1L))) A[i, j] <- 1L else A[j, i] <- 1L
    }
    if (!orc_has_cycle(A)) out[[length(out) + 1L]] <- A
  }
  out
}

## Union of the DAGs equivalent to D that satisfy the knowledge; NULL when the
## knowledge admits none of them.
orc_mpdag <- function(D, required = NULL, forbidden = NULL) {
  R <- nrow(D); ed <- orc_skeleton(D); vs <- orc_vstructures(D)
  keep <- Filter(function(A) {
    if (!identical(orc_vstructures(A), vs)) return(FALSE)
    if (!is.null(required) && nrow(required))
      for (k in seq_len(nrow(required)))
        if (A[required[k, 1], required[k, 2]] != 1L) return(FALSE)
    if (!is.null(forbidden) && nrow(forbidden))
      for (k in seq_len(nrow(forbidden)))
        if (A[forbidden[k, 1], forbidden[k, 2]] == 1L) return(FALSE)
    TRUE
  }, orc_all_dags(ed, R))
  if (!length(keep)) return(NULL)
  S <- Reduce(`+`, keep); n <- length(keep)
  M <- matrix(0L, R, R)
  for (e in seq_len(nrow(ed))) {
    i <- ed[e, 1]; j <- ed[e, 2]
    if (S[i, j] == n) M[i, j] <- 1L
    else if (S[j, i] == n) M[j, i] <- 1L
    else { M[i, j] <- M[j, i] <- 1L }
  }
  M
}

orc_random_dag <- function(R, p) {
  A <- matrix(0L, R, R)
  for (i in seq_len(R - 1L)) for (j in (i + 1L):R)
    if (stats::runif(1) < p) A[i, j] <- 1L
  ord <- sample(R)
  A[ord, ord]
}

bare_amat <- function(A) {
  A <- as.matrix(A)
  attributes(A) <- list(dim = dim(A))
  A
}

## ---- the four rules, one canonical instance each ---------------------------

is_dir <- function(A, a, b) A[a, b] == 1L && A[b, a] == 0L
is_undir <- function(A, a, b) A[a, b] == 1L && A[b, a] == 1L

test_that("R4 orients the edge that would otherwise create a collider at its tail", {
  # 1 - 2, 1 - 3, 1 - 4, 3 -> 4, 4 -> 2, with 3 and 2 nonadjacent.
  # Taking 2 -> 1 forces 4 -> 1 by R2, then 3 -> 1 by R2, and 3 -> 1 <- 2 is a
  # collider the graph did not have. So 1 -> 2 is compelled.
  A <- matrix(0L, 4, 4)
  A[1, 2] <- A[2, 1] <- 1L
  A[1, 3] <- A[3, 1] <- 1L
  A[1, 4] <- A[4, 1] <- 1L
  A[3, 4] <- 1L
  A[4, 2] <- 1L
  out <- mc_orient_meek(A)
  expect_true(is_dir(out, 1, 2))
  # R1 to R3 alone cannot do it: none of them applies to the edge 1 - 2 here
  expect_false(any(vapply(setdiff(1:4, c(1, 2)), function(c)
    is_dir(A, c, 1) && A[c, 2] == 0L && A[2, c] == 0L, logical(1))))  # no R1
  expect_false(any(vapply(setdiff(1:4, c(1, 2)), function(c)
    is_dir(A, 1, c) && is_dir(A, c, 2), logical(1))))                 # no R2
  # and it is withheld when the triple carrying that collider is ambiguous
  held <- mc_orient_meek(A, ambiguous = rbind(c(3L, 1L, 2L)))
  expect_true(is_undir(held, 1, 2))
  expect_gte(attr(held, "meek_blocked"), 1L)
})

test_that("R4 is idle when the only oriented edges come from v-structures", {
  # The CPDAG of a class is a fixed point of the four rules, so adding R4
  # cannot change any result the package produced before it existed.
  set.seed(101)
  for (i in seq_len(40)) {
    R <- sample(4:5, 1)
    D <- orc_random_dag(R, stats::runif(1, 0.3, 0.7))
    if (!nrow(orc_skeleton(D))) next
    cp <- orc_mpdag(D)
    expect_identical(bare_amat(mc_orient_meek(cp)), cp)
  }
})

## ---- background knowledge against the definition ---------------------------

test_that("required arcs give the union of the equivalent DAGs that contain them", {
  set.seed(202)
  tested <- 0L
  for (i in seq_len(60)) {
    R <- sample(4:5, 1)
    D <- orc_random_dag(R, stats::runif(1, 0.3, 0.7))
    cp <- orc_mpdag(D)
    und <- which(cp == 1L & t(cp) == 1L & upper.tri(cp), arr.ind = TRUE)
    if (!nrow(und)) next
    k <- sample(seq_len(min(3L, nrow(und))), 1)
    pick <- und[sample(nrow(und), k), , drop = FALSE]
    req <- matrix(as.integer(t(apply(pick, 1, function(e)
      if (D[e[1], e[2]] == 1L) e else rev(e)))), ncol = 2)
    want <- orc_mpdag(D, required = req)
    expect_false(is.null(want))
    got <- mc_orient_meek(cp, required = req)
    expect_identical(bare_amat(got), want)
    expect_true(attr(got, "knowledge")$acyclic)
    expect_true(attr(got, "knowledge")$required_kept)
    tested <- tested + 1L
  }
  expect_gt(tested, 30L)
})

test_that("forbidding the reverse of an arc is the same as requiring it", {
  set.seed(303)
  tested <- 0L
  for (i in seq_len(60)) {
    R <- sample(4:5, 1)
    D <- orc_random_dag(R, stats::runif(1, 0.3, 0.7))
    cp <- orc_mpdag(D)
    und <- which(cp == 1L & t(cp) == 1L & upper.tri(cp), arr.ind = TRUE)
    if (!nrow(und)) next
    e <- und[sample(nrow(und), 1L), ]
    req <- matrix(as.integer(if (D[e[1], e[2]] == 1L) e else rev(e)), ncol = 2)
    fb <- cbind(req[, 2L], req[, 1L])
    want <- orc_mpdag(D, forbidden = fb)
    got <- mc_orient_meek(cp, forbidden = fb)
    expect_identical(bare_amat(got), want)
    expect_identical(bare_amat(got), bare_amat(mc_orient_meek(cp, required = req)))
    tested <- tested + 1L
  }
  expect_gt(tested, 30L)
})

test_that("knowledge that contradicts a compelled edge is refused", {
  set.seed(404)
  tested <- 0L
  for (i in seq_len(60)) {
    R <- sample(4:5, 1)
    D <- orc_random_dag(R, stats::runif(1, 0.35, 0.75))
    cp <- orc_mpdag(D)
    comp <- which(cp == 1L & t(cp) == 0L, arr.ind = TRUE)
    if (!nrow(comp)) next
    e <- comp[sample(nrow(comp), 1L), ]
    expect_error(
      mc_orient_meek(cp, required = rbind(c(as.integer(e[2]), as.integer(e[1])))),
      "contradicts an edge already oriented")
    tested <- tested + 1L
  }
  expect_gt(tested, 20L)
})

## ---- contracts -------------------------------------------------------------

test_that("the knowledge tables are validated before anything is oriented", {
  A <- matrix(0L, 4, 4)
  A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1L
  expect_error(mc_orient_meek(A, required = rbind(c(1L, 2L), c(2L, 1L))),
               "required in both directions")
  expect_error(mc_orient_meek(A, required = rbind(c(1L, 2L)),
                              forbidden = rbind(c(1L, 2L))),
               "both required and forbidden")
  expect_error(mc_orient_meek(A, required = rbind(c(1L, 3L))),
               "the skeleton does not")
  expect_error(mc_orient_meek(A, required = rbind(c(1L, 9L))),
               "node indices in 1:4")
  expect_error(mc_orient_meek(A, required = rbind(c(2L, 2L))), "self-loops")
  expect_error(mc_orient_meek(A, required = cbind(1L)), "two columns")
  # a forbidden arc between nonadjacent responses is vacuous, not an error
  expect_silent(mc_orient_meek(A, forbidden = rbind(c(1L, 3L))))
  # a forbidden arc already oriented that way is a contradiction
  B <- A; B[2, 1] <- 0L
  expect_error(mc_orient_meek(B, forbidden = rbind(c(1L, 2L))),
               "already oriented that way")
})

test_that("the knowledge attribute records what was imposed", {
  A <- matrix(0L, 3, 3)
  A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1L
  out <- mc_orient_meek(A, required = rbind(c(1L, 2L)))
  kn <- attr(out, "knowledge")
  expect_equal(kn$n_imposed, 1L)
  expect_identical(unname(kn$required), matrix(1:2, ncol = 2))
  expect_true(kn$acyclic)
  expect_true(kn$required_kept)
  # 1 -> 2 with 2 - 3 and 1, 3 nonadjacent fires R1
  expect_true(is_dir(out, 2, 3))
  # without knowledge nothing is oriented
  plain <- mc_orient_meek(A)
  expect_true(is_undir(plain, 1, 2) && is_undir(plain, 2, 3))
  expect_equal(attr(plain, "knowledge")$n_imposed, 0L)
})

test_that("duplicated knowledge rows are harmless", {
  A <- matrix(0L, 3, 3)
  A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1L
  a <- mc_orient_meek(A, required = rbind(c(1L, 2L)))
  b <- mc_orient_meek(A, required = rbind(c(1L, 2L), c(1L, 2L)))
  expect_identical(bare_amat(a), bare_amat(b))
  expect_equal(attr(b, "knowledge")$n_imposed, 1L)
})

## ---- the pipeline ----------------------------------------------------------

make_chain_fit <- function(n = 600, seed = 3) {
  set.seed(seed)
  x <- rnorm(n)
  z1 <- rnorm(n)
  z2 <- 0.7 * z1 + rnorm(n)
  z3 <- 0.7 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

test_that("background knowledge propagates through the pipeline", {
  skip_if_not_installed("MASS")
  f <- make_chain_fit()
  plain <- mc_graph_pipeline(f)
  expect_identical(plain$label, "CPDAG")
  expect_equal(plain$n_imposed, 0L)
  # the field is always present and says that nothing was imposed
  expect_equal(plain$knowledge$n_imposed, 0L)
  expect_equal(nrow(plain$knowledge$required), 0L)
  expect_equal(nrow(plain$knowledge$forbidden), 0L)
  # the chain y1 - y2 - y3 has no collider, so nothing is oriented
  expect_equal(sum(plain$cpdag == 1L & t(plain$cpdag) == 0L), 0L)

  kn <- mc_graph_pipeline(f, required = rbind(c("y1", "y2")))
  expect_equal(kn$n_imposed, 1L)
  expect_identical(kn$label,
                   "maximally oriented graph under background knowledge")
  expect_true(kn$valid_cpdag)
  # y1 -> y2 with y1 and y3 nonadjacent fires R1, so y2 -> y3 follows
  expect_equal(unname(kn$cpdag["y1", "y2"]), 1L)
  expect_equal(unname(kn$cpdag["y2", "y1"]), 0L)
  expect_equal(unname(kn$cpdag["y2", "y3"]), 1L)
  expect_equal(unname(kn$cpdag["y3", "y2"]), 0L)
  # the skeleton and the separating sets are untouched by the knowledge
  expect_identical(kn$pc_skeleton, plain$pc_skeleton)
  expect_identical(kn$sepset, plain$sepset)
  expect_identical(names(kn$gates), names(plain$gates))
})

test_that("the knowledge can be given by index as well as by name", {
  skip_if_not_installed("MASS")
  f <- make_chain_fit()
  by_name <- mc_graph_pipeline(f, required = rbind(c("y1", "y2")))
  by_index <- mc_graph_pipeline(f, required = rbind(c(1L, 2L)))
  expect_identical(by_name$cpdag, by_index$cpdag)
  # forbidding the reverse is the same request
  by_forbid <- mc_graph_pipeline(f, forbidden = rbind(c("y2", "y1")))
  expect_identical(by_name$cpdag, by_forbid$cpdag)
})

test_that("the pipeline refuses knowledge it cannot honor", {
  skip_if_not_installed("MASS")
  f <- make_chain_fit()
  expect_error(mc_graph_pipeline(f, required = rbind(c("y1", "y3"))),
               "the skeleton does not")
  expect_error(mc_graph_pipeline(f, required = rbind(c("y1", "y9"))),
               "a response the fit does not have")
  expect_error(mc_graph_pipeline(f, orient_method = "hc",
                                 required = rbind(c("y1", "y2"))),
               "native orientation only")
  expect_error(mc_graph_pipeline(f, required = rbind(c("y1", "y2"), c("y2", "y1"))),
               "required in both directions")
})

test_that("the adjacency of the result is a plain labeled matrix", {
  skip_if_not_installed("MASS")
  f <- make_chain_fit()
  out <- mc_graph_pipeline(f, required = rbind(c("y1", "y2")))
  expect_identical(names(attributes(out$cpdag)), c("dim", "dimnames"))
  expect_identical(dimnames(out$cpdag),
                   list(f$graph$response_names, f$graph$response_names))
  txt <- capture.output(print(out))
  expect_true(any(grepl("Background knowledge: 1 arc imposed", txt)))
  expect_true(any(grepl("maximally oriented graph", txt)))
  expect_false(any(grepl("meek_conflicts", txt)))
})
