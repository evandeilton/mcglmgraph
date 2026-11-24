# Exhaustive enumeration of separating sets (mc_separating_sets and its kernel).
#
# The kernel is checked against an independent enumeration written in plain R
# on top of .ci_pcor(), the very test the skeleton uses. That comparison is the
# point of the design: the order-independent rules must differ from the standard
# one by the rule alone, never by a change of conditional-association test.

## Independent reference enumeration.
ref_sepsets <- function(Sigma, n, adj, i, j, alpha, max_ord,
                        df_residual = NULL) {
  nbi <- setdiff(which(adj[i, ] == 1L), j)
  nbj <- setdiff(which(adj[j, ] == 1L), i)
  out <- list(); seen <- character(0)
  for (nb in list(nbi, nbj)) {
    top <- min(max_ord, length(nb))
    for (s in seq.int(0L, top)) {
      combos <- if (s == 0L) list(integer(0)) else
        lapply(utils::combn(length(nb), s, simplify = FALSE),
               function(pos) nb[pos])
      for (S in combos) {
        key <- paste(sort(S), collapse = ",")
        if (key %in% seen) next
        seen <- c(seen, key)
        tst <- mcglmgraph:::.ci_pcor(Sigma, i, j, S, n,
                                     df_residual = df_residual)
        if (tst$p_value >= alpha)
          out[[length(out) + 1L]] <- list(key = key, pcor = tst$pcor,
                                          p = tst$p_value)
      }
    }
  }
  out
}

## Chain 1 -> 2 -> 3 with equal path coefficients: Sigma_13 = a^2 and the
## order-zero test accepts at this sample size, while conditioning on 2 gives
## exactly zero.
chain_sigma <- function(a = 0.4) {
  S <- matrix(c(1, a, a * a, a, 1, a, a * a, a, 1), 3, 3)
  dimnames(S) <- list(c("z1", "z2", "z3"), c("z1", "z2", "z3"))
  S
}

## i and j separated by {l} and by {m}, because Sil*Slj = Sim*Smj.
two_sepset_sigma <- function() {
  S <- matrix(c(1, 0.3, 0.6, 0.5,
                0.3, 1, 0.5, 0.6,
                0.6, 0.5, 1, 0.2,
                0.5, 0.6, 0.2, 1), 4, 4)
  dimnames(S) <- list(c("i", "j", "l", "m"), c("i", "j", "l", "m"))
  S
}

test_that("the kernel reproduces the enumeration written on top of .ci_pcor", {
  for (cfg in list(list(S = chain_sigma(), n = 100),
                   list(S = two_sepset_sigma(), n = 200))) {
    S <- cfg$S; n <- cfg$n; R <- nrow(S)
    sk <- mc_pc_skeleton(S, n)
    sep <- mc_separating_sets(S, n, sk$adjacency)
    for (e in seq_len(nrow(sep$pairs))) {
      i <- sep$pairs$i[e]; j <- sep$pairs$j[e]
      ref <- ref_sepsets(S, n, sk$adjacency, i, j, 0.05, R - 2L)
      got <- sep$sets[sep$sets$pair == e, , drop = FALSE]
      expect_equal(nrow(got), length(ref))
      keys_got <- vapply(seq_len(nrow(got)), function(r)
        paste(which(sep$members[which(sep$sets$pair == e)[r], ] == 1L),
              collapse = ","), character(1))
      expect_setequal(keys_got, vapply(ref, `[[`, character(1), "key"))
      ord <- match(keys_got, vapply(ref, `[[`, character(1), "key"))
      expect_equal(got$pcor, vapply(ref[ord], `[[`, numeric(1), "pcor"),
                   tolerance = 1e-12)
      expect_equal(got$p_value, vapply(ref[ord], `[[`, numeric(1), "p"),
                   tolerance = 1e-12)
    }
  }
})

test_that("the chain gives the two sets the theory predicts, with their values", {
  a <- 0.4
  S <- chain_sigma(a); n <- 100
  sk <- mc_pc_skeleton(S, n)
  expect_equal(unname(sk$adjacency[1, 3]), 0L)      # marginal test accepts
  sep <- mc_separating_sets(S, n, sk$adjacency)
  expect_equal(nrow(sep$sets), 2L)
  expect_setequal(sep$sets$set, c("{}", "{z2}"))
  # the empty set carries the marginal correlation a^2 and the mediator kills it
  empty <- sep$sets[sep$sets$set == "{}", ]
  med <- sep$sets[sep$sets$set == "{z2}", ]
  expect_equal(empty$pcor, a * a, tolerance = 1e-12)
  expect_equal(med$pcor, 0, tolerance = 1e-10)
  expect_equal(med$p_value, 1, tolerance = 1e-9)
  # Fisher z of the empty set, from its definition
  z <- atanh(a * a) * sqrt(n - 0 - 3)
  expect_equal(empty$p_value, 2 * pnorm(-abs(z)), tolerance = 1e-12)
})

test_that("a true collider leaves the center out of every separating set", {
  # 1 and 2 independent, 3 = 0.7*1 + 0.7*2
  S <- matrix(c(1, 0, 0.7, 0, 1, 0.7, 0.7, 0.7, 1.98), 3, 3)
  dimnames(S) <- list(c("a", "b", "c"), c("a", "b", "c"))
  n <- 600
  sk <- mc_pc_skeleton(S, n)
  sep <- mc_separating_sets(S, n, sk$adjacency)
  expect_equal(nrow(sep$sets), 1L)
  expect_equal(sep$sets$set, "{}")
  expect_equal(unname(sep$members[1, "c"]), 0L)            # center absent: collider
})

test_that("a set reachable from both neighborhoods is reported once", {
  S <- two_sepset_sigma(); n <- 200
  sk <- mc_pc_skeleton(S, n)
  sep <- mc_separating_sets(S, n, sk$adjacency)
  keys <- apply(sep$members, 1L, paste, collapse = "")
  expect_equal(anyDuplicated(keys[sep$sets$pair == 1L]), 0L)
  # both {l} and {m} separate i and j, each exactly once
  expect_setequal(sep$sets$set, c("{l}", "{m}"))
})

test_that("the degrees-of-freedom contract of the first stage is honored", {
  S <- chain_sigma(); n <- 100; nu <- 92
  sk <- mc_pc_skeleton(S, n)
  raw <- mc_separating_sets(S, n, sk$adjacency)
  adj <- mc_separating_sets(S, n, sk$adjacency, df_residual = nu)
  expect_identical(raw$diagnostics$df_contract, "n - |S| - 3")
  expect_identical(adj$diagnostics$df_contract, "df_residual - |S| - 2")
  e_raw <- raw$sets[raw$sets$set == "{}", ]
  e_adj <- adj$sets[adj$sets$set == "{}", ]
  expect_equal(e_raw$pcor, e_adj$pcor, tolerance = 1e-14)   # same estimate
  expect_equal(e_adj$p_value,
               2 * pnorm(-abs(atanh(e_adj$pcor) * sqrt(nu - 0 - 2))),
               tolerance = 1e-12)
  expect_gt(e_adj$p_value, e_raw$p_value)                   # fewer df, larger p
})

test_that("max_ord truncates the enumeration and says so", {
  S <- two_sepset_sigma(); n <- 200
  sk <- mc_pc_skeleton(S, n)
  full <- mc_separating_sets(S, n, sk$adjacency)
  none <- mc_separating_sets(S, n, sk$adjacency, max_ord = 0L)
  expect_equal(nrow(none$sets), 0L)      # the marginal correlation is 0.3
  expect_true(all(none$pairs$truncated))
  expect_false(any(full$pairs$truncated))
  expect_lt(none$diagnostics$n_tests, full$diagnostics$n_tests)
})

test_that("only the requested pairs are enumerated", {
  S <- two_sepset_sigma(); n <- 200
  sk <- mc_pc_skeleton(S, n)
  one <- mc_separating_sets(S, n, sk$adjacency, pairs = rbind(c(1L, 2L)))
  expect_equal(nrow(one$pairs), 1L)
  expect_equal(c(one$pairs$i, one$pairs$j), c(1L, 2L))
})

test_that("the contracts of the arguments are enforced", {
  S <- chain_sigma(); n <- 100
  sk <- mc_pc_skeleton(S, n)
  A <- sk$adjacency
  expect_error(mc_separating_sets(S, n, A[-1, ]), "square")
  bad <- A; bad[1, 2] <- 1L; bad[2, 1] <- 0L
  expect_error(mc_separating_sets(S, n, bad), "symmetric")
  bad2 <- A; bad2[1, 2] <- bad2[2, 1] <- 2L
  expect_error(mc_separating_sets(S, n, bad2), "0/1")
  expect_error(mc_separating_sets(S, n, A, pairs = rbind(c(1L, 2L))),
               "non-adjacent")
  expect_error(mc_separating_sets(S, n, A, alpha = 0), "strictly between")
  expect_error(mc_separating_sets(S, -1, A), "finite positive")
  expect_error(mc_separating_sets(S, n, A, max_ord = -1), "non-negative")
  expect_error(mc_separating_sets(S, n, A, df_residual = 1), "greater than 2")
  expect_error(mc_separating_sets(S[1:2, 1:2, drop = FALSE], n, A), "square")
})

test_that("the printed summary reports the collection", {
  S <- chain_sigma(); n <- 100
  sk <- mc_pc_skeleton(S, n)
  sep <- mc_separating_sets(S, n, sk$adjacency)
  expect_s3_class(sep, "mcg_sepsets")
  out <- capture.output(print(sep))
  expect_true(any(grepl("Separating sets", out)))
  expect_true(any(grepl("sets found: 2", out)))
})
