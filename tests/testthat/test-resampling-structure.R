# Resampling frequencies: of the selected support, and of each orientation.
#
# Both routines share the machinery that refits on a set of rows, so the checks
# below are about what each reports and about the two resampling conventions:
# half-samples for selection frequencies, bootstrap draws for structure ones.

chain_fit <- function(n = 250, seed = 1) {
  set.seed(seed)
  K <- diag(4) * 2
  K[1, 2] <- K[2, 1] <- -0.7
  K[2, 3] <- K[3, 2] <- -0.7
  K[3, 4] <- K[4, 3] <- -0.7
  Y <- MASS::mvrnorm(n, rep(0, 4), solve(K))
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
                  x = rnorm(n))
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
}

collider_fit <- function(n = 400, seed = 2) {
  set.seed(seed)
  x <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
}

## ---- selection frequencies with both stages refitted -----------------------

test_that("the half-samples can carry the first stage as well", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- chain_fit()
  set.seed(10); a <- mc_graph_stability(f, B = 25)
  set.seed(10); b <- mc_graph_stability(f, B = 25, refit = "full")
  expect_identical(a$refit, "residual")
  expect_identical(b$refit, "full")
  expect_equal(b$B_failed, 0L)
  # the chain is recovered either way: adjacent pairs beat the mediated ones
  for (st in list(a, b)) {
    expect_gte(st$selection_prob[1, 2], st$selection_prob[1, 3])
    expect_gte(st$selection_prob[3, 4], st$selection_prob[1, 4])
    expect_true(all(st$selection_prob >= 0 & st$selection_prob <= 1))
  }
  txt <- capture.output(print(b))
  expect_true(any(grepl("both stages refitted", txt)))
})

test_that("refitting is refused where resampling rows would not preserve the design", {
  skip_if_not_installed("MASS")
  f <- chain_fit(n = 120)
  g <- f
  M <- diag(120); M[1, 2] <- M[2, 1] <- 0.5
  g$matrix_pred[[1]] <- list(M)
  expect_error(mc_graph_stability(g, B = 5, refit = "full"),
               "needs cluster")
  expect_s3_class(mc_graph_stability(g, B = 5), "mcg_stability")
})

test_that("the half-samples can draw whole clusters", {
  skip_if_not_installed("MASS")
  f <- chain_fit(n = 200)
  cl <- rep(seq_len(100), each = 2L)
  set.seed(12)
  st <- mc_graph_stability(f, B = 10, cluster = cl)
  expect_true(st$clustered)
  expect_equal(st$B_valid, 10L)
  expect_true(all(st$selection_prob >= 0 & st$selection_prob <= 1))
  # without clusters the draw is one call per subsample to the generator
  set.seed(12); a <- mc_graph_stability(f, B = 10)
  expect_false(a$clustered)
  expect_error(mc_graph_stability(f, B = 5, cluster = 1:3),
               "one identifier per observational row")
})

test_that("the half-samples do not depend on the worker count", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- chain_fit(n = 150)
  set.seed(77); a <- mc_graph_stability(f, B = 12, refit = "full", workers = 1L)
  set.seed(77); b <- mc_graph_stability(f, B = 12, refit = "full", workers = 2L)
  expect_equal(a$selection_prob, b$selection_prob, tolerance = 1e-12)
})

## ---- structure frequencies -------------------------------------------------

test_that("the four states of a pair exhaust the resamples", {
  skip_if_not_installed("MASS")
  f <- collider_fit()
  set.seed(11)
  st <- mc_cpdag_stability(f, B = 40)
  expect_s3_class(st, "mcg_cpdag_stability")
  tot <- st$states$absent + st$states$undirected + st$states$i_to_j +
    st$states$j_to_i
  expect_equal(tot, rep(1, nrow(st$states)), tolerance = 1e-12)
  # the skeleton frequency is the complement of absence, and symmetric
  for (e in seq_len(nrow(st$states))) {
    i <- st$states$i[e]; j <- st$states$j[e]
    expect_equal(unname(st$edge_prob[i, j]), 1 - st$states$absent[e],
                 tolerance = 1e-12)
    expect_equal(unname(st$edge_prob[i, j]), unname(st$edge_prob[j, i]))
  }
  expect_equal(st$B_valid + st$B_failed, 40L)
})

test_that("a true collider is recovered in most resamples", {
  skip_if_not_installed("MASS")
  f <- collider_fit()
  set.seed(12)
  st <- mc_cpdag_stability(f, B = 40)
  # y1 -> y3 <- y2, and y1, y2 nonadjacent
  row13 <- st$states[st$states$response1 == "y1" & st$states$response2 == "y3", ]
  row23 <- st$states[st$states$response1 == "y2" & st$states$response2 == "y3", ]
  row12 <- st$states[st$states$response1 == "y1" & st$states$response2 == "y2", ]
  expect_gt(row13$i_to_j, 0.8)
  expect_gt(row23$i_to_j, 0.8)
  expect_gt(row12$absent, 0.8)
  # the consensus reads those very states back
  expect_equal(unname(st$consensus["y1", "y3"]), 1L)
  expect_equal(unname(st$consensus["y3", "y1"]), 0L)
  expect_equal(unname(st$consensus["y1", "y2"]), 0L)
  expect_equal(unname(st$consensus["y2", "y1"]), 0L)
})

test_that("a chain without a collider leaves its edges undirected", {
  skip_if_not_installed("MASS")
  f <- chain_fit(n = 400, seed = 8)
  set.seed(13)
  st <- mc_cpdag_stability(f, B = 30)
  adj <- st$states[with(st$states, (i == 1 & j == 2) | (i == 2 & j == 3) |
                          (i == 3 & j == 4)), ]
  expect_true(all(adj$undirected > 0.6))
  expect_true(all(adj$i_to_j + adj$j_to_i < 0.4))
})

test_that("the order-zero mechanism shows up as ambiguity under the conservative rule", {
  skip_if_not_installed("MASS")
  f <- chain_fit(n = 300, seed = 21)
  set.seed(14)
  std <- mc_cpdag_stability(f, B = 25, orient_rule = "standard")
  cons <- mc_cpdag_stability(f, B = 25, orient_rule = "conservative")
  expect_equal(std$mean_ambiguous, 0)
  expect_gt(cons$mean_ambiguous, 0)
  # and the conservative rule leaves at least as much undirected
  expect_gte(sum(cons$states$undirected), sum(std$states$undirected))
})

test_that("the structure frequencies can carry the first stage", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 300)
  set.seed(15)
  st <- mc_cpdag_stability(f, B = 20, refit = "full")
  expect_identical(st$refit, "full")
  expect_equal(st$B_failed, 0L)
  expect_true(all(st$edge_prob >= 0 & st$edge_prob <= 1))
  txt <- capture.output(print(st))
  expect_true(any(grepl("both stages refitted", txt)))
})

test_that("the resamples do not depend on the worker count", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 200)
  set.seed(88); a <- mc_cpdag_stability(f, B = 12, workers = 1L)
  set.seed(88); b <- mc_cpdag_stability(f, B = 12, workers = 2L)
  expect_equal(a$states, b$states, tolerance = 1e-12)
  expect_equal(a$consensus, b$consensus)
})

test_that("background knowledge reaches the resamples", {
  skip_if_not_installed("MASS")
  f <- chain_fit(n = 300, seed = 31)
  set.seed(16)
  st <- mc_cpdag_stability(f, B = 20, required = rbind(c("y1", "y2")))
  row <- st$states[st$states$i == 1 & st$states$j == 2, ]
  # the arc is imposed wherever the edge survives the skeleton search
  expect_equal(row$j_to_i, 0)
  expect_equal(row$undirected, 0)
  expect_gt(row$i_to_j, 0.8)
})

test_that("the contracts of the arguments are enforced", {
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 150)
  expect_error(mc_cpdag_stability(f, B = 1), "at least two")
  expect_error(mc_cpdag_stability(f, B = 5, cluster = 1:3),
               "one identifier per observational row")
  expect_error(mc_cpdag_stability(f, B = 5, cluster = rep(1L, 150)),
               "at least two independent clusters")
  expect_error(mc_cpdag_stability(f, B = 5, orient_rule = "greedy"), "arg")
  expect_error(mc_cpdag_stability(structure(list(), class = "lm")), "inherits")
})

test_that("clusters are resampled whole", {
  skip_if_not_installed("MASS")
  f <- collider_fit(n = 200)
  set.seed(17)
  st <- mc_cpdag_stability(f, B = 15, cluster = rep(1:100, each = 2L))
  expect_true(st$clustered)
  expect_equal(st$B_valid + st$B_failed, 15L)
})
