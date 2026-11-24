# Out-of-sample comparison of candidate supports.
#
# The held-out score is rebuilt from its definition here: the Gaussian working
# criterion of the residual layer, evaluated on the held-out residuals under the
# covariance fitted on the rest, per held-out unit.

sparse_fit <- function(n = 500, seed = 1, R = 5) {
  set.seed(seed)
  K <- diag(R) * 2
  K[1, 2] <- K[2, 1] <- -0.8
  K[3, 4] <- K[4, 3] <- -0.8
  Y <- MASS::mvrnorm(n, rep(0, R), solve(K))
  d <- as.data.frame(Y)
  names(d) <- paste0("y", seq_len(R))
  d$x <- rnorm(n)
  mcglm_graph(lapply(seq_len(R), function(i)
    stats::as.formula(paste0("y", i, " ~ x"))), data = d)
}

test_that("the held-out score is the working criterion per held-out unit", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 300)
  g <- f$graph
  set.seed(21)
  cv <- mc_graph_cv(f, supports = list(unstructured = NULL), folds = 4L)

  # rebuild fold one by hand, with the same split
  set.seed(21)
  assign_fold <- sample(rep(seq_len(4L), length.out = g$n))
  E <- as.matrix(g$residuals)
  test <- which(assign_fold == 1L)
  train <- setdiff(seq_len(g$n), test)
  spec <- mc_between_graph("ggm"); spec$R <- g$R
  Str <- crossprod(scale(E[train, ], center = TRUE, scale = FALSE)) / length(train)
  b <- mcglmgraph:::.fit_between(Str, length(train), spec)
  Ste <- crossprod(scale(E[test, ], center = TRUE, scale = FALSE)) / length(test)
  want <- -0.5 * (g$R * log(2 * pi) +
                    as.numeric(determinant(b$Sigma_b, logarithm = TRUE)$modulus) +
                    sum(solve(b$Sigma_b) * Ste))
  got <- cv$per_fold$score[cv$per_fold$fold == 1L]
  expect_equal(got, want, tolerance = 1e-9)
})

test_that("the folds partition the units", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 200)
  set.seed(31)
  cv <- mc_graph_cv(f, supports = list(unstructured = NULL), folds = 5L)
  expect_equal(nrow(cv$per_fold), 5L)
  expect_setequal(cv$per_fold$fold, 1:5)
  expect_equal(nrow(cv$failures), 0L)
})

test_that("a sparse truth is predicted better by a sparse support", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 500)
  sel <- mc_graph_stepwise(f, criterion = "bic")
  set.seed(41)
  cv <- mc_graph_cv(f, supports = list(unstructured = NULL,
                                       selected = sel$edges,
                                       empty = matrix(integer(0), 0, 2)))
  best <- cv$summary$support[1L]
  expect_identical(best, "selected")
  expect_equal(cv$summary$delta[1L], 0)
  expect_true(all(cv$summary$delta >= 0))
  # the empty support, which denies the dependence, is the worst
  expect_identical(cv$summary$support[nrow(cv$summary)], "empty")
  # and the selected support uses fewer parameters than the unstructured one
  npar <- stats::setNames(cv$summary$npar, cv$summary$support)
  expect_lt(npar[["selected"]], npar[["unstructured"]])
})

test_that("the summary is an ordered aggregation of the folds", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 300)
  set.seed(51)
  cv <- mc_graph_cv(f, supports = list(unstructured = NULL,
                                       empty = matrix(integer(0), 0, 2)),
                    folds = 4L)
  for (s in cv$summary$support) {
    z <- cv$per_fold$score[cv$per_fold$support == s]
    row <- cv$summary[cv$summary$support == s, ]
    expect_equal(row$score, mean(z), tolerance = 1e-12)
    expect_equal(row$se, sd(z) / sqrt(length(z)), tolerance = 1e-12)
    expect_equal(row$n_folds, length(z))
  }
  expect_true(!is.unsorted(rev(cv$summary$score)))
  expect_equal(cv$summary$delta, cv$summary$score[1L] - cv$summary$score,
               tolerance = 1e-12)
})

test_that("the default compares the fitted support with the complete graph", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 300)
  sub <- .refit_between(f, mc_between_graph("ggm", rbind(c(1, 2), c(3, 4))))
  set.seed(61)
  cv <- mc_graph_cv(sub, folds = 4L)
  expect_setequal(cv$supports, c("unstructured", "fitted"))
  # a fit on the complete graph has nothing sparse to compare against
  set.seed(61)
  cv2 <- mc_graph_cv(f, folds = 4L)
  expect_identical(cv2$supports, "unstructured")
})

test_that("the split does not depend on the worker count", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 200)
  set.seed(71); a <- mc_graph_cv(f, folds = 4L, workers = 1L)
  set.seed(71); b <- mc_graph_cv(f, folds = 4L, workers = 2L)
  expect_equal(a$per_fold, b$per_fold, tolerance = 1e-12)
})

test_that("the DAG layer keeps its ordering across the folds", {
  skip_if_not_installed("MASS")
  set.seed(81); n <- 300; x <- rnorm(n)
  y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + rnorm(n)
  d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                   between_graph = mc_between_graph("dag", ordering = c(3L, 1L, 2L)),
                   data = d)
  set.seed(91)
  cv <- mc_graph_cv(f, supports = list(unstructured = NULL,
                                       chain = rbind(c(1, 2), c(2, 3))),
                    folds = 4L)
  expect_equal(nrow(cv$per_fold), 8L)
  expect_true(all(is.finite(cv$per_fold$score)))
})

test_that("the contracts of the arguments are enforced", {
  skip_if_not_installed("MASS")
  f <- sparse_fit(n = 120)
  expect_error(mc_graph_cv(f, folds = 1L), "at least two")
  expect_error(mc_graph_cv(f, folds = 1000L), "cannot exceed")
  expect_error(mc_graph_cv(f, supports = list(NULL)), "named list")
  expect_error(mc_graph_cv(f, supports = "complete"), "named list")
  expect_error(mc_graph_cv(structure(list(), class = "lm")), "inherits")
})
