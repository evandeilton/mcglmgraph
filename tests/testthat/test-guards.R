# Guards on inputs and on results that are not what they would be read as.

guard_data <- function(n = 300, seed = 5) {
  set.seed(seed)
  x <- rnorm(n)
  y1 <- rnorm(n)
  y2 <- 0.6 * y1 + rnorm(n)
  y3 <- 0.5 * y1 + 0.5 * y2 + rnorm(n)
  data.frame(y1, y2, y3, x)
}

test_that("a graph layer that did not converge is returned with a warning", {
  set.seed(2); R <- 6; n <- 6
  Y <- matrix(rnorm(n * R), n, R, dimnames = list(NULL, paste0("y", seq_len(R))))
  fl <- lapply(seq_len(R), function(k) stats::as.formula(paste0("y", k, " ~ 1")))
  expect_warning(f <- mcglm_graph(fl, data = as.data.frame(Y)),
                 "did not converge.*6 units for 6 responses")
  expect_false(f$graph$converged)
  # a converged fit is silent
  expect_no_warning(mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
                                data = guard_data()))
})

test_that("a transformed response is accepted and labeled by its expression", {
  d <- guard_data()
  d$z1 <- exp(d$y1 / 4)
  f <- mcglm_graph(list(log(z1) ~ x, y2 ~ x, y3 ~ x), data = d)
  expect_identical(f$graph$response_names, c("log(z1)", "y2", "y3"))
  expect_identical(colnames(f$graph$Sigma_b), c("log(z1)", "y2", "y3"))
  ref <- mcglm_graph(list(lz ~ x, y2 ~ x, y3 ~ x),
                     data = transform(d, lz = log(z1)))
  expect_equal(unname(f$graph$par), unname(ref$graph$par))
})

test_that("a contrast without sampling variability is not tested", {
  d <- guard_data()
  f <- mcglm_graph(list(y1 ~ 1, y2 ~ 1, y3 ~ 1),
                   between_graph = mc_between_graph("dag"), data = d)
  # the first stage fixes the innovation variance of the root node
  expect_lt(mc_vcov_graph(f, "two_stage")["delta_1", "delta_1"], 1e-20)
  expect_error(mc_wald_graph(f, "delta_1", rhs = 1, type = "two_stage"),
               "singular on its range")
  expect_error(mc_wald_graph(f, c("delta_1", "gamma_21"), rhs = c(1, 0),
                             type = "two_stage"), "singular on its range")
  # the same parameter under a covariance that gives it variability, and
  # the other parameters under the corrected one, are tested
  expect_s3_class(mc_wald_graph(f, "delta_1", rhs = 1, type = "sandwich"),
                  "htest")
  # and the verdict does not depend on the scale of a row of the contrast
  L <- matrix(0, 2, length(f$graph$par),
              dimnames = list(NULL, names(f$graph$par)))
  L[1, "gamma_21"] <- 1; L[2, "gamma_31"] <- 1e6
  w6 <- mc_wald_graph(f, L, type = "two_stage")
  w1 <- mc_wald_graph(f, c("gamma_21", "gamma_31"), type = "two_stage")
  expect_equal(unname(w6$statistic), unname(w1$statistic), tolerance = 1e-8)
  w <- mc_wald_graph(f, "gamma_21", type = "two_stage")
  V <- mc_vcov_graph(f, "two_stage")
  expect_equal(unname(w$statistic),
               unname(f$graph$par["gamma_21"]^2 / V["gamma_21", "gamma_21"]))
})

test_that("a sparse identity is recognized at any size", {
  ns <- asNamespace("mcglmgraph")
  I <- methods::as(Matrix::Diagonal(4500), "generalMatrix")
  expect_true(ns$.is_identity_structure(I))
  expect_true(ns$.is_identity_structure(Matrix::Diagonal(4500)))
  J <- I
  J[1, 2] <- J[2, 1] <- 0.5
  expect_false(ns$.is_identity_structure(J))
  expect_false(ns$.is_identity_structure(2 * I))
  # an explicit zero stored off the diagonal does not change the answer
  Z <- Matrix::sparseMatrix(i = c(1:5, 1L), j = c(1:5, 2L), x = c(rep(1, 5), 0))
  expect_true(ns$.is_identity_structure(Z))
  expect_true(ns$.is_identity_structure(diag(5)))
  expect_false(ns$.is_identity_structure(matrix(1, 3, 3)))
})

test_that("the selector of the selection frequencies is checked on entry", {
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = guard_data())
  expect_error(mc_graph_stability(f, B = 4, criterion = "nonsense"),
               "should be one of")
})
