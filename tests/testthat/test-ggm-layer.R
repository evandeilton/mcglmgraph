test_that("GGM bases have the correct count and symmetry", {
  b <- mc_precision_bases(4)
  expect_equal(length(b$bases), 4 + 6)          # R diagonals + R(R-1)/2 edges
  for (Z in b$bases) expect_true(isSymmetric(Z))
})

test_that("GGM analytical derivatives match numerical differentiation", {
  b <- mc_precision_bases(4, edges = rbind(c(1, 2), c(2, 3), c(3, 4)))
  eta <- c(2, 2, 2, 2, -0.5, -0.4, -0.3)
  K <- mc_build_precision(eta, b)
  S <- mc_precision_to_covariance(K)$Sigma_b
  D <- mc_precision_derivatives(S, b)
  eps <- 1e-6; d <- 5
  ep <- eta; ep[d] <- ep[d] + eps
  em <- eta; em[d] <- em[d] - eps
  Sp <- mc_precision_to_covariance(mc_build_precision(ep, b))$Sigma_b
  Sm <- mc_precision_to_covariance(mc_build_precision(em, b))$Sigma_b
  expect_lt(max(abs(D[[d]] - (Sp - Sm) / (2 * eps))), 1e-5)
})

test_that("SPD line search keeps the precision positive definite", {
  b <- mc_precision_bases(3)
  ls <- mc_spd_line_search(c(2, 2, 2, 0, 0, 0),
                           c(-0.5, -0.5, -0.5, -0.9, -0.9, -0.9), b)
  expect_true(ls$is_spd)
})

test_that("a failed SPD line search is not reported as convergence", {
  local_mocked_bindings(
    cpp_spd_line_search = function(...) list(is_spd = FALSE),
    .package = "mcglmgraph"
  )
  ans <- mcglmgraph:::.fit_between(diag(3), n = 50,
                                    spec = mc_between_graph("ggm"))
  expect_false(ans$converged)
  expect_identical(ans$status, "spd_line_search_failed")
})
