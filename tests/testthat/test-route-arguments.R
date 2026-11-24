# Arguments of the two-stage route outside mc_vcov_graph(): every function
# that forms the two-stage variance takes `method`, refuses it where it has no
# meaning, and reports the route taken; `sensitivity` is refused or warned
# about where it has no effect.

test_that("mc_wald_graph and mc_test_adequacy refuse a route for the other covariances", {
  skip_if_not_installed("MASS")
  f <- route_fit("gaussian")
  for (tp in c("model", "sandwich")) for (m in c("shortcut", "full"))
    expect_error(mc_wald_graph(f, "eta_12", type = tp, method = m),
                 "has no meaning for type")
  expect_silent(mc_wald_graph(f, "eta_12", type = "sandwich", method = "auto"))
  g <- route_fit("gaussian",
                 between_graph = mc_between_graph("ggm", rbind(c(1, 2))))
  expect_error(mc_test_adequacy(g, type = "sandwich", method = "full"),
               "has no meaning for type")
  # with the two-stage covariance the route reaches the statistic
  ws <- mc_wald_graph(f, "eta_12", type = "two_stage", method = "shortcut")
  wf <- mc_wald_graph(f, "eta_12", type = "two_stage", method = "full")
  Vs <- mc_vcov_graph(f, "two_stage", method = "shortcut")
  Vf <- mc_vcov_graph(f, "two_stage", method = "full")
  est <- f$graph$par[["eta_12"]]
  expect_equal(unname(ws$statistic), est^2 / Vs["eta_12", "eta_12"])
  expect_equal(unname(wf$statistic), est^2 / Vf["eta_12", "eta_12"])
})

test_that("summary takes the route of the two-stage correction and reports it", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- route_fit("tweedie")
  for (m in c("auto", "shortcut", "full")) {
    s <- summary(f, variability = "two_stage", method = m)
    V <- mc_vcov_graph(f, "two_stage", method = m)
    expect_identical(s$two_stage_method, attr(V, "two_stage_method"))
    expect_identical(s$coefficients$Std.Error,
                     unname(sqrt(diag(as.matrix(V)))))
  }
  # a mean-dependent variance: auto is the full form, and the shortcut differs
  expect_identical(summary(f, variability = "two_stage")$two_stage_method,
                   "full")
  expect_false(isTRUE(all.equal(
    summary(f, variability = "two_stage", method = "shortcut")$coefficients,
    summary(f, variability = "two_stage", method = "full")$coefficients)))
  expect_identical(summary(f)$two_stage_method, NA_character_)
  expect_error(summary(f, method = "full"), "has no meaning for variability")
  expect_error(summary(f, variability = "empirical", method = "shortcut"),
               "has no meaning for variability")
  expect_output(print(summary(f, variability = "two_stage")),
                "two_stage variability, full form")
  expect_output(print(summary(f)), "(model variability)", fixed = TRUE)
})

test_that("mc_dag_coefficients takes the route of the two-stage correction and reports it", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- route_fit("tweedie", between_graph = mc_between_graph("dag"))
  E <- nrow(mc_dag_coefficients(f)$edges)
  for (m in c("auto", "shortcut", "full")) {
    cf <- mc_dag_coefficients(f, variability = "two_stage", method = m)
    V <- mc_vcov_graph(f, "two_stage", method = m)
    expect_identical(cf$two_stage_method, attr(V, "two_stage_method"))
    expect_identical(cf$edges$se,
                     unname(sqrt(diag(as.matrix(V))))[seq_len(E)])
  }
  expect_identical(
    mc_dag_coefficients(f, variability = "two_stage")$two_stage_method, "full")
  expect_identical(mc_dag_coefficients(f)$two_stage_method, NA_character_)
  expect_error(mc_dag_coefficients(f, method = "full"),
               "has no meaning for variability")
  expect_error(mc_dag_coefficients(f, variability = "empirical",
                                   method = "shortcut"),
               "has no meaning for variability")
  expect_output(print(mc_dag_coefficients(f, variability = "two_stage",
                                          method = "shortcut")),
                "two_stage variability, shortcut form")
})

test_that("sensitivity is refused or warned about where it has no effect", {
  skip_if_not_installed("MASS")
  f <- route_fit("gaussian")
  expect_error(mc_vcov_graph(f, "model", sensitivity = "expected"),
               "has no meaning for type = \"model\"")
  expect_error(mc_vcov_graph(f, "sandwich", sensitivity = "expected"),
               "has no meaning for type = \"sandwich\"")
  expect_error(mc_vcov_graph(f, "two_stage", method = "shortcut",
                             sensitivity = "expected"),
               "has no meaning for method = \"shortcut\"")
  # constant variance and an undirected layer: auto is the shortcut
  expect_warning(V <- mc_vcov_graph(f, "two_stage", sensitivity = "expected"),
                 "took the shortcut")
  expect_identical(attr(V, "two_stage_method"), "shortcut")
  # the default value is accepted everywhere, and the full form uses the other
  expect_silent(mc_vcov_graph(f, "model", sensitivity = "observed"))
  expect_silent(mc_vcov_graph(f, "two_stage", method = "shortcut",
                              sensitivity = "observed"))
  Ve <- expect_silent(mc_vcov_graph(f, "two_stage", method = "full",
                                    sensitivity = "expected"))
  expect_identical(attr(Ve, "sensitivity"), "expected")
})

test_that("the skeleton skips a conditioning set that gives no valid test", {
  # a chain 1 - 2 - 3 with the fourth response apart
  S <- diag(4)
  S[1, 2] <- S[2, 1] <- S[2, 3] <- S[3, 2] <- 0.9
  S[1, 3] <- S[3, 1] <- 0.81
  # n - |S| - 3 is positive at every order with n = 6: nothing is skipped
  sk <- expect_silent(mc_pc_skeleton(S, 6, alpha = 0.5))
  expect_identical(sk$diagnostics$n_unusable, 0L)
  # with n = 4 the factor is 1 at order zero and not positive from order one
  expect_warning(sk4 <- mc_pc_skeleton(S, 4, alpha = 0.5),
                 "could not be formed", class = "mcg_unusable_tests")
  d <- sk4$diagnostics
  expect_gt(d$n_unusable, 0L)
  expect_lt(d$n_unusable, d$n_tests)
  expect_identical(d$min_df, 1)
  # a skipped set removes no edge: the result is the one of order zero alone
  sk0 <- expect_silent(mc_pc_skeleton(S, 4, alpha = 0.5, max_ord = 0))
  expect_identical(sk4$adjacency, sk0$adjacency)
  expect_identical(sk4$sepset, sk0$sepset)
  expect_identical(sk0$diagnostics$n_unusable, 0L)
  expect_identical(d$n_tests - d$n_unusable, sk0$diagnostics$n_tests)
  A <- matrix(0L, 4, 4); A[1:3, 1:3] <- 1L; diag(A) <- 0L
  expect_identical(sk4$adjacency, A)
  # the enumeration of separating sets skips and counts the same kind of set
  ss <- mc_separating_sets(S, 4, sk4$adjacency, alpha = 0.5)
  expect_gt(ss$diagnostics$n_unusable, 0L)
  expect_true(all(ss$sets$size == 0L))
  expect_true(all(ss$pairs$truncated))
  # every test invalid: nothing is removed
  expect_warning(sk3 <- mc_pc_skeleton(S, 3), "could not be formed")
  expect_identical(sum(sk3$adjacency), 12L)
  expect_identical(sk3$diagnostics$n_unusable, sk3$diagnostics$n_tests)
  expect_identical(sk3$diagnostics$min_df, NA_real_)
})

test_that("a block with no valid partial correlation is skipped as well", {
  # |r| > 1: the block of the first two responses is not a covariance
  S <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), 3)
  ci <- getFromNamespace(".ci_pcor", "mcglmgraph")
  expect_error(ci(S, 1, 2, integer(0), 100), "valid partial correlation")
  expect_null(ci(S, 1, 2, integer(0), 100, strict = FALSE))
  expect_error(ci(diag(3), 1, 2, 3L, 4), "Insufficient Fisher")
  expect_null(ci(diag(3), 1, 2, 3L, 4, strict = FALSE))
  expect_warning(sk <- mc_pc_skeleton(S, 100), class = "mcg_unusable_tests")
  expect_gt(sk$diagnostics$n_unusable, 0L)
  expect_identical(sk$adjacency[1, 2], 1L)
  expect_null(sk$sepset[[2L]])
})

test_that("an order with no test is counted without changing the result", {
  set.seed(8); R <- 6
  A <- matrix(rnorm(R * R), R); S <- crossprod(A) / R + diag(R) * 0.3
  # df_residual = 4: orders 0 and 1 give a test, orders 2 to 4 give none
  expect_warning(sk <- mc_pc_skeleton(S, 50, alpha = 0.95, df_residual = 4),
                 class = "mcg_unusable_tests")
  ref <- expect_silent(mc_pc_skeleton(S, 50, alpha = 0.95, df_residual = 4,
                                      max_ord = 1))
  expect_identical(sk$adjacency, ref$adjacency)
  expect_identical(sk$sepset, ref$sepset)
  d <- sk$diagnostics
  expect_type(d$n_tests, "integer")
  expect_identical(d$n_tests - d$n_unusable, ref$diagnostics$n_tests)
  # the skipped orders hold every subset of the frozen neighborhoods
  adj <- ref$adjacency
  expect_identical(d$n_unusable, as.integer(sum(vapply(seq_len(R), function(i) {
    nb <- sum(adj[i, ]) - 1L
    sum(adj[i, ]) * sum(choose(nb, 2:4))
  }, numeric(1)))))
})

test_that("the pipeline does not certify an output whose skeleton skipped a test", {
  skip_if_not_installed("MASS")
  f <- route_fit("gaussian")
  ok <- expect_silent(mc_graph_pipeline(f))
  expect_identical(length(ok$gates), 8L)
  expect_true(ok$gates[["complete_separation"]])
  expect_identical(ok$separation_diagnostics$n_unusable, 0L)
  w <- 0L
  cp <- withCallingHandlers(
    mc_graph_pipeline(f, alpha = 0.9, df_residual = 3),
    mcg_unusable_tests = function(c) { w <<- w + 1L
                                       invokeRestart("muffleWarning") })
  expect_identical(w, 1L)
  expect_gt(cp$separation_diagnostics$n_unusable, 0L)
  expect_false(cp$gates[["complete_separation"]])
  expect_false(cp$valid_cpdag)
  expect_identical(cp$label, "second-order partially directed structure")
  expect_output(print(cp), "Separation tests skipped as not valid")
  expect_output(print(cp), "complete_separation")
  # a count beyond the integer range is a double and still prints
  cp$separation_diagnostics$n_unusable <- 3e10
  expect_output(print(cp), "skipped as not valid: 30000000000")
  expect_warning(vs <- mc_find_vstructures(f, alpha = 0.9, df_residual = 3),
                 class = "mcg_unusable_tests")
  expect_s3_class(vs, "data.frame")
})

test_that("the stability of the structure reports skipped tests once", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  f <- route_fit("gaussian", n = 120)
  set.seed(5)
  ok <- expect_silent(mc_cpdag_stability(f, B = 4))
  expect_identical(ok$B_incomplete, 0L)
  w <- character(0)
  set.seed(5)
  st <- withCallingHandlers(
    mc_cpdag_stability(f, B = 4, alpha = 0.9, df_residual = 3),
    warning = function(c) { w <<- c(w, conditionMessage(c))
                            invokeRestart("muffleWarning") })
  expect_length(w, 1L)
  expect_match(w, "in 4 of 4 valid resamples")
  expect_identical(st$B_incomplete, 4L)
  expect_identical(st$B_valid, 4L)
  expect_identical(st$valid_prob, 0)
  expect_output(print(st), "skipped separation tests: 4")
})
