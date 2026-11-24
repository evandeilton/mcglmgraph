# Edge scores in closed form, the single separation kernel, intervals on a
# transformed scale and the summary object.

tool_fit <- function(n = 300, seed = 7, edges = NULL, type = "ggm") {
  set.seed(seed)
  x <- rnorm(n)
  y1 <- rnorm(n)
  y2 <- 0.6 * y1 + rnorm(n)
  y3 <- 0.5 * y2 + rnorm(n)
  y4 <- 0.3 * y1 + rexp(n)
  mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x),
              between_graph = mc_between_graph(type, edges),
              data = data.frame(y1, y2, y3, y4, x))
}

# The efficient score of one absent edge from the generic assembly: the
# parametrization augmented with the edge at zero, psi and J from the traces.
generic_score <- function(f, i, j, variability) {
  ns <- asNamespace("mcglmgraph")
  g <- f$graph; R <- g$R
  cur <- g$meta$bases$edges
  spec <- mc_between_graph("ggm", rbind(cur, c(i, j)))
  bld <- ns$.make_builder(spec, R)
  aug <- attr(bld, "meta")$bases$edges
  pos <- R + which(aug[, 1L] == i & aug[, 2L] == j)
  par <- numeric(R + nrow(aug))
  par[seq_len(R)] <- g$par[seq_len(R)]
  par[R + match(paste(cur[, 1L], cur[, 2L]), paste(aug[, 1L], aug[, 2L]))] <-
    g$par[g$edge_pos]
  b <- bld(par)
  pj <- ns$cpp_pearson_between(b$Sigma, g$Shat, g$n, b$Dlist)
  lam <- setdiff(seq_along(par), pos)
  a <- numeric(length(par)); a[pos] <- 1
  a[lam] <- -solve(pj$J[lam, lam], pj$J[lam, pos])
  V <- if (variability == "model") pj$J else
    crossprod(ns$.pearson_scores(g$residuals, b$Sigma, b$Dlist))
  sum(a * pj$psi)^2 / as.numeric(t(a) %*% V %*% a)
}

test_that("the closed-form efficient score equals the generic assembly", {
  f <- tool_fit(edges = rbind(c(1, 2), c(2, 3)))
  for (v in c("model", "empirical")) {
    es <- mc_edge_scores(f, "score", variability = v)
    expect_identical(nrow(es), 4L)
    ref <- vapply(seq_len(nrow(es)), function(e)
      generic_score(f, es$i[e], es$j[e], v), numeric(1))
    expect_equal(es$stat, ref, tolerance = 1e-9, info = v)
  }
  # from the empty support as well
  f0 <- asNamespace("mcglmgraph")$.refit_between(
    f, mc_between_graph("ggm", matrix(0L, 0, 2)))
  es0 <- mc_edge_scores(f0, "score")
  ref0 <- vapply(seq_len(nrow(es0)), function(e)
    generic_score(f0, es0$i[e], es0$j[e], "model"), numeric(1))
  expect_equal(es0$stat, ref0, tolerance = 1e-9)
})

test_that("each block of the edge scores is the one of the joint call", {
  f <- tool_fit(edges = rbind(c(1, 2), c(2, 3)))
  both <- mc_edge_scores(f)
  expect_equal(mc_edge_scores(f, "wald"), both[both$present, ],
               ignore_attr = TRUE)
  expect_equal(mc_edge_scores(f, "score"), both[!both$present, ],
               ignore_attr = TRUE)
  expect_false(anyNA(both$stat))
})

test_that("the separation kernel is the partial correlation with Fisher's z", {
  set.seed(3)
  A <- matrix(rnorm(25), 5); S <- crossprod(A) / 5 + diag(5)
  ns <- asNamespace("mcglmgraph")
  for (cond in list(integer(0), 3L, c(3L, 5L))) {
    P <- solve(S[c(1, 2, cond), c(1, 2, cond)])
    pc <- -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
    r <- ns$.ci_pcor(S, 1, 2, cond, n = 80)
    expect_equal(r$pcor, pc, tolerance = 1e-12)
    expect_equal(r$p_value,
                 2 * stats::pnorm(-abs(atanh(pc)) * sqrt(80 - length(cond) - 3)),
                 tolerance = 1e-12)
    expect_identical(r$inversion_status, "exact")
    # the enumeration of separating sets reports the same test
    adj <- matrix(1L, 5, 5); diag(adj) <- 0L; adj[1, 2] <- adj[2, 1] <- 0L
    en <- ns$cpp_sepsets_pairs(S, 80, adj, matrix(c(1L, 2L), 1), 1e-12, 3L, 0)
    hit <- which(apply(en$set_members, 1L, function(m)
      identical(which(m == 1L), as.integer(sort(cond)))))
    expect_identical(en$set_pcor[hit], r$pcor)
    expect_identical(en$set_p_value[hit], r$p_value)
  }
  expect_error(ns$cpp_ci_pcor(S, 1L, 1L, integer(0), 80, 0), "invalid node")
  expect_error(ns$cpp_ci_pcor(S, 1L, 2L, 2L, 80, 0), "invalid node")
  expect_error(ns$cpp_ci_pcor(S, 1L, 9L, integer(0), 80, 0), "invalid node")
})

test_that("the Fisher interval of a partial correlation stays in its range", {
  f <- tool_fit()
  w <- mc_partial_cor_delta(f)
  z <- mc_partial_cor_delta(f, interval = "fisher")
  expect_identical(z[, c("pcor", "se", "statistic", "p_value")],
                   w[, c("pcor", "se", "statistic", "p_value")])
  q <- stats::qnorm(0.975)
  expect_equal(z$lower, tanh(atanh(w$pcor) - q * w$se / (1 - w$pcor^2)))
  expect_equal(z$upper, tanh(atanh(w$pcor) + q * w$se / (1 - w$pcor^2)))
  expect_true(all(z$lower > -1 & z$upper < 1 & z$lower < z$pcor &
                    z$pcor < z$upper))
  expect_identical(attr(z, "interval"), "fisher")
  expect_identical(attr(w, "interval"), "wald")
})

test_that("the logit interval of R squared stays in its range", {
  f <- tool_fit(type = "dag")
  w <- mc_dag_coefficients(f)$nodes
  l <- mc_dag_coefficients(f, r2_interval = "logit")$nodes
  expect_identical(l$r_squared, w$r_squared)
  # a response without parents has R squared zero, and the interval is it
  root <- w$n_parents == 0L
  expect_identical(sum(root), 1L)
  expect_identical(c(l$lower_r_squared[root], l$upper_r_squared[root]), c(0, 0))
  q <- stats::qnorm(0.975); r2 <- w$r_squared[!root]; se <- w$se_r_squared[!root]
  expect_equal(w$lower_r_squared[!root], r2 - q * se)
  expect_equal(l$lower_r_squared[!root],
               stats::plogis(stats::qlogis(r2) - q * se / (r2 * (1 - r2))))
  expect_true(all(l$lower_r_squared[!root] > 0 & l$upper_r_squared[!root] < 1))
})

test_that("summary() returns an object that prints both stages", {
  f <- tool_fit()
  expect_silent(s <- summary(f))
  expect_s3_class(s, "summary.mcglm_graph")
  expect_identical(s$variability, "model")
  expect_true(s$converged)
  expect_equal(s$coefficients$Std.Error, unname(sqrt(diag(f$graph$vcov))))
  expect_equal(s$partial_cor, mc_partial_cor(f))
  expect_true("Correlation" %in% names(s$first_stage))
  out <- capture.output(print(s))
  expect_true(any(grepl("Regression:", out, fixed = TRUE)))
  expect_true(any(grepl("Graph parameters (model variability):", out,
                        fixed = TRUE)))
  e <- summary(f, variability = "empirical")
  expect_equal(e$coefficients$Std.Error,
               unname(sqrt(diag(mc_vcov_graph(f, "sandwich")))))
  expect_true(any(grepl("empirical variability", capture.output(print(e)))))
})

test_that("the arrows of a DAG stop at the border of the nodes", {
  ns <- asNamespace("mcglmgraph")
  # a device that is not square, where user units differ between the axes
  grDevices::pdf(NULL, width = 8, height = 4); on.exit(grDevices::dev.off())
  graphics::plot(0, 0, type = "n", xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3))
  inch <- function(x, y) c(graphics::grconvertX(x, "user", "inches"),
                           graphics::grconvertY(y, "user", "inches"))
  r <- 0.1
  for (p in list(c(-1, 0, 1, 0), c(0, -1, 0, 1), c(-0.5, -0.8, 0.9, 0.4))) {
    en <- ns$.arrow_ends(p[1], p[2], p[3], p[4], r)
    # each end sits r inches from its node, on the segment that joins them
    expect_equal(sqrt(sum((inch(en[1], en[2]) - inch(p[1], p[2]))^2)), r)
    expect_equal(sqrt(sum((inch(en[3], en[4]) - inch(p[3], p[4]))^2)), r)
    full <- sqrt(sum((inch(p[3], p[4]) - inch(p[1], p[2]))^2))
    expect_equal(sqrt(sum((inch(en[3], en[4]) - inch(en[1], en[2]))^2)),
                 full - 2 * r)
  }
  # too short to cut
  expect_identical(ns$.arrow_ends(0, 0, 1e-3, 0, r), c(0, 0, 1e-3, 0))
  # the whole plot runs on a DAG and on a GGM fit
  expect_identical(sum(mc_plot_network(tool_fit(type = "dag"),
                                       threshold = 0) != 0), 6L)
  fg <- tool_fit()
  expect_silent(mc_plot_network(fg))
})
