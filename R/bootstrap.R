# bootstrap.R: resampling that refits both stages.
# Part of the mcglmgraph package. License: GPL-3.
#
# Every conditional variance the package reports treats the Pearson residuals
# as observed. The analytical correction of mc_vcov_graph(type = "two_stage")
# propagates the first stage by the shortcut or by the full form; the bootstrap
# here propagates it by construction, for any family, at the cost of refitting.

# The resampling rebuilds the first-stage call from the fitted object. With an
# identity within-response structure the rows are the units and can be drawn
# one by one. A richer matrix linear predictor describes dependence among rows,
# so the rows it connects have to be drawn together: whole clusters are drawn,
# and each drawn copy carries its own block of the structure.
.identity_matrix_pred <- function(fit) {
  mp <- fit$matrix_pred
  if (is.null(mp)) return(FALSE)
  all(vapply(mp, function(z) {
    if (!is.list(z) || length(z) != 1L) return(FALSE)
    .is_identity_structure(z[[1L]])
  }, logical(1)))
}

# Is this within-response structure the identity? The default one arrives as the
# sparse diagonal object the first stage builds, which carries the answer in its
# own slots, so the common case is decided without densifying anything: the
# dense copy of a structure with one row per unit would be the square of the
# sample size in memory.
#
# The dimension is read from the slot rather than through dim(). Matrix sits in
# Imports, so in a session that has not yet attached it the S4 method for dim()
# is not dispatched and dim() answers NULL on a ddiMatrix, which would send
# this function down the FALSE branch before it could reach the ddiMatrix case
# below. Slots are attributes, so attr() reaches them with no namespace loaded.
.is_identity_structure <- function(M) {
  d <- if (isS4(M)) attr(M, "Dim") else dim(M)
  if (is.null(d)) d <- dim(M)
  if (is.null(d) || length(d) != 2L || d[1L] != d[2L]) return(FALSE)
  if (inherits(M, "ddiMatrix"))
    return(identical(attr(M, "diag"), "U") ||
             (length(M@x) == d[1L] && all(abs(M@x - 1) < 1e-12)))
  if (is.matrix(M))
    return(isTRUE(all.equal(M, diag(d[1L]), check.attributes = FALSE)))
  # A sparse matrix is settled on its stored entries, at any size.
  if (isS4(M)) {
    ans <- tryCatch(Matrix::isDiagonal(Matrix::drop0(M)) &&
                      all(abs(Matrix::diag(M) - 1) < 1e-12),
                    error = function(e) NA)
    if (!is.na(ans)) return(isTRUE(ans))
  }
  if (d[1L] > 4000L) return(FALSE)          # too large to settle by coercion
  tryCatch(isTRUE(all.equal(suppressWarnings(as.matrix(M)), diag(d[1L]),
                            check.attributes = FALSE)),
           error = function(e) FALSE)
}

# One entry point for every resampling routine in the package: refit on a set of
# rows, either both stages or the graph layer alone.
.refit_on_rows <- function(fit, rows, stage, spec, ctrl, blocks = NULL) {
  if (identical(stage, "full") || identical(stage, "both"))
    .refit_both_stages(fit, rows, spec, ctrl, blocks = blocks) else
      .graph_layer_on_residuals(fit, rows, spec, ctrl)
}

# What a routine that refits the first stage needs from the object, checked
# once and refused with the same words wherever it is asked for. A structure
# that is not the identity needs clusters, and the clusters have to contain
# every pair of rows the structure connects, or drawing them whole would cut it.
.check_resampling_ready <- function(fit, stage, cluster = NULL) {
  if (!(identical(stage, "full") || identical(stage, "both"))) return(invisible(TRUE))
  if (is.null(fit$data) || nrow(fit$data) != fit$graph$n)
    stop("the fit does not carry the data its first stage was built on, so ",
         "the two stages cannot be refitted.", call. = FALSE)
  if (.identity_matrix_pred(fit)) return(invisible(TRUE))
  if (is.null(cluster))
    stop("refitting both stages with a within-response structure that is not ",
         "the identity needs cluster: the rows the structure connects must be ",
         "drawn together, and each drawn cluster keeps its own block of the ",
         "structure. Pass one identifier per row, grouping those rows, or ",
         "resample the residuals instead.", call. = FALSE)
  u <- .two_stage_units(fit, cluster)
  if (!.structure_nested_in(fit, u$unit))
    stop("the within-response structure connects rows of different clusters, ",
         "so drawing whole clusters would cut it; every matrix of matrix_pred ",
         "must vanish between distinct clusters.", call. = FALSE)
  invisible(TRUE)
}

# The within-response structure of a resample. Each drawn block of rows keeps
# its own submatrix, and the blocks are placed on the diagonal in the order they
# were drawn, which is the order of the rows of the resampled data. A cluster
# drawn twice gives two blocks that do not touch; Z[rows, rows] would instead
# join the two copies, as if they were the same units observed twice.
.resample_structure <- function(Z, blocks) {
  Matrix::bdiag(lapply(blocks, function(idx) Z[idx, idx, drop = FALSE]))
}

# The graph control of a fit, with the defaults filled in.
.graph_control <- function(fit) {
  ctrl <- fit$graph$control
  if (is.null(ctrl))
    ctrl <- list(max_iter = 100L, tol = 1e-8, score_tol = 1e-4,
                 rcond_tol = .Machine$double.eps^(2 / 3),
                 dag_solver = "closed")
  if (is.null(ctrl$dag_solver)) ctrl$dag_solver <- "closed"
  ctrl
}

# Refit the first stage on a set of rows, then the graph layer on its residuals.
# With a structure that is not the identity, `blocks` lists the rows of each
# drawn cluster, in the order in which `rows` concatenates them.
.refit_both_stages <- function(fit, rows, spec, ctrl, blocks = NULL) {
  dat <- fit$data[rows, , drop = FALSE]
  rownames(dat) <- NULL
  R <- length(fit$linear_pred)
  n0 <- nrow(fit$data)
  mp <- if (.identity_matrix_pred(fit)) {
    rep(list(mcglm::mc_id(dat)), R)
  } else {
    if (is.null(blocks) || !identical(as.integer(unlist(blocks)), as.integer(rows)))
      stop("a structure that is not the identity is resampled by whole ",
           "clusters only.", call. = FALSE)
    lapply(fit$matrix_pred, function(zl)
      lapply(zl, .resample_structure, blocks = blocks))
  }

  sub_list <- function(z) {
    if (is.null(z)) return(NULL)
    lapply(z, function(v) if (is.null(v) || length(v) != n0) v else v[rows])
  }
  w <- fit$weights
  if (!is.null(w)) {
    if (length(w) == n0 * R)
      w <- as.numeric(matrix(w, n0, R)[rows, , drop = FALSE])
    else if (length(w) == n0) w <- w[rows]
  }

  m <- mcglm::mcglm(
    linear_pred = fit$linear_pred,
    matrix_pred = mp,
    link = fit$link, variance = fit$variance, covariance = fit$covariance,
    offset = sub_list(fit$offset), Ntrial = sub_list(fit$Ntrial),
    power_fixed = fit$power_fixed, data = dat, weights = w,
    control_algorithm = fit$con, control_initial = fit$list_initial
  )
  .graph_layer_on(m, fit, spec, ctrl)
}

# Attach a graph layer to a first-stage fit, reusing the machinery of
# mcglm_graph() without going through it.
.graph_layer_on <- function(m, template, spec, ctrl) {
  R <- template$graph$R
  rn <- template$graph$response_names
  E <- as.matrix(stats::residuals(m, type = "pearson"))
  if (ncol(E) != R) E <- matrix(E, ncol = R)
  colnames(E) <- rn
  n <- nrow(E)
  Ec <- scale(E, center = TRUE, scale = FALSE)
  Shat <- crossprod(Ec) / n
  dimnames(Shat) <- list(rn, rn)
  spec <- .canonical_spec(spec, R)
  b <- .fit_between(Shat, n, spec, max_iter = ctrl$max_iter, tol = ctrl$tol,
                    score_tol = ctrl$score_tol, rcond_tol = ctrl$rcond_tol,
                    verbose = FALSE, dag_solver = ctrl$dag_solver)
  # The refitted first stage is the object returned, so a statistic that reads
  # the regression or dispersion parameters sees those of the resample.
  m$graph <- .graph_block(b, spec, Shat, E, n, rn, template$graph$control)
  class(m) <- c("mcglm_graph", setdiff(class(m), "mcglm_graph"))
  m
}

# Resamples of units drawn with replacement, keeping their number: rows, or
# whole clusters, whose rows are also returned as blocks so that a structure
# that is not the identity can be rebuilt one drawn cluster at a time. The calls
# to the generator are made in a fixed order, so a given seed always yields the
# same resamples.
.draw_units <- function(B, n, by_cluster = NULL) {
  lapply(seq_len(B), function(b) {
    if (is.null(by_cluster))
      return(list(rows = sample.int(n, n, replace = TRUE), blocks = NULL))
    k <- sample.int(length(by_cluster), length(by_cluster), replace = TRUE)
    blocks <- unname(by_cluster[k])
    list(rows = unlist(blocks, use.names = FALSE), blocks = blocks)
  })
}

# Default statistic: the graph parameters and the off-diagonal partial
# correlations, which are the two things the papers report.
.bootstrap_statistic <- function(f) {
  g <- f$graph
  P <- mc_partial_cor(f)
  ut <- upper.tri(P)
  nm <- outer(g$response_names, g$response_names,
              function(a, b) paste0("pcor_", a, "_", b))
  stats::setNames(c(as.numeric(g$par), P[ut]),
                  c(names(g$par), nm[ut]))
}

#' Bootstrap the graph layer with both stages refitted
#'
#' Resamples units, or clusters of units, with replacement and refits the
#' first-stage MCGLM and the graph layer in every resample, returning the
#' bootstrap distribution of a statistic of the fit (by default the graph
#' parameters and the off-diagonal partial correlations). The support of the
#' graph is held fixed at the one the supplied fit carries, so the result is an
#' inferential check on a given structure and not a post-selection estimator.
#'
#' @details
#' With \eqn{n} units (rows, or clusters when \code{cluster} is given), each
#' resample draws \eqn{n} units with replacement (Efron, 1979), refits, and
#' evaluates the statistic \eqn{T(\cdot)}. With \eqn{B_v} valid resamples
#' \eqn{T^{*}_1,\dots,T^{*}_{B_v}}, the bootstrap covariance is the sample
#' covariance of the replicates, \code{se} their standard deviations, the bias
#' estimate \eqn{\bar T^{*}-T(\hat{\mathrm F}_n)}, and the interval for level
#' \eqn{1-a} the percentile interval with limits at the quantiles \eqn{a/2} and
#' \eqn{1-a/2} of the replicates. A resample is dropped, and recorded with its
#' message in \code{failures}, when its refit raises an error, when its graph
#' layer does not converge, or when its statistic is not a finite vector of the
#' length obtained on the original sample; the summaries use the valid
#' resamples only, and at least two are required. Resampling
#' independent clusters is the scheme that remains valid when rows within a
#' cluster are dependent. Because both stages are refitted, the replicates carry
#' the variability of the first stage, which the model-based and sandwich
#' covariances of \code{\link{mc_vcov_graph}} take as known.
#'
#' @section What it is for:
#' The variances the package reports elsewhere are conditional on the estimated
#' Pearson residuals unless the analytical correction of
#' \code{\link{mc_vcov_graph}} is asked for. Refitting propagates the first
#' stage by construction, for any family, which makes this the yardstick
#' against which both routes of that correction, the shortcut and the full
#' form, are checked. With \code{stage = "residual"} only the residuals are
#' resampled and the first stage is held fixed, which is the comparison that
#' isolates how much the first stage contributes.
#'
#' @section What it requires:
#' The first-stage call is rebuilt from the fitted object. With an identity
#' within-response structure the rows are resampled, or whole clusters when
#' \code{cluster} is given. A richer matrix linear predictor describes
#' dependence among rows, so the rows it connects must be drawn together:
#' \code{cluster} is then required, every matrix of the structure must vanish
#' between distinct clusters, which is checked, and each drawn cluster keeps
#' its own block of the structure. A cluster drawn twice gives two blocks that
#' do not touch, as two independent units would; the submatrix of the rows
#' drawn would instead join the two copies. Without \code{cluster} such a fit
#' is refused.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param B number of resamples, an integer of at least two; default 200.
#' @param cluster optional vector with one identifier per row, with no missing
#'   values and at least two distinct values; whole clusters are resampled,
#'   with replacement, keeping the total number of clusters. Required when the
#'   within-response structure is not the identity and \code{stage = "both"}.
#' @param statistic optional function of a refitted \code{"mcglm_graph"} object
#'   returning a named numeric vector, of the same length at every resample. The
#'   default returns the graph parameters and the off-diagonal partial
#'   correlations.
#' @param stage \code{"both"} refits the first stage and the graph layer;
#'   \code{"residual"} resamples the stored Pearson residuals and refits the
#'   graph layer alone.
#' @param level confidence level of the percentile intervals, a number strictly
#'   between 0 and 1; default 0.95.
#' @param workers number of parallel workers, default 1 (sequential). The
#'   resamples are drawn in the parent before any dispatch, so a given seed
#'   yields the same resamples at any worker count.
#' @return An object of class \code{"mcg_bootstrap"}, a list with
#'   \code{estimate}, the named numeric vector of the statistic on the original
#'   sample; \code{replicates}, the \code{B_valid} by length-of-statistic matrix
#'   of valid resamples; \code{vcov}, their covariance matrix, and \code{se},
#'   their standard deviations; \code{bias}, the mean of the replicates minus
#'   \code{estimate}; \code{ci}, the two-column matrix (\code{lower},
#'   \code{upper}) of percentile intervals; \code{level} and \code{stage}, as
#'   given; the counts \code{B}, \code{B_valid} and \code{B_failed}; the data
#'   frame \code{failures} (\code{replicate}, \code{error}); \code{clustered}
#'   (logical) and \code{n_clusters} (\code{NA} without clusters). The
#'   \code{print} method shows the estimates, standard errors, biases and
#'   intervals and returns its argument invisibly.
#' @examples
#' \donttest{
#' set.seed(1); n <- 200; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' bt <- mc_graph_bootstrap(f, B = 40)
#' bt$se
#' }
#' @seealso \code{\link{mc_vcov_graph}} for the analytical forms and
#'   \code{\link{mc_graph_stability}} for selection frequencies.
#' @family resampling
#' @references
#' Efron, B. (1979). Bootstrap methods: another look at the jackknife.
#' \emph{The Annals of Statistics} 7(1), 1--26. \doi{10.1214/aos/1176344552}
#' @export
mc_graph_bootstrap <- function(fit, B = 200L, cluster = NULL,
                               statistic = NULL,
                               stage = c("both", "residual"),
                               level = 0.95, workers = 1L) {
  stopifnot(inherits(fit, "mcglm_graph"))
  stage <- match.arg(stage)
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != as.integer(B))
    stop("B must be an integer of at least two.", call. = FALSE)
  B <- as.integer(B)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (is.null(statistic)) statistic <- .bootstrap_statistic
  if (!is.function(statistic))
    stop("statistic must be a function of a fitted object.", call. = FALSE)

  g <- fit$graph
  n <- g$n

  if (!is.null(cluster)) {
    if (length(cluster) != n)
      stop("cluster must have one identifier per observational row.",
           call. = FALSE)
    if (anyNA(cluster))
      stop("cluster identifiers must not be missing.", call. = FALSE)
    cluster <- as.factor(cluster)
    if (nlevels(cluster) < 2L)
      stop("at least two independent clusters are required.", call. = FALSE)
    by_cluster <- split(seq_len(n), cluster)
  }
  .check_resampling_ready(fit, stage, cluster)

  # Drawn in the parent, before any dispatch, so the worker count changes who
  # does the arithmetic and never which rows each resample holds.
  draws <- .draw_units(B, n, if (is.null(cluster)) NULL else by_cluster)

  spec <- g$spec
  ctrl <- .graph_control(fit)

  t0 <- statistic(fit)
  if (!is.numeric(t0) || is.null(names(t0)))
    stop("statistic must return a named numeric vector.", call. = FALSE)

  one <- function(b) {
    err <- NULL
    val <- tryCatch({
      f <- .refit_on_rows(fit, draws[[b]]$rows, stage, spec, ctrl,
                          blocks = draws[[b]]$blocks)
      if (!isTRUE(f$graph$converged))
        stop("the graph layer did not converge: ", f$graph$status)
      statistic(f)
    }, error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(val) || length(val) != length(t0) || any(!is.finite(val)))
      return(list(ok = FALSE, error = if (is.null(err))
        "the statistic was not a finite vector of the expected length" else err))
    list(ok = TRUE, value = val)
  }

  res <- .mcg_lapply(seq_len(B), one, workers = workers)

  keep <- vector("list", B)
  failures <- data.frame(replicate = integer(), error = character(),
                         stringsAsFactors = FALSE)
  for (b in seq_len(B)) {
    z <- res[[b]]
    msg <- .mcg_error_message(z)
    if (!is.null(msg)) z <- list(ok = FALSE, error = msg)
    if (isTRUE(z$ok)) keep[[b]] <- z$value else
      failures <- rbind(failures, data.frame(replicate = b, error = z$error,
                                             stringsAsFactors = FALSE))
  }
  keep <- Filter(Negate(is.null), keep)
  B_valid <- length(keep)
  if (B_valid < 2L)
    stop("fewer than two resamples produced a usable fit; the first row of ",
         "the failure table reads: ",
         if (nrow(failures)) failures$error[1L] else "no diagnosis available",
         call. = FALSE)

  reps <- do.call(rbind, keep)
  colnames(reps) <- names(t0)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  ci <- t(apply(reps, 2L, stats::quantile, probs = probs, names = FALSE))
  colnames(ci) <- c("lower", "upper")

  structure(list(
    estimate = t0, replicates = reps,
    vcov = stats::cov(reps),
    se = apply(reps, 2L, stats::sd),
    bias = colMeans(reps) - t0,
    ci = ci, level = level, stage = stage,
    B = B, B_valid = B_valid, B_failed = B - B_valid, failures = failures,
    clustered = !is.null(cluster),
    n_clusters = if (is.null(cluster)) NA_integer_ else nlevels(cluster)
  ), class = "mcg_bootstrap")
}

# Residual-only resampling: the first stage is held fixed and its stored
# residuals are resampled, which is what isolates the second stage.
.graph_layer_on_residuals <- function(fit, rows, spec, ctrl) {
  g <- fit$graph
  E <- as.matrix(g$residuals)[rows, , drop = FALSE]
  Ec <- scale(E, center = TRUE, scale = FALSE)
  n <- nrow(E)
  Shat <- crossprod(Ec) / n
  dimnames(Shat) <- list(g$response_names, g$response_names)
  spec <- .canonical_spec(spec, g$R)
  b <- .fit_between(Shat, n, spec, max_iter = ctrl$max_iter, tol = ctrl$tol,
                    score_tol = ctrl$score_tol, rcond_tol = ctrl$rcond_tol,
                    verbose = FALSE, dag_solver = ctrl$dag_solver)
  out <- fit
  out$graph <- .graph_block(b, spec, Shat, E, n, g$response_names, g$control)
  out
}

#' @rdname mc_graph_bootstrap
#' @param x an object of class \code{"mcg_bootstrap"}.
#' @param ... ignored.
#' @export
print.mcg_bootstrap <- function(x, ...) {
  cat(sprintf("Bootstrap of the graph layer (%s, B = %d, level = %.2f)\n",
              if (x$stage == "both") "both stages refitted" else
                "residuals resampled", x$B, x$level))
  cat(sprintf("Valid resamples: %d | failed: %d%s\n", x$B_valid, x$B_failed,
              if (x$clustered)
                sprintf(" | clusters: %d", x$n_clusters) else ""))
  tab <- data.frame(estimate = x$estimate, se = x$se, bias = x$bias,
                    lower = x$ci[, "lower"], upper = x$ci[, "upper"])
  print(round(tab, 4))
  invisible(x)
}

# The four states a pair of responses can be in, read off an adjacency in the
# package convention. They are exclusive and exhaust the possibilities.
.pair_states <- function(amat) {
  R <- nrow(amat)
  pairs <- t(utils::combn(R, 2L))
  vapply(seq_len(nrow(pairs)), function(e) {
    i <- pairs[e, 1L]; j <- pairs[e, 2L]
    a <- amat[i, j] == 1L; b <- amat[j, i] == 1L
    if (!a && !b) 1L else if (a && b) 2L else if (a) 3L else 4L
  }, integer(1))
}

#' Resampling frequencies of the skeleton and of each orientation
#'
#' Repeats the constraint-based discovery of
#' \code{\link{mc_graph_pipeline}} over resamples and records, for every pair of
#' responses, how often the pair came out absent, undirected, or directed each
#' way. An orientation that the data barely support shows up here as a pair
#' whose two directions split the frequency, which a single run cannot reveal.
#'
#' @section What each resample is:
#' With \code{refit = "residual"} the resample is drawn from the stored Pearson
#' residuals and the first stage is held fixed, so the frequencies describe the
#' discovery step alone. With \code{refit = "full"} rows of the data are
#' resampled and both stages are refitted, so they carry the variability of the
#' whole procedure. Units are drawn with replacement, keeping the sample size,
#' which is the bootstrap convention (Efron, 1979);
#' \code{\link{mc_graph_stability}} instead uses half-samples, the convention
#' of the stability-selection literature (Meinshausen and Buhlmann, 2010). With
#' \eqn{B_v} valid resamples, the frequency of a state is the proportion of
#' valid resamples in which the pair is in it, so the four frequencies of a pair
#' sum to one; a resample whose refit or pipeline raises an error is dropped and
#' recorded in \code{failures}, and at least one valid resample is required.
#'
#' @section What the frequencies estimate:
#' Each frequency estimates, by resampling, the probability that the discovery
#' procedure, at the given \code{alpha} and orientation rule, returns that
#' state for the pair in a sample of the same size. It measures the stability
#' of the output of the procedure. It is not the probability that the edge or
#' the orientation belongs to the structure that generated the data, and it
#' carries no control of an error rate.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param B number of resamples, an integer of at least two; default 100.
#' @param refit \code{"residual"} (default) or \code{"full"}, as above;
#'   \code{"full"} corresponds to \code{stage = "both"} of
#'   \code{\link{mc_graph_bootstrap}}.
#' @param cluster optional vector with one identifier per row, with no missing
#'   values and at least two distinct values; whole clusters are resampled.
#'   Required with \code{refit = "full"} when the within-response structure is
#'   not the identity, as in \code{\link{mc_graph_bootstrap}}.
#' @param workers number of parallel workers, default 1. The resamples are drawn
#'   in the parent, so a given seed yields the same ones at any worker count.
#' @inheritParams mc_graph_pipeline
#' @return An object of class \code{"mcg_cpdag_stability"}: a list with
#'   \code{states}, a data frame with one row per pair of responses
#'   (\code{i}, \code{j}, \code{response1}, \code{response2}) and the
#'   frequencies of the four states (\code{absent}, \code{undirected},
#'   \code{i_to_j}, \code{j_to_i}); the symmetric matrix \code{edge_prob} of
#'   skeleton frequencies, one minus \code{absent}; \code{consensus}, the 0/1
#'   adjacency matrix that gives each pair its most frequent state, the first
#'   in the order above on a tie (an entry \code{[i, j]} equal to one alone
#'   reads \eqn{i\to j}, with \code{[j, i]} also one the edge is undirected);
#'   being formed pair by pair, it need not be the CPDAG of a DAG;
#'   \code{valid_prob}, the proportion of valid resamples whose output passed
#'   the validity gates of \code{\link{mc_graph_pipeline}};
#'   \code{mean_ambiguous}, the average number of ambiguous triples per valid
#'   resample; the counts \code{B}, \code{B_valid} and \code{B_failed};
#'   \code{B_incomplete}, the number of valid resamples in which the skeleton
#'   search skipped a test that could not be formed, as in
#'   \code{\link{mc_pc_skeleton}}, reported in one warning when it is not
#'   zero; the
#'   data frame \code{failures} (\code{replicate}, \code{error}); and the
#'   settings \code{refit}, \code{alpha}, \code{orient_method},
#'   \code{orient_rule}, \code{clustered} (logical) and
#'   \code{response_names}. The \code{print} method shows the settings, the
#'   counts and the table of frequencies, and returns its argument invisibly.
#' @examples
#' \donttest{
#' set.seed(1); n <- 400; x <- rnorm(n)
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
#' d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' mc_cpdag_stability(f, B = 50)
#' }
#' @seealso \code{\link{mc_graph_pipeline}} for a single run and
#'   \code{\link{mc_graph_bootstrap}} for the bootstrap of the graph
#'   parameters.
#' @family resampling
#' @references
#' Efron, B. (1979). Bootstrap methods: another look at the jackknife.
#' \emph{The Annals of Statistics} 7(1), 1--26. \doi{10.1214/aos/1176344552}
#'
#' Meinshausen, N. and Buhlmann, P. (2010). Stability selection. \emph{Journal
#' of the Royal Statistical Society, Series B} 72(4), 417--473.
#' \doi{10.1111/j.1467-9868.2010.00740.x}
#' @export
mc_cpdag_stability <- function(fit, B = 100L, alpha = 0.05,
                               orient_method = c("pc", "hc", "tabu"),
                               orient_rule = c("standard", "conservative",
                                               "majority"),
                               required = NULL, forbidden = NULL,
                               refit = c("residual", "full"),
                               cluster = NULL, max_ord = NULL,
                               df_residual = NULL, workers = 1L) {
  stopifnot(inherits(fit, "mcglm_graph"))
  orient_method <- match.arg(orient_method)
  orient_rule <- match.arg(orient_rule)
  refit <- match.arg(refit)
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != as.integer(B))
    stop("B must be an integer of at least two.", call. = FALSE)
  B <- as.integer(B)

  g <- fit$graph
  n <- g$n
  R <- g$R
  rn <- g$response_names

  if (!is.null(cluster)) {
    if (length(cluster) != n)
      stop("cluster must have one identifier per observational row.",
           call. = FALSE)
    if (anyNA(cluster))
      stop("cluster identifiers must not be missing.", call. = FALSE)
    cluster <- as.factor(cluster)
    if (nlevels(cluster) < 2L)
      stop("at least two independent clusters are required.", call. = FALSE)
    by_cluster <- split(seq_len(n), cluster)
  }
  .check_resampling_ready(fit, refit, cluster)
  # "auto" is resolved again inside every resample, on the fit of that resample;
  # a fit for which it is not exact is refused here, once, rather than B times.
  if (identical(df_residual, "auto")) invisible(.resolve_df_residual("auto", fit))
  draws <- .draw_units(B, n, if (is.null(cluster)) NULL else by_cluster)

  spec <- g$spec
  ctrl <- .graph_control(fit)

  one <- function(b) {
    err <- NULL
    out <- tryCatch({
      f <- .refit_on_rows(fit, draws[[b]]$rows, refit, spec, ctrl,
                          blocks = draws[[b]]$blocks)
      # skipped separation tests are reported once, for the whole run
      cp <- withCallingHandlers(
        mc_graph_pipeline(f, alpha = alpha, orient_method = orient_method,
                          orient_rule = orient_rule, required = required,
                          forbidden = forbidden, max_ord = max_ord,
                          df_residual = df_residual),
        mcg_unusable_tests = function(w) invokeRestart("muffleWarning"))
      list(state = .pair_states(cp$cpdag), valid = isTRUE(cp$valid_cpdag),
           incomplete = isTRUE(cp$separation_diagnostics$n_unusable > 0),
           amb = if (is.null(cp$n_ambiguous)) 0L else cp$n_ambiguous)
    }, error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(out)) return(list(ok = FALSE, error = err))
    c(list(ok = TRUE), out)
  }

  res <- .mcg_lapply(seq_len(B), one, workers = workers)

  pairs <- t(utils::combn(R, 2L))
  counts <- matrix(0L, nrow(pairs), 4L)
  n_valid <- 0L; n_amb <- 0; B_valid <- 0L; B_incomplete <- 0L
  failures <- data.frame(replicate = integer(), error = character(),
                         stringsAsFactors = FALSE)
  for (b in seq_len(B)) {
    z <- res[[b]]
    msg <- .mcg_error_message(z)
    if (!is.null(msg)) z <- list(ok = FALSE, error = msg)
    if (!isTRUE(z$ok)) {
      failures <- rbind(failures, data.frame(replicate = b, error = z$error,
                                             stringsAsFactors = FALSE))
      next
    }
    B_valid <- B_valid + 1L
    counts[cbind(seq_len(nrow(pairs)), z$state)] <-
      counts[cbind(seq_len(nrow(pairs)), z$state)] + 1L
    n_valid <- n_valid + as.integer(z$valid)
    n_amb <- n_amb + z$amb
    B_incomplete <- B_incomplete + as.integer(isTRUE(z$incomplete))
  }
  if (B_valid < 1L)
    stop("no resample produced a usable structure; the first failure reads: ",
         if (nrow(failures)) failures$error[1L] else "no diagnosis available",
         call. = FALSE)

  if (B_incomplete > 0L)
    warning("in ", B_incomplete, " of ", B_valid, " valid resamples some ",
            "tests of conditional association could not be formed and were ",
            "skipped; those outputs do not pass the validity gates.",
            call. = FALSE)

  prob <- counts / B_valid
  states <- data.frame(
    i = pairs[, 1L], j = pairs[, 2L],
    response1 = rn[pairs[, 1L]], response2 = rn[pairs[, 2L]],
    absent = prob[, 1L], undirected = prob[, 2L],
    i_to_j = prob[, 3L], j_to_i = prob[, 4L],
    row.names = NULL, stringsAsFactors = FALSE
  )
  edge_prob <- matrix(0, R, R, dimnames = list(rn, rn))
  edge_prob[cbind(pairs[, 1L], pairs[, 2L])] <- 1 - states$absent
  edge_prob[cbind(pairs[, 2L], pairs[, 1L])] <- 1 - states$absent

  consensus <- matrix(0L, R, R, dimnames = list(rn, rn))
  best <- max.col(prob, ties.method = "first")
  for (e in seq_len(nrow(pairs))) {
    i <- pairs[e, 1L]; j <- pairs[e, 2L]
    if (best[e] == 2L) { consensus[i, j] <- consensus[j, i] <- 1L }
    else if (best[e] == 3L) consensus[i, j] <- 1L
    else if (best[e] == 4L) consensus[j, i] <- 1L
  }

  structure(list(
    states = states, edge_prob = edge_prob, consensus = consensus,
    valid_prob = n_valid / B_valid, mean_ambiguous = n_amb / B_valid,
    B = B, B_valid = B_valid, B_failed = B - B_valid,
    B_incomplete = B_incomplete, failures = failures,
    refit = refit, alpha = alpha, orient_method = orient_method,
    orient_rule = orient_rule, clustered = !is.null(cluster),
    response_names = rn
  ), class = "mcg_cpdag_stability")
}

#' @rdname mc_cpdag_stability
#' @param x an object of class \code{"mcg_cpdag_stability"}.
#' @param ... ignored.
#' @export
print.mcg_cpdag_stability <- function(x, ...) {
  cat(sprintf("Structure frequencies (%s, %s rule, alpha = %.3g, B = %d)\n",
              x$orient_method, x$orient_rule, x$alpha, x$B))
  cat(sprintf("Resamples: %s | valid: %d | failed: %d\n",
              if (identical(x$refit, "full"))
                "rows resampled, both stages refitted" else
                  "residuals resampled, first stage held fixed",
              x$B_valid, x$B_failed))
  if (isTRUE(x$B_incomplete > 0L))
    cat(sprintf("Valid resamples with skipped separation tests: %d\n",
                x$B_incomplete))
  cat(sprintf("Outputs passing the validity gates: %.0f%%",
              100 * x$valid_prob))
  if (x$mean_ambiguous > 0)
    cat(sprintf(" | ambiguous triples per resample: %.2f", x$mean_ambiguous))
  cat("\n\n")
  print(format(x$states[, c("response1", "response2", "absent", "undirected",
                            "i_to_j", "j_to_i")], digits = 2))
  invisible(x)
}
