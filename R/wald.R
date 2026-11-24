# wald.R: Wald tests of linear hypotheses on the graph layer.
# Part of the mcglmgraph package. License: GPL-3.
#
# One statistic serves every test of the graph layer: the joint test of all
# edges, the adequacy of a selected support against the saturated graph, and any
# contrast a user writes down. They differ only in the contrast and in the fit
# the contrast is applied to, so they share this code and the covariance of
# mc_vcov_graph(), and a change of convention cannot reach one and miss another.

# A contrast given as a matrix, as one row, or as the names of the parameters it
# selects, returned as a numeric q x q_G matrix in the order of the parameters.
.contrast_matrix <- function(L, par_names) {
  D <- length(par_names)
  if (is.character(L)) {
    idx <- match(L, par_names)
    if (anyNA(idx))
      stop("L names parameters the fit does not have: ",
           paste(L[is.na(idx)], collapse = ", "), ".", call. = FALSE)
    if (anyDuplicated(idx))
      stop("L names a parameter more than once.", call. = FALSE)
    M <- matrix(0, length(idx), D)
    M[cbind(seq_along(idx), idx)] <- 1
    colnames(M) <- par_names
    return(M)
  }
  if (is.null(dim(L))) L <- matrix(L, nrow = 1L)
  L <- as.matrix(L)
  if (!is.numeric(L) || any(!is.finite(L)))
    stop("L must be a finite numeric matrix, a vector with one entry per ",
         "parameter, or a vector of parameter names.", call. = FALSE)
  if (ncol(L) != D)
    stop("L must have one column per graph parameter (", D, "), not ",
         ncol(L), ".", call. = FALSE)
  if (!is.null(colnames(L))) {
    if (!setequal(colnames(L), par_names))
      stop("the column names of L are not the names of the graph parameters.",
           call. = FALSE)
    L <- L[, par_names, drop = FALSE]
  }
  if (!nrow(L)) stop("L has no rows, so there is nothing to test.", call. = FALSE)
  storage.mode(L) <- "double"
  L
}

#' Wald test of a linear hypothesis on the graph layer
#'
#' Tests \eqn{H_0:\mathbf L\boldsymbol\theta=\mathbf c} on the parameters of the
#' graph layer by the generalized Wald statistic,
#' \deqn{W=(\mathbf L\hat{\boldsymbol\theta}-\mathbf c)^\top
#'   \left[\mathbf L\,\widehat{\mathrm{Var}}(\hat{\boldsymbol\theta})\,
#'   \mathbf L^\top\right]^{-1}(\mathbf L\hat{\boldsymbol\theta}-\mathbf c),}
#' referred to a chi-squared distribution with as many degrees of freedom as the
#' rank of \eqn{\mathbf L}. The covariance
#' \eqn{\widehat{\mathrm{Var}}(\hat{\boldsymbol\theta})} is any of the three
#' forms of \code{\link{mc_vcov_graph}}: model-based, sandwich, or sandwich
#' corrected for the first stage.
#'
#' @details
#' Let \eqn{q=\mathrm{rank}(\mathbf L)}. Under \eqn{H_0}, with the number of
#' responses fixed and the number of units \eqn{n\to\infty}, the asymptotic
#' normality of \eqn{\sqrt n(\hat{\boldsymbol\theta}-\boldsymbol\theta_0)}
#' with covariance \eqn{\mathbf J^{-1}\mathbf V\mathbf J^{-\top}} gives
#' \eqn{\sqrt n(\mathbf L\hat{\boldsymbol\theta}-\mathbf c)\xrightarrow{d}
#' \mathcal N(\mathbf 0,\boldsymbol\Xi)}, with
#' \eqn{\boldsymbol\Xi=\mathbf L\mathbf J^{-1}\mathbf V\mathbf J^{-\top}
#' \mathbf L^\top} positive definite when \eqn{\mathbf J} is nonsingular and
#' \eqn{\mathbf V} is positive definite; \eqn{W} then converges to
#' \eqn{\chi^2_q} once \eqn{\boldsymbol\Xi} is replaced by a consistent
#' estimator. For a single edge, \eqn{W=\hat\theta_e^2/
#' \widehat{\mathrm{Var}}(\hat\theta_e)} and \eqn{q=1}. The sensitivity
#' \eqn{\mathbf J} is the expected one in the GGM layer, where it coincides with
#' the curvature, and the observed curvature of the working criterion in the
#' DAG layer. With \code{cluster}, the contributions are summed within
#' independent clusters before the variability is formed, and the asymptotic
#' statement is in the number of clusters.
#'
#' @section Scale:
#' The statistic can be written with a per-unit curvature and variability,
#' \eqn{W=n(\mathbf L\hat{\boldsymbol\theta}-\mathbf c)^\top[\mathbf L
#' \mathbf H^{-1}\mathbf V\mathbf H^{-\top}\mathbf L^\top]^{-1}
#' (\mathbf L\hat{\boldsymbol\theta}-\mathbf c)}, with \eqn{\mathbf H} the
#' observed curvature of the graph block, equal to the expected sensitivity
#' \eqn{\mathbf J} in the GGM layer. The covariance returned by
#' \code{\link{mc_vcov_graph}} is already on the scale of
#' \eqn{\hat{\boldsymbol\theta}}, \eqn{\mathbf H_n^{-1}(n\mathbf V)
#' \mathbf H_n^{-1}=n^{-1}\mathbf H^{-1}\mathbf V\mathbf H^{-1}}, so the two
#' expressions are the same number.
#'
#' @section A contrast that is not of full rank:
#' Redundant rows add no restriction. The statistic is then formed on the range
#' of \eqn{\mathbf L\,\widehat{\mathrm{Var}}\,\mathbf L^\top} and the degrees of
#' freedom are the rank of \eqn{\mathbf L}; a right-hand side that the redundant
#' rows contradict is refused. The statistic is invariant to any invertible
#' recombination of the rows, applied to \eqn{\mathbf L} and \eqn{\mathbf c}
#' alike.
#'
#' @section A contrast without variability:
#' A hypothesis is refused when the covariance of the contrast is singular on
#' its range, judged against the model-based variance of the same contrast.
#' That happens with \code{type = "two_stage"} for a parameter the first stage
#' fixes exactly, such as the innovation variance of a root node of a DAG with
#' Gaussian intercept-only margins, whose corrected variance is zero.
#'
#' @section What the reference distribution requires:
#' The chi-squared reference holds for a hypothesis fixed before looking at the
#' data, under the regularity of the second-order estimating functions and, with
#' \code{type = "model"}, the working identity \eqn{\mathbf V=\mathbf J}.
#' Applied to parameters of a support chosen by a search over graphs, it is
#' descriptive. Hypotheses on the first stage, the regression and dispersion
#' parameters, are outside this function, which acts on the graph layer only;
#' the \pkg{htmcglm} package tests those.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param L the contrast: a numeric matrix with one column per graph parameter,
#'   in the order of \code{coef(fit, type = "graph")} (or with those names as
#'   column names, in any order); a numeric vector, read as one row; or a
#'   character vector of parameter names, which selects them.
#' @param rhs the right-hand side \eqn{\mathbf c}, a single number or one value
#'   per row of \code{L}.
#' @param type the covariance: \code{"sandwich"} (default), \code{"model"} or
#'   \code{"two_stage"}, as in \code{\link{mc_vcov_graph}}.
#' @inheritParams mc_vcov_graph
#' @return An object of class \code{"htest"}: the Wald \code{statistic}
#'   \eqn{W} (named \code{"X-squared"}); \code{parameter}, the degrees of
#'   freedom \code{df}, equal to the rank of \code{L}; the \code{p.value};
#'   the \code{estimate} \eqn{\mathbf L\hat{\boldsymbol\theta}} and the
#'   \code{null.value} \eqn{\mathbf c}, both named by the rows of \code{L};
#'   the \code{method} and \code{data.name} labels; the numeric matrix
#'   \code{contrast} used, with one row per row of \code{L} and one column per
#'   graph parameter; and the covariance \code{type}.
#' @seealso \code{\link{mc_test_independence}} and
#'   \code{\link{mc_test_adequacy}}, which are this test with a particular
#'   contrast, and \code{\link{mc_vcov_graph}} for the covariance.
#' @family inference on the graph layer
#' @references
#' Godambe, V. P. (1960). An optimum property of regular maximum likelihood
#' estimation. \emph{Annals of Mathematical Statistics} 31(4), 1208--1211.
#' \doi{10.1214/aoms/1177705693}
#'
#' de Freitas, L. A. C. and Bonat, W. H. (2022). Hypothesis tests for multiple
#' responses regression models in R: the htmcglm package. arXiv preprint
#' arXiv:2208.01568. \doi{10.48550/arXiv.2208.01568}
#' @examples
#' \donttest{
#' set.seed(1); n <- 300
#' K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
#' Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
#' d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' mc_wald_graph(f, c("eta_13", "eta_23"))          # two edges at once
#' mc_wald_graph(f, c(1, -1, 0, 0, 0, 0))           # eta_11 = eta_22
#' }
#' @export
mc_wald_graph <- function(fit, L, rhs = 0,
                          type = c("sandwich", "model", "two_stage"),
                          cluster = NULL, small_sample = FALSE,
                          method = c("auto", "shortcut", "full")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  type <- match.arg(type)
  method <- match.arg(method)
  .assert_graph_inference(fit, "Wald testing")
  g <- fit$graph
  par <- as.numeric(g$par)
  pn <- names(g$par)
  if (is.null(pn)) pn <- paste0("theta_", seq_along(par))
  L <- .contrast_matrix(L, pn)
  q <- nrow(L)
  if (!is.numeric(rhs) || any(!is.finite(rhs)) ||
      !(length(rhs) %in% c(1L, q)))
    stop("rhs must be one finite number or one per row of L.", call. = FALSE)
  rhs <- rep_len(as.numeric(rhs), q)

  V <- mc_vcov_graph(fit, type = type, cluster = cluster,
                     small_sample = small_sample,
                     method = if (type == "two_stage") method else "auto")
  V <- as.matrix(V)
  est <- as.numeric(L %*% par)
  d <- est - rhs
  M <- L %*% V %*% t(L)
  M <- 0.5 * (M + t(M))
  k <- qr(L)$rank
  # A contrast whose variance vanishes next to its model-based variance is a
  # quantity the covariance takes as fixed, and no statistic is defined for it.
  # Each row is put on the scale of its model-based standard deviation, so the
  # verdict does not change when a row of L is rescaled.
  ev <- eigen(M, symmetric = TRUE)
  keep <- seq_len(k)
  sref <- sqrt(pmax(diag(L %*% as.matrix(g$vcov) %*% t(L)), 0))
  evs <- if (all(is.finite(sref)) && all(sref > 0))
    eigen(M / outer(sref, sref), symmetric = TRUE, only.values = TRUE)$values else
      rep(0, q)
  if (any(!is.finite(evs)) || any(evs[keep] <= 1e-10) ||
      any(ev$values[keep] <= 0))
    stop("the covariance of the contrast is singular on its range: under the ",
         type, " covariance some tested combination has no sampling ",
         "variability, as a parameter the first stage fixes, and it cannot ",
         "be tested.", call. = FALSE)
  if (k == q) {
    W <- as.numeric(t(d) %*% solve(M) %*% d)
  } else {
    P <- ev$vectors[, keep, drop = FALSE]
    off <- d - P %*% crossprod(P, d)
    if (sqrt(sum(off^2)) > 1e-8 * max(1, sqrt(sum(d^2))))
      stop("rhs is not consistent with the redundant rows of L.", call. = FALSE)
    W <- sum(crossprod(P, d)^2 / ev$values[keep])
  }
  structure(list(
    statistic = c("X-squared" = W),
    parameter = c(df = k),
    p.value = stats::pchisq(W, k, lower.tail = FALSE),
    estimate = stats::setNames(est, rownames(L)),
    null.value = stats::setNames(rhs, rownames(L)),
    method = sprintf("Wald test of %d linear restriction%s on the graph layer (%s covariance)",
                     k, if (k == 1L) "" else "s", type),
    data.name = sprintf("%s layer with %d parameters", toupper(g$type),
                        length(par)),
    contrast = L, type = type), class = "htest")
}

#' Adequacy of a support against the saturated graph
#'
#' Tests whether the edges a support omits can all be zero at once. The
#' saturated graph is refitted on the same residuals, complete in the undirected
#' layer and complete along the same ordering in the directed one, and the
#' parameters of the omitted edges are tested jointly by
#' \code{\link{mc_wald_graph}} under the covariance of that saturated fit
#' (an adequacy check of a selected structure). With \eqn{E_{\max}} the set of
#' all pairs of responses and \eqn{\hat E} the edges of the support, the null
#' hypothesis is \eqn{H_0:\theta_e=0} for every
#' \eqn{e\in E_{\max}\setminus\hat E}, against the saturated graph, and the
#' statistic is the Wald statistic of \code{\link{mc_wald_graph}} with the
#' contrast that selects those \eqn{|E_{\max}|-|\hat E|} parameters.
#'
#' @section A diagnostic, not a test at nominal level:
#' If the omitted edges were fixed before looking at the data, the statistic
#' would be referred to a chi-squared distribution with as many degrees of
#' freedom as edges omitted. When the support was chosen by a search on the same
#' data, that reference is only the one of a fixed support: the result is
#' descriptive, with no claim of nominal level after selection. A large p-value
#' says that the zeros of the support are compatible with the residual
#' covariance; a small one says that at least one omitted edge carries partial
#' association. The returned object records this in
#' \code{descriptive_after_selection}.
#'
#' @param fit an object of class \code{"mcglm_graph"} whose support omits at
#'   least one edge; the saturated graph is refitted from its residuals and
#'   must converge.
#' @inheritParams mc_wald_graph
#' @return An object of class \code{"htest"} as in \code{\link{mc_wald_graph}},
#'   whose degrees of freedom are the number of omitted edges, with the extra
#'   components \code{omitted} (a two-column matrix \code{i}, \code{j} of the
#'   omitted pairs of responses), \code{n_saturated} (the number of edges of
#'   the saturated graph) and \code{descriptive_after_selection}, always
#'   \code{TRUE}.
#' @examples
#' \donttest{
#' set.seed(1); n <- 400
#' K <- diag(4) * 2; K[1, 2] <- K[2, 1] <- -0.7; K[3, 4] <- K[4, 3] <- -0.7
#' Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K))
#' d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
#'                 x = rnorm(n))
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
#' sel <- mc_graph_stepwise(f, criterion = "bic", trace = FALSE)
#' mc_test_adequacy(sel$fit)
#' }
#' @seealso \code{\link{mc_wald_graph}}, \code{\link{mc_test_independence}}
#'   and \code{\link{mc_graph_stepwise}}.
#' @family inference on the graph layer
#' @export
mc_test_adequacy <- function(fit, type = c("sandwich", "model", "two_stage"),
                             cluster = NULL, small_sample = FALSE,
                             method = c("auto", "shortcut", "full")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  type <- match.arg(type)
  method <- match.arg(method)
  g <- fit$graph
  R <- g$R
  if (identical(g$type, "ggm")) {
    full <- mc_between_graph("ggm")
    have <- .validate_edges(g$spec$edges, R)
    all_pairs <- .validate_edges(NULL, R)
  } else {
    ordering <- g$spec$ordering
    if (is.null(ordering)) ordering <- seq_len(R)
    full <- mc_between_graph("dag", ordering = ordering)
    have <- .validate_edges(g$spec$edges, R)
    all_pairs <- .validate_edges(NULL, R)
  }
  key <- function(e) paste(e[, 1L], e[, 2L], sep = ":")
  omitted <- all_pairs[!key(all_pairs) %in% key(have), , drop = FALSE]
  if (!nrow(omitted))
    stop("the support is already the saturated graph, so no edge is omitted ",
         "and there is nothing to test.", call. = FALSE)

  sat <- .refit_between(fit, full)
  if (!isTRUE(sat$graph$converged))
    stop("the saturated graph did not converge on these residuals (",
         sat$graph$status, "), so the omitted edges cannot be tested.",
         call. = FALSE)
  # positions of the omitted edges among the parameters of the saturated fit,
  # located by the edges themselves rather than by parameter names
  if (identical(g$type, "ggm")) {
    pos <- R + match(key(omitted), key(mc_precision_bases(R)$edges))
  } else {
    eo <- .dag_setup(full, R)$edges_original      # (parent, child), gamma order
    un <- cbind(pmin(eo[, 1L], eo[, 2L]), pmax(eo[, 1L], eo[, 2L]))
    pos <- match(key(omitted), key(un))
  }
  if (anyNA(pos))
    stop("the omitted edges could not be located among the parameters of the ",
         "saturated fit.", call. = FALSE)
  Lsel <- matrix(0, length(pos), length(sat$graph$par))
  Lsel[cbind(seq_along(pos), pos)] <- 1
  w <- mc_wald_graph(sat, Lsel, 0, type = type, cluster = cluster,
                     small_sample = small_sample, method = method)
  colnames(omitted) <- c("i", "j")
  w$method <- sprintf(paste0("Wald test of adequacy against the saturated graph ",
                             "(%s covariance; descriptive after selection)"),
                      type)
  w$data.name <- sprintf("%d of %d edges omitted from a %s layer",
                         nrow(omitted), nrow(all_pairs), toupper(g$type))
  w$omitted <- omitted
  w$n_saturated <- nrow(all_pairs)
  w$descriptive_after_selection <- TRUE
  w
}
