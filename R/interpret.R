# interpret.R: interpretation and visualization of the estimated structure.
# Part of the mcglmgraph package. License: GPL-3.

#' Partial correlations between responses
#'
#' Computes the between-response partial correlation matrix
#' \eqn{\rho_{ij\mid\mathrm{rest}} = -K_{ij}/\sqrt{K_{ii} K_{jj}}} from the
#' estimated precision matrix \eqn{K_b = \Sigma_b^{-1}} of the graph layer
#' (Lauritzen, 1996). It is invariant to the scale of the Pearson residuals,
#' and it is computed from \code{fit$graph$K_b} for both the GGM and the DAG
#' layer. A zero entry signals the absence of second-order conditional
#' association between two responses given the others, which coincides with
#' conditional independence only under joint normality. In Gaussian-copula
#' models, that interpretation applies only when the precision is estimated on
#' the latent Gaussian scale.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @return An \eqn{R \times R} matrix of partial correlations with unit
#'   diagonal, with the response names as dimnames.
#' @references
#' Lauritzen, S. L. (1996). \emph{Graphical Models}. Oxford University Press.
#' @seealso \code{\link{mc_partial_cor_delta}} for standard errors and
#'   intervals, \code{\link{mc_get_graph}}, \code{\link{mc_plot_network}}.
#' @family graph interpretation
#' @examples
#' set.seed(1); n <- 200; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_partial_cor(f)
#' @export
mc_partial_cor <- function(fit) {
  stopifnot(inherits(fit, "mcglm_graph"))
  P <- cpp_partial_cor(as.matrix(fit$graph$K_b))
  dimnames(P) <- list(fit$graph$response_names, fit$graph$response_names)
  P
}

#' Extract the estimated dependence graph
#'
#' Thresholds the fitted graph layer into an adjacency matrix and an edge list.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param threshold an edge is kept when the absolute value of its statistic
#'   exceeds \code{threshold}. The statistic depends on the layer: for a GGM fit
#'   it is the partial correlation of \code{\link{mc_partial_cor}}; for a DAG
#'   fit it is the raw structural coefficient \eqn{\gamma}, which is on the
#'   scale of the Pearson residuals and is not a correlation. Only edges of the
#'   fitted support can be kept.
#' @return A list with the \code{adjacency} matrix (\eqn{R \times R}, 0/1
#'   integer, with the response names as dimnames), the \code{edges} matrix and
#'   \code{directed}. For GGM fits the adjacency is symmetric and
#'   \code{edges} has one row \code{(i, j)}, \eqn{i < j}, per retained edge.
#'   For DAG fits the adjacency has a one in position \eqn{(i, j)} for an arrow
#'   from \eqn{i} to \eqn{j}, directed from each earlier response to its later
#'   response in the specified ordering, and \code{edges} has the columns
#'   \code{from} and \code{to}. An empty selection gives an edge matrix with
#'   no rows.
#' @seealso \code{\link{mc_partial_cor}}, \code{\link{mc_plot_network}}
#' @family graph interpretation
#' @examples
#' set.seed(1); n <- 200; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_get_graph(f)$adjacency
#' @export
mc_get_graph <- function(fit, threshold = 1e-4) {
  stopifnot(inherits(fit, "mcglm_graph"))
  R <- fit$graph$R
  rn <- fit$graph$response_names
  A <- matrix(0L, R, R, dimnames = list(rn, rn))
  if (fit$graph$type == "dag") {
    dag <- .dag_setup(fit$graph$spec, R)
    gamma <- fit$graph$par[seq_len(nrow(dag$edges_original))]
    keep <- which(abs(gamma) > threshold)
    if (length(keep)) for (e in keep) {
      i <- dag$edges_original[e, 1L]
      j <- dag$edges_original[e, 2L]
      A[i, j] <- 1L
    }
    ed <- if (length(keep)) dag$edges_original[keep, , drop = FALSE] else
      matrix(integer(0), ncol = 2L)
    colnames(ed) <- c("from", "to")
    return(list(adjacency = A, edges = ed, directed = TRUE))
  }
  P <- mc_partial_cor(fit)
  for (i in seq_len(R - 1L)) for (j in (i + 1L):R) {
    if (abs(P[i, j]) > threshold) A[i, j] <- A[j, i] <- 1L
  }
  ed <- which(upper.tri(A) & A == 1L, arr.ind = TRUE)
  colnames(ed) <- c("i", "j")
  list(adjacency = A, edges = ed, directed = FALSE)
}

#' Plot the estimated dependence network
#'
#' Draws the response variables as nodes on a circular layout. GGM edges are
#' weighted by partial correlations; DAG arrows are weighted by the structural
#' coefficients \eqn{\gamma} of the ordered modified Cholesky factorization,
#' and each arrow stops at the border of the node it points to. Edge width
#' encodes magnitude, saturating at an absolute weight of one, and color
#' encodes sign (blue positive, red negative). The graphical parameters are
#' restored on exit.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param threshold an edge is drawn when the absolute value of its weight
#'   exceeds \code{threshold}: the partial correlation for a GGM fit, the
#'   coefficient \eqn{\gamma} for a DAG fit. The coefficient is on the scale of
#'   the Pearson residuals and is not bounded by one, so the same value of
#'   \code{threshold} does not select comparably in the two layers.
#' @param edge_labels logical; if \code{TRUE}, annotate edges with the
#'   rounded weight.
#' @return Invisibly, the \eqn{R \times R} matrix of weights (partial
#'   correlations for the GGM layer; for the DAG layer, \eqn{\gamma} in the
#'   position \eqn{(i, j)} of each arrow from \eqn{i} to \eqn{j}, zero
#'   elsewhere), with the response names as dimnames. The plot is drawn on the
#'   current graphics device.
#' @seealso \code{\link{mc_get_graph}}, \code{\link{mc_partial_cor}}
#' @family graph interpretation
#' @examples
#' set.seed(1); n <- 200; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_plot_network(f)
#' @export
mc_plot_network <- function(fit, threshold = 0.05, edge_labels = FALSE) {
  stopifnot(inherits(fit, "mcglm_graph"))
  R <- fit$graph$R
  dag_mode <- fit$graph$type == "dag"
  if (dag_mode) {
    W <- matrix(0, R, R)
    dag <- .dag_setup(fit$graph$spec, R)
    gamma <- fit$graph$par[seq_len(nrow(dag$edges_original))]
    if (length(gamma)) for (e in seq_along(gamma))
      W[dag$edges_original[e, 1L], dag$edges_original[e, 2L]] <- gamma[e]
  } else W <- mc_partial_cor(fit)
  ang <- seq(0, 2 * pi, length.out = R + 1L)[-(R + 1L)] + pi / 2
  xx <- cos(ang); yy <- sin(ang)
  op <- graphics::par(mar = c(1, 1, 2, 1)); on.exit(graphics::par(op))
  graphics::plot(xx, yy, type = "n", xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
                 axes = FALSE, xlab = "", ylab = "",
                 main = "Between-response dependence network")
  pairs <- if (dag_mode) which(W != 0, arr.ind = TRUE) else
    which(upper.tri(W), arr.ind = TRUE)
  # Radius of a node disk (pch 21 at cex 4.5), in inches: arrows stop at its
  # border, so that the head, which carries the direction, stays visible.
  node_r <- 0.375 * 0.5 * graphics::par("cin")[2L] * 4.5
  if (nrow(pairs)) for (e in seq_len(nrow(pairs))) {
    i <- pairs[e, 1L]; j <- pairs[e, 2L]
    w <- W[i, j]
    if (abs(w) > threshold) {
      draw <- if (dag_mode) graphics::arrows else graphics::segments
      x0 <- xx[i]; y0 <- yy[i]; x1 <- xx[j]; y1 <- yy[j]
      if (dag_mode) {
        en <- .arrow_ends(x0, y0, x1, y1, node_r)
        x0 <- en[1L]; y0 <- en[2L]; x1 <- en[3L]; y1 <- en[4L]
      }
      # the width saturates at |w| = 1, since a raw coefficient is unbounded
      args <- list(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
                   lwd = 1 + 4 * min(abs(w), 1),
                   col = if (w > 0) "steelblue" else "firebrick")
      if (dag_mode) args$length <- 0.1
      do.call(draw, args)
      if (edge_labels) {
        graphics::text((xx[i] + xx[j]) / 2, (yy[i] + yy[j]) / 2,
                       round(w, 2), cex = 0.7)
      }
    }
  }
  graphics::points(xx, yy, pch = 21, bg = "grey90", cex = 4.5)
  graphics::text(xx, yy, fit$graph$response_names, cex = 0.9)
  dimnames(W) <- list(fit$graph$response_names, fit$graph$response_names)
  invisible(W)
}

# A segment shortened by `r` inches at both ends, in user coordinates of the
# open device. Lengths are taken in inches, so the cut does not depend on the
# aspect of the device; a segment shorter than 2r is returned as it is.
.arrow_ends <- function(x0, y0, x1, y1, r) {
  ix <- graphics::grconvertX(c(x0, x1), "user", "inches")
  iy <- graphics::grconvertY(c(y0, y1), "user", "inches")
  len <- sqrt(diff(ix)^2 + diff(iy)^2)
  if (!is.finite(len) || len <= 2 * r) return(c(x0, y0, x1, y1))
  ux <- diff(ix) / len; uy <- diff(iy) / len
  ix <- ix + c(1, -1) * r * ux; iy <- iy + c(1, -1) * r * uy
  bx <- graphics::grconvertX(ix, "inches", "user")
  by <- graphics::grconvertY(iy, "inches", "user")
  c(bx[1L], by[1L], bx[2L], by[2L])
}

#' Delta-method inference for the between-response partial correlations
#'
#' Attaches standard errors and confidence intervals to the partial
#' correlations of \code{\link{mc_partial_cor}} by the delta method. The
#' gradient of \eqn{\rho_{ij\mid\mathrm{rest}} = -K_{ij}/\sqrt{K_{ii}K_{jj}}}
#' with respect to the covariance parameters \eqn{\theta_d} is
#' \deqn{[\mathbf g_{ij}]_d = -\frac{[\dot{\mathbf K}_d]_{ij}}{\sqrt{K_{ii}K_{jj}}}
#'  - \frac{\rho_{ij\mid\mathrm{rest}}}{2}\left\{
#'    \frac{[\dot{\mathbf K}_d]_{ii}}{K_{ii}} +
#'    \frac{[\dot{\mathbf K}_d]_{jj}}{K_{jj}}\right\},}
#' and the variance follows as \eqn{\mathbf g_{ij}^\top
#' \widehat{\mathrm{Var}}(\hat{\boldsymbol\theta})\,\mathbf g_{ij}}
#' with \eqn{\dot{\mathbf K}_d = \partial\mathbf K_b/\partial\theta_d}. In
#' the GGM layer \eqn{\dot{\mathbf K}_d} is the canonical base
#' \eqn{\mathbf Z_d}; the implementation obtains it from
#' \eqn{\dot{\mathbf K}_d = -\mathbf K_b
#' \dot{\boldsymbol\Sigma}_d \mathbf K_b}, which holds for any
#' parametrization, so the same expression also covers the DAG layer, where the
#' parameters are \eqn{(\boldsymbol\gamma,\boldsymbol\delta)}.
#'
#' @section Scope of the inference:
#' With \code{variability = "model"} or \code{"empirical"} the variance refers
#' to the second stage conditionally on the estimated residuals: it propagates
#' the uncertainty of the covariance parameters, not the uncertainty inherited
#' from the first-stage fit, which \code{"two_stage"} adds. With
#' \code{variability = "model"} it also assumes the working identity
#' \eqn{\mathbf V = \mathbf J}; outside the Gaussian fourth-moment regime the
#' empirical option is the appropriate one. Applied to a support chosen by a
#' search over graphs, the interval remains descriptive, since the selection is
#' not incorporated. The function refuses fits whose numerical diagnostics
#' explicitly mark graph-layer inference as uncertified.
#'
#' @section Structural zeros:
#' In the GGM layer, a pair outside the specified precision support is zero by
#' construction: no parameter moves it, the gradient vanishes and the reported
#' standard error is zero up to rounding error, so that the statistic of such a
#' row is 0 or \code{NA} and carries no information. In the DAG layer, structural zeros of the
#' precision are instead determined by the moral graph. In particular, two
#' non-adjacent parents of a common child can have a nonzero precision entry
#' even though no arrow joins them. Therefore \code{present} denotes adjacency
#' in the specified GGM support or in the moral graph of the specified DAG.
#' Rows with \code{present = FALSE} are structural constraints and should not
#' be read as precise estimates of a null association; a \code{TRUE} value says
#' only that the entry is structurally allowed, since parameter cancellation
#' can still make it zero.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param variability \code{"model"} for the model-based covariance
#'   \eqn{\mathbf J_n^{-1}}, \code{"empirical"} for the sandwich form
#'   \eqn{\mathbf H_n^{-1}(n\hat{\mathbf V})\mathbf H_n^{-1}}, with
#'   \eqn{\mathbf H_n} the observed curvature of the graph block, or
#'   \code{"two_stage"} for the sandwich corrected for the estimation of the
#'   residuals. When the correction takes the shortcut, the two empirical forms
#'   give the partial correlations exactly the same variance, since their
#'   gradients annihilate the shortcut's cross sensitivity algebraically; with
#'   the full form they agree to first order under the conditions of the
#'   section "What the correction reaches" of \code{\link{mc_vcov_graph}}.
#'   Whether that common variance is right depends on the same conditions: a
#'   residual correlation that moves with the covariates puts both off, and
#'   only the full form, \code{\link{mc_vcov_graph}(type = "two_stage", method
#'   = "full")}, accounts for it. The model-based form agrees with the empirical
#'   ones only when, in addition, \eqn{\mathbf V=\mathbf J}, which fails
#'   outside the Gaussian fourth-moment regime.
#' @param cluster optional vector with one identifier per observational row,
#'   honored by the two empirical forms (\code{"empirical"} and
#'   \code{"two_stage"}). Contributions are summed within
#'   clusters before the variability is formed, which is what to use when the
#'   unit of dependence is not the row of the data matrix.
#' @param small_sample logical; multiply the clustered variability by
#'   \eqn{G/(G-1)}, with \eqn{G} the number of clusters (default
#'   \code{FALSE}).
#' @param level confidence level of the interval, a single number in
#'   \eqn{(0, 1)} (default 0.95).
#' @param interval \code{"wald"} (default), symmetric on the scale of the
#'   partial correlation, or \code{"fisher"}, symmetric on the scale of
#'   \eqn{\operatorname{atanh}(\rho)}, where the standard error is
#'   \eqn{\mathrm{se}/(1-\rho^2)}, and mapped back by \eqn{\tanh}, so that
#'   the bounds stay inside \eqn{(-1, 1)}. The standard error, the statistic
#'   and the p-value are the same under both.
#' @return A data frame with one row per pair of responses
#'   (\eqn{q = R(R-1)/2} rows) and columns \code{response1},
#'   \code{response2}, \code{present} (logical support indicator),
#'   \code{pcor} (the partial correlation), \code{se} (standard error),
#'   \code{statistic} (Wald statistic \code{pcor / se}, \code{NA} when the
#'   standard error is exactly zero), \code{p_value} (two-sided, normal
#'   reference) and \code{lower}, \code{upper} (interval bounds). The rows
#'   follow the pairs \eqn{(i, j)}, \eqn{i < j}, in lexicographic order. The
#'   attribute \code{"gradient"} holds the matrix of gradients
#'   \eqn{\mathbf g_{ij}}, with one row per pair and one column per graph
#'   parameter, and the attributes \code{"variability"}, \code{"level"} and
#'   \code{"interval"} record the choices used.
#' @examples
#' set.seed(1); n <- 300
#' K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
#' Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
#' d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_partial_cor_delta(f)
#' mc_partial_cor_delta(f, variability = "empirical")
#' @references
#' Lauritzen, S. L. (1996). \emph{Graphical Models}. Oxford University Press.
#' @seealso \code{\link{mc_partial_cor}} for the point estimates,
#'   \code{\link{mc_vcov_graph}} for the covariance of the graph parameters
#'   and \code{\link{mc_test_independence}} for the joint Wald test.
#' @family graph interpretation
#' @export
mc_partial_cor_delta <- function(fit,
                                 variability = c("model", "empirical",
                                                 "two_stage"),
                                 cluster = NULL, small_sample = FALSE,
                                 level = 0.95,
                                 interval = c("wald", "fisher")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  variability <- match.arg(variability)
  interval <- match.arg(interval)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a single number strictly between 0 and 1.")
  .assert_graph_inference(fit, "Partial-correlation delta inference")
  g <- fit$graph
  if (is.null(g$Dlist))
    stop("The fit does not carry the covariance derivatives.")

  # Through the exported accessor, so that the clustering arguments reach the
  # variability here as they do everywhere else.
  V <- mc_vcov_graph(fit, type = switch(variability, model = "model",
                                        empirical = "sandwich",
                                        two_stage = "two_stage"),
                     cluster = cluster, small_sample = small_sample)
  pairs <- t(utils::combn(g$R, 2L))
  storage.mode(pairs) <- "integer"

  out <- cpp_pcor_delta(as.matrix(g$K_b), g$Dlist,
                        pairs, as.matrix(V))
  se <- sqrt(out$variance)
  present <- .pcor_support_matrix(fit)[cbind(pairs[, 1L], pairs[, 2L])] == 1L
  z <- ifelse(se > 0, out$pcor / se, NA_real_)
  half <- stats::qnorm(1 - (1 - level) / 2) * se
  rho <- as.numeric(out$pcor)
  if (interval == "fisher") {
    # Wald on atanh(rho), whose standard error is se / (1 - rho^2), mapped back
    zz <- atanh(rho)
    hz <- half / (1 - rho^2)
    lower <- tanh(zz - hz); upper <- tanh(zz + hz)
    # at the boundary the transform is not defined and the interval is rho
    edge <- !is.finite(zz)
    lower[edge] <- upper[edge] <- rho[edge]
  } else {
    lower <- rho - half; upper <- rho + half
  }

  res <- data.frame(
    response1 = g$response_names[pairs[, 1L]],
    response2 = g$response_names[pairs[, 2L]],
    present = present,
    pcor = as.numeric(out$pcor),
    se = se,
    statistic = z,
    p_value = 2 * stats::pnorm(-abs(z)),
    lower = lower,
    upper = upper,
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(res, "gradient") <- out$gradient
  attr(res, "variability") <- variability
  attr(res, "level") <- level
  attr(res, "interval") <- interval
  res
}

# Structural support of the precision matrix used by mc_partial_cor_delta().
# A GGM encodes that support directly. For an ordered DAG, K = T' D^{-1} T can
# additionally connect co-parents, so its structural support is the moral graph,
# not merely the directed skeleton.
.pcor_support_matrix <- function(fit) {
  if (fit$graph$type != "dag") return(.present_matrix(fit))

  R <- fit$graph$R
  dag <- .dag_setup(fit$graph$spec, R)
  edges <- dag$edges_original
  support <- matrix(0L, R, R)
  if (nrow(edges) == 0L) return(support)

  for (e in seq_len(nrow(edges))) {
    i <- edges[e, 1L]
    j <- edges[e, 2L]
    support[i, j] <- support[j, i] <- 1L
  }
  for (child in seq_len(R)) {
    parents <- edges[edges[, 2L] == child, 1L]
    if (length(parents) < 2L) next
    co_parents <- t(utils::combn(parents, 2L))
    support[co_parents] <- 1L
    support[cbind(co_parents[, 2L], co_parents[, 1L])] <- 1L
  }
  support
}
