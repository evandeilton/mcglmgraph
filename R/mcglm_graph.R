# mcglm_graph.R: main user interface.
# Part of the mcglmgraph package. License: GPL-3.

#' Fit a graph-structured multivariate covariance generalized linear model
#'
#' Fits a multivariate covariance generalized linear model (MCGLM) and then
#' estimates a graph-structured covariance for its Pearson residuals, in a
#' two-stage residual analysis. The between-response covariance is
#' parametrized either by a Gaussian graphical model (GGM, linear in the
#' precision matrix) or by a directed acyclic graph (DAG, modified Cholesky
#' factorization).
#'
#' @details
#' \strong{First stage.} The mean and within-response dispersion parameters
#' are estimated by \code{\link[mcglm]{mcglm}} with the quasi-score and
#' Pearson estimating functions of Bonat and Jorgensen (2016)
#' \doi{10.1111/rssc.12145}, under the usual correlation-based
#' between-response structure of that package. The \eqn{n \times R} matrix
#' \eqn{E} of Pearson residuals, with rows \eqn{r_i^\top}, is then formed, and
#' \deqn{\hat S_e = \frac{1}{n}\sum_{i=1}^{n}(r_i - \bar r)(r_i - \bar r)^\top}
#' is its empirical covariance, centered at the column means \eqn{\bar r} and
#' with divisor \eqn{n}.
#'
#' \strong{Second stage.} The between-response covariance \eqn{\Sigma_b}
#' with precision \eqn{K_b = \Sigma_b^{-1}} is estimated from \eqn{\hat S_e}
#' by the second-order Pearson estimating function
#' \deqn{\psi_d(\theta) = -\frac{n}{2}\left\{
#'   \mathrm{tr}(\Sigma_b^{-1} \dot\Sigma_{b,d}) -
#'   \mathrm{tr}(\Sigma_b^{-1} \dot\Sigma_{b,d} \Sigma_b^{-1} \hat S_e)
#'   \right\},\qquad
#'   \dot\Sigma_{b,d} = \partial\Sigma_b/\partial\theta_d,}
#' whose sensitivity matrix, in the accumulated scale, has entries
#' \deqn{[J_n]_{da} = \frac{n}{2}\,\mathrm{tr}(\Sigma_b^{-1} \dot\Sigma_{b,d}
#'   \Sigma_b^{-1} \dot\Sigma_{b,a}).}
#' The vector of graph parameters is \eqn{\theta = \eta} for the GGM layer,
#' with \eqn{K_b = \sum_d \eta_d Z_d} (see \code{\link{mc_precision_bases}}),
#' and \eqn{\theta = (\gamma, \delta)} for the DAG layer, with
#' \eqn{\Sigma_b = T^{-1} D T^{-\top}} (see \code{\link{mc_between_graph}}).
#' The order of the parameters within \eqn{\theta} is given in
#' \code{\link{mc_build_sigma_between_graph}}. Both layers keep \eqn{R} free
#' scale parameters, the diagonal of \eqn{K_b} or the innovation variances
#' \eqn{\delta}: \eqn{\Sigma_b} is the covariance of the Pearson residuals and
#' is not normalized to a correlation matrix. Partial correlations and the
#' support of the precision are unaffected by that scale.
#'
#' In the GGM layer, \eqn{\dot\Sigma_{b,d} = -\Sigma_b Z_d \Sigma_b} reduces
#' the two expressions to
#' \deqn{\psi_{\eta_d} = \frac{n}{2}\,\mathrm{tr}\{Z_d(\Sigma_b - \hat S_e)\},
#'   \qquad [J_n]_{da} = \frac{n}{2}\,\mathrm{tr}(Z_d \Sigma_b Z_a \Sigma_b),}
#' so that the root matches \eqn{\Sigma_b} to \eqn{\hat S_e} on the diagonal
#' and on the edges of the graph, which are the covariance selection equations
#' of Dempster (1972). The equations are solved by Newton iterations
#' \eqn{\eta \leftarrow \eta + \alpha J_n^{-1}\psi}, started at the
#' corresponding entries of \eqn{\hat S_e^{-1}}, with the step \eqn{\alpha}
#' halved from 1, at most 30 times, until the smallest eigenvalue of \eqn{K_b}
#' exceeds \eqn{10^{-10}} (see \code{\link{mc_spd_line_search}}); when no such
#' step is found the iteration stops at the current iterate with status
#' \code{"spd_line_search_failed"}. The DAG layer has a closed-form root, node
#' by node (see \code{\link{mc_dag_closed_root}}), and can also be solved by
#' Newton iterations on \eqn{(\gamma, \log\delta)}.
#'
#' A fit is declared converged when three conditions hold: the relative change
#' of the parameters in the last iteration falls below \code{tol}, the largest
#' absolute coordinate of the estimating function, divided by \eqn{n}, is at
#' most \code{score_tol}, and the reciprocal condition number of the sensitivity
#' matrix exceeds \code{rcond_tol}. The closed-form DAG solution takes no step,
#' so there the first condition is replaced by the requirements that no
#' innovation variance is floored and no parent block is ill-conditioned.
#'
#' The graph layer is stored in \code{$graph} and does not replace the
#' covariance components of the first-stage \code{"mcglm"} object.
#' Consequently, inherited \pkg{mcglm} methods such as \code{vcov},
#' \code{gof}, \code{pAIC}, \code{pBIC} and \code{plogLik} refer to the
#' first-stage fit; graph-specific results are exposed through the
#' \code{mc_*} helpers. Because the graph layer is estimated from fitted
#' residuals, its default (model-based) covariance does not account for the
#' estimation of the first stage; see \code{\link{mc_vcov_graph}} for the
#' corrected forms.
#'
#' @param linear_pred a list of model formulas, one per response; at least two
#'   responses are required.
#' @param between_graph a specification from \code{\link{mc_between_graph}}
#'   (default: the GGM layer on the complete graph).
#' @param data a data frame containing the variables in \code{linear_pred}. It
#'   must hold complete cases: a missing or non-finite response value, a
#'   missing covariate, or a non-finite \code{offset}, \code{Ntrial} or
#'   \code{weights} entry passed through \code{...}, is refused with an error
#'   that names the response and the number of rows affected, because the
#'   graph layer needs one complete vector of Pearson residuals per unit.
#'   Removing or imputing incomplete rows is left to the user.
#' @param matrix_pred optional list of within-response covariance bases, as in
#'   \code{\link[mcglm]{mcglm}}, one element per response. When \code{NULL}, an
#'   identity structure (\code{\link[mcglm]{mc_id}}) is used for every
#'   response.
#' @param control a named list of control parameters for the between-response
#'   engine; unknown names are an error. The components are:
#'   \describe{
#'     \item{\code{max_iter}}{maximum number of Newton iterations
#'       (default 100).}
#'     \item{\code{tol}}{tolerance on the relative parameter change
#'       \eqn{\max_d |\theta_d^{(t+1)} - \theta_d^{(t)}| /
#'       (1 + |\theta_d^{(t)}|)} between consecutive iterations, computed on
#'       \eqn{\eta} or on \eqn{(\gamma, \delta)} (default \code{1e-8}).}
#'     \item{\code{score_tol}}{bound on the maximum absolute value of the
#'       estimating function divided by \eqn{n} at the solution (default
#'       \code{sqrt(tol)}).}
#'     \item{\code{rcond_tol}}{smallest admissible reciprocal condition number
#'       of the sensitivity matrix; below it a ridge is used in linear solves
#'       and the fit is not certified (default
#'       \code{.Machine$double.eps^(2/3)}).}
#'     \item{\code{verbose}}{logical; report the iterations (default
#'       \code{FALSE}).}
#'     \item{\code{dag_solver}}{how the DAG layer is solved: \code{"closed"}
#'       (default) takes the node-by-node root of
#'       \code{\link{mc_dag_closed_root}}, and \code{"newton"} runs the scoring
#'       iteration shared with the GGM layer. The two solve the same equations,
#'       the first exactly and without iterating, and agree up to the
#'       convergence tolerance of the iteration; the closed route is much
#'       faster when the number of responses is large. It has no effect on the
#'       GGM layer, which has no closed root.}
#'   }
#' @return An object of class \code{c("mcglm_graph", "mcglm")}. In addition to
#'   all the standard \code{mcglm} components, it carries a list \code{$graph}
#'   with the components
#'   \itemize{
#'     \item \code{type}, \code{spec}: the layer (\code{"ggm"} or
#'       \code{"dag"}) and the canonical \code{\link{mc_between_graph}}
#'       specification, with the number of responses recorded;
#'     \item \code{par}, \code{vcov}, \code{J}, \code{npar}, \code{edge_pos}:
#'       the named numeric vector of graph parameters (\code{eta_kk} then
#'       \code{eta_ij} for the GGM layer; \code{gamma_ji}, the coefficient of
#'       response \eqn{i} in the equation of response \eqn{j}, then
#'       \code{delta_j} for the DAG layer), their model-based covariance
#'       \eqn{J_n^{-1}}, the accumulated sensitivity matrix \eqn{J_n}, the
#'       number of parameters, and the integer positions of the edge parameters
#'       in \code{par}. The matrix \code{vcov} carries the attributes
#'       \code{"inversion_status"} and \code{"inversion_ridge"};
#'     \item \code{Sigma_b}, \code{K_b}: the fitted \eqn{R \times R}
#'       between-response covariance and precision, in the order of the
#'       responses;
#'     \item \code{Dlist}: the list of \code{npar} matrices
#'       \eqn{\dot\Sigma_{b,d}}, the derivatives of \code{Sigma_b} with respect
#'       to each element of \code{par};
#'     \item \code{Shat}, \code{residuals}, \code{n}, \code{R},
#'       \code{response_names}: the empirical covariance of the Pearson
#'       residuals \eqn{\hat S_e}, the \eqn{n \times R} matrix of residuals
#'       (not centered), the number of units and responses, and the response
#'       labels;
#'     \item \code{iter}, \code{converged}, \code{status}, \code{solver}: the
#'       number of iterations (zero for the closed-form DAG solution), the
#'       logical convergence indicator, a status string (\code{"converged"} or
#'       \code{"closed_form"} for a certified fit, otherwise the reason why it
#'       is not, such as \code{"max_iter"}, \code{"spd_line_search_failed"},
#'       \code{"linear_solve_failed"}, \code{"score_not_converged"} or
#'       \code{"sensitivity_ill_conditioned"}), and the solver used
#'       (\code{"newton"} or \code{"closed"});
#'     \item \code{score_norm}, \code{score_norm_unit}, \code{rcond_J},
#'       \code{score_tol}, \code{rcond_tol}, \code{convergence_gates},
#'       \code{stopping_rule}: the largest absolute coordinate of the
#'       estimating function at the solution (accumulated and divided by
#'       \eqn{n}), the reciprocal condition number of \eqn{J_n}, the two
#'       thresholds, the named logical vector of the three conditions
#'       (\code{parameter_change}, \code{score}, \code{sensitivity}) and a
#'       character description of the rule, so that a fit stopped by the
#'       parameter-change rule can be checked against the size of the
#'       estimating function it leaves;
#'     \item \code{inversion_status}, \code{inversion_ridge},
#'       \code{linear_solve_status}, \code{linear_solve_ridge},
#'       \code{inference_ok}: whether the inversion of \eqn{J_n} and the linear
#'       solves of the iterations were exact (\code{"exact"}, or
#'       \code{"exact_blockwise"} for the closed-form DAG solution) or
#'       \code{"ridge_regularized"}, the ridge added in each (zero when
#'       exact), and a
#'       logical that is \code{TRUE} when the inversion was exact and the three
#'       convergence conditions hold, in which case graph-layer inference is
#'       certified;
#'     \item \code{meta}, \code{init_info}, \code{control}: a list of builder
#'       metadata (layer, numbers of responses, edges and parameters, and the
#'       bases of the GGM layer or the ordering and edge matrices of the DAG
#'       layer), the diagnostics of the starting value (a list for the DAG
#'       layer solved by Newton iterations, \code{NULL} otherwise) and the
#'       control list used.
#'   }
#'   A graph layer that did not converge is returned with a warning.
#' @references
#' Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance generalized
#' linear models. \emph{Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)}, 65(5), 649-675. \doi{10.1111/rssc.12145}
#'
#' Dempster, A. P. (1972). Covariance selection. \emph{Biometrics}, 28(1),
#' 157-175. \doi{10.2307/2528966}
#'
#' Lauritzen, S. L. (1996). \emph{Graphical Models}. Oxford University Press.
#' @seealso \code{\link{mc_between_graph}} to specify the layer,
#'   \code{\link{mc_partial_cor}}, \code{\link{mc_partial_cor_delta}} and
#'   \code{\link{mc_get_graph}} to interpret it,
#'   \code{\link{mc_vcov_graph}} for the covariance of the graph parameters,
#'   and \code{\link{mc_build_sigma_between_graph}} for the map from the
#'   parameters to \eqn{\Sigma_b}.
#' @examples
#' set.seed(1); n <- 200; x <- rnorm(n)
#' y1 <- 1 + 0.5 * x + rnorm(n)
#' y2 <- 0.4 * y1 + 0.3 * x + rnorm(n)
#' y3 <- 0.2 * y2 + rnorm(n)
#' dat <- data.frame(y1, y2, y3, x)
#' fit <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
#'                    between_graph = mc_between_graph("ggm"), data = dat)
#' fit
#' coef(fit)            # inherited mcglm method
#' mc_partial_cor(fit)  # graph layer
#' @export
mcglm_graph <- function(linear_pred, between_graph = mc_between_graph("ggm"),
                        data, matrix_pred = NULL, control = list(), ...) {
  stopifnot(is.list(linear_pred), inherits(between_graph, "mcg_graph_spec"))
  if (!requireNamespace("mcglm", quietly = TRUE))
    stop("Package 'mcglm' is required. Please install it.")
  R <- length(linear_pred)
  if (R < 2L) stop("mcglm_graph requires at least two responses.")
  known_control <- c("max_iter", "tol", "score_tol", "rcond_tol", "verbose",
                     "dag_solver")
  unknown_control <- setdiff(names(control), known_control)
  if (length(unknown_control))
    stop("Unknown graph control parameter(s): ",
         paste(unknown_control, collapse = ", "), ".")
  ctrl <- utils::modifyList(
    list(max_iter = 100L, tol = 1e-8, score_tol = NULL,
         rcond_tol = .Machine$double.eps^(2 / 3), verbose = FALSE,
         dag_solver = "closed"),
    control
  )
  if (is.null(ctrl$score_tol)) ctrl$score_tol <- sqrt(ctrl$tol)
  ctrl$dag_solver <- match.arg(ctrl$dag_solver, c("closed", "newton"))

  # The graph layer needs one complete residual vector per unit, so missing
  # values are refused here rather than silently changing n.
  .check_complete_cases(linear_pred, data, list(...))

  # within-response identity structure when not supplied
  if (is.null(matrix_pred)) {
    Z0 <- mcglm::mc_id(data)
    matrix_pred <- rep(list(Z0), R)
  }

  # mean + dispersion via the genuine mcglm engine
  m <- mcglm::mcglm(linear_pred = linear_pred, matrix_pred = matrix_pred,
                    data = data, ...)

  # deparse1() gives one label per response also for a transformed left side
  resp_names <- vapply(linear_pred, function(f) deparse1(f[[2L]]),
                       character(1))
  E <- as.matrix(stats::residuals(m, type = "pearson")) # n x R Pearson residuals
  if (ncol(E) != R) E <- matrix(E, ncol = R)
  colnames(E) <- resp_names
  # Non-finite residuals would make Shat and every score undefined.
  bad <- colSums(!is.finite(E))
  if (any(bad > 0L))
    stop("the Pearson residuals are not finite for ",
         paste0(resp_names[bad > 0L], " (", bad[bad > 0L], " row",
                ifelse(bad[bad > 0L] == 1L, "", "s"), ")", collapse = ", "),
         "; the graph layer needs complete, finite residuals for every unit.",
         call. = FALSE)
  n <- nrow(E)

  # between-response covariance under the graph layer
  Ec <- scale(E, center = TRUE, scale = FALSE)
  Shat <- crossprod(Ec) / n
  dimnames(Shat) <- list(resp_names, resp_names)

  between_graph <- .canonical_spec(between_graph, R)
  bfit <- .fit_between(Shat, n, between_graph,
                       max_iter = ctrl$max_iter, tol = ctrl$tol,
                       score_tol = ctrl$score_tol,
                       rcond_tol = ctrl$rcond_tol,
                       verbose = ctrl$verbose,
                       dag_solver = ctrl$dag_solver)
  # attach the graph layer as a self-contained sub-object (no field collisions)
  m$graph <- .graph_block(bfit, between_graph, Shat, E, n, resp_names, ctrl)
  class(m) <- c("mcglm_graph", class(m))
  # A fit that is not a certified root is returned with its status, and said so.
  if (!isTRUE(bfit$converged))
    warning("the graph layer did not converge (status: ", bfit$status, "). ",
            if (n <= R) paste0("There are ", n, " units for ", R,
                               " responses, so the residual covariance is ",
                               "singular. ") else "",
            "Its estimates are not a root of the estimating equations; see ",
            "fit$graph$convergence_gates.", call. = FALSE)
  m
}

# The `$graph` block of a fit, from the result of .fit_between(). One builder
# for every route, so that no field carries the diagnostic of another fit.
.graph_block <- function(bfit, spec, Shat, residuals, n, response_names,
                         control) {
  dimnames(bfit$Sigma_b) <- dimnames(bfit$K_b) <-
    list(response_names, response_names)
  list(
    type           = spec$type,
    par            = bfit$par,
    Sigma_b        = bfit$Sigma_b,
    K_b            = bfit$K_b,
    vcov           = bfit$vcov,
    J              = bfit$J,
    Shat           = Shat,
    residuals      = residuals,
    edge_pos       = bfit$edge_pos,
    npar           = bfit$npar,
    iter           = bfit$iter,
    converged      = bfit$converged,
    status         = bfit$status,
    score_norm     = bfit$score_norm,
    score_norm_unit = bfit$score_norm_unit,
    rcond_J        = bfit$rcond_J,
    score_tol      = bfit$score_tol,
    rcond_tol      = bfit$rcond_tol,
    convergence_gates = bfit$convergence_gates,
    stopping_rule  = bfit$stopping_rule,
    solver         = bfit$solver,
    inversion_status = bfit$inversion_status,
    inversion_ridge = bfit$inversion_ridge,
    inference_ok   = bfit$inference_ok,
    linear_solve_status = bfit$linear_solve_status,
    linear_solve_ridge = bfit$linear_solve_ridge,
    meta           = bfit$meta,
    init_info      = bfit$init_info,
    spec           = spec,
    Dlist          = bfit$Dlist,
    R              = length(response_names),
    n              = n,
    response_names = response_names,
    control        = control
  )
}

# Refuse missing or non-finite responses, missing covariates and non-finite
# offset, Ntrial or weights passed through `...`, naming the response and the
# number of rows affected.
.check_complete_cases <- function(linear_pred, data, dots = list()) {
  msgs <- character(0)
  resp <- character(length(linear_pred))
  for (k in seq_along(linear_pred)) {
    f <- linear_pred[[k]]
    nm <- deparse(f[[2L]])
    resp[k] <- nm
    # A frame that cannot be built here is left to mcglm; the test of finite
    # Pearson residuals after the first stage still guards the graph layer.
    mf <- tryCatch(stats::model.frame(f, data = data, na.action = stats::na.pass),
                   error = function(e) NULL)
    if (is.null(mf)) next
    y <- as.matrix(stats::model.response(mf))
    ok_y <- if (is.numeric(y)) is.finite(y) else !is.na(y)
    n_y <- sum(rowSums(!ok_y) > 0L)
    if (n_y > 0L)
      msgs <- c(msgs, sprintf(
        "response %s has %d row%s with a missing or non-finite value",
        nm, n_y, if (n_y == 1L) "" else "s"))
    if (ncol(mf) > 1L) {
      n_x <- sum(!stats::complete.cases(mf[, -1L, drop = FALSE]))
      if (n_x > 0L)
        msgs <- c(msgs, sprintf("the covariates of %s have %d row%s with a missing value",
                                nm, n_x, if (n_x == 1L) "" else "s"))
    }
  }
  # mcglm takes these as lists with one numeric vector (or NULL) per response.
  for (arg in c("offset", "Ntrial", "weights")) {
    v <- dots[[arg]]
    if (is.null(v)) next
    if (!is.list(v)) v <- list(v)
    for (k in seq_along(v)) {
      x <- v[[k]]
      if (is.null(x) || !is.numeric(x)) next
      n_bad <- sum(!is.finite(x))
      if (n_bad > 0L)
        msgs <- c(msgs, sprintf("%s of %s has %d missing or non-finite value%s",
                                arg, if (length(v) == length(resp)) resp[k] else
                                  paste("element", k),
                                n_bad, if (n_bad == 1L) "" else "s"))
    }
  }
  if (length(msgs))
    stop(paste(msgs, collapse = "; "), ". mcglm_graph() needs complete ",
         "cases: remove or impute the incomplete rows before fitting.",
         call. = FALSE)
  invisible(TRUE)
}

# Refit only the between-response graph for a new spec (mean unchanged).
# Reuses the cached residuals/Shat. The mcglm base components are deliberately
# left intact, so inherited methods continue to describe the first-stage fit.
.refit_between <- function(fit, spec) {
  g <- fit$graph
  spec <- .canonical_spec(spec, g$R)
  ctrl <- g$control
  if (is.null(ctrl))
    ctrl <- list(max_iter = 100L, tol = 1e-8, score_tol = 1e-4,
                 rcond_tol = .Machine$double.eps^(2 / 3), verbose = FALSE,
                 dag_solver = "closed")
  if (is.null(ctrl$dag_solver)) ctrl$dag_solver <- "closed"
  bfit <- .fit_between(
    g$Shat, g$n, spec, max_iter = ctrl$max_iter, tol = ctrl$tol,
    score_tol = ctrl$score_tol, rcond_tol = ctrl$rcond_tol,
    verbose = ctrl$verbose, dag_solver = ctrl$dag_solver
  )
  # Every fitted field is rebuilt; the residuals, Shat and n are the cached ones.
  fit$graph <- .graph_block(bfit, spec, g$Shat, g$residuals, g$n,
                            g$response_names, g$control)
  fit
}

#' @rdname mcglm_graph
#' @section Two spellings of the same function:
#' \code{mcglmgraph()} and \code{mcglm_graph()} are the same object: the second
#' name is bound to the first, not wrapped around it, so there is no forwarding
#' layer and the two cannot drift apart. \code{mcglmgraph()} matches the
#' package name and is the preferred spelling for new code; \code{mcglm_graph()}
#' is kept for existing scripts. Neither is deprecated.
#'
#' The class of the returned object is \code{"mcglm_graph"} under both
#' spellings, and it is the class on which the S3 methods are registered.
#' @export
mcglmgraph <- mcglm_graph
