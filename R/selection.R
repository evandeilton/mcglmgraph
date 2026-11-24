# selection.R: graph structure selection: edge scores, stepwise search,
# L1-penalized path, penalty selection and selection-frequency aggregation.
# Part of the mcglmgraph package. License: GPL-3.

#' Edge score and Wald statistics
#'
#' Computes, for every pair of responses, a chi-squared(1) statistic that ranks
#' edges of the Gaussian graphical model (GGM) layer without refitting the whole
#' model. The precision matrix is parametrized linearly,
#' \eqn{K_b = \sum_d \eta_d Z_d}, with one diagonal parameter per response and
#' one parameter \eqn{\eta_e} per edge, so that an edge is absent exactly when
#' \eqn{\eta_e = 0}. Present edges are scored by the Wald statistic
#' \eqn{\hat\eta_e^2 / \widehat{\mathrm{Var}}(\hat\eta_e)}, with the variance
#' read from the model-based covariance \eqn{J^{-1}} stored in the fit. Absent
#' edges are scored by the \emph{efficient} Pearson score, that is, the score of
#' the candidate coordinate after projecting out the nuisance parameters
#' \eqn{\lambda} of the current model,
#' \deqn{u_e = \psi_e - J_{e\lambda} J_{\lambda\lambda}^{-1}\psi_\lambda,
#'       \qquad S_e = u_e^2 / \widehat{\mathrm{Var}}(u_e),}
#' where \eqn{\psi} is the second-order Pearson estimating function and \eqn{J}
#' its sensitivity. With \code{variability = "model"} the variance is taken
#' under the working identity \eqn{V = J} (second Bartlett identity), in which
#' case it equals the Schur complement
#' \eqn{J_{e\cdot\lambda} = J_{ee} - J_{e\lambda}J_{\lambda\lambda}^{-1}
#' J_{\lambda e}}; with \code{variability = "empirical"} it uses the empirical
#' fourth-moment variability of the Pearson contributions, which is the
#' appropriate choice when the second Bartlett identity fails.
#'
#' Dividing \eqn{\psi_e^2} by the unadjusted \eqn{J_{ee}} ignores the
#' estimation of the remaining parameters and, because
#' \eqn{J_{e\cdot\lambda}\le J_{ee}}, understates the statistic. That quantity
#' is returned as \code{stat_unadjusted}, for ranking, but it does not carry a
#' chi-squared calibration.
#'
#' The reported \code{p_value} is a working reference distribution: it is the
#' chi-squared(1) tail of \code{stat}. Its validity requires the regularity
#' conditions of the second-order estimating functions and, for
#' \code{variability = "model"}, the working identity \eqn{V = J}. Applied to
#' edges chosen by a search over graphs, it remains descriptive.
#'
#' @param fit an object of class \code{"mcglm_graph"} (GGM layer).
#' @param which one of \code{"both"} (default), \code{"wald"} (present edges
#'   only) or \code{"score"} (absent edges only).
#' @param variability \code{"model"} (default, working \eqn{V = J}) or
#'   \code{"empirical"} (fourth-moment variability) for the absent-edge scores.
#'   It does not affect the Wald statistics of present edges.
#' @details The function refuses fits whose numerical diagnostics explicitly
#'   mark graph-layer inference as uncertified. Responses are indexed
#'   \eqn{i<j} over all \eqn{R(R-1)/2} pairs.
#' @return A data frame with one row per pair of responses (restricted by
#'   \code{which}) and columns \code{i}, \code{j} (response indices with
#'   \eqn{i<j}), \code{present} (logical, edge in the support of the fit),
#'   \code{stat} (Wald statistic for present edges, efficient score for absent
#'   ones), \code{p_value} (chi-squared(1) upper tail of \code{stat}),
#'   \code{stat_unadjusted} (unadjusted score \eqn{\psi_e^2/J_{ee}}, only for
#'   absent edges, \code{NA} otherwise) and \code{variance} (estimated variance
#'   of \eqn{\hat\eta_e} for present edges, of \eqn{u_e} for absent ones). The
#'   attribute \code{"variability"} records the choice used. An edge whose
#'   estimated variance is not finite and positive has \code{stat} and
#'   \code{p_value} equal to \code{NA}, and \code{\link{mc_graph_stepwise}}
#'   does not take a present edge in that state as a removal candidate.
#' @seealso \code{\link{mc_graph_stepwise}}, which uses these statistics to
#'   rank candidate edges.
#' @family structure selection
#' @examples
#' \donttest{
#' set.seed(1); n <- 200; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_edge_scores(f)
#' }
#' @export
mc_edge_scores <- function(fit, which = c("both", "wald", "score"),
                           variability = c("model", "empirical")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  if (fit$graph$type != "ggm")
    stop("Edge scores are defined for the GGM layer.")
  which <- match.arg(which)
  variability <- match.arg(variability)
  .assert_graph_inference(fit, "Edge-score inference")
  g <- fit$graph
  R <- g$R; n <- g$n
  Sig <- g$Sigma_b; Shat <- g$Shat
  present <- .present_matrix(fit)
  pairs <- t(utils::combn(R, 2L))
  is_pres <- present[cbind(pairs[, 1L], pairs[, 2L])] == 1L
  stat <- rep(NA_real_, nrow(pairs))
  stat_raw <- rep(NA_real_, nrow(pairs))
  vare <- rep(NA_real_, nrow(pairs))
  # Only the block that is asked for is computed.
  if (any(is_pres) && which != "score") {   # present edges: Wald [cpp_wald_edges]
    pp <- pairs[is_pres, , drop = FALSE]
    pos <- vapply(seq_len(nrow(pp)),
                  function(e) .edge_param_pos(fit, pp[e, 1L], pp[e, 2L]),
                  integer(1))
    stat[is_pres] <- cpp_wald_edges(g$par, g$vcov,
                                    as.integer(pos - 1L))
    vare[is_pres] <- diag(g$vcov)[pos]
  }
  if (any(!is_pres) && which != "wald") {    # absent edges: efficient score
    ap <- matrix(as.integer(pairs[!is_pres, , drop = FALSE]), ncol = 2L)
    eff <- .efficient_edge_scores(fit, ap, variability = variability)
    stat[!is_pres] <- eff$stat
    vare[!is_pres] <- eff$variance
    stat_raw[!is_pres] <- cpp_score_edges(Sig, Shat, n, ap)   # unadjusted J_ee form
  }
  out <- data.frame(i = pairs[, 1L], j = pairs[, 2L], present = is_pres,
                    stat = stat,
                    p_value = stats::pchisq(stat, 1, lower.tail = FALSE),
                    stat_unadjusted = stat_raw, variance = vare)
  attr(out, "variability") <- variability
  switch(which, both = out, wald = out[out$present, , drop = FALSE],
         score = out[!out$present, , drop = FALSE])
}

#' Stepwise selection of the dependence graph
#'
#' Greedy search over supports of the GGM layer guided by an information
#' criterion (see \code{\link{mc_compute_criterion}}). The backward search
#' starts from the complete graph and, at each step, removes the edge whose
#' removal most reduces the criterion; the forward search starts from the empty
#' graph and adds the edge whose inclusion most reduces it; \code{"both"} starts
#' from the complete graph and considers removals and additions at every step.
#' Candidates are ranked by \code{\link{mc_edge_scores}} (smallest Wald
#' statistic for removal, largest efficient score for addition), and only the
#' \code{k_confirm} top-ranked ones are refitted per step, which is the score
#' pruning that keeps the search affordable. The search stops when no
#' candidate lowers the criterion, or after \code{max_steps} steps. Because the
#' mean structure is held fixed, only the between-response covariance is
#' refitted at each step, which is justified by the asymptotic orthogonality of
#' the first-stage and graph-layer estimating functions.
#'
#' The search is defined on the undirected layer: every candidate is a GGM
#' support, and the edge scores that rank the candidates are those of the
#' precision. A fit carrying the directed layer is refused rather than
#' converted, since the conversion would silently replace the structure the
#' caller fitted by a different one.
#'
#' Candidates whose refit fails to converge or yields a non-finite criterion
#' are skipped and recorded in the component \code{candidate_failures}.
#'
#' @param fit an object of class \code{"mcglm_graph"} (GGM layer).
#' @param direction \code{"backward"} (default), \code{"forward"} or
#'   \code{"both"}.
#' @param criterion \code{"bic"} (default), \code{"aic"}, \code{"qic"},
#'   \code{"gbic"} or its alias \code{"clbic"}, or \code{"ebic"}.
#' @param max_steps maximum number of steps; \code{NULL} (default) means the
#'   number of possible edges, \eqn{R(R-1)/2}.
#' @param k_confirm number of top-ranked candidates refitted per step
#'   (default 3).
#' @param trace logical; record the accepted moves in the component
#'   \code{trace}.
#' @param ebic_gamma the parameter of the extended BIC, used only with
#'   \code{criterion = "ebic"}; see \code{\link{mc_compute_criterion}}.
#' @param cluster optional independent-cluster identifier, one per residual
#'   row, passed to \code{\link{mc_compute_criterion}}: Pearson contributions
#'   are summed within clusters in QIC and GBIC, and the number of clusters
#'   replaces the number of rows in the logarithmic penalty of BIC, GBIC and
#'   EBIC.
#' @param small_sample logical; apply the cluster factor \eqn{G/(G-1)}, with
#'   \eqn{G} the number of clusters, to the empirical variability of the
#'   clustered criteria (default \code{FALSE}); ignored without
#'   \code{cluster}.
#' @return An object of class \code{"mcg_selection"}, a list with
#'   \code{method} (\code{"stepwise"}), \code{criterion}, \code{graph} (the
#'   adjacency matrix of the selected graph), \code{edges} (two-column matrix
#'   of selected edges), \code{n_edges}, \code{trace} (data frame with
#'   \code{step}, \code{action}, \code{edge_i}, \code{edge_j} and
#'   \code{criterion} for each accepted move, or \code{NULL}), \code{fit} (the
#'   \code{"mcglm_graph"} refitted on the selected support) and
#'   \code{candidate_failures} (data frame of the candidates that could not be
#'   evaluated). The components \code{lambda}, \code{lambda_path},
#'   \code{cpdag} and \code{vstructures} are \code{NA} or \code{NULL} here.
#' @seealso \code{\link{mc_select_lambda}} for the L1 route,
#'   \code{\link{mc_graph_stability}} for selection frequencies,
#'   \code{\link{mcg_selection-methods}}.
#' @family structure selection
#' @examples
#' \donttest{
#' set.seed(1); n <- 300
#' K <- diag(4) * 2; K[1,2]<-K[2,1]<- -0.6; K[2,3]<-K[3,2]<- -0.6; K[3,4]<-K[4,3]<- -0.6
#' Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K)); x <- rnorm(n)
#' d <- data.frame(y1=Y[,1], y2=Y[,2], y3=Y[,3], y4=Y[,4], x=x)
#' f <- mcglm_graph(list(y1~x,y2~x,y3~x,y4~x), data = d)
#' sel <- mc_graph_stepwise(f, criterion = "bic")
#' sel$n_edges
#' }
#' @export
mc_graph_stepwise <- function(fit, direction = c("backward", "forward", "both"),
                              criterion = c("bic", "aic", "qic", "gbic", "clbic",
                                            "ebic"),
                              max_steps = NULL, k_confirm = 3L, trace = TRUE,
                              ebic_gamma = 0.5, cluster = NULL,
                              small_sample = FALSE) {
  stopifnot(inherits(fit, "mcglm_graph"))
  if (!identical(fit$graph$type, "ggm"))
    stop("mc_graph_stepwise() searches over supports of the GGM layer; this fit ",
         "carries a ", toupper(fit$graph$type), " layer. Refit it with ",
         "mc_between_graph(\"ggm\") to select an undirected support.",
         call. = FALSE)
  direction <- match.arg(direction); criterion <- match.arg(criterion)
  R <- fit$graph$R
  all_edges <- t(utils::combn(R, 2L))
  if (is.null(max_steps)) max_steps <- nrow(all_edges)
  crit <- function(f) unname(mc_compute_criterion(
    f, criterion, ebic_gamma = ebic_gamma, cluster = cluster,
    small_sample = small_sample))

  cur_edges <- if (direction == "forward") all_edges[0, , drop = FALSE] else all_edges
  cur_fit <- .refit_between(fit, mc_between_graph("ggm", .nz_edges(cur_edges)))
  if (!isTRUE(cur_fit$graph$converged))
    stop("The initial stepwise model did not converge: ",
         if (is.null(cur_fit$graph$status)) "status unavailable" else
           cur_fit$graph$status, ".")
  cur_crit <- crit(cur_fit)
  if (length(cur_crit) != 1L || !is.finite(cur_crit))
    stop("The initial stepwise criterion is not finite.")
  tr <- list()
  failures <- list()
  safe_candidate <- function(edges, step, action, edge) {
    err <- NULL
    ff <- tryCatch(
      .refit_between(fit, mc_between_graph("ggm", .nz_edges(edges))),
      error = function(e) {
        err <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(ff) || !isTRUE(ff$graph$converged)) {
      failures[[length(failures) + 1L]] <<- data.frame(
        step = step, action = action, edge_i = edge[1L], edge_j = edge[2L],
        status = if (!is.null(err)) paste0("error: ", err) else
          if (is.null(ff$graph$status)) "non-converged" else
            ff$graph$status, stringsAsFactors = FALSE
      )
      return(NULL)
    }
    cv <- tryCatch(crit(ff), error = function(e) NA_real_)
    if (length(cv) != 1L || !is.finite(cv)) {
      failures[[length(failures) + 1L]] <<- data.frame(
        step = step, action = action, edge_i = edge[1L], edge_j = edge[2L],
        status = "non-finite criterion", stringsAsFactors = FALSE
      )
      return(NULL)
    }
    list(crit = cv, edges = edges, fit = ff, action = action, edge = edge)
  }

  for (step in seq_len(max_steps)) {
    improved <- FALSE; best <- NULL
    if (direction %in% c("backward", "both") && nrow(cur_edges) > 0) {
      es <- mc_edge_scores(cur_fit, "wald")
      # An edge whose Wald statistic is undefined (variance not finite and
      # positive) cannot be ranked as the weakest; it is left out and reported.
      undefined <- !is.finite(es$stat)
      if (any(undefined))
        warning("step ", step, ": the Wald statistic of edge(s) ",
                paste0("(", es$i[undefined], ",", es$j[undefined], ")",
                       collapse = ", "),
                " is undefined (variance not finite and positive); they are not ",
                "candidates for removal.", call. = FALSE)
      es <- es[!undefined, , drop = FALSE]
      es <- es[order(es$stat), , drop = FALSE]
      cand <- utils::head(es, k_confirm)
      for (cc in seq_len(nrow(cand))) {
        e_ij <- c(cand$i[cc], cand$j[cc])
        new_edges <- cur_edges[!(cur_edges[, 1] == e_ij[1] &
                                 cur_edges[, 2] == e_ij[2]), , drop = FALSE]
        candidate <- safe_candidate(new_edges, step, "remove", e_ij)
        if (!is.null(candidate) &&
            (is.null(best) || candidate$crit < best$crit))
          best <- candidate
      }
    }
    if (direction %in% c("forward", "both") && nrow(cur_edges) < nrow(all_edges)) {
      es <- mc_edge_scores(cur_fit, "score")
      es <- es[order(-es$stat), , drop = FALSE]
      cand <- utils::head(es, k_confirm)
      for (cc in seq_len(nrow(cand))) {
        e_ij <- c(cand$i[cc], cand$j[cc])
        # canonical order, so the support does not depend on insertion order
        new_edges <- .validate_edges(rbind(cur_edges, e_ij), R)
        candidate <- safe_candidate(new_edges, step, "add", e_ij)
        if (!is.null(candidate) &&
            (is.null(best) || candidate$crit < best$crit))
          best <- candidate
      }
    }
    if (!is.null(best) && best$crit < cur_crit - 1e-8) {
      cur_edges <- best$edges; cur_fit <- best$fit; cur_crit <- best$crit
      improved <- TRUE
      if (trace) tr[[length(tr) + 1L]] <- data.frame(
        step = step, action = best$action, edge_i = best$edge[1],
        edge_j = best$edge[2], criterion = cur_crit)
    }
    if (!improved) break
  }

  trace_df <- if (length(tr)) do.call(rbind, tr) else NULL
  out <- .mcg_selection("stepwise", criterion, cur_edges, cur_fit,
                        trace = trace_df)
  out$candidate_failures <- if (length(failures)) do.call(rbind, failures) else
    data.frame(step = integer(), action = character(), edge_i = integer(),
               edge_j = integer(), status = character())
  out
}

# The penalty path starts from the diagonal fit eta_kk = 1 / S_kk and takes its
# largest penalty from the KKT residual there, so every residual variance must
# be finite and positive; a response with constant Pearson residuals makes both
# undefined (condition S_kk > 0 for lambda_max).
.check_residual_variance <- function(Shat, what = "the penalty path") {
  v <- diag(Shat)
  bad <- !is.finite(v) | v <= 0
  if (any(bad)) {
    nm <- colnames(Shat)
    if (is.null(nm)) nm <- as.character(seq_along(v))
    stop(what, " needs a finite positive residual variance for every ",
         "response; response(s) ", paste(nm[bad], collapse = ", "),
         " have zero or non-finite residual variance.", call. = FALSE)
  }
  invisible(TRUE)
}

#' L1-penalized estimation of the precision graph
#'
#' Fits the GGM layer with an L1 penalty on the edge parameters,
#' \deqn{\hat{\boldsymbol\eta}^{(\lambda)} = \arg\min_{\boldsymbol\eta}
#'   \Bigl\{-\ell^p_e(\boldsymbol\eta) + \lambda\sum_{e}|\eta_e|\Bigr\},}
#' where \eqn{\ell^p_e} is the Gaussian working log-criterion of the residual
#' layer (see \code{\link{mc_graph_ploglik}}), \eqn{K_b=\sum_d\eta_d Z_d} is
#' the linearly parametrized precision matrix, and the sum runs over the edge
#' parameters only: the diagonal parameters, which set the scale of each
#' response, are not penalized. The penalty is on the accumulated scale of the
#' criterion, so the smallest value that zeroes every edge is
#' \eqn{\lambda_{\max}=n\max_{i<j}|\hat S_{e,ij}|}, with \eqn{\hat S_e} the
#' sample covariance of the Pearson residuals: for
#' \eqn{\lambda\ge\lambda_{\max}} the solution is the empty graph. The problem
#' is solved by the proximal gradient method (iterative
#' shrinkage-thresholding): a gradient step
#' \eqn{\boldsymbol\eta+s\boldsymbol\psi_{\boldsymbol\eta}}, with
#' \eqn{\boldsymbol\psi_{\boldsymbol\eta}=\partial\ell^p_e/\partial\boldsymbol\eta}
#' the Pearson estimating function, followed by soft-thresholding of the edge
#' coordinates at \eqn{s\lambda}. Backtracking halves the
#' step until the precision is positive definite and the sufficient-decrease
#' condition of Beck and Teboulle (2009) holds, under which the penalized
#' objective does not increase, so an admissible matrix alone is not sufficient
#' for accepting an iteration.
#'
#' @param fit an object of class \code{"mcglm_graph"} (GGM layer).
#' @param lambda one finite non-negative penalty parameter, on the accumulated
#'   scale of the criterion (see the objective above); \code{lambda = 0} gives
#'   the unpenalized fit of the complete graph.
#' @param max_iter maximum number of proximal iterations (default 300).
#' @param tol tolerance on the relative change of the parameter vector,
#'   \eqn{\max_d|\eta_d^{+}-\eta_d|/(1+|\eta_d|)} (default \eqn{10^{-6}}).
#' @param precondition logical; must be \code{FALSE}. \code{TRUE} is refused
#'   with an error, since the preconditioned step is not a valid proximal step.
#'   See the section "Preconditioned step".
#' @param kkt_tol tolerance for the per-unit KKT residual required to certify
#'   convergence, in addition to \code{tol} on the relative parameter change
#'   (default \eqn{10^{-4}}).
#' @param start optional starting value of the parameter vector, diagonal
#'   parameters first and then one per edge of the complete graph, as returned
#'   in \code{eta}; it must give a positive definite precision. The path of
#'   \code{\link{mc_select_lambda}} uses it for warm starts.
#' @return A list with \code{eta} (parameter vector, diagonal parameters
#'   first and then one per edge of the complete graph), \code{edges}
#'   (two-column matrix of the pairs whose coefficient is not exactly zero),
#'   \code{n_selected} (their number), \code{Sigma_b} and \code{K_b} (fitted
#'   covariance and precision), \code{iterations}, \code{converged} (logical),
#'   \code{status} (\code{"converged"}, \code{"kkt_not_certified"},
#'   \code{"maximum_iterations"} or \code{"objective_line_search_failed"}),
#'   \code{objective} (final value of the penalized objective),
#'   \code{objective_path} (its value along the iterations), \code{kkt}
#'   (KKT residual, accumulated scale), \code{kkt_unit} (per unit, that is,
#'   divided by \eqn{n}) and \code{kkt_tol} (the tolerance used).
#' @section Convergence certificate:
#' A small relative change of the parameter vector does not by itself certify a
#' stationary point of a non-smooth objective. Convergence is therefore declared
#' only when the relative change is below \code{tol} \emph{and} the per-unit KKT
#' residual is below \code{kkt_tol}: when the first holds and the second does
#' not, the iteration goes on. The optimality conditions are
#' \eqn{\hat\Sigma_{b,kk}=\hat S_{e,kk}} (diagonal),
#' \eqn{n(\hat S_{e,ij}-\hat\Sigma_{b,ij})+\lambda\,\mathrm{sgn}(\hat K_{b,ij})=0}
#' for a nonzero edge and
#' \eqn{|n(\hat S_{e,ij}-\hat\Sigma_{b,ij})|\le\lambda} for a zero edge. The KKT
#' residual \eqn{r_{\mathrm{KKT}}} is their largest violation over the
#' coordinates:
#' \eqn{\tfrac{n}{2}|\hat S_{e,kk}-\hat\Sigma_{b,kk}|} on the diagonal,
#' \eqn{|n(\hat S_{e,ij}-\hat\Sigma_{b,ij})+\lambda\,\mathrm{sgn}(\hat K_{b,ij})|}
#' on a nonzero edge and
#' \eqn{\max\{|n(\hat S_{e,ij}-\hat\Sigma_{b,ij})|-\lambda,\,0\}} on a zero
#' edge. The certificate is \eqn{r_{\mathrm{KKT}}/n\le} \code{kkt_tol}. If
#' \code{max_iter} is reached in that state the status is
#' \code{"kkt_not_certified"}, otherwise \code{"maximum_iterations"}, and in
#' both cases \code{converged} is \code{FALSE}.
#'
#' @section Step size:
#' With \eqn{f=-\ell^p_e} the smooth part of the objective, whose gradient is
#' \eqn{-\boldsymbol\psi_{\boldsymbol\eta}}, every iteration starts from the
#' step \eqn{s=\{\max_d J_{n,dd}+10^{-8}\}^{-1}} and halves it until the
#' proximal point \eqn{\boldsymbol\eta^{+}} gives a positive definite
#' precision and
#' \deqn{f(\boldsymbol\eta^{+})\le f(\boldsymbol\eta)-
#'   \boldsymbol\psi_{\boldsymbol\eta}^\top(\boldsymbol\eta^{+}-\boldsymbol\eta)
#'   +\frac{1}{2s}\|\boldsymbol\eta^{+}-\boldsymbol\eta\|^2 .}
#' The starting step is a scale and not the inverse of a Lipschitz constant;
#' the condition is what guarantees the decrease.
#'
#' @section Starting value:
#' The iteration starts from the diagonal fit that is compatible with the empty
#' graph, \eqn{\eta_{kk}=1/\widehat S_{e,kk}} and every edge at zero, which is
#' also the fit from which the largest penalty of the path is computed. Both
#' need \eqn{\widehat S_{e,kk}>0} for every response; a response whose Pearson
#' residuals have no variance is refused with an error. The
#' problem has a unique solution, so the start changes the number of
#' iterations and not the target; certified solutions from different starts
#' agree up to the tolerance.
#' @section Candidate set:
#' The path is always built over the complete graph on the \eqn{R} responses,
#' whatever support the supplied fit carries: the penalty, not the input
#' support, decides which edges survive. Passing a fit that was restricted to a
#' sparse support therefore does not restrict the candidates, and the returned
#' \code{edges} may include pairs absent from that fit. To keep a support fixed,
#' set it through the \code{edges} argument of \code{\link{mc_between_graph}}
#' and refit, rather than penalizing.
#' @section Preconditioned step:
#' A preconditioned step, which multiplies the estimating function by
#' \eqn{J^{-1}} and then applies coordinatewise soft-thresholding, is not the
#' proximal operator of the L1 penalty in the non-diagonal metric induced by
#' \eqn{J}: it can stop at a point that satisfies the parameter-change test
#' while grossly violating the KKT conditions. \code{precondition = TRUE}
#' therefore raises an error, and the fit always uses iterative
#' shrinkage-thresholding, whose certificate is the KKT residual above.
#' @references Beck, A. and Teboulle, M. (2009). A fast iterative
#'   shrinkage-thresholding algorithm for linear inverse problems. \emph{SIAM
#'   Journal on Imaging Sciences} 2(1), 183--202. \doi{10.1137/080716542}
#'
#'   Friedman, J., Hastie, T. and Tibshirani, R. (2008). Sparse inverse
#'   covariance estimation with the graphical lasso. \emph{Biostatistics}
#'   9(3), 432--441. \doi{10.1093/biostatistics/kxm045}
#'
#'   Yuan, M. and Lin, Y. (2007). Model selection and estimation in the
#'   Gaussian graphical model. \emph{Biometrika} 94(1), 19--35.
#'   \doi{10.1093/biomet/asm018}
#' @seealso \code{\link{mc_select_lambda}} to choose \code{lambda} by an
#'   information criterion.
#' @family structure selection
#' @examples
#' \donttest{
#' set.seed(1); n <- 300; x <- rnorm(n)
#' y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- rnorm(n)
#' d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' pen <- mc_graph_penalized(f, lambda = 20)
#' pen$edges
#' pen$status
#' }
#' @export
mc_graph_penalized <- function(fit, lambda, max_iter = 300L, tol = 1e-6,
                               precondition = FALSE, kkt_tol = 1e-4,
                               start = NULL) {
  stopifnot(inherits(fit, "mcglm_graph"), fit$graph$type == "ggm")
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
    stop("lambda must be one finite non-negative number.")
  if (isTRUE(precondition))
    stop("precondition = TRUE is disabled: the preconditioned step is not the ",
         "proximal operator in the metric of J and can report convergence at ",
         "a point that violates the KKT conditions. Use the default ",
         "precondition = FALSE.")
  R <- fit$graph$R; n <- fit$graph$n; Shat <- fit$graph$Shat
  .check_residual_variance(Shat, "the penalized fit")
  bases <- mc_precision_bases(R)
  E <- nrow(bases$edges)
  edge_idx <- (R + 1L):(R + E)
  eta <- c(1 / diag(Shat), rep(0, E))
  if (!is.null(start)) {
    start <- as.numeric(start)
    if (length(start) != R + E || any(!is.finite(start)))
      stop("start must hold ", R + E, " finite values: the diagonal ",
           "parameters and one per edge of the complete graph.", call. = FALSE)
    if (!isTRUE(cpp_check_spd(cpp_build_precision(start, bases$cube))$is_spd))
      stop("start does not give a positive definite precision.", call. = FALSE)
    eta <- start
  }
  conv <- FALSE
  status <- "maximum_iterations"
  # Smooth part f = -l^p of the objective, infinite outside the SPD cone, and
  # the penalty; the objective is their sum.
  smooth <- function(par) {
    K <- cpp_build_precision(par, bases$cube)
    chk <- cpp_check_spd(K)
    if (!isTRUE(chk$is_spd)) return(Inf)
    ld <- as.numeric(determinant(K, logarithm = TRUE)$modulus)
    0.5 * n * (sum(K * Shat) - ld)
  }
  penalty <- function(par) lambda * sum(abs(par[edge_idx]))
  fs <- smooth(eta)
  obj <- fs + penalty(eta)
  objective_path <- obj
  kkt <- NA_real_
  param_gate <- FALSE
  kkt_at <- function(par) {
    Kc <- cpp_build_precision(par, bases$cube)
    Sc <- cpp_precision_to_cov(Kc)
    pjc <- cpp_pearson_between(Sc, Shat, n,
                               cpp_precision_derivatives(Sc, bases$cube))
    .kkt_residual(par, pjc$psi, lambda, edge_idx)
  }
  for (it in seq_len(max_iter)) {
    K <- cpp_build_precision(eta, bases$cube)
    Sig <- cpp_precision_to_cov(K)
    D <- cpp_precision_derivatives(Sig, bases$cube)
    pj <- cpp_pearson_between(Sig, Shat, n, D)
    kkt <- .kkt_residual(eta, pj$psi, lambda, edge_idx)
    base_step <- 1 / (max(diag(pj$J), na.rm = TRUE) + 1e-8)
    accepted <- FALSE
    eta_new <- eta
    fs_new <- fs
    psi <- as.numeric(pj$psi)
    for (bt in 0:40) {
      step <- base_step * 0.5^bt
      eta_star <- eta + step * psi
      proposal <- eta_star
      proposal[edge_idx] <- cpp_soft_threshold(
        eta_star[edge_idx], step * lambda)
      # Sufficient decrease: the quadratic model of f with curvature 1/step
      # majorizes f at the proposal (grad f = -psi); slack at rounding level.
      move <- as.numeric(proposal) - eta
      candidate_fs <- smooth(proposal)
      if (is.finite(candidate_fs) &&
          candidate_fs <= fs - sum(psi * move) + sum(move^2) / (2 * step) +
            1e-15 * (1 + abs(fs))) {
        eta_new <- as.numeric(proposal)
        fs_new <- candidate_fs
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      status <- "objective_line_search_failed"
      break
    }
    rel <- max(abs(eta_new - eta) / (1 + abs(eta)))
    eta <- eta_new
    fs <- fs_new
    obj <- fs + penalty(eta)
    objective_path <- c(objective_path, obj)
    # Both gates: a small step and a small KKT residual at the accepted
    # iterate. A small step alone does not stop the iteration.
    param_gate <- rel < tol
    if (param_gate) {
      kkt <- kkt_at(eta)
      if (kkt / n <= kkt_tol) {
        conv <- TRUE
        status <- "converged"
        break
      }
    }
  }
  if (!conv && identical(status, "maximum_iterations")) {
    kkt <- kkt_at(eta)
    if (param_gate) status <- "kkt_not_certified"
  }
  # Soft-thresholding returns exact zeros, so the support is read as such.
  active <- which(eta[edge_idx] != 0)
  ed <- bases$edges[active, , drop = FALSE]
  K <- cpp_build_precision(eta, bases$cube)
  list(eta = eta, edges = ed, n_selected = length(active),
       Sigma_b = cpp_precision_to_cov(K), K_b = K,
       iterations = it, converged = conv, status = status,
       objective = obj, objective_path = objective_path,
       kkt = kkt, kkt_unit = kkt / n, kkt_tol = kkt_tol)
}

#' Select the penalty parameter along an L1 path
#'
#' Builds a geometric grid of penalty values down from the smallest value that
#' zeroes every edge, \eqn{\lambda_{\max}=n\max_{i<j}|\hat S_{e,ij}|}, to
#' \code{lambda_min_ratio} times that value. At each grid point the penalized
#' model of \code{\link{mc_graph_penalized}} is fitted, its support is refitted
#' without penalization, and the information criterion is evaluated on that
#' refit. The selected penalty is
#' \deqn{\hat\lambda=\arg\min_{\lambda}\ \mathrm{crit}\bigl(\hat{\mathcal G}_\lambda\bigr),}
#' with \eqn{\hat{\mathcal G}_\lambda} the support of the penalized solution
#' at \eqn{\lambda}; ties go to the larger penalty. The grid starts at
#' \eqn{\max(\lambda_{\max},10^{-3})}: the lower limit keeps the geometric grid
#' defined when the residual covariances are numerically zero, and it does not
#' change the first solution, which is the empty graph at any penalty of at
#' least \eqn{\lambda_{\max}}.
#'
#' @section Independent fits, warm starts and the unpenalized end:
#' By default every grid point is an independent fit from the same starting
#' value, the diagonal fit \eqn{\eta_{kk}=1/\widehat S_{e,kk}}, which also
#' gives the largest penalty of the path; both need
#' \eqn{\widehat S_{e,kk}>0} for every response, and a response whose Pearson
#' residuals have no variance is refused with an error. The points can be
#' evaluated in any order and with any number of
#' workers, and the selected penalty does not depend on either. With
#' \code{warm_start = TRUE} the grid is run from the largest penalty down, each
#' fit starting at the solution of the one before, which is the usual way of
#' computing a path; the fits are then sequential and \code{workers} is not
#' used. The problem has a unique solution, so the two routes share the target
#' and differ in the number of iterations. Certified solutions satisfy the
#' optimality conditions only up to the tolerance, so their supports can differ
#' in pairs whose zero-edge optimality condition, the scaled gap between the
#' residual covariance and the fitted covariance, is close to \eqn{\lambda},
#' within a margin proportional to the tolerance and to the conditioning of the
#' curvature. With
#' \code{include_zero = TRUE} the unpenalized fit of the complete graph is added
#' to the candidates as the end of the path, \eqn{\lambda=0}; the geometric grid
#' never reaches it. By default the grid is independent and does not include
#' that end.
#'
#' @param fit an object of class \code{"mcglm_graph"} (GGM layer).
#' @param n_lambda number of grid points (default 30).
#' @param lambda_min_ratio ratio of the smallest to the largest penalty
#'   (default \eqn{10^{-3}}).
#' @param criterion \code{"bic"} (default), \code{"aic"}, \code{"qic"},
#'   \code{"gbic"} or its alias \code{"clbic"}, or \code{"ebic"}.
#' @param refit logical; if \code{TRUE} (default) the returned fit is the
#'   unpenalized refit of the selected support. If \code{FALSE}, the penalized
#'   solution is attached to the fit instead; its covariance is kept for
#'   diagnostics but it is flagged as not certified for edge-wise inference.
#' @param workers number of parallel workers for the penalty grid. The grid
#'   points are independent fits, so the selected penalty does not depend on
#'   this value. Sequential by default; see \code{\link{mc_graph_stability}}
#'   for why it is opt-in. Not used with \code{warm_start = TRUE}.
#' @param warm_start logical; run the grid in sequence from the largest
#'   penalty, each fit starting at the previous solution.
#' @param include_zero logical; add the unpenalized complete graph as the
#'   candidate at \eqn{\lambda=0}.
#' @param ebic_gamma the parameter of the extended BIC, used only with
#'   \code{criterion = "ebic"}.
#' @param cluster optional independent-cluster identifier, one per residual
#'   row, passed to \code{\link{mc_compute_criterion}}: Pearson contributions
#'   are summed within clusters in QIC and GBIC, and the number of clusters
#'   replaces the number of rows in the logarithmic penalty of BIC, GBIC and
#'   EBIC.
#' @param small_sample logical; apply the cluster factor \eqn{G/(G-1)}, with
#'   \eqn{G} the number of clusters, to the empirical variability of the
#'   clustered criteria (default \code{FALSE}); ignored without
#'   \code{cluster}.
#' @param max_iter,tol,kkt_tol passed to \code{\link{mc_graph_penalized}} at
#'   every grid point.
#' @return An object of class \code{"mcg_selection"} (see
#'   \code{\link{mc_graph_stepwise}} for the common components) with
#'   \code{method = "penalized"}, the selected \code{lambda}, a logical
#'   \code{warm_start}, and the table \code{lambda_path}, a data frame with one
#'   row per grid point (plus one for \eqn{\lambda=0} when
#'   \code{include_zero = TRUE}) and columns \code{lambda}, \code{n_edges},
#'   \code{criterion}, \code{penalized_converged}, \code{penalized_status},
#'   \code{refit_converged} and \code{refit_status}.
#' @section Which grid points compete:
#' The criterion is evaluated on the unpenalized refit of the support of every
#' grid point whose penalized fit was certified, by the parameter change and
#' by the KKT residual, and whose refit converged. A point without the
#' certificate stays in \code{lambda_path} with its status and does not enter
#' the comparison. An error is raised when no grid point qualifies.
#' @seealso \code{\link{mc_graph_penalized}}, \code{\link{mc_compute_criterion}}.
#' @family structure selection
#' @examples
#' \donttest{
#' set.seed(1); n <- 300; x <- rnorm(n)
#' y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- rnorm(n)
#' d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' sel <- mc_select_lambda(f, n_lambda = 10)
#' sel$edges
#' sel$lambda_path[, c("lambda", "n_edges", "criterion", "penalized_status")]
#' }
#' @export
mc_select_lambda <- function(fit, n_lambda = 30L, lambda_min_ratio = 1e-3,
                             criterion = c("bic", "aic", "qic", "gbic", "clbic",
                                           "ebic"),
                             refit = TRUE, workers = 1L, warm_start = FALSE,
                             include_zero = FALSE, ebic_gamma = 0.5,
                             cluster = NULL, small_sample = FALSE,
                             max_iter = 300L, tol = 1e-6, kkt_tol = 1e-4) {
  stopifnot(inherits(fit, "mcglm_graph"), fit$graph$type == "ggm")
  criterion <- match.arg(criterion)
  for (a in c("max_iter", "tol", "kkt_tol")) {
    v <- get(a)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0 ||
        (a == "max_iter" && (v < 1 || v != round(v))))
      stop(a, if (a == "max_iter") " must be one positive integer." else
        " must be one finite non-negative number.", call. = FALSE)
  }
  for (a in c("warm_start", "include_zero"))
    if (!is.logical(get(a)) || length(get(a)) != 1L || is.na(get(a)))
      stop(a, " must be a single TRUE or FALSE.", call. = FALSE)
  R <- fit$graph$R; n <- fit$graph$n; Shat <- fit$graph$Shat
  .check_residual_variance(Shat, "the penalty path")
  diag_fit <- .refit_between(fit, mc_between_graph("ggm", matrix(0L, 0, 2)))
  if (!isTRUE(diag_fit$graph$converged))
    stop("The diagonal starting model did not converge: ",
         if (is.null(diag_fit$graph$status)) "status unavailable" else
           diag_fit$graph$status, ".")
  Sd <- diag_fit$graph$Sigma_b
  lam_max <- 0
  for (i in seq_len(R - 1L)) for (j in (i + 1L):R)
    lam_max <- max(lam_max, abs(n * (Sd[i, j] - Shat[i, j])))
  lam_max <- max(lam_max, 1e-3)
  lambdas <- exp(seq(log(lam_max), log(lam_max * lambda_min_ratio),
                     length.out = n_lambda))
  path <- data.frame(
    lambda = lambdas, n_edges = NA_real_, criterion = NA_real_,
    penalized_converged = FALSE, penalized_status = NA_character_,
    refit_converged = FALSE, refit_status = NA_character_,
    stringsAsFactors = FALSE
  )
  best <- list(crit = Inf, edges = NULL, lambda = NA, penalized = NULL)
  crit_of <- function(ff) tryCatch(
    unname(mc_compute_criterion(
      ff, criterion, ebic_gamma = ebic_gamma, cluster = cluster,
      small_sample = small_sample)),
    error = function(e) NA_real_)

  # A grid point: the penalized fit, then the unpenalized refit of its support
  # on which the criterion is evaluated.
  one <- function(k, start = NULL) {
    pen <- tryCatch(mc_graph_penalized(fit, lambdas[k], start = start,
                                       max_iter = max_iter, tol = tol,
                                       kkt_tol = kkt_tol),
                    error = function(e) {
                      list(converged = FALSE,
                           status = paste0("error: ", conditionMessage(e)))
                    })
    out <- list(penalized_converged = isTRUE(pen$converged),
                penalized_status = if (is.null(pen$status))
                  "status unavailable" else pen$status,
                n_edges = NA_real_, refit_converged = FALSE,
                refit_status = NA_character_, criterion = NA_real_,
                pen = pen)
    if (!isTRUE(pen$converged)) return(out)
    out$n_edges <- pen$n_selected
    refit_error <- NULL
    ff <- tryCatch(
      .refit_between(fit, mc_between_graph("ggm", .nz_edges(pen$edges))),
      error = function(e) {
        refit_error <<- paste0("error: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(ff)) {
      out$refit_status <- refit_error
      return(out)
    }
    out$refit_converged <- isTRUE(ff$graph$converged)
    out$refit_status <- if (is.null(ff$graph$status))
      "status unavailable" else ff$graph$status
    if (!isTRUE(ff$graph$converged)) return(out)
    cv <- crit_of(ff)
    if (length(cv) == 1L && is.finite(cv)) out$criterion <- cv
    out
  }

  if (isTRUE(warm_start)) {
    if (!identical(as.integer(workers), 1L))
      warning("warm_start runs the grid in sequence; workers is not used.",
              call. = FALSE)
    # From the largest penalty down, each fit starting at the last solution
    # that was certified; an uncertified one is not propagated.
    res <- vector("list", length(lambdas))
    start <- NULL
    for (k in seq_along(lambdas)) {
      res[[k]] <- one(k, start)
      if (isTRUE(res[[k]]$penalized_converged)) start <- res[[k]]$pen$eta
    }
  } else {
    # The grid points are independent fits: no warm starts, no shared state,
    # no random draws. They can therefore be evaluated in any order, and with
    # any worker count, without changing the selected penalty. The reduction
    # below is sequential so that ties resolve to the smallest index.
    res <- .mcg_lapply(seq_along(lambdas), one, workers = workers)
  }

  for (k in seq_along(lambdas)) {
    z <- res[[k]]
    msg <- .mcg_error_message(z)
    if (!is.null(msg)) {
      path$penalized_status[k] <- paste0("error: ", msg)
      next
    }
    path$penalized_converged[k] <- z$penalized_converged
    path$penalized_status[k] <- z$penalized_status
    path$n_edges[k] <- z$n_edges
    path$refit_converged[k] <- z$refit_converged
    path$refit_status[k] <- z$refit_status
    path$criterion[k] <- z$criterion
    if (is.finite(z$criterion) && z$criterion < best$crit)
      best <- list(crit = z$criterion, edges = z$pen$edges,
                   lambda = lambdas[k], penalized = z$pen)
  }

  if (isTRUE(include_zero)) {
    # The end of the path is the unpenalized complete graph: no penalized fit
    # is needed to reach it, since its solution is the one of the saturated
    # estimating equations.
    full_edges <- mc_precision_bases(R)$edges
    ff <- tryCatch(.refit_between(fit, mc_between_graph("ggm")),
                   error = function(e) NULL)
    ok0 <- !is.null(ff) && isTRUE(ff$graph$converged)
    cv0 <- if (ok0) crit_of(ff) else NA_real_
    path <- rbind(path, data.frame(
      lambda = 0, n_edges = nrow(full_edges), criterion = cv0,
      penalized_converged = NA, penalized_status = "unpenalized",
      refit_converged = ok0,
      refit_status = if (is.null(ff)) "error" else ff$graph$status,
      stringsAsFactors = FALSE))
    if (is.finite(cv0) && cv0 < best$crit)
      best <- list(crit = cv0, edges = full_edges, lambda = 0,
                   penalized = NULL)
  }

  if (!is.finite(best$crit))
    stop("No L1 candidate produced both a converged penalized fit and a ",
         "converged unpenalized support refit.")
  spec <- mc_between_graph("ggm", .nz_edges(best$edges))
  final <- if (isTRUE(refit) || is.null(best$penalized))
    .refit_between(fit, spec) else
      .attach_penalized_fit(fit, spec, best$penalized)
  out <- .mcg_selection("penalized", criterion, best$edges, final)
  out$lambda <- best$lambda; out$lambda_path <- path
  out$warm_start <- isTRUE(warm_start)
  out
}

.attach_penalized_fit <- function(fit, spec, pen) {
  R <- fit$graph$R
  spec <- .canonical_spec(spec, R)
  full_edges <- mc_precision_bases(R)$edges
  selected <- .nz_edges(pen$edges)
  idx <- if (nrow(selected) == 0L) integer(0) else
    match(paste(selected[, 1L], selected[, 2L], sep = ":"),
          paste(full_edges[, 1L], full_edges[, 2L], sep = ":"))
  par <- c(pen$eta[seq_len(R)], pen$eta[R + idx])
  bld <- .make_builder(spec, R)
  meta <- attr(bld, "meta")
  b <- bld(par)
  pj <- cpp_pearson_between(b$Sigma, fit$graph$Shat, fit$graph$n, b$Dlist)
  names(par) <- .par_names(spec, meta)
  fit$graph$type <- "ggm"
  fit$graph$par <- par
  fit$graph$Sigma_b <- b$Sigma
  fit$graph$K_b <- b$K
  fit$graph$J <- pj$J
  rcond_tol <- fit$graph$rcond_tol
  if (is.null(rcond_tol)) rcond_tol <- .Machine$double.eps^(2 / 3)
  fit$graph$vcov <- .safe_solve(pj$J, rcond_tol = rcond_tol)
  fit$graph$edge_pos <- meta$edge_pos
  fit$graph$npar <- meta$npar
  fit$graph$iter <- pen$iterations
  fit$graph$converged <- pen$converged
  fit$graph$status <- pen$status
  fit$graph$score_norm <- max(abs(as.numeric(pj$psi)))
  fit$graph$score_norm_unit <- fit$graph$score_norm / fit$graph$n
  fit$graph$rcond_J <- tryCatch(
    rcond(0.5 * (pj$J + t(pj$J))), error = function(e) NA_real_
  )
  fit$graph$rcond_tol <- rcond_tol
  fit$graph$convergence_gates <- c(
    kkt = isTRUE(pen$converged) && is.finite(pen$kkt_unit) &&
      pen$kkt_unit <= pen$kkt_tol,
    sensitivity = is.finite(fit$graph$rcond_J) &&
      fit$graph$rcond_J > rcond_tol
  )
  fit$graph$stopping_rule <- paste0(
    "penalized relative change and kkt_unit <= ", format(pen$kkt_tol)
  )
  fit$graph$inversion_status <- attr(
    fit$graph$vcov, "inversion_status", exact = TRUE
  )
  fit$graph$inversion_ridge <- attr(
    fit$graph$vcov, "inversion_ridge", exact = TRUE
  )
  # A covariance is retained for diagnostics, but an unpenalized support refit
  # is required before conventional edge-wise inference is certified.
  fit$graph$inference_ok <- FALSE
  fit$graph$linear_solve_status <- "not_applicable_penalized"
  fit$graph$linear_solve_ridge <- 0
  fit$graph$meta <- meta
  fit$graph$init_info <- NULL
  fit$graph$spec <- spec
  fit$graph$Dlist <- b$Dlist
  fit
}

# ---- internal helpers -------------------------------------------------------

.present_matrix <- function(fit) {
  R <- fit$graph$R
  ed <- fit$graph$spec$edges
  A <- matrix(0L, R, R)
  if (is.null(ed)) {                              # complete graph
    A[upper.tri(A)] <- 1L; A <- A + t(A)
  } else if (nrow(ed) > 0) {
    for (e in seq_len(nrow(ed))) { i <- ed[e, 1]; j <- ed[e, 2]; A[i, j] <- A[j, i] <- 1L }
  }
  A
}

# Edges of a GGM fit in the order of its edge parameters. They are read from
# the bases the fit was built with, so they match par[edge_pos] by definition.
.ggm_param_edges <- function(fit) {
  g <- fit$graph
  ed <- g$meta$bases$edges
  if (is.null(ed)) ed <- .validate_edges(g$spec$edges, g$R)
  ed
}

# Position of the parameter of edge (i, j) in a GGM fit, matched by the pair.
.edge_param_pos <- function(fit, i, j) {
  ed <- .ggm_param_edges(fit)
  e <- which(ed[, 1L] == min(i, j) & ed[, 2L] == max(i, j))
  if (length(e) != 1L)
    stop("edge (", i, ", ", j, ") is not in the support of the fit.",
         call. = FALSE)
  fit$graph$R + e
}

.nz_edges <- function(edges) {
  if (is.null(edges) || nrow(edges) == 0) return(matrix(0L, 0, 2))
  matrix(as.integer(edges), ncol = 2)
}

# Efficient (nuisance-adjusted) score for candidate edges of the GGM layer.
#
# The current parametrization is augmented with the candidate edge parameter
# held at zero, which leaves Sigma unchanged. Writing e for the candidate
# coordinate and lambda for the remaining ones,
#   u_e      = psi_e - J_{e,lam} J_{lam,lam}^{-1} psi_lam = a' psi,
#   var(u_e) = a' V a,   a = (-J_{lam,lam}^{-1} J_{lam,e}, 1),
# with V = J (working Bartlett identity) or V the empirical fourth-moment
# variability. Under V = J the variance collapses to the Schur complement
# J_{e.lam} <= J_ee, which is why dividing by J_ee alone understates the
# statistic.
#
# Every piece has a closed form in the linear-precision layer. For a parameter
# on the pair (k, l), with weight w = 1/2 on the diagonal and 1 on an edge,
#   psi         = n w (Sigma_kl - S_kl),
#   J[(k,l),(i,j)] = n w w' (Sigma_ki Sigma_lj + Sigma_kj Sigma_li),
#   psi_unit    = w (Sigma_kl - e_k e_l),
# and J_{lam,lam} is the sensitivity of the fit, whose inverse the fit stores.
.efficient_edge_scores <- function(fit, edges, variability = "model") {
  g <- fit$graph
  R <- g$R; n <- g$n
  S <- as.matrix(g$Sigma_b); Sh <- as.matrix(g$Shat)
  # pairs and weights of the current parameters, in the order of g$par
  cur <- .nz_edges(.ggm_param_edges(fit))
  pk <- c(seq_len(R), cur[, 1L]); pl <- c(seq_len(R), cur[, 2L])
  w <- c(rep(0.5, R), rep(1, nrow(cur)))
  Vll <- as.matrix(g$vcov)
  psi_l <- n * w * (S[cbind(pk, pl)] - Sh[cbind(pk, pl)])
  Vpsi <- as.numeric(Vll %*% psi_l)
  emp <- !identical(variability, "model")
  if (emp) {
    Ec <- scale(as.matrix(g$residuals), center = TRUE, scale = FALSE)
    Psi_l <- .pearson_scores(g$residuals, g$Sigma_b, g$Dlist)
  }
  q <- nrow(edges)
  stat <- numeric(q); vv <- numeric(q)
  for (e in seq_len(q)) {
    i <- min(edges[e, ]); j <- max(edges[e, ])
    Jle <- n * w * (S[pk, i] * S[pl, j] + S[pk, j] * S[pl, i])
    Jee <- n * (S[i, i] * S[j, j] + S[i, j]^2)
    u <- n * (S[i, j] - Sh[i, j]) - sum(Jle * Vpsi)
    cvec <- as.numeric(Vll %*% Jle)
    var_u <- if (emp)
      sum(((S[i, j] - Ec[, i] * Ec[, j]) - as.numeric(Psi_l %*% cvec))^2) else
        Jee - sum(Jle * cvec)
    stat[e] <- if (is.finite(var_u) && var_u > 0) u * u / var_u else NA_real_
    vv[e] <- var_u
  }
  list(stat = stat, variance = vv)
}

# KKT / subdifferential residual of the L1-penalized Gaussian working problem
#   min f(eta) + lambda * sum |eta_E|,   f = 0.5 n [tr(K Shat) - log|K|],
# using grad f = -psi. For a free coordinate the condition is |grad| = 0; for a
# penalized coordinate it is |grad + lambda sign(eta)| = 0 when eta is nonzero
# and |grad| <= lambda when it is zero. Returned on the accumulated scale of
# psi; divide by n for the per-unit scale used by the tolerance.
.kkt_residual <- function(eta, psi, lambda, edge_idx, zero_tol = 0) {
  gr <- -as.numeric(psi)
  viol <- abs(gr)
  for (p in edge_idx) {
    viol[p] <- if (abs(eta[p]) > zero_tol)
      abs(gr[p] + lambda * sign(eta[p])) else max(abs(gr[p]) - lambda, 0)
  }
  max(viol)
}

.mcg_selection <- function(method, criterion, edges, fit, trace = NULL) {
  g <- mc_get_graph(fit)
  structure(list(method = method, criterion = criterion,
                 graph = g$adjacency, edges = .nz_edges(edges),
                 n_edges = if (is.null(edges)) 0L else nrow(.nz_edges(edges)),
                 trace = trace, lambda = NA, lambda_path = NULL,
                 cpdag = NULL, vstructures = NULL, fit = fit),
            class = "mcg_selection")
}

#' Methods for graph structure selection results
#'
#' Print, summarize and plot the result of \code{\link{mc_graph_stepwise}} or
#' \code{\link{mc_select_lambda}}. \code{print} shows the selector, the
#' criterion, the number and list of selected edges and, for the L1 route, the
#' selected penalty. \code{summary} adds the search trace (stepwise), the head
#' of the penalty path (L1), and the v-structures and CPDAG when the object
#' carries them. \code{plot} draws the network of the selected fit with
#' \code{\link{mc_plot_network}}.
#'
#' @param x,object an object of class \code{"mcg_selection"}.
#' @param ... for \code{plot}, further arguments passed to
#'   \code{\link{mc_plot_network}}; ignored by \code{print} and
#'   \code{summary}.
#' @return \code{print} and \code{summary} return the object invisibly;
#'   \code{plot} returns invisibly the value of \code{mc_plot_network} applied
#'   to the selected fit (or the object itself when it carries no fit).
#' @seealso \code{\link{mc_graph_stepwise}}, \code{\link{mc_select_lambda}}.
#' @examples
#' \donttest{
#' set.seed(1); n <- 200; x <- rnorm(n)
#' y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- rnorm(n)
#' d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' sel <- mc_graph_stepwise(f)
#' print(sel)
#' summary(sel)
#' }
#' @name mcg_selection-methods
#' @export
print.mcg_selection <- function(x, ...) {
  cat(sprintf("Graph structure selection (%s, criterion = %s)\n",
              x$method, x$criterion))
  cat(sprintf("Selected edges: %d\n", x$n_edges))
  if (!is.null(x$lambda) && !is.na(x$lambda))
    cat(sprintf("Selected lambda: %.4f\n", x$lambda))
  if (!is.null(x$edges) && x$n_edges > 0) { cat("Edges:\n"); print(x$edges) }
  invisible(x)
}

#' @rdname mcg_selection-methods
#' @export
summary.mcg_selection <- function(object, ...) {
  print(object)
  if (!is.null(object$trace)) { cat("\nSearch trace:\n"); print(object$trace) }
  if (!is.null(object$lambda_path)) {
    cat("\nLambda path (head):\n"); print(utils::head(object$lambda_path))
  }
  if (!is.null(object$vstructures) && nrow(object$vstructures) > 0) {
    cat("\nV-structures:\n"); print(object$vstructures)
  }
  if (!is.null(object$cpdag)) { cat("\nCPDAG:\n"); print(object$cpdag) }
  invisible(object)
}

#' @rdname mcg_selection-methods
#' @export
plot.mcg_selection <- function(x, ...) {
  if (inherits(x$fit, "mcglm_graph")) return(invisible(mc_plot_network(x$fit, ...)))
  invisible(x)
}

#' Selection-frequency aggregation of the between-response graph structure
#'
#' Repeats graph selection in \code{B} half-samples and records the empirical
#' selection frequency \eqn{\hat\Pi_e} of every edge, the proportion of valid
#' subsamples in which the selector retains it. Edges with
#' \eqn{\hat\Pi_e\ge\pi_{\mathrm{thr}}} are returned as recurrent. The
#' frequency estimates the probability that the selector retains the edge in a
#' random half-sample; it is not the probability that the edge belongs to the
#' graph. With a selector tuned by an information criterion, the default, the
#' procedure is a diagnostic aggregation scheme and not the error-controlled
#' stability selection of Meinshausen and Buehlmann (2010), because the number
#' of edges each subsample selects is not controlled: the recurrent set then
#' carries no formal bound on its false selections. Passing \code{q} gives the
#' version for which such a bound is reported, under the conditions stated
#' below.
#'
#' @section Stability selection with error control:
#' With \code{q}, every subsample selects the first \code{q} edges to enter
#' along the L1 path of \code{\link{mc_graph_penalized}}, run with warm starts
#' over a geometric grid of 50 penalties from the value that zeroes every edge
#' in that subsample, \eqn{\lambda_{\max}} (or \eqn{10^{-3}} if larger), down
#' to \eqn{10^{-3}} times it; edges that enter at the same grid point are
#' ranked by the size of their coefficient there. Grid points whose fit is not
#' certified are skipped, and fewer than \code{q} edges are selected when
#' fewer enter along the grid, so \eqn{q} is an upper limit on the number of
#' edges per subsample. The
#' expected number of falsely selected edges among those whose frequency
#' reaches \code{pi_thr} is then bounded by Theorem 1 of Meinshausen and
#' Buehlmann (2010),
#' \deqn{E(V)\le\frac{q^2}{(2\pi_{\mathrm{thr}}-1)\,|E_{\max}|},\qquad
#'   \tfrac12<\pi_{\mathrm{thr}}\le1,}
#' where \eqn{V} is the number of falsely selected edges, \eqn{q} the number of
#' edges each subsample selects and \eqn{|E_{\max}|=R(R-1)/2} the number of
#' pairs, reported as \code{error_bound}. The bound
#' rests on two conditions the package cannot check: the selection indicators
#' of the null edges are exchangeable, and the selector is no worse than random
#' guessing. It is a bound on an expected count, not on a probability, and it
#' refers to frequencies computed over all the subsamples; when some
#' subsamples fail, \code{selection_prob_conservative} is the frequency that
#' counts a failed subsample as selecting no edge, whereas the recurrent set is
#' read from the frequencies among the valid ones. With
#' \code{sampling = "complementary"} the subsamples come in
#' \code{B} pairs of disjoint halves, the complementary pairs of Shah and
#' Samworth (2013). What the bound controls is the edges of the residual
#' layer as the second stage sees them; with \code{refit = "residual"} the
#' first stage is held fixed.
#'
#' @section What each half-sample is:
#' With \code{refit = "residual"} the half-sample is drawn from the stored
#' Pearson residuals and the first stage is held fixed, so the frequencies
#' describe the second stage alone. With \code{refit = "full"} the half-sample
#' is drawn from the rows of the data, the first-stage MCGLM (multivariate
#' covariance generalized linear model) is refitted on it, and the selector
#' runs on the residuals of that refit, so the frequencies carry the
#' variability of both stages. The second form accounts for the estimation of
#' the first stage and costs an MCGLM fit per half-sample; the first is the
#' default and is cheaper.
#'
#' A half-sample holds \eqn{\lfloor n/2\rfloor} of the \eqn{n} rows, drawn
#' without replacement, and at least two. With
#' \code{sampling = "complementary"} one random permutation of the rows is
#' split into two disjoint halves of \eqn{\lfloor n/2\rfloor} rows each (at
#' least one), so that with \eqn{n} odd one row is left out of the pair.
#'
#' With \code{cluster}, half of the clusters are drawn instead of half of the
#' rows, under the same rounding and the same minimum sizes counted in
#' clusters, and every row of a drawn cluster enters the half-sample. That is the
#' form to use when the unit of dependence is a cluster of rows, and it is
#' required with \code{refit = "full"} when the within-response structure is
#' not the identity: each drawn cluster then keeps its own block of the
#' structure, as in \code{\link{mc_graph_bootstrap}}.
#'
#' @param fit an object of class \code{"mcglm_graph"} (GGM layer).
#' @param B number of subsamples (default 100); with
#'   \code{sampling = "complementary"}, the number of pairs, so \eqn{2B}
#'   subsamples are drawn.
#' @param pi_thr selection-frequency threshold in \eqn{[0,1]} (default 0.6).
#' @param method per-subsample selector: \code{"stepwise"}
#'   (\code{\link{mc_graph_stepwise}}, default) or \code{"penalized"}
#'   (\code{\link{mc_select_lambda}}). It is set to \code{"penalized"} when
#'   \code{q} is given.
#' @param criterion information criterion passed to the selector
#'   (\code{"bic"}, the default, \code{"aic"}, \code{"qic"}, \code{"gbic"},
#'   its alias \code{"clbic"}, or \code{"ebic"}); not used when \code{q} is
#'   given.
#' @param refit \code{"residual"} (default) resamples the stored Pearson
#'   residuals with the first stage held fixed; \code{"full"} resamples rows
#'   of the data and
#'   refits both stages in every half-sample.
#' @param workers number of parallel workers for the resampling loop. The
#'   subsamples are drawn in the parent before any dispatch, so a given seed
#'   yields the same subsamples at any worker count and the result does not
#'   depend on it. Sequential by default, deliberately: user code often
#'   parallelizes at an outer level, and nesting a second level of workers
#'   underneath oversubscribes the machine. Forking is used where available; on
#'   Windows the call falls back to sequential evaluation. Inside a parallel
#'   worker the argument is ignored, with a warning.
#' @param cluster optional vector with one identifier per observational row;
#'   half of the clusters, rounded down, are drawn without replacement in
#'   every half-sample (at least two with \code{sampling = "half"}, at least
#'   one per half with \code{sampling = "complementary"}).
#' @param q optional integer between 1 and \eqn{R(R-1)/2}: the number of edges
#'   each subsample selects along the L1 path, which turns on the
#'   error-controlled procedure described above and requires
#'   \code{method = "penalized"} when \code{method} is given.
#' @param sampling \code{"half"} (default), \code{B} independent half-samples
#'   drawn without replacement, or \code{"complementary"}, \code{B} pairs of
#'   disjoint halves.
#' @param ... further arguments passed to the per-subsample selector
#'   (\code{\link{mc_graph_stepwise}} or \code{\link{mc_select_lambda}});
#'   ignored when \code{q} is given. For the \code{print} method, ignored.
#' @return An object of class \code{"mcg_stability"} with the empirical
#'   selection-frequency matrix among valid subsamples
#'   (\code{selection_prob}), its conservative counterpart divided by all
#'   attempted subsamples (\code{selection_prob_conservative}), the numbers of
#'   valid and failed subsamples (\code{B_valid}, \code{B_failed}), and the
#'   recurrent edge set (\code{stable_edges}, a two-column matrix of the pairs
#'   whose frequency among valid subsamples is at least \code{pi_thr}). With
#'   \code{q} it also carries \code{q} and the
#'   \code{error_bound}, \code{NA} when \code{pi_thr} is not above one half.
#'   The list also records \code{B}, \code{pi_thr}, \code{method},
#'   \code{criterion}, \code{refit}, \code{sampling}, \code{clustered},
#'   \code{n_subsamples} and the data frame \code{failures} (subsample index
#'   and error message of each failed subsample). The matrices have the
#'   response names as dimnames and hold the frequency of each pair in both
#'   triangles. \code{print} returns its argument invisibly.
#' @references Meinshausen, N. and Buehlmann, P. (2010). Stability selection.
#'   \emph{Journal of the Royal Statistical Society: Series B} 72(4),
#'   417--473. \doi{10.1111/j.1467-9868.2010.00740.x}
#'
#'   Shah, R. D. and Samworth, R. J. (2013). Variable selection with error
#'   control: another look at stability selection. \emph{Journal of the
#'   Royal Statistical Society: Series B} 75(1), 55--80.
#'   \doi{10.1111/j.1467-9868.2011.01034.x}
#' @seealso \code{\link{mc_graph_stepwise}}, \code{\link{mc_select_lambda}},
#'   \code{\link{mc_graph_penalized}}.
#' @family structure selection
#' @examples
#' \donttest{
#' set.seed(1); n <- 300
#' K <- diag(4)*2; K[1,2]<-K[2,1]<- -0.6; K[2,3]<-K[3,2]<- -0.6; K[3,4]<-K[4,3]<- -0.6
#' Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K)); x <- rnorm(n)
#' d <- data.frame(y1=Y[,1], y2=Y[,2], y3=Y[,3], y4=Y[,4], x=x)
#' f <- mcglm_graph(list(y1~x,y2~x,y3~x,y4~x), data = d)
#' mc_graph_stability(f, B = 40)
#' }
#' @export
mc_graph_stability <- function(fit, B = 100L, pi_thr = 0.6,
                               method = c("stepwise", "penalized"),
                               criterion = "bic",
                               refit = c("residual", "full"),
                               workers = 1L, cluster = NULL, q = NULL,
                               sampling = c("half", "complementary"), ...) {
  stopifnot(inherits(fit, "mcglm_graph"), fit$graph$type == "ggm")
  method_given <- !missing(method)
  method <- match.arg(method)
  refit <- match.arg(refit)
  sampling <- match.arg(sampling)
  E_max <- choose(fit$graph$R, 2L)
  if (!is.null(q)) {
    if (length(q) != 1L || !is.finite(q) || q < 1 || q != as.integer(q) ||
        q > E_max)
      stop("q must be one integer between 1 and the number of pairs, ", E_max,
           ".", call. = FALSE)
    q <- as.integer(q)
    if (method_given && method != "penalized")
      stop("q selects the first edges to enter along the L1 path, so it ",
           "applies to method = \"penalized\" only.", call. = FALSE)
    method <- "penalized"
  } else {
    # refused here, once, and not inside every subsample
    criterion <- match.arg(criterion, c("bic", "aic", "qic", "gbic", "clbic",
                                        "ebic"))
  }
  by_cluster <- NULL
  if (!is.null(cluster)) {
    u <- .two_stage_units(fit, cluster)
    by_cluster <- unname(split(seq_len(fit$graph$n), u$unit))
  }
  .check_resampling_ready(fit, refit, cluster)
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")
  if (length(pi_thr) != 1L || !is.finite(pi_thr) ||
      pi_thr < 0 || pi_thr > 1)
    stop("pi_thr must be one finite number in [0, 1].")
  B <- as.integer(B)
  g <- fit$graph
  E <- as.matrix(g$residuals); n <- nrow(E); R <- g$R
  half <- max(2L, floor(n / 2))
  count <- matrix(0, R, R)
  failures <- data.frame(subsample = integer(), error = character(),
                         stringsAsFactors = FALSE)
  B_valid <- 0L

  # Every subsample is drawn here, in the parent, before any dispatch. That is
  # what makes workers = 4 return exactly what workers = 1 returns for a given
  # seed: the worker count changes who does the arithmetic, never which rows
  # each subsample contains.
  # Without clusters the draw is a single call per subsample, so a given seed
  # always gives the same half-samples.
  draws <- if (identical(sampling, "half")) {
    if (is.null(by_cluster)) {
      lapply(seq_len(B), function(b) {
        rows <- sample.int(n, half)
        list(rows = rows, blocks = NULL)
      })
    } else {
      halfc <- max(2L, floor(length(by_cluster) / 2))
      lapply(seq_len(B), function(b) {
        blocks <- by_cluster[sample.int(length(by_cluster), halfc)]
        list(rows = unlist(blocks, use.names = FALSE), blocks = blocks)
      })
    }
  } else {
    # complementary pairs: one permutation per pair, split into two disjoint
    # halves of the rows, or of the clusters
    units <- if (is.null(by_cluster)) as.list(seq_len(n)) else by_cluster
    hu <- max(1L, floor(length(units) / 2))
    unlist(lapply(seq_len(B), function(b) {
      p <- sample.int(length(units))
      lapply(list(p[seq_len(hu)], p[hu + seq_len(hu)]), function(k) {
        blocks <- units[k]
        list(rows = unlist(blocks, use.names = FALSE),
             blocks = if (is.null(by_cluster)) NULL else blocks)
      })
    }), recursive = FALSE)
  }
  n_sub <- length(draws)

  spec <- fit$graph$spec
  ctrl <- .graph_control(fit)

  # The residual route swaps the residual block and leaves the rest of the
  # object alone, which is what the selector expects: it refits the graph
  # layer itself, from the complete graph, so refitting it here first would be
  # work thrown away. The full route
  # has to refit the first stage, and its graph layer comes along with it.
  make_fb <- function(rows, blocks) {
    if (identical(refit, "full"))
      return(.refit_both_stages(fit, rows, spec, ctrl, blocks = blocks))
    Eb <- E[rows, , drop = FALSE]
    Ec <- scale(Eb, center = TRUE, scale = FALSE)
    fb <- fit
    fb$graph$residuals <- Eb
    fb$graph$Shat <- crossprod(Ec) / nrow(Eb)
    fb$graph$n <- nrow(Eb)
    fb
  }

  one <- function(b) {
    error <- NULL
    fb <- tryCatch(make_fb(draws[[b]]$rows, draws[[b]]$blocks),
                   error = function(e) {
                     error <<- conditionMessage(e)
                     NULL
                   })
    if (is.null(fb))
      return(list(ok = FALSE, error = paste0("refit failed: ", error)))
    if (!is.null(q)) {
      ed <- tryCatch(.first_q_edges(fb, q), error = function(e) {
        error <<- conditionMessage(e)
        NULL
      })
      if (is.null(ed)) return(list(ok = FALSE, error = error))
      return(list(ok = TRUE, edges = ed))
    }
    sel <- tryCatch(
      if (method == "penalized")
        mc_select_lambda(fb, criterion = criterion, ...)
      else
        mc_graph_stepwise(fb, criterion = criterion, ...),
      error = function(e) {
        error <<- conditionMessage(e)
        NULL
      })
    valid <- inherits(sel, "mcg_selection") &&
      inherits(sel$fit, "mcglm_graph") &&
      isTRUE(sel$fit$graph$converged) &&
      !is.null(sel$edges)
    if (!valid) {
      if (is.null(error))
        error <- "selector returned an invalid or non-converged fit"
      return(list(ok = FALSE, error = error))
    }
    list(ok = TRUE, edges = sel$edges)
  }

  res <- .mcg_lapply(seq_len(n_sub), one, workers = workers)

  for (b in seq_len(n_sub)) {
    z <- res[[b]]
    # A forked worker reports a failure as a try-error instead of raising it.
    msg <- .mcg_error_message(z)
    if (!is.null(msg)) z <- list(ok = FALSE, error = msg)
    if (!isTRUE(z$ok)) {
      failures <- rbind(
        failures,
        data.frame(subsample = b, error = z$error, stringsAsFactors = FALSE)
      )
      next
    }
    B_valid <- B_valid + 1L
    if (nrow(z$edges) == 0L) next
    ed <- z$edges
    for (e in seq_len(nrow(ed))) {
      i <- min(ed[e, ]); j <- max(ed[e, ])
      count[i, j] <- count[i, j] + 1
    }
  }
  prob_conservative <- (count + t(count)) / n_sub
  prob <- if (B_valid > 0L) (count + t(count)) / B_valid else
    matrix(NA_real_, R, R)
  dimnames(prob) <- list(g$response_names, g$response_names)
  dimnames(prob_conservative) <- dimnames(prob)
  ut <- if (B_valid > 0L)
    which(upper.tri(prob) & prob >= pi_thr, arr.ind = TRUE) else
      matrix(integer(), ncol = 2L)
  stable <- if (nrow(ut)) unname(cbind(ut[, 1], ut[, 2])) else matrix(integer(0), 0, 2)
  out <- list(selection_prob = prob,
              selection_prob_conservative = prob_conservative,
              stable_edges = stable, B = B, B_valid = B_valid,
              B_failed = n_sub - B_valid, failures = failures,
              pi_thr = pi_thr, method = method, criterion = criterion,
              refit = refit, clustered = !is.null(cluster),
              sampling = sampling, n_subsamples = n_sub)
  if (!is.null(q)) {
    out$q <- q
    out$error_bound <- .stability_bound(q, pi_thr, E_max)
    out$criterion <- NA_character_
  }
  structure(out, class = "mcg_stability")
}

# Theorem 1 of Meinshausen and Buehlmann (2010): E(V) <= q^2 / ((2 pi - 1) p),
# defined for a threshold above one half.
.stability_bound <- function(q, pi_thr, p) {
  if (pi_thr <= 0.5) return(NA_real_)
  q^2 / ((2 * pi_thr - 1) * p)
}

# The first q edges to enter along the L1 path of a subsample, from the penalty
# that zeroes every edge downwards, each fit starting at the last certified
# solution. Edges entering at the same grid point are ranked by the size of
# their coefficient there.
.first_q_edges <- function(fb, q, n_lambda = 50L, lambda_min_ratio = 1e-3) {
  g <- fb$graph
  R <- g$R; n <- g$n; S <- g$Shat
  .check_residual_variance(S, "the subsample path")
  full <- mc_precision_bases(R)$edges
  E <- nrow(full)
  lam_max <- max(n * max(abs(S[upper.tri(S)])), 1e-3)
  lambdas <- exp(seq(log(lam_max), log(lam_max * lambda_min_ratio),
                     length.out = n_lambda))
  entry <- rep(NA_integer_, E)
  size <- rep(0, E)
  start <- NULL
  certified <- 0L
  for (k in seq_along(lambdas)) {
    pen <- mc_graph_penalized(fb, lambdas[k], start = start)
    if (!isTRUE(pen$converged)) next
    certified <- certified + 1L
    start <- pen$eta
    coef <- abs(pen$eta[R + seq_len(E)])
    new <- coef > 0 & is.na(entry)
    entry[new] <- k
    size[new] <- coef[new]
    if (sum(!is.na(entry)) >= q) break
  }
  if (!certified)
    stop("no point of the L1 path was certified for this subsample.")
  entered <- which(!is.na(entry))
  pick <- entered[order(entry[entered], -size[entered])]
  full[utils::head(pick, q), , drop = FALSE]
}

#' @rdname mc_graph_stability
#' @param x an object of class \code{"mcg_stability"}.
#' @export
print.mcg_stability <- function(x, ...) {
  if (!is.null(x$q)) {
    cat(sprintf("Stability selection (first %d edges on the L1 path, B = %d%s, threshold = %.2f)\n",
                x$q, x$B, if (identical(x$sampling, "complementary"))
                  " complementary pairs" else "", x$pi_thr))
    cat(sprintf("Bound on the expected number of false edges: %s\n",
                if (is.na(x$error_bound)) "none, the threshold is not above 1/2"
                else format(signif(x$error_bound, 3))))
  } else
  cat(sprintf("Selection-frequency aggregation (%s, criterion = %s, B = %d, threshold = %.2f)\n",
              x$method, x$criterion, x$B, x$pi_thr))
  cat(sprintf("Half-samples: %s\n",
              if (identical(x$refit, "full"))
                "rows resampled, both stages refitted" else
                  "residuals resampled, first stage held fixed"))
  cat(sprintf("Valid subsamples: %d | failed: %d\n",
              if (is.null(x$B_valid)) x$B else x$B_valid,
              if (is.null(x$B_failed)) 0L else x$B_failed))
  cat(sprintf("Recurrent edges: %d\n", nrow(x$stable_edges)))
  if (nrow(x$stable_edges) > 0) { cat("Edges:\n"); print(x$stable_edges) }
  cat("\nSelection probabilities:\n"); print(round(x$selection_prob, 2))
  invisible(x)
}
