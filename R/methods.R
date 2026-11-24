# methods.R: S3 methods for the mcglm_graph class.
# Part of the mcglmgraph package. License: GPL-3.
#
# An mcglm_graph object inherits from "mcglm", so the mcglm methods (vcov,
# fitted, residuals, anova, confint, gof, pAIC, pBIC, plogLik) are reused
# natively. Only print, summary and coef are customized to surface the graph
# layer; the graph-specific Wald test is provided separately so that the
# inherited anova.mcglm remains available.

#' @describeIn mcglm_graph Print the fit, highlighting the graph layer: the
#'   layer, the number of responses, units and edges, the convergence of the
#'   between-response engine and the between-response covariance
#'   \eqn{\boldsymbol\Sigma_b}. Returns \code{x} invisibly.
#' @param x an object of class \code{"mcglm_graph"}.
#' @param ... for \code{mcglm_graph()} and \code{mcglmgraph()}: further
#'   arguments passed to \code{\link[mcglm]{mcglm}} (such as \code{link},
#'   \code{variance}, \code{power_fixed}, \code{offset}, \code{Ntrial},
#'   \code{weights}). For the methods: ignored by the \code{print} methods and
#'   by \code{plot} with \code{type = "residuals"}; passed on to the inherited
#'   mcglm method by \code{summary} and \code{coef}; passed to
#'   \code{\link{mc_plot_network}} by \code{plot} with
#'   \code{type = "network"}.
#' @export
print.mcglm_graph <- function(x, ...) {
  g <- x$graph
  cat("Multivariate Covariance GLM with graph-structured between-response covariance\n")
  cat(sprintf("Layer: %s | responses: %d | n: %d | edges: %d\n",
              toupper(g$type), g$R, g$n, length(g$edge_pos)))
  cat(sprintf("Between-response engine converged: %s (%d iterations)\n",
              g$converged, g$iter))
  cat("\nBetween-response covariance (Sigma_b):\n")
  print(round(g$Sigma_b, 4))
  cat("\nThis is also a valid 'mcglm' object: use coef(), vcov(), gof(), ...\n")
  invisible(x)
}

#' @describeIn mcglm_graph Summary combining the mcglm summary with the graph
#'   layer (graph parameters, partial correlations and information criteria).
#'   Returns an object of class \code{"summary.mcglm_graph"}, a list with the
#'   first-stage summary in \code{first_stage} and its printed lines in
#'   \code{first_stage_text}; \code{layer}, \code{R}, \code{n} and
#'   \code{n_edges}; \code{converged}, \code{status} and
#'   \code{inference_ok}; the \code{variability} used; the data frame
#'   \code{coefficients} of the graph parameters, with columns
#'   \code{Estimate}, \code{Std.Error}, \code{z} (their ratio) and
#'   \code{p.value} (two-sided, standard normal reference, for the null value
#'   zero of each parameter taken alone), the last three \code{NA} when the
#'   fit is not certified for inference; \code{partial_cor}, the matrix of
#'   partial correlations; and \code{criteria}, the AIC, BIC and QIC of the
#'   graph layer. Its print method shows the first stage and then the graph
#'   layer, and returns its argument invisibly.
#' @param object an object of class \code{"mcglm_graph"}.
#' @param variability for \code{summary}: the covariance of the graph
#'   parameters behind the standard errors, \code{"model"} (default),
#'   \code{"empirical"} for the sandwich or \code{"two_stage"}, as in
#'   \code{\link{mc_vcov_graph}}.
#' @param cluster,small_sample for \code{summary}: passed to
#'   \code{\link{mc_vcov_graph}}; they apply to \code{variability =
#'   "empirical"} and \code{"two_stage"} only and are refused with
#'   \code{"model"}.
#' @export
summary.mcglm_graph <- function(object,
                                variability = c("model", "empirical",
                                                "two_stage"),
                                cluster = NULL, small_sample = FALSE, ...) {
  variability <- match.arg(variability)
  g <- object$graph
  # The mcglm summary prints as it computes; its lines are kept for the print
  # method. A failure there is reported, never hidden.
  txt <- character(0)
  con <- textConnection("txt", "w", local = TRUE)
  sink(con)
  first <- tryCatch(NextMethod(), error = function(e) e,
                    finally = { sink(); close(con) })
  if (inherits(first, "error")) {
    warning("the first-stage summary could not be produced: ",
            conditionMessage(first), call. = FALSE)
    first <- NULL
    txt <- character(0)
  }
  pp <- g$par
  if (isFALSE(g$inference_ok)) {
    tab <- data.frame(Estimate = pp, Std.Error = NA_real_, z = NA_real_,
                      p.value = NA_real_)
  } else {
    V <- mc_vcov_graph(object, type = switch(variability, model = "model",
                                             empirical = "sandwich",
                                             two_stage = "two_stage"),
                       cluster = cluster, small_sample = small_sample)
    se <- sqrt(pmax(diag(as.matrix(V)), 0))
    z  <- ifelse(se > 0, pp / se, NA_real_)
    tab <- data.frame(Estimate = pp, Std.Error = se, z = z,
                      p.value = 2 * stats::pnorm(-abs(z)))
  }
  structure(list(
    first_stage = first, first_stage_text = txt,
    layer = g$type, R = g$R, n = g$n, n_edges = length(g$edge_pos),
    converged = isTRUE(g$converged), status = g$status,
    inference_ok = !isFALSE(g$inference_ok), variability = variability,
    coefficients = tab,
    partial_cor = mc_partial_cor(object),
    criteria = mc_compute_criterion(object, c("aic", "bic", "qic"))
  ), class = "summary.mcglm_graph")
}

#' @describeIn mcglm_graph Print the summary: the first stage as mcglm prints
#'   it, then the graph layer.
#' @export
print.summary.mcglm_graph <- function(x, ...) {
  if (length(x$first_stage_text)) cat(x$first_stage_text, sep = "\n")
  cat("\n=== Graph layer (", toupper(x$layer), ") ===\n", sep = "")
  if (!x$converged)
    cat("The graph layer did not converge (status: ", x$status, ").\n",
        sep = "")
  if (!x$inference_ok) {
    cat("Inferential diagnostics are not certified for this fit; ",
        "standard errors and Wald p-values are suppressed.\n", sep = "")
    cat("Graph parameters:\n")
  } else {
    cat("Graph parameters (", x$variability, " variability):\n", sep = "")
  }
  print(round(x$coefficients, 4))
  cat("\nPartial correlations:\n"); print(round(x$partial_cor, 3))
  cat("\nInformation criteria (graph):\n")
  print(round(x$criteria, 2))
  invisible(x)
}

#' @describeIn mcglm_graph Extract coefficients. With no \code{type}, or with
#'   a \code{type} of \code{\link[mcglm]{coef.mcglm}} (\code{"beta"},
#'   \code{"tau"}, \code{"power"}, \code{"correlation"}), the inherited mcglm
#'   method answers, with its arguments \code{std.error} and \code{response}
#'   passed through \code{...}, and returns its data frame; otherwise
#'   \code{type}, which may be abbreviated, selects the regression, graph or
#'   dispersion coefficients, returned as plain numeric vectors.
#' @export
coef.mcglm_graph <- function(object, type = NULL, ...) {
  # The inherited methods (summary, anova) ask for the mcglm types, so those
  # are answered by coef.mcglm and not refused here.
  mcglm_types <- c("beta", "tau", "power", "correlation")
  if (is.null(type) || (is.character(type) && length(type) &&
                        all(type %in% mcglm_types)))
    return(NextMethod())
  own <- c("regression", "graph", "dispersion")
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      is.na(pmatch(type, own)))
    stop("type must be one of ", paste0("\"", own, "\"", collapse = ", "),
         ", or a type of the mcglm method (",
         paste0("\"", mcglm_types, "\"", collapse = ", "), ").", call. = FALSE)
  switch(own[pmatch(type, own)],
         regression = object$Regression,
         graph      = object$graph$par,
         dispersion = object$Covariance)
}

#' @describeIn mcglm_graph Diagnostic and network plots. With
#'   \code{type = "network"} it draws the dependence network with
#'   \code{\link{mc_plot_network}} and returns invisibly the matrix of edge
#'   weights that function returns; with \code{type = "residuals"} it draws one
#'   histogram of the Pearson residuals per response, restores the graphical
#'   parameters on exit and returns \code{x} invisibly.
#' @param type for \code{coef}: optional, one of \code{"regression"},
#'   \code{"graph"} or \code{"dispersion"}, or a \code{type} accepted by the
#'   mcglm method. For \code{plot}: \code{"network"} for the dependence network
#'   or \code{"residuals"} for residual histograms.
#' @export
plot.mcglm_graph <- function(x, type = c("network", "residuals"), ...) {
  type <- match.arg(type)
  if (type == "network") return(invisible(mc_plot_network(x, ...)))
  g <- x$graph
  op <- graphics::par(mfrow = c(1, g$R)); on.exit(graphics::par(op))
  for (r in seq_len(g$R))
    graphics::hist(g$residuals[, r], main = g$response_names[r],
                   xlab = "Pearson residual", col = "lightblue", breaks = 15)
  invisible(x)
}

#' Covariance matrix of the graph-layer parameters
#'
#' Returns the estimated covariance matrix of the graph-layer parameter vector
#' \eqn{\hat{\boldsymbol\theta}} of a fitted \code{mcglm_graph} object in the
#' three forms the package supports: model-based, empirical sandwich, and
#' sandwich corrected for the estimation of the first stage.
#'
#' @details
#' The graph parameters solve second-order (Pearson) estimating functions
#' \eqn{\mathbf U_n(\boldsymbol\theta)=\sum_i\boldsymbol\psi_i(\boldsymbol\theta)}
#' built on the Pearson residuals of the first-stage MCGLM (Bonat and Jorgensen,
#' 2016). Write \eqn{\mathbf J=-E(\partial\boldsymbol\psi_i/\partial
#' \boldsymbol\theta^\top)} for the sensitivity and
#' \eqn{\mathbf V=\lim_n n^{-1}\sum_i\mathrm{Var}(\boldsymbol\psi_i)} for the
#' variability, both per unit. The Godambe information is
#' \eqn{\mathbf G=\mathbf J^\top\mathbf V^{-1}\mathbf J}. With the number
#' of responses fixed and the number of units \eqn{n\to\infty}, under
#' regularity conditions that include estimating functions unbiased at
#' \eqn{\boldsymbol\theta_0} for every unit or cluster, with independent
#' contributions over units or clusters, a root \eqn{\boldsymbol\theta_0} that
#' is interior and locally unique, a nonsingular \eqn{\mathbf J}, a finite and
#' positive definite \eqn{\mathbf V} with a Lindeberg condition, and uniform
#' control of the fourth moments of the residuals,
#' \deqn{\sqrt n(\hat{\boldsymbol\theta}-\boldsymbol\theta_0)\xrightarrow{d}
#'   \mathcal N\left(\mathbf 0,\;\mathbf G^{-1}\right),\qquad
#'   \mathbf G^{-1}=\mathbf J^{-1}\mathbf V\mathbf J^{-\top}}
#' (Godambe, 1960; Jorgensen and Knudsen, 2004). The matrix returned estimates
#' \eqn{n^{-1}\mathbf G^{-1}}, the covariance of \eqn{\hat{\boldsymbol\theta}}
#' itself. With independent clusters the same statement holds for the cluster
#' sums of the contributions, with the number of clusters in the place of
#' \eqn{n}. The robustness of the sandwich is to the fourth-order structure
#' around the root identified by the estimating functions; it does not correct
#' bias from a misspecified mean or second-order covariance, and it treats the
#' graph as fixed.
#'
#' With \code{type = "model"} the covariance is the model-based
#' \eqn{\mathbf J_n^{-1}}, valid under the working identity
#' \eqn{\mathbf V = \mathbf J} when the relevant second moment is correctly
#' specified and Gaussian fourth-moment identities hold. With
#' \code{type = "sandwich"}, also called the naive sandwich because it takes the
#' Pearson residuals as observed, it is \eqn{\mathbf H_n^{-1}\mathbf V_n\mathbf
#' H_n^{-1}}, with \eqn{\mathbf H_n} the curvature of the section "Scale
#' convention" and \eqn{\mathbf V_n=n\widehat{\mathbf V}} the sum of the outer
#' products of the Pearson contributions, an empirical fourth-moment
#' variability that does not assume that identity. In the DAG layer its block
#' for the structural coefficients of a node is the
#' heteroskedasticity-consistent estimator of White (1980) for the least-squares
#' regression of the centered residual of that node on those of its parents.
#' The unit contribution to the estimating function of parameter
#' \eqn{\theta_d} is
#' \deqn{\psi_{d,i}=-\tfrac12\left\{\operatorname{tr}(\boldsymbol\Sigma_b^{-1}
#'   \dot{\boldsymbol\Sigma}_{b,d})-\mathbf r_i^\top\boldsymbol\Sigma_b^{-1}
#'   \dot{\boldsymbol\Sigma}_{b,d}\boldsymbol\Sigma_b^{-1}\mathbf r_i\right\},}
#' and \eqn{\widehat{\mathbf V}=n^{-1}\sum_i\boldsymbol\psi_i\boldsymbol\psi_i^\top},
#' where \eqn{\mathbf r_i} is the vector of centered Pearson residuals of unit
#' \eqn{i}.
#' This empirical variability requires scores centered in expectation for each
#' independent unit or cluster, or a vanishing aggregate of their squared means.
#' A zero total score at the fitted root does not establish that condition under
#' heterogeneous second-moment misspecification. With \code{cluster}, the scores
#' are summed by cluster before taking outer products; independence of those
#' clusters remains an assumption about the data.
#'
#' @section Scale convention:
#' The model-based sensitivity stored in \code{fit$graph$J} is on the
#' accumulated scale, \eqn{\mathbf J_n=n\mathbf J}, while the variability
#' above is per unit. For a sparse DAG, the empirical sandwich and two-stage
#' correction use the observed negative Hessian \eqn{\mathbf H_n} of the working
#' criterion. Its parent blocks use the residual sample covariance; it can
#' differ from \eqn{\mathbf J_n} if the DAG misses second-moment structure.
#' For the linear-precision GGM they coincide. The empirical matrix returned
#' here is \eqn{\mathbf H_n^{-1}(n\widehat{\mathbf V})\mathbf H_n^{-1}},
#' which is the covariance of \eqn{\hat{\boldsymbol\theta}} and not of
#' \eqn{\sqrt n\,\hat{\boldsymbol\theta}}, so its square roots are directly
#' comparable with the estimates returned by \code{coef(fit, type = "graph")}.
#'
#' @section Propagating the first stage:
#' The first two forms refer to the second stage conditionally on the estimated
#' Pearson residuals: they carry the sampling variability of the estimating
#' functions of the graph layer and not the uncertainty inherited from the
#' first-stage fit. \code{type = "two_stage"} propagates it, with the number of
#' responses fixed, \eqn{n\to\infty} and the regularity conditions of the
#' Details holding for the stacked system. Stacking the two sets of estimating
#' functions gives a block triangular sensitivity, and block inversion gives the
#' influence function of the second stage,
#' \deqn{\mathrm{IF}_i=\mathbf H^{-1}\left(\boldsymbol\psi_i-\mathbf B
#'   \mathbf A^{-1}\boldsymbol\psi_{\alpha,i}\right),\qquad
#'   \mathbf B=-E\!\left(\frac{\partial\boldsymbol\psi_i}
#'   {\partial\boldsymbol\alpha^\top}\right).}
#' Here \eqn{\mathbf H=-E(\partial\boldsymbol\psi_i/\partial
#' \boldsymbol\theta^\top)} is the curvature of the graph block under the
#' distribution of the data, estimated by the observed curvature
#' \eqn{\mathbf H_n/n}, and \eqn{\mathbf A=-E(\partial
#' \boldsymbol\psi_{\alpha,i}/\partial\boldsymbol\alpha^\top)} is the
#' sensitivity of the first stage, whose parameters are
#' \eqn{\boldsymbol\alpha}; both are taken nonsingular. The covariance of
#' \eqn{\hat{\boldsymbol\theta}} is
#' \eqn{n^{-1}\mathbf H^{-1}\widetilde{\mathbf V}\mathbf H^{-\top}} with
#' \deqn{\widetilde{\mathbf V}=\mathbf V_{\theta\theta}
#'   -\mathbf B\mathbf A^{-1}\mathbf V_{\alpha\theta}
#'   -\mathbf V_{\theta\alpha}\mathbf A^{-\top}\mathbf B^\top
#'   +\mathbf B\mathbf A^{-1}\mathbf V_{\alpha\alpha}\mathbf A^{-\top}
#'   \mathbf B^\top,}
#' where \eqn{\mathbf V_{\alpha\alpha}}, \eqn{\mathbf V_{\alpha\theta}} and
#' \eqn{\mathbf V_{\theta\theta}} are the blocks of the variability of the
#' stacked contributions. The first two types keep only
#' \eqn{\mathbf V_{\theta\theta}}; the difference between the corrected and
#' the naive covariance has no determined sign. This is the correction of Murphy
#' and Topel (1985) and, as a stacked system, the two-step result of Newey and
#' McFadden (1994, Theorem 6.1).
#'
#' When every response carries a single dispersion parameter that multiplies
#' its variance function (a Poisson-Tweedie variance only with its power fixed
#' at one), the within-response structure is the identity and the prior weights
#' are constant, the first-stage term can be written, to first order, with
#' second-stage quantities only. This is the \emph{reduced form}, and these
#' three conditions are those of the shortcut below. It also uses that the
#' first stage estimates the correlation between responses without restriction,
#' as \pkg{mcglm} does, and that the first stage acts on the residual
#' covariance as a diagonal rescaling, for which conditions (a) to (c) of the
#' section "What the correction reaches" are sufficient:
#' \deqn{\mathbf B\mathbf A^{-1}\boldsymbol\psi_{\alpha,i}=\mathbf Q
#'   (\mathbf r_i\circ\mathbf r_i-\mathbf m),\qquad
#'   Q_{dr}=\frac{1}{2m_r}\bigl[\mathbf K_b\dot{\boldsymbol\Sigma}_{b,d}
#'   \mathbf K_b\widehat{\mathbf S}_e\bigr]_{rr},\qquad
#'   \dot{\boldsymbol\Sigma}_{b,d}=\partial\boldsymbol\Sigma_b/\partial\theta_d,}
#' with \eqn{\circ} the elementwise product, \eqn{\mathbf K_b=
#' \boldsymbol\Sigma_b^{-1}}, \eqn{\widehat{\mathbf S}_e} the sample covariance
#' of the Pearson residuals and \eqn{m_r} the average of the squared centered
#' Pearson residuals of response \eqn{r}, which, with constant prior weights and
#' residuals that average zero, is the value at which the first-stage Pearson
#' equations fix it. Since \eqn{r_{ri}^2} scales as
#' \eqn{1/\tau_r}, unit \eqn{i} moves \eqn{\log\hat\tau_r} by
#' \eqn{(r_{ri}^2-m_r)/(n m_r)}, which is where the factor \eqn{1/m_r} comes
#' from. Without the bias correction of \pkg{mcglm}
#' (\code{control_algorithm = list(correct = FALSE)}) \eqn{m_r=1}. With it, its
#' default, and unit weights, \eqn{m_r=(n-p_r)/n} when the responses share
#' the design (or there is a single response), \eqn{p_r} the number of
#' regression parameters: exactly for the constant variance function with the
#' identity link and for power two with the log link, and closely for the
#' other multiplicative variance functions. A constant
#' weight \eqn{w_0} turns \eqn{p_r/n} into \eqn{p_r/(nw_0^2)}, and with designs
#' that differ \eqn{m_r} is the observed average, slightly below one. Centering
#' on \eqn{m_r} makes the correction term sum to zero in the sample, and a
#' parameter the first stage fixes exactly, such as the innovation variance of
#' a root node of a DAG with Gaussian intercept-only margins, gets a zero
#' standard error. The \emph{shortcut}, which is what \code{method = "shortcut"}
#' computes, is the reduced form evaluated at the fitted covariance,
#' \eqn{Q_{dr}=[\mathbf K_b\dot{\boldsymbol\Sigma}_{b,d}]_{rr}/(2m_r)}. The
#' two coincide in every row of \eqn{\mathbf Q} whose entries the fitted and
#' the sample covariances agree on: at every support of the GGM layer, and for
#' the parameters of every node of a DAG without immoralities. For the
#' parameters of a collider with non-adjacent parents, and of its descendants,
#' the shortcut is a first-order approximation, across every response column
#' of the parent block, the parents included, that also needs the DAG to be
#' correctly specified; for the parameters of every other node it equals the
#' reduced form. The \emph{full form}
#' evaluates \eqn{\mathbf B} at the data from the first-stage contributions and
#' needs none of these conditions.
#'
#' @section When the shortcut is right:
#' The statements of this section are of first order, with the number of
#' responses fixed, \eqn{n\to\infty} and expectations read as limits of
#' averages over independent units. The reduced form rests on two facts of
#' different kinds. One is algebraic: under the three conditions of the
#' previous section, and with the correlation between responses of the first
#' stage unrestricted, the squared Pearson residual of each response, minus its
#' average over the units, is a linear combination of the centered first-stage
#' contributions of the same unit, whatever the data. That is what
#' \code{\link{mc_two_stage_available}} checks. The other is statistical. Let
#' \eqn{\mathbf M(\boldsymbol\alpha)=E\{\mathbf r_i(\boldsymbol\alpha)
#' \mathbf r_i(\boldsymbol\alpha)^\top\}} be the second moment of the Pearson
#' residuals as a function of the first-stage parameters and \eqn{\mathbf M_0}
#' its value at their limit, which has unit diagonal. Let the Pearson residuals
#' have mean zero at that limit, so that centering them, as the graph layer
#' does, changes neither \eqn{\mathbf M_0} nor its derivatives. The reduced
#' form holds when every first-stage parameter \eqn{\alpha_m} moves
#' \eqn{\mathbf M} as a diagonal rescaling, that is, when for some diagonal
#' matrix \eqn{\boldsymbol\Lambda_m}
#' \deqn{\frac{\partial\mathbf M}{\partial\alpha_m}=
#'   -(\boldsymbol\Lambda_m\mathbf M_0+\mathbf M_0\boldsymbol\Lambda_m)}
#' on the diagonal and in the entries that the matrices \eqn{\mathbf K_b
#' \dot{\boldsymbol\Sigma}_{b,d}\mathbf K_b} read: the edges in the GGM layer,
#' the block of each node and its parents in the DAG layer. Conditions (b) and
#' (c) of the section "What the correction reaches", the latter with
#' \eqn{\mathbf M_0} in the place of \eqn{\boldsymbol\Sigma_b}, give it, with
#' \eqn{\boldsymbol\Lambda_m=\operatorname{Diag}\{E(\partial\log\varsigma_{ri}/
#' \partial\alpha_m)\}}. They are sufficient and can be weakened in two ways:
#' \itemize{
#'   \item (c) can be replaced by \eqn{E(\mathbf r_i\mathbf r_i^\top\mid
#'     \mathbf x_i)=c(\mathbf x_i)\mathbf M_0}, with \eqn{c} a positive function
#'     common to the responses, and \eqn{\boldsymbol\Lambda_m} then carries
#'     \eqn{E\{c(\mathbf x_i)\,\partial\log\varsigma_{ri}/\partial\alpha_m\}}.
#'     When every response has the constant variance function (c) is not
#'     needed in any form, whatever the link, because
#'     \eqn{\partial\log\varsigma_{ri}/\partial\alpha_m} is then the same for
#'     every unit;
#'   \item (b) is used only to give the Pearson residuals mean zero and to make
#'     \eqn{E(r_{si}\,\varsigma_{ri}^{-1}\partial\mu_{ri}/\partial
#'     \boldsymbol\beta_r)=\mathbf 0} for every pair of responses \eqn{(r,s)}.
#'     When all responses have the same design, with an intercept, and the
#'     variance function and the link make that vector the row of the design
#'     times a factor that is the same for every unit, as the constant variance
#'     function with the identity link and the Tweedie variance function of
#'     power fixed at two with the log link do, the quasi-score equations of
#'     the first stage make the Pearson residuals orthogonal to the design in
#'     the sample, and (b) is not needed. In the first of these two cases the
#'     reduced form then holds with no condition on the mean or on the
#'     conditional second moments. In the second, (c), or its weaker form, is
#'     required of the Pearson residuals about the limit of the fitted mean;
#'     when it holds for those about the true mean, a misspecified mean
#'     violates it only at second order in the relative difference between
#'     the two means.
#' }
#' Outside these cases the shortcut can omit a term of first order, whose sign
#' is not determined. That is so for a mean that misses a covariate when the
#' designs of the responses differ and, with a variance function that depends
#' on the mean, for a residual correlation that moves with the covariates; the
#' variances of the partial correlations can then be affected as well. The
#' full form needs none of these conditions: with the observed sensitivity it
#' estimates the influence function of Newey and McFadden (1994, Theorem 6.1)
#' whenever the contributions of each unit, row or cluster, to the estimating
#' functions of both stages have zero mean at their limits and the units are
#' independent. Where both routes apply and the conditions hold they agree to
#' first order; where the conditions fail, the full form is the one to use, and
#' \code{method = "full"} asks for it. A bootstrap that resamples independent
#' units or clusters and refits both stages (\code{\link{mc_graph_bootstrap}})
#' propagates the first stage by construction and is the reference against
#' which both routes can be checked.
#'
#' @section The full form:
#' Outside the shortcut the three pieces are formed from the first stage itself.
#' The fitted MCGLM is rebuilt at its estimate with the routines of
#' \pkg{mcglm}, and split into contributions \eqn{\boldsymbol\psi_{\alpha,u}}
#' of its quasi-score and Pearson functions over units \eqn{u}; \eqn{\mathbf A}
#' is its sensitivity; and \eqn{\mathbf B} is formed by the chain rule on the
#' Pearson residual \eqn{r_{ri}=(y_{ri}-\mu_{ri})/\varsigma_{ri}},
#' \eqn{\varsigma_{ri}^2=C_{(ri),(ri)}},
#' \deqn{-\frac{\partial r_{ri}}{\partial\boldsymbol\alpha}=
#'   \frac{1}{\varsigma_{ri}}\frac{\partial\mu_{ri}}{\partial\boldsymbol\alpha}
#'   +r_{ri}\frac{\partial\log\varsigma_{ri}}{\partial\boldsymbol\alpha},}
#' with no condition on the mean or on the residual covariance. The
#' correlations between responses of the first stage do not enter
#' \eqn{\varsigma_{ri}}, so their block of \eqn{\mathbf B} is zero.
#'
#' The sensitivity \eqn{\mathbf A} of the first stage is where the full form
#' departs from what \pkg{mcglm} reports. By default, \code{sensitivity =
#' "observed"}, it is the Jacobian at the estimate of the first-stage
#' estimating equations as \pkg{mcglm} solved them, its bias correction of the
#' Pearson equations included when the fit used it (\code{fit$con$correct}),
#' formed by central differences of their totals, each one rebuilt with the
#' routines of \pkg{mcglm}. Divided by the number of units it estimates the
#' matrix the influence function asks for whether or not the mean and the
#' variance function of the first stage are correctly specified.
#' With \code{sensitivity = "expected"} \eqn{\mathbf A} takes its model-based
#' form, with blocks \eqn{\mathbf D^\top\mathbf C^{-1}\mathbf W\mathbf D} for
#' the regression parameters,
#' \eqn{\operatorname{tr}(\mathbf C^{-1}\dot{\mathbf C}_\lambda\mathbf W
#' \mathbf C^{-1}\dot{\mathbf C}_{\lambda'})} for the covariance parameters
#' \eqn{\lambda} and the same trace with \eqn{\partial\mathbf C/\partial
#' \beta} in the place of \eqn{\dot{\mathbf C}_{\lambda'}} for the covariance
#' parameters on the regression parameters, where \eqn{\mathbf D} is the
#' derivative of the stacked means, \eqn{\mathbf C} the first-stage covariance
#' and \eqn{\mathbf W} the diagonal matrix of prior weights. It is cheaper and
#' estimates the same limit when the mean and the covariance of the first
#' stage are correctly specified; the corrected covariances of the two forms
#' then differ by a term of order at most \eqn{n^{-1/2}} relative to the naive
#' sandwich. For a parameter whose corrected variance vanishes to first order,
#' as that of the innovation variance of a root node of a DAG does when the
#' first stage fixes the average of the squared residuals of that response and
#' their mean is zero at the limit, the variance reported under either form is
#' of lower order than the naive one, and the two forms can differ by as much
#' as their own size. When the mean or the variance function of the first stage
#' is misspecified the expected form in general no longer
#' estimates \eqn{\mathbf A}, the difference between the two corrected
#' covariances can be of first order, and the observed sensitivity is the one
#' to use. The observed form inherits the conditioning of the first stage:
#' when that stage is weakly identified, as a Poisson-Tweedie fit whose
#' estimated power and dispersion trade off against each other can be, its
#' Jacobian is close to singular and the propagated variance is numerically
#' unstable. The attribute
#' \code{"sensitivity_rcond"} reports the reciprocal condition number of the
#' sensitivity used, and comparing the two forms on the same fit is a direct
#' check. Either way the
#' block of the covariance parameters on the regression parameters is kept:
#' the \code{joint_inv_sensitivity} of \pkg{mcglm} has it at zero, it is not
#' zero when the variance depends on the mean, and without it the regression
#' columns of \eqn{\mathbf B\mathbf A^{-1}} do not vanish under the conditions
#' of the reduced form, where they must. The unit
#' contributions of the first stage are centered: the bias correction is one
#' term for the whole sample rather than a sum over units, so the contributions
#' without it sum to minus it at the estimate, and a term common to all units
#' does not enter the variability.
#'
#' The first stage is a sum over units only when its covariance is block
#' diagonal over them. With an identity within-response structure the units
#' are the rows, or the clusters when \code{cluster} is given. With any other
#' structure \code{cluster} is required and has to nest it: every matrix of the
#' matrix linear predictor must vanish between distinct clusters, which is
#' checked.
#'
#' \code{method = "auto"} takes the shortcut when it applies and the full form
#' otherwise, with one exception in the DAG layer. The shortcut evaluates the
#' reduced form at the fitted covariance, which is exact for the parameters of
#' every node of a DAG without immoralities but, for the parameters of a
#' collider with non-adjacent parents and of its descendants, only a
#' first-order approximation under correct specification; so when the fitted
#' DAG has an immorality \code{"auto"} takes the full form if
#' it is available, and otherwise keeps the shortcut with a warning. An explicit
#' \code{method = "shortcut"} is honored on such a DAG without a warning. On
#' the fits the shortcut covers the two agree to first order when the
#' conditions of the section "When the shortcut is right" hold.
#' Use \code{\link{mc_two_stage_available}} to see which route applies; a fit
#' neither covers raises an error rather than returning the uncorrected matrix
#' under the corrected name.
#' @section What the correction reaches:
#' How much of the propagation reaches a given quantity is settled, for the
#' undirected layer, by a diagonal-rescaling argument, under three conditions on
#' the first stage, which are sufficient and, as the section "When the shortcut
#' is right" says, not always necessary:
#' \describe{
#'   \item{(a)}{the Pearson residual has the form
#'     \eqn{r_{ri}=(y_{ri}-\mu_{ri})/\varsigma_{ri}}, with the mean and the
#'     scale depending only on the first-stage parameters and the covariates,
#'     which is what \pkg{mcglm} returns;}
#'   \item{(b)}{the mean is correctly specified,
#'     \eqn{E(r_{ri}\mid\mathbf x_i)=0};}
#'   \item{(c)}{the conditional residual covariance does not depend on the
#'     covariates, \eqn{E(\mathbf r_i\mathbf r_i^\top\mid\mathbf x_i)=
#'     \boldsymbol\Sigma_b}.}
#' }
#' Under the three, the propagation term is the coordinate vector of
#' \eqn{\mathbf K_b\boldsymbol\Lambda_i+\boldsymbol\Lambda_i\mathbf K_b}, with
#' \eqn{\boldsymbol\Lambda_i} diagonal: the first stage acts on the residual
#' covariance, to first order, as a diagonal rescaling. The propagation vanishes
#' on every edge whose precision entry is zero and is annihilated by the gradient
#' of every partial correlation, so for those quantities the sandwich and the
#' corrected form agree to first order. The diagonal coefficients and the fitted
#' \eqn{\boldsymbol\Sigma_b} are not invariant and do need the correction.
#' The DAG layer has the same reading, for any fixed DAG, with the second
#' moment of the residuals in the place of \eqn{\boldsymbol\Sigma_b} in (c):
#' the propagation vanishes on every structural coefficient whose target value
#' is zero and on every function of the parameters that is invariant to a
#' diagonal rescaling of the residuals, such as the partial correlations, while
#' the nonzero structural coefficients and the innovation variances need the
#' correction.
#'
#' For the shortcut the agreement on the partial correlations is exact and
#' holds whatever the data: the gradients of the partial correlations
#' annihilate \eqn{\mathbf H^{-1}\mathbf Q} algebraically, since
#' \eqn{2\mathbf H^{-1}\mathbf Q_{\cdot r}} is the direction of a rescaling, so
#' \code{type = "sandwich"} and the shortcut give them the same variance. That
#' says nothing about whether the variance is right. When the first stage does
#' not act as a rescaling, as can happen with a variance function that depends
#' on the mean and a residual correlation that moves with the covariates, both
#' can miss the same first-order term on the partial correlations, and only
#' the full form accounts for it.
#'
#' The model-based form \eqn{\mathbf J_n^{-1}} agrees with the empirical
#' forms only when, in addition, \eqn{\mathbf H=\mathbf J} and
#' \eqn{\mathbf V=\mathbf J}; the latter is the correctly specified Gaussian
#' fourth-moment regime. Outside it the difference is not a matter of the first
#' stage, and (a) to (c) do not remove it. The package checks none of (a) to
#' (c): they are statements about the process that generated the data. Applied
#' to a support chosen by a search over graphs, every form remains descriptive,
#' since the selection is not incorporated.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param type \code{"model"} for \eqn{\mathbf J_n^{-1}}, \code{"sandwich"} for
#'   the empirical sandwich, or \code{"two_stage"} for the sandwich corrected
#'   for the estimation of the residuals.
#' @param cluster optional vector with one identifier per observational row.
#'   With \code{type = "sandwich"} or \code{"two_stage"}, the contributions are
#'   summed within clusters before the variability is formed, which is the form
#'   to use when the unit of dependence is not the row of the data matrix. It
#'   must have no missing value and at least two distinct values; the default
#'   \code{NULL} takes the rows as the units. It is refused with
#'   \code{type = "model"}.
#' @param small_sample logical, default \code{FALSE}; multiply the clustered
#'   variability by \eqn{G/(G-1)}, with \eqn{G} the number of clusters. It has
#'   no effect without \code{cluster} and is refused with
#'   \code{type = "model"}.
#' @param method for \code{type = "two_stage"}, \code{"auto"} (the shortcut
#'   when it applies, the full form otherwise; on a DAG with an immorality, the
#'   full form whenever it is available), \code{"shortcut"} or \code{"full"}.
#'   It has no meaning for the other two types.
#' @param sensitivity for the full form, \code{"observed"} (default) or
#'   \code{"expected"}, the sensitivity of the first stage described in the
#'   section "The full form".
#' @return A \eqn{q_{\mathcal G} \times q_{\mathcal G}} symmetric matrix, with
#'   \eqn{q_{\mathcal G}} the number of graph parameters and their names in its
#'   dimnames, carrying the attributes \code{"type"}, \code{"variability"},
#'   \code{"edge_pos"} (positions of the edge parameters), \code{"n_clusters"},
#'   \code{"small_sample"}, which is \code{TRUE} exactly when the factor
#'   \eqn{G/(G-1)} was applied, and the inversion diagnostics
#'   \code{"inversion_status"} and \code{"inversion_ridge"} of the sensitivity
#'   solve (with \code{type = "model"} also \code{"inversion_rcond"}, the
#'   reciprocal condition number met by that solve). \code{"n_clusters"} is
#'   \code{NA} without \code{cluster}. With \code{type = "two_stage"} it also
#'   carries
#'   \code{"two_stage_method"}, the route taken, and either the factor
#'   \code{"Q"} the shortcut applies, \eqn{[\mathbf K_b\dot{\boldsymbol
#'   \Sigma}_{b,d}]_{rr}/(2m_r)}, with the averages \eqn{m_r} in
#'   \code{"center"}, or, for the full form, the cross sensitivity
#'   \code{"B"}, the first-stage sensitivity \code{"A"}, both on the
#'   accumulated scale, which \code{"sensitivity"} was used and its reciprocal
#'   condition number \code{"sensitivity_rcond"}. The function refuses fits
#'   whose numerical diagnostics mark graph-layer inference as uncertified.
#' @seealso \code{\link{mc_test_independence}} for the joint Wald test, which
#'   consumes this matrix, \code{\link{mc_wald_graph}} for general linear
#'   hypotheses, \code{\link{mc_two_stage_available}} for the route the
#'   correction takes, \code{\link{mc_graph_bootstrap}} for resampling, and
#'   \code{\link{mc_partial_cor_delta}} for the delta-method inference on the
#'   partial correlations.
#' @references
#' Bonat, W. H. and Jorgensen, B. (2016). Multivariate covariance generalized
#' linear models. \emph{Journal of the Royal Statistical Society, Series C}
#' 65(5), 649--675. \doi{10.1111/rssc.12145}
#'
#' Godambe, V. P. (1960). An optimum property of regular maximum likelihood
#' estimation. \emph{Annals of Mathematical Statistics} 31(4), 1208--1211.
#' \doi{10.1214/aoms/1177705693}
#'
#' Jorgensen, B. and Knudsen, S. J. (2004). Parameter orthogonality and bias
#' adjustment for estimating functions. \emph{Scandinavian Journal of
#' Statistics} 31(1), 93--114. \doi{10.1111/j.1467-9469.2004.00375.x}
#'
#' Murphy, K. M. and Topel, R. H. (1985). Estimation and inference in two-step
#' econometric models. \emph{Journal of Business and Economic Statistics} 3(4),
#' 370--379. \doi{10.1080/07350015.1985.10509471}
#'
#' Newey, W. K. and McFadden, D. (1994). Large sample estimation and hypothesis
#' testing. In Engle, R. F. and McFadden, D. L. (eds), \emph{Handbook of
#' Econometrics}, vol. 4, 2111--2245. Elsevier.
#' \doi{10.1016/S1573-4412(05)80005-4}
#'
#' White, H. (1980). A heteroskedasticity-consistent covariance matrix
#' estimator and a direct test for heteroskedasticity. \emph{Econometrica}
#' 48(4), 817--838. \doi{10.2307/1912934}
#' @family inference on the graph layer
#' @examples
#' \donttest{
#' set.seed(1); n <- 300
#' K <- diag(3) * 2; K[1, 2] <- K[2, 1] <- -0.7
#' Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
#' d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = rnorm(n))
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' V <- mc_vcov_graph(f, type = "sandwich")
#' sqrt(diag(V))[attr(V, "edge_pos")]      # standard errors of the edges
#' }
#' @export
mc_vcov_graph <- function(fit, type = c("model", "sandwich", "two_stage"),
                          cluster = NULL, small_sample = FALSE,
                          method = c("auto", "shortcut", "full"),
                          sensitivity = c("observed", "expected")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  type <- match.arg(type)
  method <- match.arg(method)
  sensitivity <- match.arg(sensitivity)
  if (type != "two_stage" && method != "auto")
    stop("method chooses the route of type = \"two_stage\" and has no ",
         "meaning for type = \"", type, "\".", call. = FALSE)
  if (length(small_sample) != 1L || !is.logical(small_sample) ||
      is.na(small_sample))
    stop("small_sample must be a single TRUE or FALSE.", call. = FALSE)
  if (type == "model") {
    if (!is.null(cluster))
      stop("cluster applies to an empirical variability only; the model-based ",
           "covariance has no variability term to aggregate.", call. = FALSE)
    if (isTRUE(small_sample))
      stop("small_sample applies to an empirical variability only.",
           call. = FALSE)
  }
  .assert_graph_inference(fit, "Graph-layer variance")
  g <- fit$graph
  V <- switch(type,
    model = g$vcov,
    sandwich = .robust_vcov(fit, cluster = cluster,
                            small_sample = small_sample),
    two_stage = .two_stage_vcov(fit, cluster = cluster,
                                small_sample = small_sample, method = method,
                                sensitivity = sensitivity))
  nms <- names(g$par)
  if (length(nms) == nrow(V)) dimnames(V) <- list(nms, nms)
  attr(V, "type") <- type
  attr(V, "variability") <- if (type == "model") "model" else "empirical"
  attr(V, "edge_pos") <- g$edge_pos
  attr(V, "n_clusters") <- if (is.null(cluster)) NA_integer_ else
    nlevels(as.factor(cluster))
  # Both empirical forms apply G/(G-1) to the clustered variability, and
  # neither applies it without clusters, so the flag records exactly that.
  attr(V, "small_sample") <- isTRUE(small_sample) && !is.null(cluster) &&
    type %in% c("sandwich", "two_stage")
  V
}

#' Wald test of between-response conditional association
#'
#' Tests the null hypothesis that all edge parameters of the fitted graph layer
#' are zero, under which the between-response covariance
#' \eqn{\boldsymbol\Sigma_b} of the Pearson residuals is diagonal: the absence
#' of second-order conditional association among the responses given the
#' covariates. Under joint normality of the analyzed residuals this coincides
#' with mutual conditional independence; in Gaussian-copula models the same
#' reading requires a precision matrix fitted on the latent Gaussian scale.
#'
#' @details
#' It is \code{\link{mc_wald_graph}} with the contrast \eqn{\mathbf L} that
#' selects the \eqn{q} edge parameters \eqn{\boldsymbol\theta_E} (the
#' off-diagonal precision parameters of the GGM layer, the structural
#' coefficients of the DAG layer) and \eqn{\mathbf c=\mathbf 0}:
#' \deqn{W=\hat{\boldsymbol\theta}_E^\top\,\widehat{\mathrm{Var}}(
#'   \hat{\boldsymbol\theta}_E)^{-1}\,\hat{\boldsymbol\theta}_E,}
#' with the covariance of \code{\link{mc_vcov_graph}} selected by \code{type}.
#' Under the null, with the number of responses fixed and
#' \eqn{n\to\infty}, \eqn{W} is asymptotically \eqn{\chi^2_q}, \eqn{q} being the
#' number of edges; the p-value is the upper tail of that distribution. The
#' default covariance is the model-based one, for which that reference also
#' requires the working identity \eqn{\mathbf V=\mathbf J} (a correctly
#' specified second moment and Gaussian fourth moments); \code{type =
#' "sandwich"} does not use it, and \code{type = "two_stage"} also accounts for
#' the estimation of the first stage. The edges tested are those of the support
#' of the fit, and the
#' graph is treated as fixed: for a support chosen by a search over graphs, the
#' result is descriptive. It is a standalone function so that the inherited
#' \code{anova.mcglm} method remains available on the object. The function
#' refuses fits whose numerical diagnostics mark graph-layer inference as
#' uncertified.
#'
#' @param fit an object of class \code{"mcglm_graph"} with at least one edge.
#' @param robust logical, default \code{FALSE}; if \code{TRUE}, the test uses
#'   the robust sandwich covariance \eqn{\mathbf H_n^{-1}\mathbf V_n\mathbf
#'   H_n^{-1}}, with \eqn{\mathbf H_n} the observed curvature of the graph
#'   block, built from the empirical fourth-moment variability, instead of the
#'   model-based \eqn{\mathbf J_n^{-1}} (which assumes
#'   \eqn{\mathbf V=\mathbf J}).
#'   It is equivalent to \code{type = "sandwich"}, which also offers
#'   \code{"two_stage"}.
#' @param cluster optional vector with one identifier per observational row,
#'   identifying independent clusters. Requires an empirical covariance
#'   (\code{robust = TRUE} or an empirical \code{type}): Pearson contributions
#'   are then summed within clusters before the empirical variability is formed.
#' @param small_sample logical, default \code{FALSE}; multiply the cluster
#'   variability by \eqn{G/(G-1)}, where \eqn{G} is the number of clusters.
#'   Requires an empirical covariance and has no effect without \code{cluster}.
#' @param type optional covariance, \code{"model"}, \code{"sandwich"} or
#'   \code{"two_stage"}, as in \code{\link{mc_vcov_graph}}. When given it takes
#'   precedence over \code{robust}, which must then not contradict it. The
#'   default \code{NULL} uses \code{"sandwich"} if \code{robust} is
#'   \code{TRUE} and \code{"model"} otherwise.
#' @return An object of class \code{"htest"}: the Wald \code{statistic}, the
#'   degrees of freedom in \code{parameter} (the number of edges), the
#'   \code{p.value}, and the \code{method} and \code{data.name} labels.
#' @seealso \code{\link{mc_wald_graph}} for other linear hypotheses,
#'   \code{\link{mc_test_adequacy}} for the adequacy of a selected support, and
#'   \code{\link{mc_vcov_graph}} for the covariance forms.
#' @family inference on the graph layer
#' @references
#' Godambe, V. P. (1960). An optimum property of regular maximum likelihood
#' estimation. \emph{Annals of Mathematical Statistics} 31(4), 1208--1211.
#' \doi{10.1214/aoms/1177705693}
#' @examples
#' \donttest{
#' set.seed(1); n <- 300
#' K <- diag(3) * 2; K[1,2] <- K[2,1] <- -0.7
#' Y <- matrix(rnorm(n * 3), n) %*% chol(solve(K))
#' d <- data.frame(y1 = Y[,1], y2 = Y[,2], y3 = Y[,3], x = rnorm(n))
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_test_independence(f)
#' mc_test_independence(f, robust = TRUE)
#' mc_test_independence(f, type = "two_stage")
#' }
#' @export
mc_test_independence <- function(fit, robust = FALSE, cluster = NULL,
                                 small_sample = FALSE, type = NULL) {
  stopifnot(inherits(fit, "mcglm_graph"))
  .assert_graph_inference(fit, "Wald testing")
  g <- fit$graph
  ep <- g$edge_pos
  if (length(ep) == 0L) stop("The model has no edge parameters to test.")
  if (is.null(type)) {
    type <- if (isTRUE(robust)) "sandwich" else "model"
  } else {
    type <- match.arg(type, c("model", "sandwich", "two_stage"))
    if (!missing(robust) && !identical(isTRUE(robust), type != "model"))
      stop("robust and type disagree; give one of them.", call. = FALSE)
  }
  if (type == "model" && !is.null(cluster))
    stop("cluster is available only for robust = TRUE or an empirical type.")
  if (type == "model" && isTRUE(small_sample))
    stop("small_sample is available only for robust = TRUE or an empirical type.")
  # One source of truth: the joint test is the general Wald test with the
  # contrast that selects the edges, so the two cannot drift apart. The
  # selector copies the entries of the estimate and of its covariance exactly,
  # so the statistic is the plain quadratic form in the edges, to the last bit.
  q <- length(ep)
  L <- matrix(0, q, length(g$par))
  L[cbind(seq_len(q), ep)] <- 1
  w <- mc_wald_graph(fit, L, 0, type = type, cluster = cluster,
                     small_sample = small_sample)
  structure(list(statistic = w$statistic,
                 parameter = c(df = q),
                 p.value = w$p.value,
                 method = "Wald test of between-response conditional association",
                 data.name = sprintf("%s layer with %d edges",
                                     toupper(g$type), q)),
            class = "htest")
}

.assert_graph_inference <- function(fit, operation = "Graph inference") {
  g <- fit$graph
  if (isFALSE(g$inference_ok)) {
    detail <- c(
      if (!is.null(g$status)) paste0("fit status=", g$status),
      if (!is.null(g$inversion_status))
        paste0("inversion=", g$inversion_status),
      if (!is.null(g$rcond_J) && is.finite(g$rcond_J))
        paste0("rcond_J=", format(g$rcond_J, digits = 4))
    )
    stop(operation, " is not certified for this fit",
         if (length(detail)) paste0(" (", paste(detail, collapse = "; "), ")") else "",
         ". Refit an identifiable, converged unpenalized support before inference.",
         call. = FALSE)
  }
  invisible(TRUE)
}
