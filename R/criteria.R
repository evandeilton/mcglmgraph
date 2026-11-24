# criteria.R: Gaussian working criterion and information criteria.
# Part of the mcglmgraph package. License: GPL-3.

#' Gaussian working log-criterion of the between-response fit
#'
#' Evaluates the Gaussian working log-criterion of the Pearson residuals under
#' the estimated between-response covariance,
#' \deqn{\ell^p_e = -\tfrac{n}{2}\bigl[R\log(2\pi) + \log|\hat\Sigma_b| +
#' \mathrm{tr}(\hat\Sigma_b^{-1}\hat S_e)\bigr],}
#' where \eqn{n} is the number of units, \eqn{R} the number of responses,
#' \eqn{\hat\Sigma_b} the fitted between-response covariance and
#' \eqn{\hat S_e} the sample covariance (centered, divided by \eqn{n}) of the
#' Pearson residuals of the first stage. This scalar criterion generates the
#' Pearson equations used by the graph layer, and it is the quantity from which
#' \code{\link{mc_compute_criterion}} builds the information criteria. It is
#' not presented as a composite likelihood or as a fully specified likelihood
#' for the marginal MCGLM (multivariate covariance generalized linear model)
#' responses.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @return A single number, the Gaussian working log-criterion.
#' @seealso \code{\link{mc_compute_criterion}}, \code{\link{mc_graph_cv}}.
#' @family information criteria
#' @examples
#' \donttest{
#' set.seed(1); n <- 150; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_graph_ploglik(f)
#' }
#' @export
mc_graph_ploglik <- function(fit) {
  stopifnot(inherits(fit, "mcglm_graph"))
  g <- fit$graph
  Sinv <- solve(g$Sigma_b)
  ld <- determinant(g$Sigma_b, logarithm = TRUE)$modulus
  -0.5 * g$n * (g$R * log(2 * pi) + as.numeric(ld) + sum(Sinv * g$Shat))
}

#' Information criteria for a graph-structured fit
#'
#' Computes working AIC and BIC, the quasi-likelihood information criterion
#' (QIC), a Godambe-adjusted working BIC (GBIC) and the extended BIC (EBIC), all
#' built on the Gaussian working log-criterion \eqn{\ell^p_e} of
#' \code{\link{mc_graph_ploglik}}. With \eqn{p} the number of free parameters
#' of the between-response covariance (diagonal and edge parameters) and
#' \eqn{n} the number of units (rows of the residual matrix),
#' \deqn{\mathrm{AIC} = -2\ell^p_e + 2p,\qquad
#'       \mathrm{BIC} = -2\ell^p_e + \log(n)\,p,}
#' \deqn{\mathrm{QIC} = -2\ell^p_e + 2\,d_{\mathrm{eff}},\qquad
#'       \mathrm{GBIC} = -2\ell^p_e + \log(n)\,d_{\mathrm{eff}},}
#' with the effective dimension \eqn{d_{\mathrm{eff}}=\mathrm{tr}(H^{-1}V)}.
#' Here \eqn{H} is the curvature (sensitivity) and \eqn{V} the empirical
#' variability of the Pearson estimating functions, both per unit, so that the
#' trace is invariant to the scale. Under the Gaussian working model, when
#' \eqn{H = V = J}, the effective dimension reduces to the parameter count and
#' QIC and GBIC reduce to AIC and BIC. \code{"clbic"} is an exact alias of
#' \code{"gbic"}; it does not assert a composite-likelihood factorization.
#' Smaller values indicate a better working model.
#'
#' With \code{cluster}, the logarithmic multiplier \eqn{\log(n)} of BIC, GBIC
#' and EBIC becomes \eqn{\log(G)}, with \eqn{G} the number of independent
#' clusters, and the variability \eqn{V} is formed from the contributions
#' summed within clusters. The criterion \eqn{\ell^p_e}, which remains a sum
#' over the \eqn{n} rows, the parameter count \eqn{p}, AIC and the multiplier 2
#' of QIC are not altered.
#'
#' @section The extended BIC:
#' \code{"ebic"} is the form of Foygel and Drton (2010) for Gaussian graphical
#' models, \eqn{\mathrm{EBIC}_\gamma=\mathrm{BIC}+4\gamma|E|\log R}, with
#' \eqn{|E|} the number of edge parameters of the fit and \eqn{R} the number of
#' responses. The diagonal parameters are common to every candidate support, so
#' they enter through the BIC alone and do not alter differences between
#' supports. \eqn{\gamma=0} gives the BIC back. It is offered as a working
#' criterion: its consistency for the graph is established for Gaussian
#' likelihoods, and it does not carry over, as such, to the working criterion of
#' the residual layer.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param which character vector among \code{"aic"}, \code{"bic"},
#'   \code{"qic"}, \code{"gbic"}, \code{"clbic"} and \code{"ebic"}.
#'   \code{"clbic"} is an alias of \code{"gbic"}. When left out, the first
#'   five are returned; \code{"ebic"} is computed only when requested.
#' @param ebic_gamma the parameter \eqn{\gamma\ge0} of the extended BIC
#'   (default 0.5), used only when \code{"ebic"} is requested.
#' @param cluster optional independent-cluster identifier, one per residual row
#'   (at least two distinct values, no missing values). For QIC and GBIC,
#'   Pearson contributions are summed within clusters before forming their
#'   empirical variability. BIC, GBIC and EBIC use the logarithm of the number
#'   of independent clusters, instead of the number of rows, in their penalty.
#'   Supplying the identifier does not by itself establish independence of the
#'   clusters or centering of their scores.
#' @param small_sample logical; multiply the clustered variability by
#'   \eqn{G/(G-1)}, with \eqn{G} the number of clusters; ignored without
#'   \code{cluster}.
#' @section Scope of the Godambe penalties:
#' For a sparse DAG, \eqn{H} is the observed curvature of the residual working
#' criterion. It can differ from the model-based Pearson sensitivity when the
#' fitted covariance is misspecified. For the linear-precision GGM they agree.
#' The QIC optimism calculation assumes genuinely known first-stage marginal
#' parameters. With marginals estimated on the same data, QIC remains a working
#' criterion; the reported penalty is not an exact two-stage optimism correction.
#' The empirical variability is
#' \eqn{V=n^{-1}\sum_i\psi_i\psi_i^\top}, with \eqn{\psi_i} the Pearson
#' contribution of unit \eqn{i} (summed within clusters when \code{cluster} is
#' given, then divided by \eqn{n}). It estimates the variance of the scores
#' only when their means vanish per independent unit or cluster (or the
#' aggregate contribution of squared means vanishes). A pooled zero score is
#' insufficient.
#' @return A named numeric vector of the requested criteria, carrying the
#'   attributes \code{"criterion_status"} (\code{"godambe"},
#'   \code{"godambe_pseudoinverse"}, \code{"parameter_count_fallback"} or
#'   \code{"not_required"}) and \code{"effective_dimension"} (the value of
#'   \eqn{d_{\mathrm{eff}}}, \code{NA} when no Godambe-adjusted criterion was
#'   requested). A status other than \code{"godambe"} means the effective
#'   dimension could not be computed from a well-conditioned sensitivity (a
#'   spectral pseudo-inverse or, last, the parameter count was used instead)
#'   and must be reported as such.
#' @references Foygel, R. and Drton, M. (2010). Extended Bayesian information
#'   criteria for Gaussian graphical models. \emph{Advances in Neural
#'   Information Processing Systems} 23, 604--612.
#'
#'   Pan, W. (2001). Akaike's information criterion in generalized estimating
#'   equations. \emph{Biometrics} 57(1), 120--125.
#'   \doi{10.1111/j.0006-341X.2001.00120.x}
#'
#'   Gao, X. and Song, P. X.-K. (2010). Composite likelihood Bayesian
#'   information criteria for model selection in high-dimensional data.
#'   \emph{Journal of the American Statistical Association} 105(492),
#'   1531--1540. \doi{10.1198/jasa.2010.tm09414}
#'
#'   Varin, C. and Vidoni, P. (2005). A note on composite likelihood inference
#'   and model selection. \emph{Biometrika} 92(3), 519--528.
#'   \doi{10.1093/biomet/92.3.519}
#' @seealso \code{\link{mc_graph_ploglik}}, \code{\link{mc_graph_stepwise}},
#'   \code{\link{mc_select_lambda}}, \code{\link{mc_graph_cv}}.
#' @family information criteria
#' @examples
#' \donttest{
#' set.seed(1); n <- 150; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_compute_criterion(f, c("aic", "bic", "qic", "gbic"))
#' mc_compute_criterion(f, "ebic", ebic_gamma = 0.5)
#' }
#' @export
mc_compute_criterion <- function(
    fit, which = c("aic", "bic", "qic", "gbic", "clbic", "ebic"),
    ebic_gamma = 0.5, cluster = NULL, small_sample = FALSE) {
  which <- if (missing(which)) c("aic", "bic", "qic", "gbic", "clbic") else
    match.arg(which, several.ok = TRUE)
  if ("ebic" %in% which &&
      (length(ebic_gamma) != 1L || !is.finite(ebic_gamma) || ebic_gamma < 0))
    stop("ebic_gamma must be one finite number of at least zero.",
         call. = FALSE)
  ll <- mc_graph_ploglik(fit)
  g  <- fit$graph
  n_independent <- if (is.null(cluster)) g$n else
    .two_stage_units(fit, cluster)$n_units
  p  <- g$npar
  out <- numeric(0)
  if ("aic" %in% which) out["aic"] <- -2 * ll + 2 * p
  if ("bic" %in% which) out["bic"] <- -2 * ll + log(n_independent) * p
  if ("ebic" %in% which)
    out["ebic"] <- -2 * ll + log(n_independent) * p +
      4 * ebic_gamma * length(g$edge_pos) * log(g$R)
  status <- "not_required"
  pen <- NA_real_
  if (any(c("qic", "gbic", "clbic") %in% which)) {
    Vemp <- .empirical_variability(fit, cluster = cluster,
                                  small_sample = small_sample)
    Hunit <- .graph_inference_bread(fit) / g$n
    eff <- .effective_dimension(Hunit, Vemp, p)
    pen <- eff$value
    status <- eff$status
    if ("qic"   %in% which) out["qic"]   <- -2 * ll + 2 * pen
    gbic <- -2 * ll + log(n_independent) * pen
    if ("gbic"  %in% which) out["gbic"]  <- gbic
    if ("clbic" %in% which) out["clbic"] <- gbic
  }
  attr(out, "criterion_status") <- status
  attr(out, "effective_dimension") <- pen
  out
}

# Effective dimension tr(H^{-1} V) of the Godambe-adjusted criteria.
#
# A direct solve of H can fail or be numerically meaningless when the
# curvature is ill-conditioned. Substituting the raw parameter count silently
# would turn GBIC into BIC without a trace, so the fallback is explicit and
# graded: an exact solve when H is well conditioned, otherwise a
# spectral pseudo-inverse that discards eigenvalues below
# max(eigen) * .Machine$double.eps^(1/2), and only if that also fails does the
# parameter count stand in. The status travels with the result.
.effective_dimension <- function(Junit, Vemp, p) {
  sym <- 0.5 * (Junit + t(Junit))
  rc <- tryCatch(rcond(sym), error = function(e) 0)
  if (is.finite(rc) && rc > .Machine$double.eps^(2/3)) {
    val <- tryCatch(sum(diag(solve(sym, Vemp))), error = function(e) NA_real_)
    if (is.finite(val)) return(list(value = val, status = "godambe"))
  }
  ev <- tryCatch(eigen(sym, symmetric = TRUE), error = function(e) NULL)
  if (!is.null(ev)) {
    keep <- ev$values > max(ev$values) * sqrt(.Machine$double.eps)
    if (any(keep)) {
      Jplus <- ev$vectors[, keep, drop = FALSE] %*%
        diag(1 / ev$values[keep], sum(keep)) %*%
        t(ev$vectors[, keep, drop = FALSE])
      val <- sum(diag(Jplus %*% Vemp))
      if (is.finite(val))
        return(list(value = val, status = "godambe_pseudoinverse"))
    }
  }
  list(value = p, status = "parameter_count_fallback")
}

# The DAG criterion separates by node. Its observed negative Hessian has
# H_{gamma_j,gamma_j} = n S[P_j,P_j] / delta_j, with S the residual sample
# covariance, not the covariance implied by a sparse fitted DAG. The remaining
# entries below retain the exact finite-sample curvature away from the root.
# In the linear-precision GGM the Hessian is data-independent and equals g$J.
.graph_inference_bread <- function(fit) {
  g <- fit$graph
  if (g$type != "dag") return(g$J)
  meta <- g$meta
  R <- g$R; E <- meta$E; n <- g$n
  S <- g$Shat[meta$ordering, meta$ordering, drop = FALSE]
  gamma <- as.numeric(g$par[seq_len(E)])
  delta <- as.numeric(g$par[E + seq_len(R)])
  H <- matrix(0, E + R, E + R)
  for (j in seq_len(R)) {
    idx <- which(meta$edges[, 2L] == j)
    q <- S[j, j]
    if (length(idx)) {
      P <- meta$edges[idx, 1L]
      A <- S[P, P, drop = FALSE]
      b <- S[P, j]
      gj <- gamma[idx]
      residual <- as.numeric(b - A %*% gj)
      q <- q - 2 * sum(gj * b) + drop(crossprod(gj, A %*% gj))
      H[idx, idx] <- n * A / delta[j]
      H[idx, E + j] <- n * residual / delta[j]^2
      H[E + j, idx] <- H[idx, E + j]
    }
    H[E + j, E + j] <- n * (q / delta[j]^3 - 1 / (2 * delta[j]^2))
  }
  H
}

# Robust covariance of the graph parameters using the empirical
# fourth-moment variability: Var = H^{-1} (n V) H^{-1}. It reduces to the
# model-based J^{-1} only when H = V = J.
.robust_vcov <- function(fit, cluster = NULL, small_sample = FALSE) {
  g <- fit$graph
  Vemp <- .empirical_variability(
    fit, cluster = cluster, small_sample = small_sample
  )
  # The fit carries the conditioning tolerance it was solved under. Using the
  # default here instead would let the sandwich accept a sensitivity that the
  # fitting engine itself had refused, or refuse one it had accepted.
  rct <- g$rcond_tol
  if (is.null(rct) || !is.finite(rct)) rct <- .Machine$double.eps^(2 / 3)
  Jinv <- .safe_solve(.graph_inference_bread(fit), rcond_tol = rct)
  M <- Jinv %*% (g$n * Vemp) %*% Jinv
  M <- 0.5 * (M + t(M))
  attr(M, "inversion_status") <- attr(
    Jinv, "inversion_status", exact = TRUE
  )
  attr(M, "inversion_ridge") <- attr(
    Jinv, "inversion_ridge", exact = TRUE
  )
  M
}

# Empirical (robust) per-unit variability of the Pearson estimating functions.
.empirical_variability <- function(fit, cluster = NULL,
                                   small_sample = FALSE) {
  g <- fit$graph
  .empirical_variability_D(g$residuals, g$Sigma_b, g$Dlist,
                           cluster = cluster, small_sample = small_sample)
}

# The per-unit Pearson contributions themselves, n by P. Split out of
# .empirical_variability_D() so that the two-stage correction, which needs the
# contributions and not only their outer product, uses the same arithmetic
# rather than a second copy of it.
.pearson_scores <- function(E, Sigma, D) {
  Ec <- scale(E, center = TRUE, scale = FALSE)
  Sinv <- solve(Sigma)
  P <- length(D)
  n <- nrow(Ec)
  SinvD <- lapply(D, function(Dd) Sinv %*% Dd)
  tr_term <- vapply(SinvD, function(M) sum(diag(M)), numeric(1))
  Psi <- vapply(seq_len(P), function(d) {
    A <- SinvD[[d]] %*% Sinv
    -0.5 * (tr_term[d] - rowSums((Ec %*% A) * Ec))
  }, numeric(n))
  dim(Psi) <- c(n, P)
  Psi
}

# Same quantity for an explicit derivative list, so that a parametrization
# augmented with a candidate edge can be scored without refitting.
#
# The per-unit contribution is
#   psi_id = -(1/2) [ tr(Sinv D_d) - e_i' Sinv D_d Sinv e_i ].
# The quadratic form is written as a single matrix product per parameter,
# rowSums((Ec %*% A_d) * Ec), which gives every e_i' A_d e_i at once in one
# BLAS call instead of one call per row of data; it is the same arithmetic,
# reassociated, so it agrees with the row-by-row form to rounding. This function sits on the
# critical path of every Godambe-adjusted criterion, so it is evaluated once per
# candidate edge in the stepwise search and once per penalty on the L1 path.
.empirical_variability_D <- function(E, Sigma, D, cluster = NULL,
                                     small_sample = FALSE) {
  Psi <- .pearson_scores(E, Sigma, D)
  n <- nrow(Psi)
  if (is.null(cluster)) return(crossprod(Psi) / n)
  if (length(cluster) != n)
    stop("cluster must have one identifier per observational row.")
  if (anyNA(cluster))
    stop("cluster identifiers must not be missing.")
  cluster <- as.factor(cluster)
  G <- nlevels(cluster)
  if (G < 2L) stop("at least two independent clusters are required.")
  Psi_cluster <- rowsum(Psi, cluster, reorder = FALSE)
  correction <- if (isTRUE(small_sample)) G / (G - 1) else 1
  correction * crossprod(Psi_cluster) / n
}
