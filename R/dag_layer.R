# dag_layer.R: Directed Acyclic Graph layer (R wrappers over C++ kernels)
# Part of the mcglmgraph package. License: GPL-3.

#' Build the structural matrix T of the LDL parametrization
#'
#' Assembles the unit lower-triangular matrix \eqn{T} of the modified Cholesky
#' (LDL) parametrization of a directed acyclic graph (DAG) on \eqn{R}
#' responses, in which the between-response covariance is
#' \eqn{\Sigma_b = T^{-1} D (T^{-1})^\top} with \eqn{D} diagonal and positive
#' (Pourahmadi, 1999). Each edge \eqn{i \to j} contributes one structural
#' coefficient \eqn{\gamma_{ji}}, the regression coefficient of response
#' \eqn{j} on its parent \eqn{i}, and the responses are listed in a
#' topological order of the graph (parents before children).
#'
#' @param gamma finite numeric vector of structural coefficients, exactly one
#'   per edge.
#' @param R a single positive integer, the number of responses.
#' @param edges optional integer matrix with columns \eqn{(i, j)}, \eqn{i < j},
#'   the directed edges \eqn{i \to j}. If \code{NULL}, the complete DAG on the
#'   natural ordering \eqn{1 < 2 < \dots < R} is used. Each edge must occur
#'   exactly once; duplicated rows, including reversed duplicates, are rejected
#'   because they do not define distinct structural coordinates. A row given as
#'   \eqn{(j, i)} with \eqn{j > i} is read as the same edge \eqn{i \to j}: the
#'   direction always follows the natural order of the responses. The rows are
#'   sorted lexicographically before use, and \code{gamma} is read in that
#'   sorted order, not in the order the rows were supplied. Sort the matrix
#'   yourself, or read the order back from the \code{edges} component of the
#'   result, whenever the correspondence matters.
#' @return The \eqn{R \times R} unit lower-triangular matrix \eqn{T} with
#'   \eqn{T_{ji} = -\gamma_{ji}} for each edge \eqn{i \to j} and zero
#'   elsewhere off the diagonal. Its determinant is one.
#' @references Pourahmadi, M. (1999). Joint mean-covariance models with
#'   applications to longitudinal data: unconstrained parameterisation.
#'   \emph{Biometrika} 86(3), 677-690. \doi{10.1093/biomet/86.3.677}
#' @seealso \code{\link{mc_build_D_matrix}},
#'   \code{\link{mc_cholesky_to_covariance}}
#' @family directed layer
#' @examples
#' mc_build_T_matrix(c(0.4, 0.2, 0.3), R = 3)
#' @export
mc_build_T_matrix <- function(gamma, R, edges = NULL) {
  R <- .dag_dimension(R)
  eidx <- .dag_edges(edges, R)
  gamma <- .dag_gamma(gamma, nrow(eidx))
  cpp_build_T(gamma, R, eidx)
}

#' Build the diagonal innovation-variance matrix D
#'
#' Returns the diagonal matrix \eqn{D = \mathrm{diag}(\delta_1, \dots,
#' \delta_R)} of innovation variances of the LDL parametrization,
#' \eqn{\Sigma_b = T^{-1} D (T^{-1})^\top}. The innovation variance
#' \eqn{\delta_j} is the variance of response \eqn{j} that remains after
#' regression on its parents.
#'
#' @param delta numeric vector of exactly \eqn{R} finite, strictly positive
#'   innovation variances.
#' @param R a single positive integer, the number of responses.
#' @return The \eqn{R \times R} diagonal matrix \eqn{D}.
#' @seealso \code{\link{mc_build_T_matrix}}
#' @family directed layer
#' @examples
#' mc_build_D_matrix(c(1, 0.8, 0.6), 3)
#' @export
mc_build_D_matrix <- function(delta, R) {
  R <- .dag_dimension(R)
  delta <- .dag_delta(delta, R)
  diag(delta, nrow = R, ncol = R)
}

#' Covariance from the modified Cholesky (LDL) parametrization
#'
#' Computes the between-response covariance
#' \deqn{\Sigma_b = T^{-1} D (T^{-1})^\top,}
#' where \eqn{T} is the unit lower-triangular matrix of structural
#' coefficients, \eqn{T_{ji} = -\gamma_{ji}}, and \eqn{D = \mathrm{diag}(\delta)}
#' (Pourahmadi, 1999). Because \eqn{\det T = 1}, \eqn{\Sigma_b} is positive
#' definite for any \code{gamma} and any \code{delta > 0}, so the
#' parametrization is unconstrained apart from the positivity of \code{delta}.
#'
#' @inheritParams mc_build_T_matrix
#' @param delta numeric vector of exactly \eqn{R} finite, strictly positive
#'   innovation variances.
#' @return A list with \code{Sigma_b}, the \eqn{R \times R} covariance
#'   (symmetrized), \code{T_inv}, the inverse of \eqn{T}, and \code{T}.
#' @inherit mc_build_T_matrix references
#' @family directed layer
#' @examples
#' mc_cholesky_to_covariance(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)$Sigma_b
#' @export
mc_cholesky_to_covariance <- function(gamma, delta, R, edges = NULL) {
  R <- .dag_dimension(R)
  eidx <- .dag_edges(edges, R)
  gamma <- .dag_gamma(gamma, nrow(eidx))
  delta <- .dag_delta(delta, R)
  cpp_cholesky_to_cov(gamma, delta, R, eidx)
}

#' Precision from the modified Cholesky (LDL) parametrization
#'
#' Computes the between-response precision matrix
#' \eqn{K_b = \Sigma_b^{-1} = T^\top D^{-1} T} of the LDL parametrization
#' directly from \code{gamma} and \code{delta}, without inverting
#' \eqn{\Sigma_b}. A zero entry of \eqn{K_b} is a missing edge of the
#' concentration graph. For \eqn{i < k},
#' \deqn{[K_b]_{ik}=-\frac{\gamma_{ki}}{\delta_k}+
#'   \sum_{j>k}\frac{\gamma_{ji}\gamma_{jk}}{\delta_j},}
#' so the concentration graph is contained in the moral graph of the directed
#' acyclic graph (DAG), in which the parents of a common child are joined, and
#' coincides with it except where the two terms cancel, a set of coefficients
#' of Lebesgue measure zero.
#'
#' @inheritParams mc_cholesky_to_covariance
#' @return The \eqn{R \times R} precision matrix \eqn{K_b = T^\top D^{-1} T}
#'   (symmetrized).
#' @family directed layer
#' @examples
#' mc_cholesky_to_precision(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)
#' @export
mc_cholesky_to_precision <- function(gamma, delta, R, edges = NULL) {
  R <- .dag_dimension(R)
  eidx <- .dag_edges(edges, R)
  gamma <- .dag_gamma(gamma, nrow(eidx))
  delta <- .dag_delta(delta, R)
  cpp_cholesky_to_precision(gamma, delta, R, eidx)
}

#' Analytical LDL derivatives of the covariance
#'
#' Returns the derivatives of \eqn{\Sigma_b = T^{-1} D (T^{-1})^\top} with
#' respect to the innovation variances \eqn{\delta_k} and the structural
#' coefficients \eqn{\gamma_{ji}}:
#' \deqn{\frac{\partial\Sigma_b}{\partial\delta_k}=
#'   T^{-1}E_{kk}(T^{-1})^\top,\qquad
#'   \frac{\partial\Sigma_b}{\partial\gamma_{ji}}=
#'   T^{-1}E_{ji}\Sigma_b+\Sigma_b E_{ij}(T^{-1})^\top,}
#' where \eqn{E_{ab}} is the matrix with a one in position \eqn{(a, b)} and
#' zeros elsewhere. The second form is symmetric, its second term being the
#' transpose of the first. These are the matrices that enter the estimating
#' functions of the directed layer.
#'
#' @inheritParams mc_cholesky_to_covariance
#' @return A list with \code{D_delta}, a list of \eqn{R} symmetric matrices
#'   (one per innovation variance, in the order of \code{delta}), and
#'   \code{D_gamma}, a list with one symmetric matrix per edge, in the sorted
#'   order of the edges (see \code{edges} in \code{\link{mc_build_T_matrix}}).
#' @family directed layer
#' @examples
#' d <- mc_cholesky_derivatives(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)
#' length(d$D_gamma)
#' @export
mc_cholesky_derivatives <- function(gamma, delta, R, edges = NULL) {
  R <- .dag_dimension(R)
  eidx <- .dag_edges(edges, R)
  gamma <- .dag_gamma(gamma, nrow(eidx))
  delta <- .dag_delta(delta, R)
  cpp_cholesky_derivatives(gamma, delta, R, eidx)
}

# Topological order of a set of arcs, by Kahn's algorithm. Returns NULL when the
# arcs contain a directed cycle, which is what makes the request ill posed
# rather than merely awkward.
.topological_order <- function(from, to, R) {
  indeg <- tabulate(to, nbins = R)
  order_out <- integer(0)
  live <- rep(TRUE, R)
  repeat {
    src <- which(live & indeg == 0L)
    if (!length(src)) break
    src <- src[1L]                     # smallest index: a deterministic choice
    order_out <- c(order_out, src)
    live[src] <- FALSE
    out <- which(from == src)
    if (length(out)) {
      heads <- to[out]
      indeg[heads] <- indeg[heads] - 1L
    }
    from[out] <- NA_integer_; to[out] <- NA_integer_
  }
  if (length(order_out) < R) NULL else order_out
}

# ---- closed-form layer ------------------------------------------------------

#' Closed-form root of the directed layer
#'
#' Solves the Pearson estimating equations of the modified Cholesky (LDL) layer
#' of a directed acyclic graph (DAG) directly, from the covariance of the Pearson
#' residuals. Because the structural matrix has unit determinant, the Gaussian
#' working criterion of the residual layer separates by node,
#' \deqn{\ell^p_e=-\frac n2\sum_{j=1}^R\left\{\log\delta_j+
#'   \frac{q_j(\boldsymbol\gamma_j)}{\delta_j}\right\},\qquad
#'   q_j=\widehat S_{e,jj}-2\boldsymbol\gamma_j^\top\widehat{\mathbf S}_{e,P_j,j}
#'   +\boldsymbol\gamma_j^\top\widehat{\mathbf S}_{e,P_jP_j}\boldsymbol\gamma_j,}
#' so the system is \eqn{R} independent least-squares problems and its root is
#' \deqn{\hat{\boldsymbol\gamma}_j=\widehat{\mathbf S}_{e,P_jP_j}^{-1}
#'   \widehat{\mathbf S}_{e,P_j,j},\qquad
#'   \hat\delta_j=\widehat S_{e,jj}-\widehat{\mathbf S}_{e,j,P_j}
#'   \widehat{\mathbf S}_{e,P_jP_j}^{-1}\widehat{\mathbf S}_{e,P_j,j}.}
#' For a fixed directed support the layer is therefore a sequence of
#' least-squares regressions of the Pearson residuals of each response on those
#' of its parents. A response with no parent keeps \eqn{\hat\delta_j=\widehat
#' S_{e,jj}} and contributes no coefficient. Here \eqn{P_j} is the set of
#' parents of response \eqn{j}, \eqn{\boldsymbol\gamma_j} the vector of its
#' coefficients \eqn{\gamma_{jk}}, \eqn{k \in P_j}, and \eqn{n} the number of
#' units; the criterion is written up to an additive constant.
#'
#' @section Relation to the iterative engine:
#' \code{\link{mcglm_graph}} takes this route by default. Passing
#' \code{control = list(dag_solver = "newton")} runs instead the scoring
#' iteration on \eqn{(\boldsymbol\gamma,\log\boldsymbol\delta)} that the
#' directed and undirected layers share. It solves the same system and serves
#' as an independent check on the closed form.
#'
#' @section Giving the arcs directly:
#' With \code{directed = FALSE} the edges are unordered adjacencies and their
#' direction comes from \code{ordering}, as elsewhere in the package. That is
#' awkward when the support to fit is an estimated DAG, whose arcs need not run
#' along the natural order of the responses. With \code{directed = TRUE} the two
#' columns of \code{edges} are read as parent and child, any acyclic set of arcs
#' is accepted, and a compatible ordering is derived internally by Kahn's
#' algorithm. A cyclic set is refused, since it has no modified Cholesky
#' factorization at all.
#'
#' @param Shat the \eqn{R \times R} covariance of the Pearson residuals, in the
#'   order of the responses.
#' @param ordering optional integer permutation of \code{1:R} giving the order
#'   of the responses in the factorization, parents before children, so that
#'   \code{ordering[1]} has no parent; \code{NULL} uses \code{1:R}. With
#'   \code{directed = TRUE} it is derived from the arcs when left \code{NULL},
#'   and checked against them when supplied.
#' @param directed logical; \code{TRUE} reads the columns of \code{edges} as
#'   parent and child rather than as an unordered pair.
#' @param R a single positive integer, the number of responses.
#' @param edges optional two-column integer matrix of response indices. With
#'   \code{directed = FALSE} each row is an unordered adjacency and its direction
#'   follows \code{ordering} (the earlier response is the parent). With
#'   \code{directed = TRUE} each row is an arc (parent, child); each pair of
#'   responses may appear at most once, and self-loops are refused. If
#'   \code{NULL}, the complete DAG on \code{ordering} is used.
#' @return A list with \code{gamma}, \code{delta}, the implied \code{Sigma_b}
#'   and \code{K_b} in the order of the responses, the \code{edges} actually
#'   used as (parent, child) pairs in response labels with \code{gamma} aligned
#'   to their rows, the \code{ordering}, the \code{parents} of each response,
#'   and \code{status}, which records how many nodes needed a ridge on their
#'   normal equations, how many innovation variances were lifted to the
#'   positivity floor (\eqn{10^{-10}\widehat S_{e,jj}}), and
#'   \code{ill_conditioned}, the responses whose parent block, rescaled to a
#'   unit diagonal, has a reciprocal condition number at or below
#'   \code{.Machine$double.eps^(2/3)}; their coefficients are not identified
#'   and come from a solve with a ridge of \eqn{10^{-8}} on that unit
#'   diagonal. Both tests are invariant to a rescaling of the responses, and
#'   so is the root: rescaling response \eqn{j} by \eqn{a_j} multiplies
#'   \eqn{\gamma_{jk}} by \eqn{a_j/a_k} and \eqn{\delta_j} by \eqn{a_j^2}. A
#'   nonzero status means the residual covariance is close to singular on some
#'   parent set, and \code{\link{mcglm_graph}} then does not certify the fit
#'   for inference. \code{delta} is indexed by response, like the diagonal of
#'   \code{Sigma_b}.
#' @references Pourahmadi, M. (1999). Joint mean-covariance models with
#'   applications to longitudinal data: unconstrained parameterisation.
#'   \emph{Biometrika} 86(3), 677-690. \doi{10.1093/biomet/86.3.677}
#' @seealso \code{\link{mc_dag_sensitivity}}, \code{\link{mc_dag_coefficients}}
#' @family directed layer
#' @examples
#' S <- mc_cholesky_to_covariance(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49),
#'                                R = 3)$Sigma_b
#' fit <- mc_dag_closed_root(S, R = 3)
#' round(fit$gamma, 12)     # recovers the generating coefficients
#' round(fit$delta, 12)
#' @export
mc_dag_closed_root <- function(Shat, R, edges = NULL, ordering = NULL,
                               directed = FALSE) {
  R <- .dag_dimension(R)
  if (length(directed) != 1L || !is.logical(directed) || is.na(directed))
    stop("directed must be a single TRUE or FALSE.", call. = FALSE)
  Shat <- as.matrix(Shat)
  if (nrow(Shat) != R || ncol(Shat) != R || !is.numeric(Shat) ||
      any(!is.finite(Shat)))
    stop("Shat must be a finite numeric ", R, " by ", R, " matrix.",
         call. = FALSE)
  if (!isTRUE(all.equal(Shat, t(Shat), tolerance = 1e-10)))
    stop("Shat must be symmetric.", call. = FALSE)
  Shat <- 0.5 * (Shat + t(Shat))

  if (isTRUE(directed) && !is.null(edges)) {
    arcs <- .validate_arcs(edges, R)
    topo <- .topological_order(arcs[, 1L], arcs[, 2L], R)
    if (is.null(topo))
      stop("the arcs contain a directed cycle, so they have no modified ",
           "Cholesky factorization.", call. = FALSE)
    if (is.null(ordering)) {
      ordering <- topo
    } else {
      ordering <- .check_ordering(ordering, R)
      rank_of <- match(seq_len(R), ordering)
      if (any(rank_of[arcs[, 1L]] >= rank_of[arcs[, 2L]]))
        stop("the supplied ordering is not compatible with the arcs: some ",
             "parent comes after its child.", call. = FALSE)
    }
    edges <- arcs                       # unordered pairs are formed below
  } else if (!is.null(ordering)) {
    ordering <- .check_ordering(ordering, R)
  }
  if (is.null(ordering)) ordering <- seq_len(R)

  spec <- structure(list(type = "dag", edges = edges, ordering = ordering),
                    class = "mcg_graph_spec")
  dag <- .dag_setup(spec, R)
  eidx <- dag$edges_pos                 # positional, i < j, lexicographic
  inv <- dag$inverse

  Shat_pos <- Shat[ordering, ordering, drop = FALSE]
  res <- cpp_dag_closed_root(Shat_pos, R, eidx)
  gamma <- as.numeric(res$gamma)
  delta_pos <- as.numeric(res$delta)
  cc <- cpp_cholesky_to_cov(gamma, delta_pos, R, eidx)
  K_pos <- cpp_cholesky_to_precision(gamma, delta_pos, R, eidx)

  eo <- dag$edges_original              # (parent, child) in response labels
  colnames(eo) <- c("parent", "child")
  delta <- numeric(R)
  delta[ordering] <- delta_pos          # index by response, like diag(Sigma_b)
  parents <- lapply(seq_len(R), function(j)
    if (nrow(eo)) as.integer(eo[eo[, 2L] == j, 1L]) else integer(0))

  list(gamma = gamma, delta = delta,
       Sigma_b = cc$Sigma_b[inv, inv, drop = FALSE],
       K_b = K_pos[inv, inv, drop = FALSE],
       parents = parents, edges = eo, ordering = ordering,
       status = list(n_ridge = as.integer(res$n_ridge),
                     n_floored = as.integer(res$n_floored),
                     ill_conditioned = as.integer(
                       ordering[as.integer(res$ill_conditioned)])))
}

# Arcs given as (parent, child), with no constraint that the parent index be the
# smaller one. Reversed duplicates are a contradiction here, not a repetition.
.validate_arcs <- function(arcs, R, arg = "edges") {
  arcs <- as.matrix(arcs)
  if (length(dim(arcs)) != 2L || ncol(arcs) != 2L)
    stop(arg, " must have exactly two columns.", call. = FALSE)
  if (nrow(arcs) == 0L) return(matrix(integer(0), nrow = 0L, ncol = 2L))
  if (!is.numeric(arcs) || any(!is.finite(arcs)) || any(arcs != round(arcs)))
    stop(arg, " must contain only finite integer indices.", call. = FALSE)
  storage.mode(arcs) <- "integer"
  if (any(arcs[, 1L] == arcs[, 2L]))
    stop("self-loops are not allowed.", call. = FALSE)
  if (any(arcs < 1L) || any(arcs > R))
    stop("edge indices out of range: expected values in 1:", R, ".",
         call. = FALSE)
  key <- paste(pmin(arcs[, 1L], arcs[, 2L]), pmax(arcs[, 1L], arcs[, 2L]),
               sep = ":")
  if (anyDuplicated(key))
    stop("a pair of responses appears twice among the arcs; each pair carries ",
         "at most one direction.", call. = FALSE)
  arcs
}

.check_ordering <- function(ordering, R) {
  if (!is.numeric(ordering) || length(ordering) != R ||
      any(!is.finite(ordering)) || any(ordering != round(ordering)) ||
      !setequal(as.integer(ordering), seq_len(R)))
    stop("ordering must be a permutation of seq_len(R).", call. = FALSE)
  as.integer(ordering)
}

#' Block-diagonal sensitivity of the directed layer
#'
#' Assembles the model-based Pearson sensitivity \eqn{\mathbf J_n} of the
#' directed layer from its closed form. The empirical sandwich and two-stage
#' covariance functions use the observed Hessian of the working criterion
#' instead, since a sparse directed acyclic graph (DAG) can omit second-moment
#' structure. On the accumulated scale, with \eqn{P_j} the parents of response
#' \eqn{j} and \eqn{\boldsymbol\gamma_j} their coefficients,
#' \deqn{\mathbf J_{n,\boldsymbol\gamma_j\boldsymbol\gamma_j}=
#'   \frac{n}{\delta_j}\boldsymbol\Sigma_{b,P_jP_j},\qquad
#'   J_{n,\delta_j\delta_j}=\frac{n}{2\delta_j^2},}
#' and every other block vanishes: between the coefficients of distinct
#' responses, between coefficients and innovation variances, and between
#' innovation variances of distinct responses. The generic assembly used by the
#' fitting engine forms \eqn{P^2} traces of \eqn{R \times R} products; this form
#' costs one block per node.
#'
#' @param Sigma_b the fitted \eqn{R \times R} covariance of the layer.
#' @param delta the \eqn{R} innovation variances.
#' @param n number of experimental units.
#' @inheritParams mc_build_T_matrix
#' @return The \eqn{(|E| + R) \times (|E| + R)} symmetric, block-diagonal
#'   sensitivity in the parameter order \code{c(gamma, delta)}, where
#'   \eqn{|E|} is the number of edges and the \code{gamma} entries follow the
#'   sorted order of the edges (see \code{edges} in
#'   \code{\link{mc_build_T_matrix}}).
#' @seealso \code{\link{mc_dag_closed_root}}
#' @family directed layer
#' @examples
#' cc <- mc_cholesky_to_covariance(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49), R = 3)
#' J <- mc_dag_sensitivity(cc$Sigma_b, c(1, 0.64, 0.49), n = 100, R = 3)
#' round(J, 6)
#' @export
mc_dag_sensitivity <- function(Sigma_b, delta, n, R, edges = NULL) {
  R <- .dag_dimension(R)
  eidx <- .dag_edges(edges, R)
  delta <- .dag_delta(delta, R)
  Sigma_b <- as.matrix(Sigma_b)
  if (nrow(Sigma_b) != R || ncol(Sigma_b) != R || !is.numeric(Sigma_b) ||
      any(!is.finite(Sigma_b)))
    stop("Sigma_b must be a finite numeric ", R, " by ", R, " matrix.",
         call. = FALSE)
  if (length(n) != 1L || !is.finite(n) || n <= 0)
    stop("n must be one finite positive number.", call. = FALSE)
  cpp_dag_sensitivity(0.5 * (Sigma_b + t(Sigma_b)), delta, as.numeric(n), eidx)
}

#' Scale-free coefficients of a fitted directed layer
#'
#' Reports, with standard errors and intervals, the quantities of the directed
#' layer that do not depend on the scale
#' on which the residual dependence is measured. Under a diagonal rescaling
#' \eqn{\boldsymbol\Sigma_b\mapsto\mathbf A\boldsymbol\Sigma_b\mathbf A} the raw
#' coefficients move as \eqn{\gamma_{jk}\mapsto\gamma_{jk}a_j/a_k} and the
#' innovation variances as \eqn{\delta_j\mapsto a_j^2\delta_j}, whereas
#' \deqn{\gamma^*_{jk}=\gamma_{jk}\left(\frac{\Sigma_{b,kk}}{\Sigma_{b,jj}}
#'   \right)^{1/2},\qquad R^2_{j\cdot P_j}=1-\frac{\delta_j}{\Sigma_{b,jj}}}
#' are invariant. These are the quantities to report: they play for the
#' directed layer the role that partial correlations play for the undirected
#' one.
#'
#' @section What the coefficients are:
#' \eqn{\gamma_{jk}} is the coefficient of response \eqn{k} in the linear
#' projection of the residual component of response \eqn{j} on those of its
#' parents, under the ordering of the fit, and \eqn{\delta_j} is the variance
#' left by that projection. The factorization does not by itself define a
#' causal model. Reading a coefficient as the effect of an intervention
#' requires a linear structural model with independent exogenous innovations,
#' the correct ordering, no unobserved confounding and equations that are
#' invariant under intervention. The fit checks none of these.
#'
#' @section When the first stage is such a rescaling:
#' The first stage acts on the residual covariance, to first order, as a
#' diagonal rescaling under conditions (b) and (c) of the section "What the
#' correction reaches" of \code{\link{mc_vcov_graph}}: a correctly specified
#' mean, and a conditional residual covariance that does not depend on the
#' covariates. Under them the direction in which the first stage moves
#' \eqn{\boldsymbol\Sigma_b} is \eqn{\boldsymbol\Lambda\boldsymbol\Sigma_b+
#' \boldsymbol\Sigma_b\boldsymbol\Lambda} with \eqn{\boldsymbol\Lambda}
#' diagonal, whatever the parametrization, and the quantities above do not see
#' it. The two conditions are sufficient: with the constant variance function
#' (c) is not needed, and with a variance function that depends on the mean a
#' residual correlation that changes with the covariates moves
#' \eqn{\boldsymbol\Sigma_b} in a direction the standardized coefficients do
#' see. The finer characterization of the section "When the shortcut is right"
#' of \code{\link{mc_vcov_graph}} was established for the undirected layer; for
#' this one the sufficient conditions are the ones to rely on. The package
#' checks none of them.
#'
#' @section What the standard errors cover:
#' They come from \code{\link{mc_vcov_graph}} and refer to the second stage
#' conditionally on the estimated residuals, unless \code{variability =
#' "two_stage"}. The empirical forms use the observed Hessian of the working
#' criterion;
#' the model form uses the Pearson sensitivity and requires correct
#' second-moment specification and its Gaussian working variability identity.
#' The Wald statistic of \eqn{\gamma_{jk}=0} is invariant under
#' the diagonal rescaling, because the coefficient and its standard error are
#' both multiplied by \eqn{a_j/a_k}, so under (b) and (c) the test is adequate to
#' first order where the raw coefficient is not.
#'
#' @section Intervals for the scale-free quantities:
#' Both \eqn{\gamma^*_{jk}} and \eqn{R^2_{j\cdot P_j}} are functions of
#' \eqn{(\boldsymbol\gamma,\boldsymbol\delta)}, through the coefficient itself
#' and through \eqn{\boldsymbol\Sigma_b=\mathbf T^{-1}\mathbf D\mathbf T^{-\top}},
#' so their variance follows by the delta method with the analytical gradients
#' \deqn{\frac{\partial\gamma^*_{jk}}{\partial\theta_d}=
#'   \frac{\partial\gamma_{jk}}{\partial\theta_d}\,f_{jk}+
#'   \frac{\gamma^*_{jk}}{2}\left(\frac{\dot\Sigma_{b,d,kk}}{\Sigma_{b,kk}}-
#'   \frac{\dot\Sigma_{b,d,jj}}{\Sigma_{b,jj}}\right),\qquad
#'   \frac{\partial R^2_{j\cdot P_j}}{\partial\theta_d}=
#'   -\frac{1}{\Sigma_{b,jj}}\frac{\partial\delta_j}{\partial\theta_d}+
#'   \frac{\delta_j\,\dot\Sigma_{b,d,jj}}{\Sigma_{b,jj}^2},}
#' with \eqn{f_{jk}=(\Sigma_{b,kk}/\Sigma_{b,jj})^{1/2}} and
#' \eqn{\dot{\boldsymbol\Sigma}_{b,d}=\partial\boldsymbol\Sigma_b/
#' \partial\theta_d}. That is \code{interval = "delta"}, the default. With
#' \code{interval = "fixed_factor"} the interval for \eqn{\gamma^*} is the one
#' for \eqn{\gamma} multiplied by \eqn{f_{jk}} as if the factor were known,
#' which leaves out the uncertainty of the variances in it. The interval for
#' \eqn{R^2} is always by the delta method, since it has no fixed-factor
#' analogue.
#'
#' @param fit an object of class \code{"mcglm_graph"} fitted with the DAG layer.
#' @param variability \code{"model"} for \eqn{\mathbf J_n^{-1}},
#'   \code{"empirical"} for the sandwich, or \code{"two_stage"} for the sandwich
#'   corrected for the estimation of the residuals.
#' @param level confidence level of the intervals.
#' @param interval \code{"delta"} (default) or \code{"fixed_factor"}, the
#'   interval for \eqn{\gamma^*} as described above.
#' @param cluster,small_sample passed to \code{\link{mc_vcov_graph}} for the
#'   empirical forms.
#' @param r2_interval \code{"wald"} (default), symmetric on the scale of
#'   \eqn{R^2}, or \code{"logit"}, symmetric on the scale of
#'   \eqn{\log\{R^2/(1-R^2)\}} and mapped back, so that the bounds stay inside
#'   \eqn{(0, 1)}. A response without parents has \eqn{R^2=0} by construction,
#'   and its interval is that point.
#' @return An object of class \code{"mcg_dag_coef"}: a list with \code{edges},
#'   one row per directed edge with the raw and standardized coefficient, the
#'   standard error of each, the Wald statistic of the raw coefficient, its
#'   p-value and the two intervals, and \code{nodes}, one row per response with
#'   the innovation variance and its standard error, and
#'   \eqn{R^2_{j\cdot P_j}} with its standard error and interval. The gradients
#'   of \eqn{\gamma^*} and of \eqn{R^2} travel in the attributes
#'   \code{"gradient_std"} and \code{"gradient_r2"}. The list also records
#'   \code{variability}, \code{level}, \code{interval} and
#'   \code{r2_interval}. The \code{print} method shows the main columns and
#'   returns \code{x} invisibly.
#' @seealso \code{\link{mc_vcov_graph}}, \code{\link{mc_dag_unit_diagonal}}
#' @family directed layer
#' @examples
#' \donttest{
#' set.seed(4); n <- 400; x <- rnorm(n)
#' y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- 0.4 * y2 + rnorm(n)
#' d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x),
#'                  between_graph = mc_between_graph("dag"), data = d)
#' mc_dag_coefficients(f)
#' }
#' @export
mc_dag_coefficients <- function(fit,
                                variability = c("model", "empirical",
                                                "two_stage"),
                                level = 0.95,
                                interval = c("delta", "fixed_factor"),
                                cluster = NULL, small_sample = FALSE,
                                r2_interval = c("wald", "logit")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  variability <- match.arg(variability)
  interval <- match.arg(interval)
  r2_interval <- match.arg(r2_interval)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a single number strictly between 0 and 1.",
         call. = FALSE)
  g <- fit$graph
  if (!identical(g$type, "dag"))
    stop("mc_dag_coefficients() applies to a DAG layer; this fit carries a ",
         toupper(g$type), " layer.", call. = FALSE)

  V <- as.matrix(mc_vcov_graph(
    fit, type = switch(variability, model = "model", empirical = "sandwich",
                       two_stage = "two_stage"),
    cluster = cluster, small_sample = small_sample))
  se_all <- sqrt(pmax(diag(V), 0))
  R <- g$R
  rn <- g$response_names
  dag <- .dag_setup(g$spec, R)
  eo <- dag$edges_original                # parent, child in response labels
  E <- nrow(eo)
  P <- length(g$par)
  Sd <- diag(g$Sigma_b)                   # in response order
  delta <- g$par[E + seq_len(R)]          # in positional order
  resp_of_pos <- dag$ordering             # position -> response
  q <- stats::qnorm(1 - (1 - level) / 2)
  # d Sigma_bb / d theta_d for every response b, one row per parameter
  dS <- t(vapply(g$Dlist, function(Dd) diag(as.matrix(Dd)), numeric(R)))
  dS <- matrix(dS, nrow = P)
  delta_se <- function(G) sqrt(pmax(rowSums((G %*% V) * G), 0))

  if (E) {
    parent <- eo[, 1L]; child <- eo[, 2L]
    gam <- as.numeric(g$par[seq_len(E)])
    fac <- sqrt(Sd[parent] / Sd[child])
    gstd <- gam * fac
    se <- se_all[seq_len(E)]
    z <- ifelse(se > 0, gam / se, NA_real_)
    half <- q * se
    Gs <- matrix(0, E, P)
    for (e in seq_len(E)) {
      Gs[e, e] <- fac[e]
      Gs[e, ] <- Gs[e, ] + 0.5 * gstd[e] *
        (dS[, parent[e]] / Sd[parent[e]] - dS[, child[e]] / Sd[child[e]])
    }
    se_std <- delta_se(Gs)
    if (interval == "delta") {
      lo_std <- gstd - q * se_std; up_std <- gstd + q * se_std
    } else {
      lo_std <- (gam - half) * fac; up_std <- (gam + half) * fac
    }
    edges_df <- data.frame(
      from = as.integer(parent), to = as.integer(child),
      parent = rn[parent], child = rn[child],
      gamma = gam, gamma_std = gstd, se = se, se_std = se_std,
      statistic = z, p_value = 2 * stats::pnorm(-abs(z)),
      lower = gam - half, upper = gam + half,
      lower_std = lo_std, upper_std = up_std,
      row.names = NULL, stringsAsFactors = FALSE
    )
  } else {
    Gs <- matrix(0, 0, P)
    edges_df <- data.frame(
      from = integer(0), to = integer(0), parent = character(0),
      child = character(0), gamma = numeric(0), gamma_std = numeric(0),
      se = numeric(0), se_std = numeric(0), statistic = numeric(0),
      p_value = numeric(0), lower = numeric(0), upper = numeric(0),
      lower_std = numeric(0), upper_std = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  n_parents <- vapply(seq_len(R), function(pos)
    if (E) sum(dag$edges_pos[, 2L] == pos) else 0L, integer(1))
  Sj <- Sd[resp_of_pos]
  r2 <- 1 - as.numeric(delta) / Sj
  Gr <- matrix(0, R, P)
  for (p in seq_len(R)) {
    Gr[p, E + p] <- -1 / Sj[p]
    Gr[p, ] <- Gr[p, ] + as.numeric(delta[p]) * dS[, resp_of_pos[p]] / Sj[p]^2
  }
  se_r2 <- delta_se(Gr)
  if (r2_interval == "logit") {
    # Wald on log{R^2 / (1 - R^2)}, with standard error se / {R^2 (1 - R^2)}
    hl <- q * se_r2 / (r2 * (1 - r2))
    lo_r2 <- stats::plogis(stats::qlogis(r2) - hl)
    up_r2 <- stats::plogis(stats::qlogis(r2) + hl)
    # R^2 = 0 by construction without parents, and the logit is not defined
    # at the boundary: the interval is the point
    edge <- n_parents == 0L | !is.finite(stats::qlogis(r2))
    lo_r2[edge] <- up_r2[edge] <- r2[edge]
  } else {
    lo_r2 <- r2 - q * se_r2; up_r2 <- r2 + q * se_r2
  }
  nodes_df <- data.frame(
    position = seq_len(R),
    response = rn[resp_of_pos],
    n_parents = n_parents,
    delta = as.numeric(delta),
    se_delta = se_all[E + seq_len(R)],
    r_squared = r2,
    se_r_squared = se_r2,
    lower_r_squared = lo_r2,
    upper_r_squared = up_r2,
    row.names = NULL, stringsAsFactors = FALSE
  )

  structure(list(edges = edges_df, nodes = nodes_df,
                 variability = variability, level = level,
                 interval = interval, r2_interval = r2_interval),
            gradient_std = Gs, gradient_r2 = Gr,
            class = "mcg_dag_coef")
}

#' @rdname mc_dag_coefficients
#' @param x an object of class \code{"mcg_dag_coef"}.
#' @param ... ignored.
#' @export
print.mcg_dag_coef <- function(x, ...) {
  cat(sprintf("Directed layer (%s variability, level %.2f)\n",
              x$variability, x$level))
  cat("\nEdges (parent -> child):\n")
  if (nrow(x$edges)) {
    show <- x$edges[, c("parent", "child", "gamma", "gamma_std", "se",
                        "se_std", "statistic", "p_value")]
    print(format(show, digits = 4))
  } else cat("  none\n")
  cat("\nNodes:\n")
  print(format(x$nodes[, c("response", "n_parents", "delta", "se_delta",
                           "r_squared", "se_r_squared")], digits = 4))
  invisible(x)
}

#' Innovation variances that force a unit diagonal
#'
#' Returns the innovation variances under which the covariance of the directed
#' layer has unit diagonal, that is under which it is a correlation matrix, and
#' says whether the coefficients supplied admit them.
#'
#' @section Why this is a constraint and not a reparametrization:
#' The undirected layer can fix the diagonal at one at the cost of nothing. Here
#' the requirement
#' \deqn{\delta_j=1-\boldsymbol\gamma_j^\top\boldsymbol\Sigma_{b,P_jP_j}
#'   \boldsymbol\gamma_j,\qquad\text{equivalently}\qquad
#'   \delta_1=1,\quad\delta_j=1-\sum_{k<j}L_{jk}^2\delta_k,}
#' with \eqn{\mathbf L=\mathbf T^{-1}}, leaves \eqn{|E_D|} free parameters but
#' confines them to the open set where every \eqn{\delta_j} stays positive, a
#' set defined by \eqn{R} nonlinear inequalities. It contains the origin and is
#' a proper subset of the coefficient space whenever there is an edge. That is
#' exactly the freedom the layer buys by leaving the diagonal free, which is why
#' the fitting routines do not impose it: the constraint is useful for reading a
#' covariance as a correlation and for reproducing a generator, not for
#' estimating.
#'
#' @section One representative per orbit:
#' Rescaling by \eqn{\mathbf A} sends \eqn{\delta_j} to \eqn{a_j^2\delta_j}, so
#' every orbit \eqn{\{\mathbf A\boldsymbol\Sigma_b\mathbf A\}} contains exactly
#' one matrix with unit innovation variances, reached by
#' \eqn{a_j=\delta_j^{-1/2}}. The component \code{to_unit_innovations} carries
#' that rescaling.
#'
#' @inheritParams mc_build_T_matrix
#' @return A list with \code{delta}, the innovation variances implied by the
#'   coefficients; \code{admissible}, true when all of them are positive;
#'   \code{min_delta}; \code{Sigma_b}, the covariance they produce, which has
#'   unit diagonal when admissible; \code{r_squared}, the variance of each
#'   response accounted for by its parents; and
#'   \code{to_unit_innovations}, the diagonal rescaling that sends the
#'   covariance to the representative of its orbit with unit innovations.
#' @examples
#' # a chain: the second and third responses keep less and less of their own
#' mc_dag_unit_diagonal(c(0.6, 0.6), R = 3, edges = rbind(c(1, 2), c(2, 3)))$delta
#' # coefficients too large for a correlation matrix are flagged
#' mc_dag_unit_diagonal(c(0.9, 0.9), R = 3,
#'                      edges = rbind(c(1, 3), c(2, 3)))$admissible
#' @seealso \code{\link{mc_dag_coefficients}} for the scale-free quantities.
#' @family directed layer
#' @export
mc_dag_unit_diagonal <- function(gamma, R, edges = NULL) {
  R <- .dag_dimension(R)
  eidx <- .dag_edges(edges, R)
  gamma <- .dag_gamma(gamma, nrow(eidx))
  Tm <- cpp_build_T(gamma, R, eidx)
  L <- solve(Tm)                       # unit lower triangular

  delta <- numeric(R)
  for (j in seq_len(R)) {
    prev <- if (j > 1L) sum(L[j, seq_len(j - 1L)]^2 * delta[seq_len(j - 1L)]) else 0
    delta[j] <- 1 - prev
  }
  admissible <- all(delta > 0)
  Sigma <- L %*% diag(delta, R) %*% t(L)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  list(delta = delta, admissible = admissible, min_delta = min(delta),
       Sigma_b = Sigma,
       r_squared = 1 - delta / diag(Sigma),
       to_unit_innovations = if (admissible) 1 / sqrt(delta) else
         rep(NA_real_, R))
}

# ---- internal helpers -------------------------------------------------------

# Single edge validator for the whole package. Every entry point (the
# low-level LDL wrappers, .dag_setup(), .check_edges() and mc_between_graph())
# goes through it, so the contract a caller sees does not depend on the entry
# point: a non-integer index is refused, never truncated to another edge.
#
# Returns an integer matrix with i < j, ordered lexicographically. R = NULL
# skips the upper range check, when R is not yet known.
.validate_edges <- function(edges, R, arg = "edges") {
  if (!is.null(R) || is.null(edges)) R <- .dag_dimension(R)
  if (is.null(edges)) {
    return(if (R >= 2L) t(utils::combn(R, 2L)) else
      matrix(integer(0), nrow = 0L, ncol = 2L))
  }
  edges <- as.matrix(edges)
  if (length(dim(edges)) != 2L || ncol(edges) != 2L)
    stop(arg, " must have exactly two columns.", call. = FALSE)
  if (nrow(edges) == 0L)
    return(matrix(integer(0), nrow = 0L, ncol = 2L))
  if (!is.numeric(edges) || any(!is.finite(edges)) || any(edges != round(edges)))
    stop(arg, " must contain only finite integer indices.", call. = FALSE)
  storage.mode(edges) <- "integer"
  edges <- cbind(pmin(edges[, 1L], edges[, 2L]),
                 pmax(edges[, 1L], edges[, 2L]))
  if (any(edges[, 1L] == edges[, 2L]))
    stop("self-loops are not allowed.", call. = FALSE)
  if (is.null(R)) {
    if (any(edges < 1L))
      stop("edge indices out of range: expected positive values.",
           call. = FALSE)
  } else if (any(edges < 1L) || any(edges > R)) {
    stop("edge indices out of range: expected values in 1:", R, ".",
         call. = FALSE)
  }
  if (anyDuplicated(paste(edges[, 1L], edges[, 2L], sep = ":")))
    stop("duplicated edges are not allowed; each pair defines one parameter.",
         call. = FALSE)
  edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
}

.dag_edges <- function(edges, R) .validate_edges(edges, R)

.dag_dimension <- function(R) {
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) ||
      R < 1 || R != floor(R) || R > .Machine$integer.max)
    stop("R must be a single positive integer.")
  as.integer(R)
}

.dag_gamma <- function(gamma, n_edges) {
  if (!is.numeric(gamma) || length(gamma) != n_edges ||
      any(!is.finite(gamma)))
    stop(sprintf("gamma must contain exactly %d finite value%s, one per edge.",
                 n_edges, if (n_edges == 1L) "" else "s"))
  as.numeric(gamma)
}

.dag_delta <- function(delta, R) {
  if (!is.numeric(delta) || length(delta) != R || any(!is.finite(delta)))
    stop(sprintf("delta must contain exactly %d finite positive value%s.",
                 R, if (R == 1L) "" else "s"))
  if (any(delta <= 0)) stop("delta must be strictly positive.")
  as.numeric(delta)
}
