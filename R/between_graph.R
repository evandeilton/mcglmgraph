# between_graph.R: unified between-response covariance specification and the
# second-order estimating-function engine for Sigma_b.
# Part of the mcglmgraph package. License: GPL-3.

#' Specify the between-response graph layer
#'
#' Builds the specification object that controls how the between-response
#' covariance \eqn{\Sigma_b} of the Pearson residuals is parametrized in
#' \code{\link{mcglm_graph}}. Two layers are available.
#' \itemize{
#'   \item \code{"ggm"}, the Gaussian graphical model (GGM): the precision
#'     matrix is linear in the parameters,
#'     \eqn{K_b = \Sigma_b^{-1} = \sum_d \eta_d Z_d}, with one diagonal
#'     parameter per response and one parameter per edge (see
#'     \code{\link{mc_precision_bases}}). A missing edge sets the
#'     corresponding precision entry, and hence the partial correlation, to
#'     zero.
#'   \item \code{"dag"}, the directed acyclic graph (DAG) layer, through the
#'     modified Cholesky (LDL) factorization
#'     \eqn{\Sigma_b = T^{-1} D T^{-\top}}, with \eqn{T = I - \Gamma} unit
#'     lower triangular in the given ordering, \eqn{\Gamma} holding the
#'     structural coefficients \eqn{\gamma_{ji}} on the edges, and
#'     \eqn{D = \mathrm{diag}(\delta_1, \dots, \delta_R)} the innovation
#'     variances. Each response is regressed on its parents.
#' }
#'
#' @param type character, the layer: \code{"ggm"} for the Gaussian graphical
#'   model (precision-based) or \code{"dag"} for the directed acyclic graph
#'   (modified Cholesky).
#' @param edges optional integer matrix with columns \eqn{(i, j)} giving the
#'   edges as response positions; \code{NULL} means the complete graph. A
#'   matrix with no rows is the empty graph. The edges are stored in canonical
#'   form, each pair with \eqn{i < j} and the rows in lexicographic order,
#'   which is the order of the edge parameters, so the order in which they are
#'   supplied does not matter. In the DAG layer an edge is an adjacency and its
#'   direction is set by \code{ordering}.
#' @param ordering optional integer permutation of \code{1:R} defining the
#'   ordering of the DAG layer: an edge points from the response that comes
#'   first in \code{ordering} to the one that comes later. The default is the
#'   natural order. It has no meaning for \code{type = "ggm"}, whose precision
#'   support is unordered, and supplying it there is an error rather than a
#'   silently ignored argument.
#' @return An object of class \code{"mcg_graph_spec"}: a list with components
#'   \code{type}, \code{edges} (canonical edge matrix or \code{NULL}) and
#'   \code{ordering} (integer vector or \code{NULL}).
#' @references
#' Dempster, A. P. (1972). Covariance selection. \emph{Biometrics}, 28(1),
#' 157-175. \doi{10.2307/2528966}
#'
#' Lauritzen, S. L. (1996). \emph{Graphical Models}. Oxford University Press.
#'
#' Pourahmadi, M. (1999). Joint mean-covariance models with applications to
#' longitudinal data: unconstrained parameterisation. \emph{Biometrika}, 86(3),
#' 677-690. \doi{10.1093/biomet/86.3.677}
#' @seealso \code{\link{mcglm_graph}}, \code{\link{mc_build_sigma_between_graph}}
#' @examples
#' mc_between_graph("ggm")
#' mc_between_graph("dag", edges = rbind(c(1, 2), c(2, 3)))
#' @export
mc_between_graph <- function(type = c("ggm", "dag"), edges = NULL,
                             ordering = NULL) {
  type <- match.arg(type)
  if (type == "ggm" && !is.null(ordering))
    stop("ordering applies to the DAG layer only; the GGM precision support ",
         "is unordered. Drop the argument or use type = \"dag\".",
         call. = FALSE)
  if (!is.null(ordering)) {
    if (!is.numeric(ordering) || any(!is.finite(ordering)) ||
        any(ordering != round(ordering)) || anyDuplicated(ordering) ||
        !identical(sort(as.integer(ordering)), seq_along(ordering)))
      stop("ordering must be a permutation of seq_len(R).", call. = FALSE)
    ordering <- as.integer(ordering)
  }
  # Canonical edges: every consumer then reads them in the parameter order.
  if (!is.null(edges))
    edges <- .validate_edges(edges, if (is.null(ordering)) NULL else
      length(ordering))
  structure(list(type = type, edges = edges, ordering = ordering),
            class = "mcg_graph_spec")
}

# Canonical copy of a specification for R responses: validated edges in the
# parameter order and R recorded, whatever way the specification was built.
.canonical_spec <- function(spec, R) {
  if (!is.null(spec$edges)) spec$edges <- .validate_edges(spec$edges, R)
  spec$R <- as.integer(R)
  spec
}

#' Build the between-response covariance and its derivatives under a graph
#'
#' Maps a vector of graph parameters to the between-response covariance
#' \eqn{\Sigma_b} and to the list of its derivatives with respect to each
#' parameter, in the list form consumed by the estimating-function engine of
#' \code{\link{mcglm_graph}}. It takes the place that the correlation-based
#' between-response builder has in a multivariate covariance generalized linear
#' model (MCGLM). For the Gaussian graphical model (GGM) layer,
#' \eqn{\Sigma_b = (\sum_d \eta_d Z_d)^{-1}} and
#' \eqn{\partial\Sigma_b/\partial\eta_d = -\Sigma_b Z_d \Sigma_b}; for the
#' directed acyclic graph (DAG) layer, \eqn{\Sigma_b = T^{-1} D T^{-\top}},
#' with \eqn{T = I - \Gamma} and \eqn{D = \mathrm{diag}(\delta)}, and the
#' derivatives are taken with respect to \eqn{\gamma} and \eqn{\delta}, the
#' innovation variances themselves and not their logarithms. Matrices are
#' returned in the order of the responses, whatever the DAG ordering.
#'
#' @param theta_G numeric vector of graph parameters. For \code{"ggm"}: the
#'   \eqn{\eta} vector (diagonal parameters, then edges in canonical order,
#'   length \eqn{R + |E|}). For \code{"dag"}: \code{c(gamma, delta)}, with the
#'   \eqn{|E|} coefficients \eqn{\gamma}, one per edge, sorted
#'   lexicographically by the positions of its endpoints in the ordering (the
#'   canonical edge order under the natural ordering), followed by the \eqn{R}
#'   innovation variances \eqn{\delta}, indexed by position in the ordering.
#'   The GGM vector must give a positive-definite precision, otherwise an error
#'   is raised; the DAG vector is expected to have \eqn{\delta > 0}, which is
#'   not checked.
#' @param n_resp integer, number of responses \eqn{R}.
#' @param spec an object of class \code{"mcg_graph_spec"} from
#'   \code{\link{mc_between_graph}}.
#' @param inverse logical; if \code{TRUE}, also return the precision matrix.
#' @return A list with \code{Sigmab} (the \eqn{R \times R} covariance),
#'   \code{D_Sigmab} (list of derivatives, one per element of
#'   \code{theta_G}), \code{edge_pos} (positions of the edge parameters in
#'   \code{theta_G}) and, when \code{inverse = TRUE}, \code{Kb} (the precision
#'   matrix).
#' @seealso \code{\link{mc_between_graph}}, \code{\link{mc_precision_bases}}
#' @examples
#' sp <- mc_between_graph("ggm")
#' out <- mc_build_sigma_between_graph(c(2, 2, 2, -0.4, -0.3, 0), 3, sp)
#' out$Sigmab
#' @export
mc_build_sigma_between_graph <- function(theta_G, n_resp, spec,
                                         inverse = FALSE) {
  bld <- .make_builder(spec, n_resp)
  b <- bld(as.numeric(theta_G))
  out <- list(Sigmab = b$Sigma, D_Sigmab = b$Dlist, edge_pos = b$edge_pos)
  if (inverse) out$Kb <- b$K
  out
}

# ---- builder factory --------------------------------------------------------
# Returns a closure par -> list(Sigma, K, Dlist, edge_pos, R, npar, layer info)

.make_builder <- function(spec, R) {
  R <- as.integer(R)
  if (spec$type == "ggm") {
    bases <- mc_precision_bases(R, spec$edges)
    E <- nrow(bases$edges)
    edge_pos <- if (E > 0) (R + 1L):(R + E) else integer(0)
    fun <- function(par) {
      K <- cpp_build_precision(par, bases$cube)
      Sigma <- cpp_precision_to_cov(K)
      Dlist <- cpp_precision_derivatives(Sigma, bases$cube)
      list(Sigma = Sigma, K = K, Dlist = Dlist, edge_pos = edge_pos)
    }
    attr(fun, "meta") <- list(type = "ggm", bases = bases, R = R, E = E,
                              npar = R + E, edge_pos = edge_pos)
  } else {
    dag <- .dag_setup(spec, R)
    eidx <- dag$edges_pos
    E <- nrow(eidx)
    edge_pos <- seq_len(E)
    fun <- function(par) {
      gamma <- par[seq_len(E)]
      delta <- par[E + seq_len(R)]
      cc <- cpp_cholesky_to_cov(gamma, delta, R, eidx)
      K_pos <- cpp_cholesky_to_precision(gamma, delta, R, eidx)
      der <- cpp_cholesky_derivatives(gamma, delta, R, eidx)
      inv <- dag$inverse
      Dlist <- lapply(c(der$D_gamma, der$D_delta),
                      function(M) M[inv, inv, drop = FALSE])
      list(Sigma = cc$Sigma_b[inv, inv, drop = FALSE],
           K = K_pos[inv, inv, drop = FALSE], Dlist = Dlist,
           edge_pos = edge_pos)
    }
    attr(fun, "meta") <- list(type = "dag", edges = eidx,
                              edges_original = dag$edges_original,
                              ordering = dag$ordering, R = R, E = E,
                              npar = E + R, edge_pos = edge_pos)
  }
  fun
}

# ---- initial values ---------------------------------------------------------

.init_par <- function(spec, R, Shat) {
  Kfull <- .safe_solve(Shat)
  if (spec$type == "ggm") {
    bases <- mc_precision_bases(R, spec$edges)
    par <- numeric(length(bases$names))
    for (k in seq_len(R)) par[k] <- Kfull[k, k]
    E <- nrow(bases$edges)
    for (e in seq_len(E)) {
      i <- bases$edges[e, 1L]; j <- bases$edges[e, 2L]
      par[R + e] <- Kfull[i, j]
    }
    # ensure SPD start
    if (!cpp_check_spd(cpp_build_precision(par, bases$cube))$is_spd) {
      par <- c(diag(Kfull), rep(0, E))
    }
    par
  } else {
    dag <- .dag_setup(spec, R)
    eidx <- dag$edges_pos
    Shat_pos <- Shat[dag$ordering, dag$ordering, drop = FALSE]
    fac <- .spd_factor(Shat_pos)
    Lc <- t(fac$chol)                         # Shat_pos(reg.) = Lc Lc^T
    dL <- diag(Lc)
    L0 <- Lc %*% diag(1 / dL, R)              # unit lower
    T0 <- solve(L0)                           # T = L0^{-1}
    gamma <- numeric(nrow(eidx))
    for (e in seq_len(nrow(eidx))) {
      i <- eidx[e, 1L]; j <- eidx[e, 2L]
      gamma[e] <- -T0[j, i]
    }
    delta <- dL^2
    out <- c(gamma, delta)
    attr(out, "init_info") <- list(
      spd_adjustment = fac$adjustment,
      min_eigen_before = fac$min_eigen_before,
      eigen_floor = fac$eigen_floor
    )
    out
  }
}

# ---- Newton engine for the between-response covariance ----------------------
# Solves the second-order Pearson estimating equations for Sigma_b given the
# sample covariance Shat of the Pearson residuals (graph block of the two-stage
# fit). Returns parameters, Sigma_b, K_b, sensitivity and vcov.

.fit_between <- function(Shat, n, spec, par0 = NULL,
                         max_iter = 100L, tol = 1e-8,
                         score_tol = sqrt(tol),
                         rcond_tol = .Machine$double.eps^(2 / 3),
                         verbose = FALSE,
                         dag_solver = c("closed", "newton")) {
  dag_solver <- match.arg(dag_solver)
  R <- nrow(Shat)
  bld <- .make_builder(spec, R)
  meta <- attr(bld, "meta")
  # The closed route needs no starting value, so the spectral regularization
  # that produces one is not run on that path.
  closed_dag <- identical(meta$type, "dag") && identical(dag_solver, "closed") &&
    is.null(par0)
  par <- if (!is.null(par0)) par0 else
    if (closed_dag) numeric(meta$npar) else .init_par(spec, R, Shat)
  init_info <- attr(par, "init_info", exact = TRUE)
  attr(par, "init_info") <- NULL

  conv <- FALSE
  status <- "max_iter"
  linear_solve_status <- "not_attempted"
  linear_solve_ridge <- 0

  # The directed layer has a closed-form root: the estimating equations
  # separate by node and each block is a least-squares problem. Taking it
  # removes the iteration, the step damping and the line search from the path,
  # at the cost of nothing, since it solves the same system exactly. It is the
  # default route; dag_solver = "newton" selects the iteration below.
  if (meta$type == "dag" && dag_solver == "closed") {
    eidx <- meta$edges
    Shat_pos <- Shat[meta$ordering, meta$ordering, drop = FALSE]
    root <- cpp_dag_closed_root(Shat_pos, R, eidx, rcond_tol = rcond_tol)
    par <- c(as.numeric(root$gamma), as.numeric(root$delta))
    it <- 0L
    conv <- TRUE
    status <- "closed_form"
    linear_solve_status <- if (as.integer(root$n_ridge) > 0L)
      "ridge_regularized" else "exact"
    if (as.integer(root$n_floored) > 0L) {
      conv <- FALSE
      status <- "innovation_variance_floored"
    }
    # A near-singular parent block leaves the coefficients of that node
    # unidentified: the ridge root is kept, but it is not a certified root.
    ill_nodes <- as.integer(root$ill_conditioned)
    if (length(ill_nodes)) {
      conv <- FALSE
      lab <- colnames(Shat)
      if (is.null(lab)) lab <- as.character(seq_len(R))
      status <- paste0("parent_block_ill_conditioned: ",
                       paste(lab[meta$ordering[ill_nodes]], collapse = ", "),
                       if (as.integer(root$n_floored) > 0L)
                         "; innovation_variance_floored" else "")
    }
  } else
  for (it in seq_len(max_iter)) {
    b <- bld(par)
    pj <- cpp_pearson_between(b$Sigma, Shat, n, b$Dlist)
    if (meta$type == "ggm") {
      Jinv <- tryCatch(.safe_solve(pj$J, rcond_tol = rcond_tol),
                       error = function(e) NULL)
      step <- if (is.null(Jinv)) NULL else as.numeric(Jinv %*% pj$psi)
      if (is.null(step) || any(!is.finite(step))) {
        status <- "linear_solve_failed"
        if (verbose) message("GGM sensitivity solve failed; stopping at the last admissible iterate.")
        break
      }
      linear_solve_status <- attr(Jinv, "inversion_status", exact = TRUE)
      iter_ridge <- attr(Jinv, "inversion_ridge", exact = TRUE)
      if (is.null(iter_ridge)) iter_ridge <- 0
      linear_solve_ridge <- max(
        linear_solve_ridge, iter_ridge
      )
      ls <- cpp_spd_line_search(par, step, meta$bases$cube, 30L)
      if (!isTRUE(ls$is_spd)) {
        status <- "spd_line_search_failed"
        if (verbose)
          message("SPD line search failed; stopping at the last admissible iterate.")
        break
      }
      par_new <- ls$eta_new
    } else {
      # DAG: Newton step in the unconstrained log-innovation parametrization
      # xi = log(delta), so delta = exp(xi) > 0 and Sigma_b stays positive
      # definite by construction for any (gamma, xi), so no boundary projection
      # on delta is needed. The estimating function and sensitivity are mapped
      # to the xi scale by the diagonal Jacobian d(delta)/d(xi) = diag(delta).
      delta <- par[meta$E + seq_len(R)]
      scl <- c(rep(1, meta$E), delta)
      psi_xi <- pj$psi * scl
      J_xi <- pj$J * outer(scl, scl)
      Jinv <- tryCatch(.safe_solve(J_xi, rcond_tol = rcond_tol),
                       error = function(e) NULL)
      step_xi <- if (is.null(Jinv)) NULL else as.numeric(Jinv %*% psi_xi)
      if (is.null(step_xi) || any(!is.finite(step_xi))) {
        status <- "linear_solve_failed"
        if (verbose) message("DAG sensitivity solve failed; stopping at the last admissible iterate.")
        break
      }
      linear_solve_status <- attr(Jinv, "inversion_status", exact = TRUE)
      iter_ridge <- attr(Jinv, "inversion_ridge", exact = TRUE)
      if (is.null(iter_ridge)) iter_ridge <- 0
      linear_solve_ridge <- max(
        linear_solve_ridge, iter_ridge
      )
      cur_xi <- c(par[seq_len(meta$E)], log(delta))
      smax <- max(abs(step_xi))
      damp <- if (smax > 2) 2 / smax else 1     # cap |d xi| for stability
      new_xi <- cur_xi + damp * step_xi
      par_new <- c(new_xi[seq_len(meta$E)], exp(new_xi[meta$E + seq_len(R)]))
    }
    rel <- max(abs(par_new - par) / (1 + abs(par)))
    par <- par_new
    if (verbose) message(sprintf("iter %d: rel=%.2e", it, rel))
    if (rel < tol) {
      conv <- TRUE
      status <- "converged"
      break
    }
  }

  # The line search returns a one-column matrix and the closed root a plain
  # vector; both layers store a plain vector, so coef(fit, "graph") has one
  # shape.
  par <- as.numeric(par)
  b <- bld(par)
  # On the closed route the estimating function and the sensitivity have block
  # forms of their own, one block per node, so neither needs the generic
  # assembly of P^2 traces of R x R products. The two agree to rounding, which
  # the test suite checks; what changes is the cost, which stops growing with
  # the square of the number of parameters.
  if (meta$type == "dag" && dag_solver == "closed") {
    eidx <- meta$edges
    ordp <- meta$ordering
    Shat_pos <- Shat[ordp, ordp, drop = FALSE]
    Sigma_pos <- b$Sigma[ordp, ordp, drop = FALSE]
    delta_hat <- par[meta$E + seq_len(R)]
    pj <- list(
      psi = cpp_dag_pearson_closed(Shat_pos, par[seq_len(meta$E)],
                                   delta_hat, n, eidx),
      J = cpp_dag_sensitivity(Sigma_pos, delta_hat, n, eidx)
    )
  } else {
    pj <- cpp_pearson_between(b$Sigma, Shat, n, b$Dlist)
  }
  J <- pj$J
  # The same block structure inverts block by block. The general solve is the
  # fallback, so a node whose parent block is singular goes through
  # .safe_solve() rather than silently producing a zero block.
  vcov <- NULL
  if (meta$type == "dag" && dag_solver == "closed") {
    binv <- cpp_dag_sensitivity_inverse(Sigma_pos, delta_hat, n, meta$edges)
    if (as.integer(binv$n_singular) == 0L && all(is.finite(binv$inverse))) {
      vcov <- binv$inverse
      attr(vcov, "inversion_status") <- "exact_blockwise"
      attr(vcov, "inversion_ridge") <- 0
    }
  }
  if (is.null(vcov)) vcov <- .safe_solve(J, rcond_tol = rcond_tol)
  inversion_status <- attr(vcov, "inversion_status", exact = TRUE)
  inversion_ridge <- attr(vcov, "inversion_ridge", exact = TRUE)
  names(par) <- .par_names(spec, meta)
  score_norm <- max(abs(as.numeric(pj$psi)))
  score_norm_unit <- score_norm / n
  rcond_J <- tryCatch(rcond(0.5 * (J + t(J))), error = function(e) NA_real_)
  # The parameter-change gate asks whether the last Newton step was small. A
  # closed-form solve took no step, so the gate is satisfied by construction and
  # the score gate below is what actually certifies the root.
  convergence_gates <- c(
    parameter_change = isTRUE(conv),
    score = is.finite(score_norm_unit) && score_norm_unit <= score_tol,
    sensitivity = is.finite(rcond_J) && rcond_J > rcond_tol
  )
  if (isTRUE(conv) && !all(convergence_gates)) {
    conv <- FALSE
    status <- if (!convergence_gates[["score"]]) "score_not_converged" else
      "sensitivity_ill_conditioned"
  }
  list(par = par, Sigma_b = b$Sigma, K_b = b$K, J = J, vcov = vcov,
       edge_pos = meta$edge_pos, npar = meta$npar, iter = it,
       converged = conv, status = status, meta = meta,
       init_info = init_info, Dlist = b$Dlist,
       score_norm = score_norm, score_norm_unit = score_norm_unit,
       rcond_J = rcond_J, score_tol = score_tol, rcond_tol = rcond_tol,
       convergence_gates = convergence_gates,
       stopping_rule = paste0(
         "relative parameter change < tol; score_norm_unit <= ",
         format(score_tol), "; rcond_J > ", format(rcond_tol)
       ),
       inversion_status = inversion_status,
       inversion_ridge = inversion_ridge,
       inference_ok = inversion_status %in% c("exact", "exact_blockwise") &&
         all(convergence_gates),
       linear_solve_status = linear_solve_status,
       linear_solve_ridge = linear_solve_ridge,
       solver = if (meta$type == "dag") dag_solver else "newton")
}

.par_names <- function(spec, meta) {
  if (meta$type == "ggm") return(meta$bases$names)
  # paste0() recycles a zero-length index into one name, so the empty
  # support needs its own branch.
  eo <- meta$edges_original
  gam <- if (nrow(eo) > 0L) paste0("gamma_", eo[, 2L], eo[, 1L]) else
    character(0)
  c(gam, paste0("delta_", meta$ordering))
}

# Convert response labels to positions in the requested DAG ordering. Edges are
# treated as unordered adjacencies and oriented from the earlier to the later
# response in that ordering. The numerical LDL kernels continue to work on the
# natural positional order; fitted matrices are mapped back to response order.
.dag_setup <- function(spec, R) {
  ordering <- spec$ordering
  if (is.null(ordering)) ordering <- seq_len(R)
  ordering <- as.integer(ordering)
  if (length(ordering) != R || !setequal(ordering, seq_len(R)))
    stop("ordering must be a permutation of seq_len(R).")

  # Same validator as the low-level LDL API: the contract must not depend on
  # whether the caller came through mcglm_graph() or through mc_build_T_matrix().
  edges <- .validate_edges(spec$edges, R)
  if (nrow(edges) > 0L) {
    pos <- match(seq_len(R), ordering)
    edges_pos <- cbind(pos[edges[, 1L]], pos[edges[, 2L]])
    edges_pos <- t(apply(edges_pos, 1L, sort))
    edges_pos <- edges_pos[order(edges_pos[, 1L], edges_pos[, 2L]),
                           , drop = FALSE]
  } else {
    edges_pos <- matrix(integer(0), ncol = 2L)
  }
  edges_original <- if (nrow(edges_pos) > 0L)
    cbind(ordering[edges_pos[, 1L]], ordering[edges_pos[, 2L]]) else
      matrix(integer(0), ncol = 2L)
  list(ordering = ordering, inverse = match(seq_len(R), ordering),
       edges_pos = edges_pos, edges_original = edges_original)
}

.safe_solve <- function(M, rcond_tol = .Machine$double.eps^(2 / 3),
                        ridge_rel = 1e-8) {
  M <- 0.5 * (as.matrix(M) + t(as.matrix(M)))
  if (!nrow(M) || nrow(M) != ncol(M) || any(!is.finite(M)))
    stop("Cannot invert a non-finite or non-square matrix.")
  rc <- tryCatch(rcond(M), error = function(e) NA_real_)
  out <- if (is.finite(rc) && rc > rcond_tol)
    tryCatch(solve(M), error = function(e) NULL) else NULL
  status <- "exact"
  ridge <- 0
  if (is.null(out) || any(!is.finite(out))) {
    scale <- max(abs(diag(M)), max(abs(M)), 1)
    ridge <- ridge_rel * scale
    out <- tryCatch(solve(M + diag(ridge, nrow(M))),
                    error = function(e) NULL)
    status <- "ridge_regularized"
  }
  if (is.null(out) || any(!is.finite(out)))
    stop("Matrix inversion failed even after scale-relative regularization.")
  out <- 0.5 * (out + t(out))
  attr(out, "inversion_status") <- status
  attr(out, "inversion_rcond") <- rc
  attr(out, "inversion_ridge") <- ridge
  out
}

# Deterministic spectral regularization used only to initialize the DAG layer.
# It preserves an already well-conditioned covariance and lifts eigenvalues
# below a scale-relative floor before the Cholesky factorization.
.spd_factor <- function(M, rel_floor = 1e-8) {
  M <- 0.5 * (as.matrix(M) + t(as.matrix(M)))
  if (!nrow(M) || nrow(M) != ncol(M) || any(!is.finite(M)))
    stop("The DAG initialization requires a finite square covariance matrix.")
  ee <- eigen(M, symmetric = TRUE)
  scale <- max(abs(ee$values))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  floor <- rel_floor * scale
  values <- pmax(ee$values, floor)
  adjusted <- ee$vectors %*% diag(values, nrow = length(values)) %*%
    t(ee$vectors)
  adjusted <- 0.5 * (adjusted + t(adjusted))
  list(
    chol = chol(adjusted),
    adjustment = max(values - ee$values),
    min_eigen_before = min(ee$values),
    eigen_floor = floor
  )
}
