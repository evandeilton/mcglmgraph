# ggm_layer.R: Gaussian Graphical Model layer (R wrappers over C++ kernels)
# Part of the mcglmgraph package. License: GPL-3.

#' Build the base matrices of the GGM precision predictor
#'
#' Constructs the base matrices \eqn{Z_d} of the linear-in-precision predictor
#' of the Gaussian graphical model (GGM) layer,
#' \deqn{K_b = \sum_{d=1}^{D} \eta_d Z_d,\qquad \Sigma_b = K_b^{-1},}
#' where \eqn{K_b} is the between-response precision matrix and
#' \eqn{D = R + |E|}. There are \eqn{R} diagonal bases
#' \eqn{Z_{kk} = e_k e_k^\top} and one symmetric base
#' \eqn{Z_{ij} = e_i e_j^\top + e_j e_i^\top} per edge \eqn{(i, j)}, so that
#' \eqn{\eta_{kk} = [K_b]_{kk}} and \eqn{\eta_{ij} = [K_b]_{ij}}. The pair
#' \eqn{(i, j)} is conditionally independent given the other responses, in the
#' Gaussian sense, exactly when it is outside the edge set (Dempster, 1972;
#' Lauritzen, 1996).
#'
#' @param R integer, number of responses (\eqn{R \ge 2}).
#' @param edges optional integer matrix with two columns giving the edges
#'   \eqn{(i, j)}, \eqn{i < j}, as response positions in \code{1:R}. If
#'   \code{NULL}, the complete graph is used; a matrix with no rows is the
#'   empty graph. The edges are stored sorted in lexicographic order, which is
#'   the order of the edge parameters.
#' @return An object of class \code{"mcg_bases"}: a list with components
#'   \code{bases} (list of \eqn{D} symmetric \eqn{R \times R} matrices),
#'   \code{cube} (the same bases as an \eqn{R \times R \times D} array, the
#'   form consumed by the C++ kernels), \code{names} (parameter names
#'   \code{eta_kk} and \code{eta_ij}), \code{type} (\code{"diag"} or
#'   \code{"edge"} for each base), \code{edges} (the canonical edge matrix) and
#'   \code{R}.
#' @references
#' Dempster, A. P. (1972). Covariance selection. \emph{Biometrics}, 28(1),
#' 157-175. \doi{10.2307/2528966}
#'
#' Lauritzen, S. L. (1996). \emph{Graphical Models}. Oxford University Press.
#' @seealso \code{\link{mc_build_precision}}, \code{\link{mc_between_graph}}
#' @family GGM layer
#' @examples
#' b <- mc_precision_bases(3)
#' length(b$bases)            # 3 diagonal + 3 edges = 6
#' mc_precision_bases(3, edges = rbind(c(1, 2), c(2, 3)))$names
#' @export
mc_precision_bases <- function(R, edges = NULL) {
  R <- .dag_dimension(R)
  if (R < 2L) stop("R must be at least 2.", call. = FALSE)
  edges <- .check_edges(edges, R)
  E <- nrow(edges)
  D <- R + E
  cube <- array(0, dim = c(R, R, D))
  nms <- character(D)
  typ <- character(D)
  for (k in seq_len(R)) {
    cube[k, k, k] <- 1
    nms[k] <- paste0("eta_", k, k)
    typ[k] <- "diag"
  }
  for (e in seq_len(E)) {
    i <- edges[e, 1L]; j <- edges[e, 2L]
    cube[i, j, R + e] <- 1
    cube[j, i, R + e] <- 1
    nms[R + e] <- paste0("eta_", i, j)
    typ[R + e] <- "edge"
  }
  bases <- lapply(seq_len(D), function(d) cube[, , d, drop = TRUE])
  structure(list(bases = bases, cube = cube, names = nms, type = typ,
                 edges = edges, R = R),
            class = "mcg_bases")
}

#' Assemble the between-response precision matrix
#'
#' Computes \eqn{K_b = \sum_d \eta_d Z_d} for a parameter vector \eqn{\eta}
#' and the bases \eqn{Z_d} of \code{\link{mc_precision_bases}}. The result is
#' symmetrized; positive definiteness is not checked here (see
#' \code{\link{mc_check_spd}}).
#'
#' @param eta numeric vector of length \eqn{D = R + |E|}: the diagonal
#'   parameters first, then the edge parameters, in the order of \code{bases}.
#' @param bases an object of class \code{"mcg_bases"} from
#'   \code{\link{mc_precision_bases}}.
#' @return The \eqn{R \times R} symmetric precision matrix \eqn{K_b}.
#' @family GGM layer
#' @examples
#' b <- mc_precision_bases(3)
#' mc_build_precision(c(2, 2, 2, -0.4, -0.3, 0), b)
#' @export
mc_build_precision <- function(eta, bases) {
  stopifnot(inherits(bases, "mcg_bases"))
  cpp_build_precision(as.numeric(eta), bases$cube)
}

#' Convert a precision matrix into the between-response covariance
#'
#' Inverts a symmetric positive-definite precision matrix,
#' \eqn{\Sigma_b = K_b^{-1}}, and reports its log-determinant and condition
#' number. An error is raised when \code{K} is not positive definite in the
#' sense of \code{\link{mc_check_spd}}, that is, when its smallest eigenvalue
#' does not exceed \code{1e-10}.
#'
#' @param K symmetric positive-definite precision matrix.
#' @return A list with \code{Sigma_b} (the covariance \eqn{K^{-1}}),
#'   \code{logdet} (log-determinant of \eqn{K}) and \code{cond} (condition
#'   number of \eqn{K}, the ratio of its largest to its smallest eigenvalue).
#' @family GGM layer
#' @examples
#' K <- mc_build_precision(c(2, 2, 2, -0.4, -0.3, 0), mc_precision_bases(3))
#' mc_precision_to_covariance(K)$Sigma_b
#' @export
mc_precision_to_covariance <- function(K) {
  spd <- cpp_check_spd(K)
  if (!spd$is_spd) stop("Precision matrix is not positive definite.")
  Sigma <- cpp_precision_to_cov(K)
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  list(Sigma_b = Sigma, logdet = sum(log(ev)), cond = max(ev) / min(ev))
}

#' Analytical derivatives of the covariance with respect to eta
#'
#' Returns the derivatives of \eqn{\Sigma_b = K_b^{-1}} with respect to the
#' precision parameters,
#' \deqn{\frac{\partial \Sigma_b}{\partial \eta_d} = -\Sigma_b Z_d \Sigma_b,
#'   \qquad d = 1, \dots, D,}
#' which follow from \eqn{\partial K_b^{-1} = -K_b^{-1}(\partial K_b)K_b^{-1}}
#' and \eqn{\partial K_b/\partial\eta_d = Z_d}. Each matrix is symmetrized.
#'
#' @param Sigma_b the between-response covariance matrix \eqn{\Sigma_b}.
#' @param bases an object of class \code{"mcg_bases"}.
#' @return A list of \eqn{D} symmetric \eqn{R \times R} matrices in the order
#'   of \code{bases}.
#' @family GGM layer
#' @examples
#' b <- mc_precision_bases(3)
#' S <- mc_precision_to_covariance(mc_build_precision(c(2,2,2,-.4,-.3,0), b))$Sigma_b
#' length(mc_precision_derivatives(S, b))
#' @export
mc_precision_derivatives <- function(Sigma_b, bases) {
  stopifnot(inherits(bases, "mcg_bases"))
  cpp_precision_derivatives(as.matrix(Sigma_b), bases$cube)
}

#' Check positive-definiteness of a matrix
#'
#' A symmetric matrix is declared positive definite when the smallest
#' eigenvalue of its symmetric part \eqn{(K + K^\top)/2} exceeds \code{tol}.
#'
#' @param K square numeric matrix; only its symmetric part is used.
#' @param tol tolerance on the smallest eigenvalue (default \code{1e-10}).
#' @return A list with \code{is_spd} (logical, \code{TRUE} when the smallest
#'   eigenvalue exceeds \code{tol}) and \code{min_eig} (the smallest
#'   eigenvalue).
#' @family GGM layer
#' @examples
#' mc_check_spd(diag(3))
#' @export
mc_check_spd <- function(K, tol = 1e-10) {
  cpp_check_spd(as.matrix(K), tol)
}

#' Backtracking line search preserving positive-definiteness
#'
#' Given a Newton direction \eqn{d}, finds the largest step
#' \deqn{\alpha = \max\{2^{-k}: k = 0, 1, \dots, k_{\max},\;
#'   \lambda_{\min}(K_b(\eta + 2^{-k} d)) > 10^{-10}\},}
#' with \eqn{\lambda_{\min}} the smallest eigenvalue and \eqn{k_{\max}} given
#' by \code{max_halvings}. The positivity test uses the numerical threshold
#' \eqn{10^{-10}} in place of zero. Only positive definiteness is checked, not
#' a decrease of a merit function. By Weyl's inequality,
#' \eqn{\lambda_{\min}(K_b + \alpha\,\Delta K) \ge \lambda_{\min}(K_b) -
#' \alpha\,\|\Delta K\|_2}, with \eqn{\Delta K = \sum_d d_d Z_d}, so an
#' admissible step exists for small enough \eqn{\alpha} whenever
#' \eqn{K_b(\eta)} passes the test; it may still lie below
#' \eqn{2^{-k_{\max}}}, in which case the search fails. The fitting engine of
#' \code{\link{mcglm_graph}} calls the same search with 30 halvings.
#'
#' @param eta_old current parameter vector \eqn{\eta}, with \eqn{K_b(\eta)}
#'   positive definite.
#' @param delta_eta the Newton direction \eqn{d}, of the same length.
#' @param bases an object of class \code{"mcg_bases"} from
#'   \code{\link{mc_precision_bases}}.
#' @param max_halvings maximum number \eqn{k_{\max}} of step halvings
#'   (default 20), so that at most \code{max_halvings + 1} steps are tried.
#' @return A list with components
#'   \describe{
#'     \item{\code{eta_new}}{a one-column matrix with the last point tried,
#'       \eqn{\eta + \alpha d}: the accepted point when the search succeeds.}
#'     \item{\code{step_size}}{the accepted step \eqn{\alpha} when the search
#'       succeeds.}
#'     \item{\code{halvings}}{the number of halvings performed.}
#'     \item{\code{is_spd}}{logical, \code{TRUE} when an admissible step was
#'       found.}
#'   }
#'   When no admissible step is found, \code{is_spd} is \code{FALSE},
#'   \code{eta_new} is the point tried with the smallest step,
#'   \eqn{2^{-k_{\max}}}, which fails the test and must not be used,
#'   \code{halvings} is \code{max_halvings + 1} and \code{step_size} is
#'   \eqn{2^{-(k_{\max}+1)}}, half the step that produced \code{eta_new}.
#' @references
#' Horn, R. A. and Johnson, C. R. (2013). \emph{Matrix Analysis}, 2nd ed.
#' Cambridge University Press. \doi{10.1017/CBO9781139020411}
#' @seealso \code{\link{mc_check_spd}}
#' @family GGM layer
#' @examples
#' b <- mc_precision_bases(3)
#' mc_spd_line_search(c(2,2,2,0,0,0), c(-.5,-.5,-.5,-.4,-.4,-.4), b)$halvings
#' @export
mc_spd_line_search <- function(eta_old, delta_eta, bases, max_halvings = 20L) {
  stopifnot(inherits(bases, "mcg_bases"))
  cpp_spd_line_search(as.numeric(eta_old), as.numeric(delta_eta),
                      bases$cube, as.integer(max_halvings))
}

# ---- internal helpers -------------------------------------------------------

# The GGM layer shares the single validator of .validate_edges(), so the same
# malformed input is refused the same way in both layers.
.check_edges <- function(edges, R) .validate_edges(edges, R)
