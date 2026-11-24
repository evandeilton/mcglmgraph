# metrics.R: comparing an estimated structure with a known one.
# Part of the mcglmgraph package. License: GPL-3.
#
# Standard confusion-count definitions of structure recovery, over the pairs of
# responses.

# Edge indicators over the pairs in upper-triangular order, from whatever the
# caller has at hand: an adjacency, a logical vector already in that order, or a
# two-column matrix of edges.
.edge_indicator <- function(x, R = NULL, arg = "x") {
  if (is.logical(x) && is.null(dim(x))) return(as.integer(x))
  if (is.numeric(x) && is.null(dim(x)) && all(x %in% c(0, 1)))
    return(as.integer(x))
  x <- as.matrix(x)
  if (ncol(x) == 2L && (is.null(R) || nrow(x) != ncol(x) || nrow(x) != R)) {
    if (is.null(R))
      stop("R is needed to turn an edge list into an indicator.", call. = FALSE)
    ed <- .validate_edges(x, R, arg = arg)
    pr <- t(utils::combn(R, 2L))
    key <- paste(pr[, 1L], pr[, 2L])
    return(as.integer(key %in% paste(ed[, 1L], ed[, 2L])))
  }
  if (nrow(x) != ncol(x))
    stop(arg, " must be a square adjacency, an edge list, or an indicator.",
         call. = FALSE)
  A <- (x != 0) | (t(x) != 0)
  as.integer(A[upper.tri(A)])
}

.check_amat <- function(x, arg) {
  x <- as.matrix(x)
  if (!nrow(x) || nrow(x) != ncol(x))
    stop(arg, " must be a square adjacency matrix.", call. = FALSE)
  if (!is.numeric(x) || any(!x %in% c(0, 1)) || any(diag(x) != 0))
    stop(arg, " must hold zeros and ones with a zero diagonal.", call. = FALSE)
  storage.mode(x) <- "integer"
  x
}

#' Recovery of an undirected support
#'
#' Compares an estimated set of edges with a known one over the
#' \eqn{R(R-1)/2} pairs of responses and returns the usual counts and rates.
#' With true positives \eqn{TP}, false positives \eqn{FP}, false negatives
#' \eqn{FN} and true negatives \eqn{TN},
#' \deqn{\mathrm{TPR}=\frac{TP}{TP+FN},\quad
#'       \mathrm{FPR}=\frac{FP}{FP+TN},\quad
#'       \mathrm{precision}=\frac{TP}{TP+FP},\quad
#'       F_1=\frac{2TP}{2TP+FP+FN},}
#' and the structural Hamming distance is \eqn{\mathrm{SHD}=FP+FN}. The
#' comparison ignores direction: an adjacency is reduced to its skeleton first
#' (a pair is an edge when either entry is nonzero).
#'
#' @param estimated,truth an adjacency matrix, a logical (or 0/1) indicator
#'   over the pairs in upper-triangular order, or a two-column matrix of edges.
#' @param R number of responses, needed only when an argument is an edge list;
#'   when omitted it is taken from \code{truth} if that is a square matrix.
#' @return A named numeric vector with \code{tpr}, \code{fpr},
#'   \code{precision}, \code{f1}, the structural Hamming distance \code{shd},
#'   and the counts \code{tp}, \code{fp}, \code{fn} and \code{tn}. A rate whose
#'   denominator is empty is \code{NA}: precision is undefined when nothing was
#'   selected, the true positive rate when there was nothing to find, and the
#'   false positive rate when the true graph is complete. \eqn{F_1} is
#'   \code{NA} only when both edge sets are empty; it is zero when nothing was
#'   selected and the truth has edges.
#' @examples
#' truth <- rbind(c(1, 2), c(2, 3))
#' est <- rbind(c(1, 2), c(1, 3))
#' mc_structure_metrics(est, truth, R = 3)
#' @seealso \code{\link{mc_cpdag_metrics}} for partially directed structures.
#' @family structure recovery
#' @export
mc_structure_metrics <- function(estimated, truth, R = NULL) {
  if (is.null(R)) {
    d <- dim(as.matrix(truth))
    if (!is.null(d) && length(d) == 2L && d[1L] == d[2L] && d[1L] > 1L)
      R <- d[1L]
  }
  e <- .edge_indicator(estimated, R, "estimated")
  t <- .edge_indicator(truth, R, "truth")
  if (length(e) != length(t))
    stop("the two structures describe a different number of pairs.",
         call. = FALSE)
  cnt <- cpp_confusion(as.integer(e), as.integer(t))
  TP <- cnt[1L]; FP <- cnt[2L]; FN <- cnt[3L]; TN <- cnt[4L]
  P <- TP + FN; N <- FP + TN
  c(tpr = if (P) TP / P else NA_real_,
    fpr = if (N) FP / N else NA_real_,
    precision = if (TP + FP) TP / (TP + FP) else NA_real_,
    f1 = if (2 * TP + FP + FN) 2 * TP / (2 * TP + FP + FN) else NA_real_,
    shd = FP + FN, tp = TP, fp = FP, fn = FN, tn = TN)
}

#' Recovery of a partially directed structure
#'
#' Compares an estimated completed partially directed acyclic graph (CPDAG)
#' with a known one, separating what the skeleton recovers from what the
#' orientation does.
#'
#' @section What each measure counts:
#' Every pair of responses is in one of four states: absent, undirected, or
#' directed either way. \code{cpdag_shd} counts the pairs whose states differ.
#' The skeleton measures (\code{cpdag_skeleton_f1}, \code{cpdag_skeleton_shd})
#' are those of \code{\link{mc_structure_metrics}} and ignore direction.
#' \code{cpdag_orient_recall} is the proportion of the pairs the truth directs
#' that the estimate carries in the same state, and \code{cpdag_overorient} the
#' proportion of the pairs the truth leaves undirected while joining them that
#' the estimate directs, which is the rate at which a direction is invented. A
#' v-structure is a triple \eqn{i \to k \leftarrow j} with \eqn{i} and \eqn{j}
#' nonadjacent. V-structure recall is the proportion of the true v-structures
#' that the estimate declares, and precision the proportion of the declared ones
#' that are true; recall is conditional on the truth containing one and
#' precision on some v-structure having been declared, and either is \code{NA}
#' when its condition is empty.
#'
#' @param estimated,truth square adjacency matrices of the same dimension with
#'   zeros and ones and a zero diagonal, in the package convention: entries
#'   \eqn{(i,j)} and \eqn{(j,i)} both set mean an undirected edge, and only
#'   \eqn{(i,j)} set means the arc \eqn{i \to j}.
#' @return A named numeric vector with \code{cpdag_shd},
#'   \code{cpdag_skeleton_f1}, \code{cpdag_skeleton_shd},
#'   \code{cpdag_vstructure_recall}, \code{cpdag_vstructure_precision},
#'   \code{cpdag_orient_recall} and \code{cpdag_overorient}.
#' @examples
#' truth <- matrix(0L, 3, 3)
#' truth[1, 3] <- truth[2, 3] <- 1L          # 1 -> 3 <- 2
#' est <- truth
#' mc_cpdag_metrics(est, truth)
#' @seealso \code{\link{mc_true_cpdag}} to build the true CPDAG of a DAG,
#'   \code{\link{mc_structure_metrics}}.
#' @family structure recovery
#' @export
mc_cpdag_metrics <- function(estimated, truth) {
  E <- .check_amat(estimated, "estimated")
  T <- .check_amat(truth, "truth")
  if (nrow(E) != nrow(T))
    stop("the two structures describe a different number of responses.",
         call. = FALSE)
  es <- cpp_pair_states(E)
  ts <- cpp_pair_states(T)
  sm <- mc_structure_metrics(as.integer(es != 1L), as.integer(ts != 1L))

  true_edge <- ts != 1L
  true_dir <- ts %in% c(3L, 4L)
  est_dir <- es %in% c(3L, 4L)

  key <- function(M) {
    v <- cpp_vstructure_set(M)
    if (!nrow(v)) return(character(0))
    unique(paste(v[, 1L], v[, 2L], v[, 3L], sep = "->"))
  }
  vt <- key(T); ve <- key(E)

  c(cpdag_shd = sum(es != ts),
    cpdag_skeleton_f1 = unname(sm[["f1"]]),
    cpdag_skeleton_shd = unname(sm[["shd"]]),
    cpdag_vstructure_recall = if (length(vt)) mean(vt %in% ve) else NA_real_,
    cpdag_vstructure_precision = if (length(ve)) mean(ve %in% vt) else NA_real_,
    cpdag_orient_recall = if (any(true_dir))
      mean(es[true_dir] == ts[true_dir]) else NA_real_,
    cpdag_overorient = if (any(!true_dir & true_edge))
      mean(est_dir[!true_dir & true_edge]) else NA_real_)
}

#' CPDAG of a known directed structure
#'
#' Builds the completed partially directed acyclic graph (CPDAG) of the
#' equivalence class of a directed acyclic graph (DAG) given through the
#' structural matrix of the
#' modified Cholesky parametrization, where a nonzero \eqn{T_{ji}} codes the
#' arc \eqn{i \to j}. The equivalence class is the set of DAGs with the same
#' skeleton and the same v-structures (Verma and Pearl, 1990; Andersson,
#' Madigan and Perlman, 1997). The v-structures are written first and Meek's
#' rules complete the orientation, which gives the object an estimated CPDAG
#' should be compared against.
#'
#' @section Arcs anywhere, not only below the diagonal:
#' A structural matrix of the modified Cholesky parametrization is lower
#' triangular, so its arcs run from a lower index to a higher one and reading
#' only the lower triangle would suffice. A DAG that arrives from elsewhere need
#' not be numbered in a topological order, and its arcs can sit anywhere, so
#' every off-diagonal entry is read. For a triangular input the two readings
#' agree. A pair carrying an arc in both directions, or a set of arcs with a
#' cycle, is refused.
#'
#' @param T a finite numeric square matrix whose entry \eqn{(j, i)} is nonzero
#'   when \eqn{i \to j}; the structural matrix of the directed layer is one
#'   such, and the transpose of a directed adjacency is another.
#' @param tol non-negative number; entries not above it in absolute value count
#'   as zero (default 0).
#' @return The adjacency matrix of the CPDAG, in the package convention (see
#'   \code{\link{mc_cpdag_metrics}}), with the dimnames of \code{T}.
#' @references Verma, T. and Pearl, J. (1990). Equivalence and synthesis of
#'   causal models. \emph{Proceedings of the Sixth Conference on Uncertainty in
#'   Artificial Intelligence}, 220--227.
#'
#'   Andersson, S. A., Madigan, D. and Perlman, M. D. (1997). A
#'   characterization of Markov equivalence classes for acyclic digraphs.
#'   \emph{The Annals of Statistics} 25(2), 505--541.
#'   \doi{10.1214/aos/1031833662}
#' @examples
#' T <- diag(3)
#' T[3, 1] <- -0.5; T[3, 2] <- -0.4        # 1 -> 3 <- 2, a collider
#' mc_true_cpdag(T)
#' @seealso \code{\link{mc_cpdag_metrics}}.
#' @family structure recovery
#' @export
mc_true_cpdag <- function(T, tol = 0) {
  T <- as.matrix(T)
  if (!nrow(T) || nrow(T) != ncol(T) || !is.numeric(T) || any(!is.finite(T)))
    stop("T must be a finite numeric square matrix.", call. = FALSE)
  if (length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("tol must be one finite non-negative number.", call. = FALSE)
  R <- nrow(T)
  dir <- matrix(0L, R, R)
  for (j in seq_len(R)) for (i in seq_len(R))
    if (i != j && abs(T[j, i]) > tol) dir[i, j] <- 1L
  if (any(dir == 1L & t(dir) == 1L))
    stop("T carries an arc in both directions between two responses.",
         call. = FALSE)
  if (.has_directed_cycle(dir))
    stop("the arcs of T contain a directed cycle.", call. = FALSE)
  skel <- ((dir + t(dir)) > 0) * 1L
  amat <- skel
  for (k in seq_len(R)) {
    pa <- which(dir[, k] == 1L)
    if (length(pa) < 2L) next
    cmb <- utils::combn(pa, 2L)
    for (q in seq_len(ncol(cmb))) {
      i <- cmb[1L, q]; j <- cmb[2L, q]
      if (skel[i, j] == 0L) {
        amat[i, k] <- 1L; amat[k, i] <- 0L
        amat[j, k] <- 1L; amat[k, j] <- 0L
      }
    }
  }
  out <- mc_orient_meek(amat)
  attributes(out) <- list(dim = dim(out), dimnames = dimnames(T))
  out
}
