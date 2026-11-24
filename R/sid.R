# sid.R: structural intervention distance.
# Part of the mcglmgraph package. License: GPL-3.

.check_dag <- function(A, arg) {
  A <- .check_amat(A, arg)
  if (any(A == 1L & t(A) == 1L))
    stop(arg, " must be a DAG: it carries an undirected edge.", call. = FALSE)
  if (.has_directed_cycle(A))
    stop(arg, " must be acyclic.", call. = FALSE)
  A
}

#' Structural intervention distance
#'
#' Counts the ordered pairs of responses for which the estimated graph would
#' send one to the wrong adjustment set (Peters and Buehlmann, 2015). For the
#' pair \eqn{(i, j)} the estimate prescribes adjusting for the parents of
#' \eqn{i} that it carries; the pair is counted when that set is not a valid
#' adjustment set for the effect of \eqn{i} on \eqn{j} in the true graph. The
#' distance is zero when the estimate gives a valid adjustment set for every
#' pair, and it is not symmetric in its two arguments.
#'
#' @section Why it is not the Hamming distance:
#' The structural Hamming distance counts the edges two graphs disagree about
#' and says nothing about what the disagreement costs. Two estimates at the same
#' Hamming distance from the truth can be worth very different things, one
#' still giving the right adjustment set for every pair and the other for none.
#' This distance asks that question instead.
#'
#' @section Which adjustment sets count as valid:
#' Validity is decided by the adjustment criterion of Shpitser, VanderWeele and
#' Robins (2010), in the form given by Perkovic et al. (2018), which is
#' necessary and sufficient, and not by the back-door criterion of Pearl
#' (2009). The back-door criterion is sufficient without being necessary:
#' it refuses any set containing a descendant of \eqn{i}, while a descendant
#' that lies on no causal path from \eqn{i} to \eqn{j} does no harm, so using
#' it would count pairs the adjustment in fact recovers. Writing
#' \eqn{\mathrm{cn}(i,j)} for the nodes other than \eqn{i} on the directed
#' paths from \eqn{i} to \eqn{j} in the true graph, and
#' \eqn{\mathrm{forb}(i,j)} for their descendants (each node counted among its
#' own) together with \eqn{i}, the parent set \eqn{S} is valid when it does not
#' meet \eqn{\mathrm{forb}(i,j)} and when it d-separates \eqn{i} from \eqn{j}
#' in the graph with the arcs from \eqn{i} into \eqn{\mathrm{cn}(i,j)} removed. With
#' no causal path from \eqn{i} to \eqn{j} the second condition asks for
#' separation in the whole graph, which is what recovers a null effect.
#'
#' One case needs care. When \eqn{j} itself belongs to the adjustment set the
#' formula collapses to the marginal of \eqn{j}, since conditioning on it and
#' averaging over it returns it unchanged. That is right exactly when the effect
#' is null, and wrong otherwise. Reading it as always wrong would make the
#' distance of a graph to itself nonzero, because a parent of \eqn{i} always
#' sits in its adjustment set and is never its descendant.
#'
#' @section A partially directed estimate gives an interval:
#' A completed partially directed acyclic graph (CPDAG) does not prescribe one
#' adjustment set, so the distance is not one number but the range over the
#' members of the class. The members are enumerated with
#' \code{\link{mc_class_members}} and the result reports the lower and upper end,
#' the member attaining each, and whether the enumeration had to be capped.
#' When it was capped, the interval covers the enumerated members only.
#'
#' @section What the distance measures:
#' The reading in terms of interventions, in which a valid adjustment set
#' recovers the effect of one response on another, holds when the true graph
#' is a causal DAG: a structural model with independent exogenous innovations,
#' no unobserved confounding and no selection bias, whose equations are
#' invariant under intervention. Without those assumptions the distance remains
#' a purely graphical measure of disagreement between the two structures,
#' computed from the parent sets of the estimate. It compares graphs and uses
#' no data.
#'
#' @param truth the true DAG, as a square 0/1 adjacency matrix where a one in
#'   position \eqn{(a, b)} means \eqn{a \to b}. It must be acyclic and have no
#'   undirected edge.
#' @param estimated the estimated graph with the same number of nodes and the
#'   same convention: a DAG, a CPDAG (undirected edges carry a one in both
#'   positions), or an object of class \code{"mcg_cpdag"}.
#' @param max_members cap on the number of members enumerated for a partially
#'   directed estimate.
#' @return For a directed estimate, a single number, the count of offending
#'   pairs, with the \eqn{R \times R} 0/1 matrix of the offending ordered
#'   pairs attached as the attribute \code{"pairs"} (a one in \eqn{(i, j)}
#'   marks the pair \eqn{(i, j)}). For a partially directed one, an object of
#'   class \code{"mcg_sid"}, a list with \code{lower} and \code{upper}, the
#'   smallest and largest distance over the members; \code{values}, the
#'   distance of each member; \code{n_members}; \code{which_lower} and
#'   \code{which_upper}, the positions in the member list (as returned by
#'   \code{\link{mc_class_members}}) attaining each end; and \code{truncated}.
#'   The \code{print} method returns \code{x} invisibly.
#' @references Peters, J. and Buehlmann, P. (2015). Structural intervention
#'   distance for evaluating causal graphs. \emph{Neural Computation} 27(3),
#'   771-799. \doi{10.1162/NECO_a_00708}
#'
#'   Shpitser, I., VanderWeele, T. and Robins, J. M. (2010). On the validity of
#'   covariate adjustment for estimating causal effects. \emph{Proceedings of
#'   the 26th Conference on Uncertainty in Artificial Intelligence}, 527-536.
#'
#'   Perkovic, E., Textor, J., Kalisch, M. and Maathuis, M. H. (2018). Complete
#'   graphical characterization and construction of adjustment sets in Markov
#'   equivalence classes of ancestral graphs. \emph{Journal of Machine Learning
#'   Research} 18(220), 1-62.
#'
#'   Pearl, J. (2009). \emph{Causality: Models, Reasoning and Inference}, 2nd
#'   edition. Cambridge University Press.
#' @examples
#' # truth 1 -> 2 -> 3; an estimate that reverses the first arc
#' G <- matrix(0L, 3, 3); G[1, 2] <- G[2, 3] <- 1L
#' H <- matrix(0L, 3, 3); H[2, 1] <- H[2, 3] <- 1L
#' mc_sid(G, H)
#' mc_sid(G, G)          # zero against itself
#' @seealso \code{\link{mc_cpdag_metrics}} for the Hamming-type measures,
#'   \code{\link{mc_class_members}} for the members of a class.
#' @family equivalence class
#' @export
mc_sid <- function(truth, estimated, max_members = 1000L) {
  G <- .check_dag(truth, "truth")
  if (inherits(estimated, "mcg_cpdag")) estimated <- estimated$cpdag
  H <- .check_amat(estimated, "estimated")
  if (nrow(G) != nrow(H))
    stop("the two graphs describe a different number of responses.",
         call. = FALSE)

  undirected <- any(H == 1L & t(H) == 1L)
  if (!undirected) {
    H <- .check_dag(H, "estimated")
    res <- cpp_sid(G, H)
    out <- as.numeric(res$distance)
    pairs <- matrix(as.integer(res$pairs), nrow(G), nrow(G))
    dimnames(pairs) <- dimnames(G)
    attr(out, "pairs") <- pairs
    return(out)
  }

  members <- mc_class_members(H, max_members = max_members)
  if (!length(members))
    stop("the estimated class has no consistent extension.", call. = FALSE)
  vals <- vapply(members, function(A)
    as.numeric(cpp_sid(G, A)$distance), numeric(1))
  structure(list(lower = min(vals), upper = max(vals), values = vals,
                 n_members = length(members),
                 which_lower = which.min(vals), which_upper = which.max(vals),
                 truncated = isTRUE(attr(members, "truncated"))),
            class = "mcg_sid")
}

#' @rdname mc_sid
#' @param x an object of class \code{"mcg_sid"}.
#' @param ... ignored.
#' @export
print.mcg_sid <- function(x, ...) {
  cat(sprintf("Structural intervention distance over %d member%s%s\n",
              x$n_members, if (x$n_members == 1L) "" else "s",
              if (x$truncated) " (enumeration truncated)" else ""))
  cat(sprintf("Interval: [%g, %g]\n", x$lower, x$upper))
  invisible(x)
}
