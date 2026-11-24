# ges.R: greedy equivalence search over Markov equivalence classes.
# Part of the mcglmgraph package. License: GPL-3.

# The CPDAG of a directed graph: write the v-structures, then propagate.
.dag_to_cpdag <- function(dag) {
  R <- nrow(dag)
  skel <- ((dag + t(dag)) > 0) * 1L
  amat <- skel
  for (k in seq_len(R)) {
    pa <- which(dag[, k] == 1L & dag[k, ] == 0L)
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
  attributes(out) <- list(dim = dim(out))
  out
}

# Back to a CPDAG after an operator has left a partially directed graph.
.pdag_to_cpdag <- function(amat) {
  storage.mode(amat) <- "integer"
  res <- cpp_pdag_to_dag(amat)
  if (!isTRUE(res$ok)) return(NULL)
  .dag_to_cpdag(matrix(as.integer(res$dag), nrow(amat), ncol(amat)))
}

.apply_insert <- function(amat, x, y, set) {
  amat[x, y] <- 1L; amat[y, x] <- 0L
  for (t in set) { amat[t, y] <- 1L; amat[y, t] <- 0L }
  amat
}

.apply_delete <- function(amat, x, y, set) {
  amat[x, y] <- 0L; amat[y, x] <- 0L
  for (h in set) {
    amat[y, h] <- 1L; amat[h, y] <- 0L
    amat[x, h] <- 1L; amat[h, x] <- 0L
  }
  amat
}

# Total penalized score of a CPDAG, through any member of its class.
.class_score <- function(S, n, amat, penalty) {
  res <- cpp_pdag_to_dag(amat)
  if (!isTRUE(res$ok)) return(-Inf)
  dag <- matrix(as.integer(res$dag), nrow(amat), ncol(amat))
  sum(vapply(seq_len(nrow(dag)), function(j)
    cpp_ges_local_score(S, n, j, which(dag[, j] == 1L & dag[j, ] == 0L),
                        penalty), numeric(1)))
}

#' Greedy equivalence search on the residual covariance
#'
#' Searches over Markov equivalence classes of directed acyclic graphs (DAGs)
#' by the two-phase greedy equivalence search (GES) of Chickering (2002): a
#' forward phase that starts from the empty class (or from \code{start}) and
#' applies, at each step, the insertion that raises the score most, and a
#' backward phase that starts from what the forward phase left and applies the
#' best deletion. Each phase stops when no operator improves the score. The
#' score of a class, computed through any of its members, is the penalized
#' Gaussian working criterion
#' \deqn{\sum_{j=1}^{R}\left\{-\frac{n}{2}\log\hat\delta_j-
#'   \frac{\lambda}{2}\,(|P_j|+1)\right\},}
#' where \eqn{\hat\delta_j} is the residual variance of the regression of
#' response \eqn{j} on its parents \eqn{P_j} in the residual covariance (the
#' closed-form innovation variance of the directed layer), \eqn{|P_j|+1} counts
#' the coefficients and the innovation variance of the node, and \eqn{\lambda}
#' is the \code{penalty}. Larger is better. Up to an additive constant the
#' score is \eqn{-1/2} times the criterion \eqn{-2\ell^p_e+\lambda q}, with
#' \eqn{\ell^p_e} the maximized working criterion and \eqn{q=R+|E|} the number
#' of parameters of the member, so that maximizing it with
#' \eqn{\lambda=\log n} minimizes the Bayesian information criterion (BIC) of
#' the working criterion. A parent set whose block of the
#' residual covariance is numerically singular, or that leaves a nonpositive
#' residual variance, scores \code{-Inf} and is never chosen.
#'
#' @section Why a search over classes is well posed here:
#' The working criterion of the residual layer depends on the data through the
#' residual covariance alone, so two members of an equivalence class fit the
#' same covariance with the same number of parameters and score alike. That is
#' score equivalence, and without it a search over classes would not be defined.
#' It also means the output is a class and not a graph: what the second-order
#' separations identify is what the score can distinguish. The score is also
#' decomposable, a sum of terms that each depend on one node and its parents,
#' which is what lets an operator be scored locally. The large-sample
#' optimality of the search shown by Chickering (2002) requires in addition a
#' locally consistent score, a property that is not claimed for the working
#' criterion on Pearson residuals.
#'
#' @section How it differs from the constraint-based route:
#' \code{\link{mc_graph_pipeline}} tests separations one at a time at a fixed
#' level and then orients; this one trades every edge against a penalty on the
#' number of parameters and never tests a hypothesis. The two answer the same
#' question by different means and need not agree, which is the point of having
#' both. Neither is run through an external library.
#'
#' @param fit an object of class \code{"mcglm_graph"}, whose residual covariance
#'   the search runs on.
#' @param penalty single finite non-negative number, the multiplier
#'   \eqn{\lambda} on the parameter count. \code{NULL} (default) uses
#'   \code{log(n)}, with \code{n} the number of units of the fit, which is
#'   the Bayesian information criterion; \code{2} gives Akaike's.
#' @param start optional square 0/1 adjacency matrix on the responses, in the
#'   package convention, to start the forward phase from; it should be a CPDAG.
#'   The default is the empty class.
#' @param max_set largest number of candidate neighbors enumerated by an
#'   operator (a node pair whose candidate set is larger is skipped), a guard
#'   against the exponential growth in dense neighborhoods. A single non-negative
#'   integer.
#' @param max_steps cap on the number of steps of each phase, a single positive
#'   integer.
#' @return An object of class \code{"mcg_ges"}: a list with the \code{cpdag}
#'   adjacency matrix (labeled by response); its \code{score}; the
#'   \code{penalty} and \code{n} used; \code{steps}, a named integer vector
#'   with the number of \code{forward} and \code{backward} steps; and
#'   \code{trace}, a data frame with one row per operator applied (\code{phase},
#'   \code{step}, the nodes \code{x} and \code{y}, the set \code{set} of
#'   neighbors involved, comma separated, and the score \code{gain}).
#'   \code{response_names} is also recorded. The \code{print} method returns
#'   \code{x} invisibly.
#' @references Chickering, D. M. (2002). Optimal structure identification with
#'   greedy search. \emph{Journal of Machine Learning Research} 3, 507-554.
#' @family structure discovery
#' @examples
#' \donttest{
#' set.seed(1); n <- 600; x <- rnorm(n)
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.8 * z1 + 0.8 * z2 + rnorm(n)
#' d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' mc_ges(f)
#' }
#' @seealso \code{\link{mc_graph_pipeline}} for the constraint-based route,
#'   \code{\link{mc_class_members}} to enumerate the DAGs of the class found.
#' @export
mc_ges <- function(fit, penalty = NULL, start = NULL, max_set = 12L,
                   max_steps = 200L) {
  stopifnot(inherits(fit, "mcglm_graph"))
  g <- fit$graph
  S <- as.matrix(g$Shat)
  n <- g$n
  R <- g$R
  if (is.null(penalty)) penalty <- log(n)
  if (length(penalty) != 1L || !is.finite(penalty) || penalty < 0)
    stop("penalty must be one finite non-negative number.", call. = FALSE)
  if (length(max_set) != 1L || !is.finite(max_set) || max_set < 0)
    stop("max_set must be one non-negative integer.", call. = FALSE)
  if (length(max_steps) != 1L || !is.finite(max_steps) || max_steps < 1)
    stop("max_steps must be one positive integer.", call. = FALSE)
  max_set <- as.integer(max_set); max_steps <- as.integer(max_steps)

  amat <- if (is.null(start)) matrix(0L, R, R) else {
    a <- as.matrix(start)
    if (!nrow(a) || nrow(a) != R || ncol(a) != R || any(!a %in% c(0, 1)) ||
        any(diag(a) != 0))
      stop("start must be a square 0/1 adjacency on the responses.",
           call. = FALSE)
    storage.mode(a) <- "integer"
    a
  }

  trace <- list()
  steps <- c(forward = 0L, backward = 0L)

  for (phase in c("forward", "backward")) {
    repeat {
      if (steps[[phase]] >= max_steps) break
      best <- if (phase == "forward")
        cpp_ges_best_insert(S, n, amat, penalty, max_set) else
          cpp_ges_best_delete(S, n, amat, penalty, max_set)
      if (!isTRUE(best$found)) break
      x <- as.integer(best$x); y <- as.integer(best$y)
      set <- as.integer(best$set)
      cand <- if (phase == "forward") .apply_insert(amat, x, y, set) else
        .apply_delete(amat, x, y, set)
      nxt <- .pdag_to_cpdag(cand)
      if (is.null(nxt)) break              # the operator left no extension
      amat <- nxt
      steps[[phase]] <- steps[[phase]] + 1L
      trace[[length(trace) + 1L]] <- data.frame(
        phase = phase, step = steps[[phase]], x = x, y = y,
        set = paste(set, collapse = ","), gain = as.numeric(best$gain),
        stringsAsFactors = FALSE)
    }
  }

  dimnames(amat) <- list(g$response_names, g$response_names)
  structure(list(
    cpdag = amat,
    score = .class_score(S, n, amat, penalty),
    penalty = penalty, steps = steps,
    trace = if (length(trace)) do.call(rbind, trace) else
      data.frame(phase = character(), step = integer(), x = integer(),
                 y = integer(), set = character(), gain = numeric()),
    response_names = g$response_names, n = n), class = "mcg_ges")
}

#' @rdname mc_ges
#' @param x an object of class \code{"mcg_ges"}.
#' @param ... ignored.
#' @export
print.mcg_ges <- function(x, ...) {
  cat(sprintf("Greedy equivalence search (penalty = %.3g, n = %d)\n",
              x$penalty, x$n))
  cat(sprintf("Steps: %d forward, %d backward | score: %.4f\n",
              x$steps[["forward"]], x$steps[["backward"]], x$score))
  ed <- sum(x$cpdag == 1L | t(x$cpdag) == 1L) / 2
  cat(sprintf("Skeleton edges: %d\n", ed))
  cat("Adjacency:\n")
  print(x$cpdag)
  invisible(x)
}
