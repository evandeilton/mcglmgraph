# cpdag.R: constraint-based structure learning: PC-style skeleton with
# separating sets, v-structure and Meek orientation, CPDAG pipeline.
# License: GPL-3.
#
# Why constraint-based and not the concentration (GGM) graph: the GGM precision
# selection conditions on *all* other responses at once, so it returns the moral
# graph, in which the parents of a collider are made adjacent. On a moral
# skeleton the unshielded triple that identifies a collider no longer exists, so
# colliders become unrecoverable. The PC skeleton instead removes an edge as
# soon as *some* subset separates the pair, and records that separating set,
# which is exactly the information the v-structure rule needs. Everything below
# uses only the estimated between-response covariance and the sample size, so
# the discovery is self-contained inside mcglmgraph; bnlearn is offered solely
# as an optional score-based alternative for the orientation step.
#
# Adjacency convention (amat):
#   amat[i, j] == 1 && amat[j, i] == 1  -> undirected edge i - j
#   amat[i, j] == 1 && amat[j, i] == 0  -> directed edge   i -> j
#   amat[i, j] == 0 && amat[j, i] == 0  -> no edge

# ---- second-order conditional-association test ------------------------------
# Partial correlation of i and j given the set S, read off the inverse of the
# (i, j, S) sub-block of the covariance, with Fisher's z transform. Returns the
# partial correlation and the two-sided p-value of H0: rho_{ij . S} = 0. Fisher's
# z is a working approximation: for dependent responses fitted by estimating
# functions its calibration is only approximate.
.ci_pcor <- function(Sigma, i, j, S, n, df_residual = NULL) {
  # With raw observations the Fisher factor is sqrt(n - |S| - 3). If the caller
  # supplies nu_e, the residual degrees of freedom after the first-stage fit,
  # it is sqrt(nu_e - |S| - 2).
  df <- if (is.null(df_residual)) n - length(S) - 3 else
    df_residual - length(S) - 2
  if (!is.finite(df) || df <= 0)
    stop("Insufficient Fisher degrees of freedom: need a positive value, got ",
         format(df), " for conditioning order ", length(S), ".")
  # One kernel serves the skeleton and the enumeration of separating sets.
  r <- cpp_ci_pcor(as.matrix(Sigma), as.integer(i), as.integer(j),
                   as.integer(S), as.numeric(n),
                   if (is.null(df_residual)) 0 else as.numeric(df_residual))
  if (!isTRUE(r$usable))
    stop("Conditional-association block does not define a valid partial ",
         "correlation.")
  list(pcor = r$pcor, p_value = r$p_value, df = r$df,
       inversion_status = if (isTRUE(r$regularized)) "ridge_regularized" else
         "exact",
       inversion_ridge = r$ridge)
}

# ---- shared validation of the separation inputs -----------------------------
# mc_pc_skeleton() and mc_separating_sets() read the same covariance under the
# same contract, so they must refuse the same input with the same words. A
# caller should not learn a different contract depending on which door it came
# through.
.check_separation_inputs <- function(Sigma, n, alpha, df_residual) {
  Sigma <- as.matrix(Sigma)
  if (!is.numeric(Sigma) || !nrow(Sigma) || nrow(Sigma) != ncol(Sigma) ||
      any(!is.finite(Sigma)))
    stop("Sigma must be a finite numeric square matrix.")
  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-10)))
    stop("Sigma must be symmetric.")
  if (nrow(Sigma) < 2L)
    stop("Sigma must describe at least two responses.")
  if (length(n) != 1L || !is.finite(n) || n <= 0)
    stop("n must be one finite positive number.")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must lie strictly between zero and one.")
  if (!is.null(df_residual) &&
      (length(df_residual) != 1L || !is.finite(df_residual) ||
       df_residual <= 2))
    stop("df_residual must be one finite residual degrees-of-freedom value greater than 2.")
  Sigma
}

# ---- residual degrees of freedom --------------------------------------------

#' Residual degrees of freedom of the first stage
#'
#' Returns \eqn{\nu_e}, the residual degrees of freedom that the Fisher factor
#' of the conditional-association test uses in the form \eqn{\nu_e-|S|-2},
#' when the first stage fixes it exactly, and says why it does not otherwise.
#'
#' @section When the correction is exact:
#' The Pearson residuals of a Gaussian response with identity link and constant
#' variance are the least-squares residuals of its design, scaled. When every
#' response has that form and the same column space of design, of rank
#' \eqn{h}, the residual covariance is the sample covariance after projecting
#' out a common \eqn{h}-dimensional space. For Gaussian responses its partial
#' correlations are then distributed as those of a sample covariance with
#' \eqn{n-h} degrees of freedom, that of \eqn{n-h+1} observations with an
#' estimated mean, so \eqn{\nu_e=n-h}. With the intercept alone, \eqn{h=1} and
#' \eqn{\nu_e-|S|-2=n-|S|-3}, the factor used when no correction is supplied.
#' Outside that case no finite-sample correction is established and none is
#' offered: the function returns \code{exact = FALSE} with the reasons, and the
#' approximation by \eqn{n} (\code{df_residual = NULL} in the separation
#' routines) is the one to use. The check requires constant variance, the
#' identity link, an identity within-response structure, unit first-stage
#' weights, and designs of equal rank spanning the same column space for all
#' responses.
#'
#' @param fit an object of class \code{"mcglm_graph"} or \code{"mcglm"}.
#' @return An object of class \code{"mcg_residual_df"}: a list with \code{df},
#'   the value of \eqn{\nu_e} or \code{NA}; \code{exact}, a logical; \code{n},
#'   the number of observations;
#'   \code{rank}, the common rank \eqn{h} of the designs or \code{NA}; and
#'   \code{why}, a character vector with the reasons the correction is not
#'   exact, empty when it is. The \code{print} method returns \code{x}
#'   invisibly.
#' @seealso \code{\link{mc_pc_skeleton}}, whose \code{df_residual} accepts it,
#'   and \code{\link{mc_graph_pipeline}}, where \code{df_residual = "auto"}
#'   calls it.
#' @family structure discovery
#' @examples
#' \donttest{
#' set.seed(1); n <- 60; x <- rnorm(n); g <- gl(3, 20)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x, g = g)
#' f <- mcglm_graph(list(y1 ~ x + g, y2 ~ x + g, y3 ~ x + g), data = d)
#' mc_residual_df(f)          # 60 - 4
#' }
#' @export
mc_residual_df <- function(fit) {
  if (!inherits(fit, "mcglm"))
    stop("fit must be an object of class \"mcglm_graph\" or \"mcglm\".",
         call. = FALSE)
  n <- as.integer(fit$n_obs)
  why <- character(0)
  vr <- tryCatch(unlist(fit$variance), error = function(e) NULL)
  lk <- tryCatch(unlist(fit$link), error = function(e) NULL)
  X <- fit$list_X
  if (is.null(vr) || is.null(lk) || is.null(X) || !length(n) || is.na(n))
    return(structure(list(df = NA_real_, exact = FALSE, n = n,
                          rank = NA_integer_,
                          why = "the first stage of the fit could not be read"),
                     class = "mcg_residual_df"))
  if (any(vr != "constant"))
    why <- c(why, "not every response has constant variance")
  if (any(lk != "identity"))
    why <- c(why, "not every response has the identity link")
  if (!.identity_matrix_pred(fit))
    why <- c(why, "the within-response structure is not the identity")
  w <- tryCatch(.first_stage_weights(fit, n * length(X)),
                error = function(e) NULL)
  if (is.null(w) || any(abs(w - 1) > 1e-12))
    why <- c(why, "the first stage is weighted")
  ranks <- vapply(X, function(M) qr(as.matrix(M))$rank, integer(1))
  h <- qr(do.call(cbind, lapply(X, as.matrix)))$rank
  if (any(ranks != h))
    why <- c(why, "the designs of the responses do not span the same column space")
  exact <- !length(why)
  structure(list(df = if (exact) as.numeric(n - h) else NA_real_,
                 exact = exact, n = n,
                 rank = if (exact) as.integer(h) else NA_integer_, why = why),
            class = "mcg_residual_df")
}

#' @rdname mc_residual_df
#' @param x an object of class \code{"mcg_residual_df"}.
#' @param ... ignored.
#' @export
print.mcg_residual_df <- function(x, ...) {
  if (isTRUE(x$exact)) {
    cat(sprintf("Residual degrees of freedom: %g (n = %d, rank of the design = %d)\n",
                x$df, x$n, x$rank))
  } else {
    cat("Residual degrees of freedom: not exact for this fit\n")
    for (w in x$why) cat("  - ", w, "\n", sep = "")
    cat("The approximation by n (df_residual = NULL) is the one to use.\n")
  }
  invisible(x)
}

# df_residual as the separation routines take it: NULL for the approximation by
# n, a number, an "mcg_residual_df" object, or "auto" when a fit is at hand. A
# correction that is not exact is refused, never approximated.
.resolve_df_residual <- function(df_residual, fit = NULL) {
  if (is.null(df_residual) || is.numeric(df_residual)) return(df_residual)
  if (identical(df_residual, "auto")) {
    if (is.null(fit))
      stop("df_residual = \"auto\" needs the fitted model, which this function ",
           "does not receive; pass mc_residual_df(fit) instead.", call. = FALSE)
    df_residual <- mc_residual_df(fit)
  }
  if (inherits(df_residual, "mcg_residual_df")) {
    if (!isTRUE(df_residual$exact))
      stop("the residual degrees of freedom are not exact for this fit (",
           paste(df_residual$why, collapse = "; "), "). Leave df_residual = ",
           "NULL for the approximation by n, or supply a value.", call. = FALSE)
    return(df_residual$df)
  }
  stop("df_residual must be NULL, a number, \"auto\" or the result of ",
       "mc_residual_df().", call. = FALSE)
}

# ---- PC-stable skeleton -----------------------------------------------------

#' PC-style skeleton with separating sets
#'
#' Constraint-based skeleton discovery from a between-response covariance
#' matrix. Starting from the complete graph, an edge \eqn{i - j} is deleted as
#' soon as a conditioning set \eqn{S} of adjacent nodes makes the partial
#' correlation \eqn{\rho_{ij \cdot S}} non-significant, and that set is stored as
#' the separating set of the pair. The neighborhoods are frozen within each
#' conditioning order and both endpoint neighborhoods are examined, as in the
#' order-independent skeleton stage of the PC-stable algorithm (Colombo and
#' Maathuis, 2014), a variant of the PC algorithm (Spirtes, Glymour and
#' Scheines, 2000; Kalisch and Buehlmann, 2007).
#'
#' @section The test:
#' The partial correlation \eqn{\hat\rho_{ij\cdot S}} is read off the inverse
#' of the \eqn{(i, j, S)} block of \code{Sigma} and clamped to
#' \eqn{[-0.999999, 0.999999]} before Fisher's transform
#' \eqn{z=\operatorname{atanh}(\hat\rho_{ij\cdot S})\sqrt{n-|S|-3}} (with
#' \eqn{\nu_e-|S|-2} in place of \eqn{n-|S|-3} when \code{df_residual} is
#' given), so that a perfectly collinear block gives a large but finite
#' statistic. A pair is declared conditionally uncorrelated, and its edge
#' removed, when the two-sided normal p-value of \eqn{z} is at least
#' \code{alpha}. A block that is nearly singular is inverted after a ridge of
#' \eqn{10^{-8}} times its scale and the event is counted in
#' \code{diagnostics}; a conditioning order whose Fisher degrees of freedom are
#' not positive, or a block that gives no valid partial correlation even after
#' the ridge, is an error. Because \code{Sigma} is usually computed from the
#' estimated residuals of responses that need not be Gaussian, the normal
#' reference for \eqn{z} is a working asymptotic calibration and not an exact
#' one. The factor takes
#' \code{n} as the number of independent units, the rows of the data for a fit
#' with an identity within-response structure. There is no \code{cluster}
#' argument: with dependent rows, \code{n} overstates the information and the
#' test is anti-conservative.
#'
#' @param Sigma an \eqn{R \times R} covariance matrix of the responses, e.g.
#'   \code{fit$graph$Shat}.
#' @param n number of experimental units, taken as independent.
#' @param alpha significance level of the Fisher z test of a zero partial
#'   correlation, a single number strictly between zero and one.
#' @param max_ord maximum size of the conditioning set; \code{NULL} uses
#'   \eqn{R - 2}.
#' @param df_residual optional residual degrees of freedom \eqn{\nu_e} from the
#'   first-stage fit, as a number or as the result of
#'   \code{\link{mc_residual_df}}. When supplied, the Fisher factor uses
#'   \eqn{\nu_e - |S| - 2}; otherwise it uses \eqn{n - |S| - 3}.
#' @return A list with the symmetric integer \code{adjacency} matrix of the
#'   skeleton (one for an edge, zero otherwise); the \code{sepset} list of
#'   length \eqn{R^2}, whose elements \code{(i - 1) * R + j} and
#'   \code{(j - 1) * R + i} hold the integer vector of nodes that separated
#'   \eqn{i} from \eqn{j}, \code{integer(0)} for the empty set, and
#'   \code{NULL} when the pair was never separated; the \code{df_residual}
#'   used; and \code{diagnostics}, a list with the number of tests
#'   (\code{n_tests}), the number of ridge-regularized blocks
#'   (\code{n_regularized}), the largest ridge (\code{max_ridge}), the smallest
#'   Fisher degrees of freedom (\code{min_df}) and the contract used
#'   (\code{df_contract}).
#' @references Spirtes, P., Glymour, C. and Scheines, R. (2000).
#'   \emph{Causation, Prediction, and Search}, 2nd edition. MIT Press.
#'
#'   Kalisch, M. and Buehlmann, P. (2007). Estimating high-dimensional directed
#'   acyclic graphs with the PC-algorithm. \emph{Journal of Machine Learning
#'   Research} 8, 613-636.
#'
#'   Colombo, D. and Maathuis, M. H. (2014). Order-independent constraint-based
#'   causal structure learning. \emph{Journal of Machine Learning Research} 15,
#'   3921-3962.
#' @seealso \code{\link{mc_separating_sets}}, \code{\link{mc_find_vstructures}},
#'   \code{\link{mc_graph_pipeline}}
#' @family structure discovery
#' @examples
#' set.seed(1); n <- 600
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
#' S <- cov(cbind(z1, z2, z3))
#' mc_pc_skeleton(S, n)$adjacency
#' @export
mc_pc_skeleton <- function(Sigma, n, alpha = 0.05, max_ord = NULL,
                           df_residual = NULL) {
  df_residual <- .resolve_df_residual(df_residual)
  Sigma <- .check_separation_inputs(Sigma, n, alpha, df_residual)
  R <- nrow(Sigma)
  adj <- matrix(1L, R, R); diag(adj) <- 0L
  sepset <- vector("list", R * R)
  if (is.null(max_ord)) max_ord <- R - 2L
  if (length(max_ord) != 1L || !is.finite(max_ord) || max_ord < 0)
    stop("max_ord must be one non-negative integer.")
  max_ord <- min(as.integer(max_ord), R - 2L)
  n_tests <- 0L
  n_regularized <- 0L
  max_ridge <- 0
  min_df <- Inf
  ord <- 0L
  repeat {
    adj0 <- adj                       # frozen neighborhoods (PC-stable)
    has_candidate <- FALSE
    for (i in seq_len(R)) for (j in seq_len(R)) {
      if (i == j || adj[i, j] == 0L) next
      nb <- setdiff(which(adj0[i, ] == 1L), j)
      if (length(nb) < ord) next
      has_candidate <- TRUE
      # Enumerate conditioning sets by position and map back into nb. Indexing
      # positions avoids utils::combn() treating a length-one nb as a scalar n
      # (which would silently condition on seq_len(nb) instead of nb itself).
      combos <- if (ord == 0L) list(integer(0)) else
        lapply(utils::combn(length(nb), ord, simplify = FALSE),
               function(pos) nb[pos])
      for (S in combos) {
        tst <- .ci_pcor(Sigma, i, j, S, n, df_residual = df_residual)
        n_tests <- n_tests + 1L
        n_regularized <- n_regularized +
          as.integer(tst$inversion_status != "exact")
        max_ridge <- max(max_ridge, tst$inversion_ridge)
        min_df <- min(min_df, tst$df)
        if (tst$p_value >= alpha) {                # i _||_ j | S
          adj[i, j] <- adj[j, i] <- 0L
          sepset[[(i - 1L) * R + j]] <- S
          sepset[[(j - 1L) * R + i]] <- S
          break
        }
      }
    }
    ord <- ord + 1L
    if (!has_candidate || ord > max_ord) break
  }
  list(
    adjacency = adj, sepset = sepset, df_residual = df_residual,
    diagnostics = list(
      n_tests = n_tests, n_regularized = n_regularized,
      max_ridge = max_ridge,
      min_df = if (is.finite(min_df)) min_df else NA_real_,
      df_contract = if (is.null(df_residual)) "n - |S| - 3" else
        "df_residual - |S| - 2"
    )
  )
}

# ---- exhaustive separating sets ---------------------------------------------

#' All separating sets of a pair of responses
#'
#' Enumerates every conditioning set that makes a pair of responses
#' conditionally uncorrelated at the second order, given a skeleton. The
#' skeleton search of \code{\link{mc_pc_skeleton}} stops at the first such set
#' and records only that one, which is all the standard v-structure rule needs.
#' The order-independent orientation rules need the whole collection, and they
#' need it on the final skeleton rather than on the neighborhoods that happened
#' to be live when the edge was deleted, so it is rebuilt here from scratch: for
#' the pair \eqn{(i, j)} every subset of the neighbors of \eqn{i} other than
#' \eqn{j}, and every subset of the neighbors of \eqn{j} other than \eqn{i}, of
#' size at most \code{max_ord}, is tested, and the separating ones are returned.
#' A set reachable from both neighborhoods is reported once.
#'
#' @section Why the collection is informative on its own:
#' Two sets can separate the same pair with very different strength, and the
#' skeleton keeps whichever the enumeration reached first. Reporting the
#' collection shows how much of the orientation rests on that arbitrary choice.
#'
#' @param Sigma an \eqn{R \times R} covariance matrix of the responses.
#' @param n number of experimental units, taken as independent; the test is
#'   the one of \code{\link{mc_pc_skeleton}}, with the same clamp of the
#'   partial correlation and no \code{cluster} argument.
#' @param adjacency the symmetric 0/1 skeleton, for instance the
#'   \code{adjacency} component of \code{\link{mc_pc_skeleton}}.
#' @param pairs optional integer matrix with two columns giving the pairs to
#'   examine; \code{NULL} takes every non-adjacent pair. The pairs must be
#'   non-adjacent in \code{adjacency}, since a separating set is not defined for
#'   an edge of the skeleton.
#' @param alpha significance level of the conditional-association test, on the
#'   same contract as \code{\link{mc_pc_skeleton}}.
#' @param max_ord largest conditioning set enumerated; \code{NULL} uses
#'   \eqn{R - 2}.
#' @param df_residual optional residual degrees of freedom \eqn{\nu_e} of the
#'   first-stage fit, a number or the result of \code{\link{mc_residual_df}},
#'   passed to the same Fisher contract used by the skeleton.
#' @return An object of class \code{"mcg_sepsets"}: a list with \code{pairs}
#'   (one row per pair examined, with the number of separating sets found and a
#'   \code{truncated} flag), \code{sets} (one row per separating set, with its
#'   label, size, partial correlation and p-value), the 0/1 matrix
#'   \code{members} of set membership, and \code{diagnostics}. The
#'   \code{truncated} flag of a pair is \code{TRUE} when its collection may be
#'   incomplete: a neighborhood of the pair is larger than \code{max_ord}, or
#'   some conditioning set gave no valid test, because its Fisher degrees of
#'   freedom are not positive or its block gives no valid partial correlation.
#'   Such a set is skipped and counted in \code{n_unusable}; unlike in
#'   \code{\link{mc_pc_skeleton}}, it does not stop the enumeration.
#'   The object also records \code{alpha}, \code{max_ord}, \code{df_residual}
#'   and \code{response_names}, and \code{diagnostics} with the number of tests
#'   (\code{n_tests}), of ridge-regularized blocks (\code{n_regularized}), of
#'   unusable blocks (\code{n_unusable}), the largest ridge (\code{max_ridge})
#'   and the Fisher contract (\code{df_contract}). The \code{print} method
#'   returns \code{x} invisibly.
#' @seealso \code{\link{mc_graph_pipeline}}, whose \code{orient_rule} argument
#'   consumes this object.
#' @references Ramsey, J., Zhang, J. and Spirtes, P. (2006).
#'   Adjacency-faithfulness and conservative causal inference. \emph{Proceedings
#'   of the 22nd Conference on Uncertainty in Artificial Intelligence}, 401-408.
#'
#'   Colombo, D. and Maathuis, M. H. (2014). Order-independent constraint-based
#'   causal structure learning. \emph{Journal of Machine Learning Research} 15,
#'   3921-3962.
#' @family structure discovery
#' @examples
#' set.seed(1); n <- 600
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
#' S <- cov(cbind(z1, z2, z3))
#' sk <- mc_pc_skeleton(S, n)
#' mc_separating_sets(S, n, sk$adjacency)$sets
#' @export
mc_separating_sets <- function(Sigma, n, adjacency, pairs = NULL,
                               alpha = 0.05, max_ord = NULL,
                               df_residual = NULL) {
  df_residual <- .resolve_df_residual(df_residual)
  Sigma <- .check_separation_inputs(Sigma, n, alpha, df_residual)
  R <- nrow(Sigma)
  adjacency <- as.matrix(adjacency)
  if (!nrow(adjacency) || nrow(adjacency) != R || ncol(adjacency) != R)
    stop("adjacency must be square with one row per response.")
  if (!is.numeric(adjacency) || any(!adjacency %in% c(0, 1)) ||
      any(diag(adjacency) != 0))
    stop("adjacency must be a 0/1 matrix with a zero diagonal.")
  if (!identical(unname(adjacency), unname(t(adjacency))))
    stop("adjacency must be symmetric; separating sets are defined on the skeleton.")
  storage.mode(adjacency) <- "integer"

  if (is.null(max_ord)) max_ord <- R - 2L
  if (length(max_ord) != 1L || !is.finite(max_ord) || max_ord < 0)
    stop("max_ord must be one non-negative integer.")
  max_ord <- min(as.integer(max_ord), R - 2L)

  if (is.null(pairs)) {
    idx <- which(upper.tri(adjacency) & adjacency == 0L, arr.ind = TRUE)
    pairs <- matrix(as.integer(idx), ncol = 2L)
  } else {
    pairs <- .validate_edges(pairs, R, arg = "pairs")
    if (nrow(pairs) && any(adjacency[pairs] == 1L))
      stop("pairs must hold non-adjacent responses; a separating set is not ",
           "defined for an edge of the skeleton.", call. = FALSE)
  }
  storage.mode(pairs) <- "integer"

  nm <- colnames(Sigma)
  if (is.null(nm) || anyNA(nm)) nm <- paste0("V", seq_len(R))

  res <- cpp_sepsets_pairs(
    Sigma, as.numeric(n), adjacency, pairs, as.numeric(alpha),
    as.integer(max_ord),
    if (is.null(df_residual)) 0 else as.numeric(df_residual)
  )

  set_pair <- as.integer(res$set_pair)
  n_sets <- length(set_pair)
  members <- matrix(as.integer(res$set_members), nrow = n_sets, ncol = R)
  colnames(members) <- nm
  label <- if (n_sets) vapply(seq_len(n_sets), function(s) {
    ix <- which(members[s, ] == 1L)
    if (!length(ix)) "{}" else paste0("{", paste(nm[ix], collapse = ","), "}")
  }, character(1)) else character(0)

  sets <- data.frame(
    pair = set_pair,
    i = pairs[set_pair, 1L], j = pairs[set_pair, 2L],
    response1 = nm[pairs[set_pair, 1L]], response2 = nm[pairs[set_pair, 2L]],
    set = label,
    size = as.integer(res$set_size),
    pcor = as.numeric(res$set_pcor),
    p_value = as.numeric(res$set_p_value),
    row.names = NULL, stringsAsFactors = FALSE
  )
  pairs_df <- data.frame(
    i = pairs[, 1L], j = pairs[, 2L],
    response1 = nm[pairs[, 1L]], response2 = nm[pairs[, 2L]],
    n_sets = tabulate(set_pair, nbins = nrow(pairs)),
    truncated = as.logical(as.integer(res$truncated)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  structure(list(
    pairs = pairs_df, sets = sets, members = members,
    alpha = alpha, max_ord = max_ord, df_residual = df_residual,
    response_names = nm,
    diagnostics = list(
      n_tests = as.integer(res$n_tests),
      n_regularized = as.integer(res$n_regularized),
      n_unusable = as.integer(res$n_unusable),
      max_ridge = as.numeric(res$max_ridge),
      df_contract = if (is.null(df_residual)) "n - |S| - 3" else
        "df_residual - |S| - 2"
    )
  ), class = "mcg_sepsets")
}

#' @rdname mc_separating_sets
#' @param x an object of class \code{"mcg_sepsets"}.
#' @param ... ignored.
#' @export
print.mcg_sepsets <- function(x, ...) {
  cat(sprintf("Separating sets (alpha = %.3g, max order = %d, df = %s)\n",
              x$alpha, x$max_ord, x$diagnostics$df_contract))
  cat(sprintf("Pairs examined: %d | sets found: %d | tests: %d\n",
              nrow(x$pairs), nrow(x$sets), x$diagnostics$n_tests))
  if (any(x$pairs$truncated))
    cat(sprintf("Pairs with a truncated enumeration: %d\n",
                sum(x$pairs$truncated)))
  if (x$diagnostics$n_regularized > 0L)
    cat(sprintf("Regularized blocks: %d (max ridge %.3g)\n",
                x$diagnostics$n_regularized, x$diagnostics$max_ridge))
  if (nrow(x$sets)) {
    cat("\n")
    print(x$sets[, c("response1", "response2", "set", "pcor", "p_value")])
  }
  invisible(x)
}

# ---- orientation of the unshielded triples ----------------------------------
# For every unshielded triple i - k - j of the skeleton, k is a collider
# i -> k <- j exactly when it does not belong to a set that separates i and j:
# when k is a collider the pair i, j is dependent only after conditioning on k,
# so k cannot belong to any such set. Which sets are consulted is what the three
# rules below disagree about.
# Decide each unshielded triple under an order-independent rule.
#
# For the triple (i, k, j) with center k, let SS be the collection of sets that
# separate i and j on the final skeleton. The conservative rule of Ramsey, Zhang
# and Spirtes (2006) declares a collider when k belongs to none of them, a
# non-collider when it belongs to all of them, and leaves the triple ambiguous
# otherwise. The majority rule of Colombo and Maathuis (2014) decides by the
# proportion of sets containing k, with a tie left ambiguous. An empty
# collection is ambiguous: there is no evidence either way, which happens when
# the enumeration was truncated.
.decide_triples <- function(triples, sep, rule) {
  nt <- nrow(triples)
  out <- data.frame(is_collider = logical(nt), ambiguous = logical(nt),
                    n_sets = integer(nt), prop_center = rep(NA_real_, nt))
  if (nt == 0L) return(out)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = ":")
  pid <- match(key(triples[, 1], triples[, 3]),
               key(sep$pairs$i, sep$pairs$j))
  for (t in seq_len(nt)) {
    k <- triples[t, 2]
    rows <- if (is.na(pid[t])) integer(0) else which(sep$sets$pair == pid[t])
    out$n_sets[t] <- length(rows)
    if (!length(rows)) { out$ambiguous[t] <- TRUE; next }
    inc <- sep$members[rows, k] == 1L
    prop <- mean(inc)
    out$prop_center[t] <- prop
    if (identical(rule, "conservative")) {
      if (!any(inc)) out$is_collider[t] <- TRUE
      else if (!all(inc)) out$ambiguous[t] <- TRUE
    } else {                                     # majority
      if (prop < 0.5) out$is_collider[t] <- TRUE
      else if (prop == 0.5) out$ambiguous[t] <- TRUE
    }
  }
  out
}

# ---- v-structure orientation ------------------------------------------------

# Orient the colliders of the skeleton under one of three rules.
#
# "standard" reads the single separating set the skeleton recorded, which is the
# rule of Spirtes, Glymour and Scheines and the default. "conservative" and
# "majority" consult the whole collection returned
# by mc_separating_sets(), and can leave a triple ambiguous, in which case no
# orientation is written for it and the triple is passed on so that Meek's rules
# do not settle by propagation what the separations left open.
#
# Two colliders can require opposite directions on the same edge, for instance
# on the skeleton 1-2-3-4 when the routine records both 1->2<-3 and 2->3<-4.
# Writing each collider as it is found lets the later write silently overturn
# the earlier one, and the returned object then contradicts a v-structure it
# reports. Orientations are therefore collected first, conflicting edges are
# left undirected, and the conflicts are returned so that the caller can refuse
# the CPDAG label.
.orient_vstructures <- function(skel, rule = c("standard", "conservative",
                                               "majority"), sep = NULL) {
  rule <- match.arg(rule)
  adj <- skel$adjacency; R <- nrow(adj); sepset <- skel$sepset
  storage.mode(adj) <- "integer"
  triples <- cpp_unshielded_triples(adj)          # columns (i, k, j), k center
  amat <- adj
  empty_vs <- data.frame(i = integer(0), k = integer(0), j = integer(0),
                         is_collider = logical(0), applied = logical(0),
                         ambiguous = logical(0), n_sets = integer(0),
                         prop_center = numeric(0))
  conflicts <- data.frame(from = integer(0), to = integer(0))
  if (nrow(triples) == 0)
    return(list(amat = amat, vstructures = empty_vs, conflicts = conflicts,
                ambiguous = triples[0, , drop = FALSE]))

  if (rule == "standard") {
    coll <- vapply(seq_len(nrow(triples)), function(t) {
      i <- triples[t, 1]; k <- triples[t, 2]; j <- triples[t, 3]
      !(k %in% sepset[[(i - 1L) * R + j]])
    }, logical(1))
    amb <- rep(FALSE, nrow(triples))
    n_sets <- rep(NA_integer_, nrow(triples))
    prop <- rep(NA_real_, nrow(triples))
  } else {
    if (is.null(sep))
      stop("The ", rule, " rule needs the separating sets of the skeleton.",
           call. = FALSE)
    dec <- .decide_triples(triples, sep, rule)
    coll <- dec$is_collider; amb <- dec$ambiguous
    n_sets <- dec$n_sets; prop <- dec$prop_center
  }

  req <- matrix(FALSE, R, R)                      # req[a, b]: some collider needs a -> b
  for (t in which(coll)) {
    i <- triples[t, 1]; k <- triples[t, 2]; j <- triples[t, 3]
    req[i, k] <- TRUE; req[j, k] <- TRUE
  }
  clash <- req & t(req)                           # both directions requested
  if (any(clash)) {
    idx <- which(clash & upper.tri(clash), arr.ind = TRUE)
    conflicts <- data.frame(from = as.integer(idx[, 1]),
                            to = as.integer(idx[, 2]))
  }
  # A collider is applied only when neither of its two arcs is contested, so a
  # conflicting triple contributes no orientation at all rather than half of
  # one: writing a single leg would assert a collider the data did not settle.
  applied <- vapply(seq_len(nrow(triples)), function(t) {
    if (!coll[t]) return(FALSE)
    i <- triples[t, 1]; k <- triples[t, 2]; j <- triples[t, 3]
    !clash[i, k] && !clash[j, k]
  }, logical(1))
  for (t in which(applied)) {
    i <- triples[t, 1]; k <- triples[t, 2]; j <- triples[t, 3]
    amat[i, k] <- 1L; amat[k, i] <- 0L
    amat[j, k] <- 1L; amat[k, j] <- 0L
  }
  vs <- data.frame(i = triples[, 1], k = triples[, 2], j = triples[, 3],
                   is_collider = coll, applied = applied, ambiguous = amb,
                   n_sets = n_sets, prop_center = prop)
  list(amat = amat, vstructures = vs, conflicts = conflicts,
       ambiguous = triples[amb, , drop = FALSE])
}

# ---- validity gates for the partially directed output -----------------------

# The two gates below decide combinatorial properties of the adjacency and are
# evaluated once per resample by the structure-frequency routines, so they live
# in C++. The R functions stay as the names the rest of the package calls.
.has_directed_cycle <- function(amat) {
  storage.mode(amat) <- "integer"
  cpp_has_directed_cycle(amat)
}

# Existence of a consistent DAG extension (Dor and Tarsi, 1992): repeatedly
# remove a node with no outgoing directed edge whose undirected neighbors are
# adjacent to all its other neighbors; a graph that empties out admits an
# extension that preserves the directed edges and the v-structures.
.has_consistent_extension <- function(amat) {
  storage.mode(amat) <- "integer"
  cpp_has_consistent_extension(amat)
}

# Every collider that was actually written must survive Meek propagation.
.vstructures_preserved <- function(amat, vs) {
  if (is.null(vs) || nrow(vs) == 0) return(TRUE)
  keep <- vs$is_collider & vs$applied
  if (!any(keep)) return(TRUE)
  sub <- vs[keep, , drop = FALSE]
  all(vapply(seq_len(nrow(sub)), function(t) {
    i <- sub$i[t]; k <- sub$k[t]; j <- sub$j[t]
    amat[i, k] == 1L && amat[k, i] == 0L &&
      amat[j, k] == 1L && amat[k, j] == 0L
  }, logical(1)))
}

# Separating sets needed by an order-independent rule: only the endpoint pairs
# of the unshielded triples of the skeleton are enumerated, since those are the
# only pairs whose collection the rules consult.
.sepsets_for_triples <- function(Sigma, n, skel, alpha, max_ord, df_residual) {
  adj <- skel$adjacency
  storage.mode(adj) <- "integer"
  triples <- cpp_unshielded_triples(adj)
  prs <- if (nrow(triples))
    unique(cbind(pmin(triples[, 1], triples[, 3]),
                 pmax(triples[, 1], triples[, 3]))) else
      matrix(integer(0), nrow = 0L, ncol = 2L)
  mc_separating_sets(Sigma, n, adj, pairs = prs, alpha = alpha,
                     max_ord = max_ord, df_residual = df_residual)
}

#' Detect v-structures in the constraint-based skeleton
#'
#' Runs the PC skeleton search on the estimated between-response covariance and,
#' for each unshielded triple \eqn{(i, k, j)} (\eqn{i - k - j} with \eqn{i} and
#' \eqn{j} nonadjacent), decides whether \eqn{k} is a collider
#' \eqn{i \to k \leftarrow j}, under one of three rules. A collider is
#' identified by the absence of \eqn{k} from the sets that separate \eqn{i}
#' and \eqn{j}, which is why the search conditions on subsets and not on all
#' other responses at once.
#'
#' @section The three rules:
#' \code{"standard"} consults the single set the skeleton recorded for the pair
#' and declares a collider when the center is absent from it. It is the rule of
#' Spirtes, Glymour and Scheines and it depends on which separating set the
#' enumeration happened to reach first. \code{"conservative"} enumerates every
#' set that separates the endpoints on the final skeleton and declares a
#' collider only when the center is in none of them, a non-collider when it is
#' in all of them, and leaves the triple ambiguous otherwise
#' (Ramsey, Zhang and Spirtes, 2006). \code{"majority"} decides by the
#' proportion of separating sets containing the center, with a tie left
#' ambiguous (Colombo and Maathuis, 2014). Combined with the order-independent
#' skeleton, the last two make the procedure independent of the order of the
#' responses.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param alpha significance level of the second-order separation test.
#' @param max_ord maximum conditioning order passed to
#'   \code{\link{mc_pc_skeleton}} and to the enumeration of separating sets.
#' @param df_residual optional residual degrees of freedom passed to
#'   \code{\link{mc_pc_skeleton}} and to the reported separating-set tests:
#'   \code{NULL} for the approximation by \eqn{n}, a number, or \code{"auto"},
#'   which takes \eqn{\nu_e} from \code{\link{mc_residual_df}} and refuses a fit
#'   for which it is not exact.
#' @param orient_rule \code{"standard"} (default), \code{"conservative"} or
#'   \code{"majority"}.
#' @return A data frame with columns \code{i}, \code{k}, \code{j},
#'   \code{is_collider}, \code{applied}, \code{ambiguous}, \code{n_sets},
#'   \code{prop_center}, \code{pcor} and \code{p_value}. The last two refer to
#'   the test of \eqn{i} and \eqn{j} given the set the skeleton recorded;
#'   \code{applied} is \code{FALSE} for a collider whose orientation was
#'   withheld because another collider required the opposite direction on one of
#'   its edges. Under \code{"standard"} the columns \code{n_sets} and
#'   \code{prop_center} are \code{NA}, since no collection is formed.
#' @references Spirtes, P., Glymour, C. and Scheines, R. (2000).
#'   \emph{Causation, Prediction, and Search}, 2nd edition. MIT Press.
#'
#'   Ramsey, J., Zhang, J. and Spirtes, P. (2006). Adjacency-faithfulness and
#'   conservative causal inference. \emph{Proceedings of the 22nd Conference on
#'   Uncertainty in Artificial Intelligence}, 401-408.
#'
#'   Colombo, D. and Maathuis, M. H. (2014). Order-independent constraint-based
#'   causal structure learning. \emph{Journal of Machine Learning Research} 15,
#'   3921-3962.
#' @seealso \code{\link{mc_pc_skeleton}}, \code{\link{mc_separating_sets}},
#'   \code{\link{mc_graph_pipeline}}
#' @family structure discovery
#' @examples
#' \donttest{
#' set.seed(1); n <- 600
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7*z1 + 0.7*z2 + rnorm(n)  # collider
#' d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_find_vstructures(f)
#' }
#' @export
mc_find_vstructures <- function(fit, alpha = 0.05, max_ord = NULL,
                                df_residual = NULL,
                                orient_rule = c("standard", "conservative",
                                                "majority")) {
  stopifnot(inherits(fit, "mcglm_graph"))
  orient_rule <- match.arg(orient_rule)
  df_residual <- .resolve_df_residual(df_residual, fit)
  Sigma <- fit$graph$Shat; n <- fit$graph$n; R <- nrow(Sigma)
  skel <- mc_pc_skeleton(Sigma, n, alpha = alpha, max_ord = max_ord,
                         df_residual = df_residual)
  sep <- if (orient_rule == "standard") NULL else
    .sepsets_for_triples(Sigma, n, skel, alpha, max_ord, df_residual)
  vs <- .orient_vstructures(skel, rule = orient_rule, sep = sep)$vstructures
  if (nrow(vs) == 0)
    return(data.frame(i = integer(0), k = integer(0), j = integer(0),
                      is_collider = logical(0), applied = logical(0),
                      ambiguous = logical(0), n_sets = integer(0),
                      prop_center = numeric(0),
                      pcor = numeric(0), p_value = numeric(0)))
  rep_pcor <- numeric(nrow(vs)); rep_p <- numeric(nrow(vs))
  for (t in seq_len(nrow(vs))) {
    S <- skel$sepset[[(vs$i[t] - 1L) * R + vs$j[t]]]
    tst <- .ci_pcor(Sigma, vs$i[t], vs$j[t], S, n,
                    df_residual = df_residual)
    rep_pcor[t] <- tst$pcor; rep_p[t] <- tst$p_value
  }
  data.frame(vs, pcor = rep_pcor, p_value = rep_p)
}

# ---- Meek orientation rules -------------------------------------------------

#' Apply Meek's orientation rules
#'
#' Orients every edge whose direction is implied by the edges already oriented,
#' by repeated application of Meek's four rules until no further edge can be
#' oriented.
#'
#' @section The four rules:
#' \deqn{\begin{array}{ll}
#'   \mathrm{R1} & c \to a,\ a - b,\ c\ \mathrm{and}\ b\ \mathrm{nonadjacent}
#'     \Rightarrow a \to b\\
#'   \mathrm{R2} & a \to c \to b,\ a - b \Rightarrow a \to b\\
#'   \mathrm{R3} & a - b,\ a - c,\ a - d,\ c \to b,\ d \to b,\ c\
#'     \mathrm{and}\ d\ \mathrm{nonadjacent} \Rightarrow a \to b\\
#'   \mathrm{R4} & a - b,\ a - c,\ a - d,\ c \to d,\ d \to b,\ c\
#'     \mathrm{and}\ b\ \mathrm{nonadjacent} \Rightarrow a \to b
#' \end{array}}
#' R1 to R3 are complete when the only oriented edges are those of the
#' v-structures, so the completed partially directed acyclic graph (CPDAG) of a
#' Markov equivalence class needs no more than those three (Meek, 1995). R4 is
#' idle in that case and becomes necessary as soon as
#' edges are oriented from outside, through \code{required} or \code{forbidden},
#' because the graph is then no longer the CPDAG of any class and the first
#' three rules stop being complete. For an input that admits a consistent
#' extension, the four together return the maximally
#' oriented partially directed graph: an edge is left undirected exactly when
#' the background knowledge and the v-structures admit an acyclic extension with
#' each of its two directions.
#'
#' @section Background knowledge:
#' \code{required} lists arcs that must be present with the stated direction and
#' \code{forbidden} arcs that must not. An undirected edge whose reverse is
#' forbidden is oriented, since only one direction is left. Both are applied
#' before propagation and are checked for internal consistency: an arc required
#' in both directions, an arc both required and forbidden, an arc required
#' between responses the skeleton does not join, and an arc required against an
#' edge already oriented the other way are all refused. An incompatibility that
#' only surfaces during propagation, such as a cycle, is reported in the
#' returned attributes rather than raised, because the partially oriented graph
#' is still informative about where the incompatibility lies.
#'
#' @section Ambiguous triples:
#' Under the conservative and majority rules an unshielded triple can be left
#' undecided, and propagation must not settle by implication what the
#' separations left open. R1, R3 and R4 all orient an edge in order to avoid
#' creating a collider at its tail, so each is withheld when the triple carrying
#' that collider is ambiguous. R2 only forbids a cycle and involves no
#' unshielded triple, so it is applied unconditionally. Leaving
#' \code{ambiguous} at \code{NULL} reproduces the unrestricted rules exactly.
#'
#' @param amat a square 0/1 adjacency matrix with zero diagonal, in the package
#'   convention: \code{amat[i, j] = amat[j, i] = 1} is an undirected edge
#'   \eqn{i - j}, \code{amat[i, j] = 1} with \code{amat[j, i] = 0} is an arc
#'   \eqn{i \to j}, and both zero is no edge.
#' @param ambiguous optional integer matrix with three columns \eqn{(i, k, j)},
#'   one row per unshielded triple whose center could not be decided. Both
#'   orders of the endpoints are treated as the same triple.
#' @param required optional integer matrix with two columns \eqn{(from, to)},
#'   arcs the orientation must contain.
#' @param forbidden optional integer matrix with two columns \eqn{(from, to)},
#'   arcs the orientation must not contain.
#' @return The oriented adjacency matrix, with attributes
#'   \code{"meek_conflicts"} (edges two rules wanted in opposite directions,
#'   left undirected), \code{"meek_iterations"}, \code{"meek_blocked"} (the
#'   number of orientations withheld by an ambiguous triple) and
#'   \code{"knowledge"}, a list recording the arcs imposed and whether the
#'   result stayed acyclic and kept every required arc.
#' @references Meek, C. (1995). Causal inference and causal explanation with
#'   background knowledge. \emph{Proceedings of the 11th Conference on
#'   Uncertainty in Artificial Intelligence}, 403-410.
#'
#'   Verma, T. and Pearl, J. (1990). Equivalence and synthesis of causal models.
#'   \emph{Proceedings of the 6th Conference on Uncertainty in Artificial
#'   Intelligence}, 220-227.
#' @family structure discovery
#' @examples
#' # 1 -> 3 <- 2 with 3 - 4: R1 orients 3 -> 4
#' a <- matrix(0L, 4, 4)
#' a[1, 3] <- a[2, 3] <- 1L
#' a[3, 4] <- a[4, 3] <- 1L
#' mc_orient_meek(a)
#'
#' # background knowledge on a fully undirected chain 1 - 2 - 3
#' b <- matrix(0L, 3, 3)
#' b[1, 2] <- b[2, 1] <- b[2, 3] <- b[3, 2] <- 1L
#' mc_orient_meek(b, required = rbind(c(1L, 2L)))   # 1 -> 2 forces 2 -> 3
#' @export
mc_orient_meek <- function(amat, ambiguous = NULL, required = NULL,
                           forbidden = NULL) {
  amat <- as.matrix(amat)
  if (!nrow(amat) || nrow(amat) != ncol(amat) || any(!amat %in% 0:1) ||
      any(diag(amat) != 0))
    stop("amat must be a square 0/1 adjacency matrix with a zero diagonal.")
  storage.mode(amat) <- "integer"
  R <- nrow(amat)

  amb_keys <- integer(0)
  if (!is.null(ambiguous)) {
    ambiguous <- as.matrix(ambiguous)
    if (ncol(ambiguous) != 3L)
      stop("ambiguous must have exactly three columns (i, k, j).")
    if (nrow(ambiguous)) {
      if (!is.numeric(ambiguous) || any(!is.finite(ambiguous)) ||
          any(ambiguous != round(ambiguous)) ||
          any(ambiguous < 1L) || any(ambiguous > R))
        stop("ambiguous must contain node indices in 1:", R, ".")
      storage.mode(ambiguous) <- "integer"
      # keys are 0-based inside the kernel, with both endpoint orders present
      amb_keys <- unique(c(
        (ambiguous[, 1] - 1L) * R * R + (ambiguous[, 2] - 1L) * R +
          (ambiguous[, 3] - 1L),
        (ambiguous[, 3] - 1L) * R * R + (ambiguous[, 2] - 1L) * R +
          (ambiguous[, 1] - 1L)
      ))
    }
  }

  kn <- .apply_knowledge(amat, required, forbidden)
  amat <- kn$amat

  res <- cpp_meek_rules(amat, as.integer(amb_keys))
  out <- matrix(as.integer(res$amat), R, R, dimnames = dimnames(amat))
  conf <- matrix(as.logical(res$conflicts), R, R)

  ci <- which(conf & upper.tri(conf), arr.ind = TRUE)
  attr(out, "meek_conflicts") <- if (nrow(ci))
    data.frame(from = as.integer(ci[, 1L]), to = as.integer(ci[, 2L])) else
      data.frame(from = integer(0), to = integer(0))
  attr(out, "meek_iterations") <- as.integer(res$iterations)
  attr(out, "meek_blocked") <- as.integer(res$blocked)
  attr(out, "knowledge") <- list(
    required = kn$required, forbidden = kn$forbidden,
    n_imposed = kn$n_imposed,
    acyclic = !.has_directed_cycle(out),
    required_kept = .arcs_present(out, kn$required)
  )
  out
}

# Impose the background knowledge on the adjacency, refusing the requests that
# contradict one another or the edges already oriented. Returns the adjacency
# with the arcs written and the validated tables.
.apply_knowledge <- function(amat, required, forbidden) {
  R <- nrow(amat)
  required <- .validate_arcs_or_null(required, R, "required")
  forbidden <- .validate_arcs_or_null(forbidden, R, "forbidden")
  n_imposed <- 0L

  if (nrow(required) && nrow(forbidden)) {
    both <- intersect(paste(required[, 1L], required[, 2L], sep = ":"),
                      paste(forbidden[, 1L], forbidden[, 2L], sep = ":"))
    if (length(both))
      stop("the same arc is both required and forbidden: ",
           paste(sub(":", " -> ", both), collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(required)) {
    rev <- paste(required[, 2L], required[, 1L], sep = ":")
    if (any(rev %in% paste(required[, 1L], required[, 2L], sep = ":")))
      stop("an arc is required in both directions.", call. = FALSE)
    for (k in seq_len(nrow(required))) {
      a <- required[k, 1L]; b <- required[k, 2L]
      if (amat[a, b] == 0L && amat[b, a] == 0L)
        stop("a required arc joins responses that the skeleton does not: ",
             a, " -> ", b, ".", call. = FALSE)
      if (amat[a, b] == 0L && amat[b, a] == 1L)
        stop("a required arc contradicts an edge already oriented the other ",
             "way: ", a, " -> ", b, ".", call. = FALSE)
      if (amat[a, b] == 1L && amat[b, a] == 1L) {
        amat[b, a] <- 0L
        n_imposed <- n_imposed + 1L
      }
    }
  }
  if (nrow(forbidden)) {
    for (k in seq_len(nrow(forbidden))) {
      a <- forbidden[k, 1L]; b <- forbidden[k, 2L]
      if (amat[a, b] == 1L && amat[b, a] == 0L)
        stop("a forbidden arc is already oriented that way: ", a, " -> ", b,
             ".", call. = FALSE)
      # only one direction is left for the edge, so it is no longer a choice
      if (amat[a, b] == 1L && amat[b, a] == 1L) {
        amat[a, b] <- 0L
        n_imposed <- n_imposed + 1L
      }
    }
  }
  list(amat = amat, required = required, forbidden = forbidden,
       n_imposed = n_imposed)
}

.validate_arcs_or_null <- function(arcs, R, arg) {
  if (is.null(arcs)) return(matrix(integer(0), nrow = 0L, ncol = 2L))
  arcs <- as.matrix(arcs)
  if (length(dim(arcs)) != 2L || ncol(arcs) != 2L)
    stop(arg, " must have exactly two columns (from, to).", call. = FALSE)
  if (nrow(arcs) == 0L) return(matrix(integer(0), nrow = 0L, ncol = 2L))
  if (!is.numeric(arcs) || any(!is.finite(arcs)) || any(arcs != round(arcs)) ||
      any(arcs < 1L) || any(arcs > R))
    stop(arg, " must contain node indices in 1:", R, ".", call. = FALSE)
  if (any(arcs[, 1L] == arcs[, 2L]))
    stop(arg, " must not contain self-loops.", call. = FALSE)
  storage.mode(arcs) <- "integer"
  unique(arcs)
}

# Are all these arcs present and oriented that way in the adjacency?
.arcs_present <- function(amat, arcs) {
  if (!nrow(arcs)) return(TRUE)
  all(vapply(seq_len(nrow(arcs)), function(k)
    amat[arcs[k, 1L], arcs[k, 2L]] == 1L &&
      amat[arcs[k, 2L], arcs[k, 1L]] == 0L, logical(1)))
}

#' Convert a skeleton to the set of allowed bnlearn arcs
#'
#' Lists both directions of every edge of a skeleton as named arcs, in the
#' \code{from}/\code{to} format of \pkg{bnlearn}. The complement of this set
#' within all ordered pairs of nodes is a blacklist that confines a score-based
#' search to the skeleton.
#'
#' @param edges two-column integer matrix of undirected skeleton edges, as
#'   indices into \code{node_names}; \code{NULL} or a matrix with no rows gives
#'   an empty result.
#' @param node_names character vector of node labels.
#' @return A data frame with character columns \code{from} and \code{to},
#'   containing both directions of every skeleton edge. Its complement can be
#'   supplied to \pkg{bnlearn} as a blacklist.
#' @family structure discovery
#' @examples
#' mc_skeleton_to_bnlearn(rbind(c(1, 2), c(2, 3)), c("y1", "y2", "y3"))
#' @export
mc_skeleton_to_bnlearn <- function(edges, node_names) {
  if (is.null(edges) || nrow(edges) == 0)
    return(data.frame(from = character(0), to = character(0)))
  arcs <- rbind(
    data.frame(from = node_names[edges[, 1]], to = node_names[edges[, 2]]),
    data.frame(from = node_names[edges[, 2]], to = node_names[edges[, 1]]))
  arcs
}

# ---- full pipeline ----------------------------------------------------------

# Background-knowledge arcs given either as response labels or as indices.
# Labels are what a user of the fitted object has at hand; indices are what the
# adjacency is written in.
.knowledge_to_index <- function(arcs, response_names, arg) {
  if (is.null(arcs)) return(NULL)
  if (is.data.frame(arcs)) arcs <- as.matrix(arcs)
  arcs <- as.matrix(arcs)
  if (length(dim(arcs)) != 2L || ncol(arcs) != 2L)
    stop(arg, " must have exactly two columns (from, to).", call. = FALSE)
  if (nrow(arcs) == 0L) return(NULL)
  if (is.character(arcs)) {
    idx <- match(arcs, response_names)
    if (anyNA(idx))
      stop(arg, " names a response the fit does not have: ",
           paste(unique(arcs[is.na(idx)]), collapse = ", "), ".", call. = FALSE)
    arcs <- matrix(idx, ncol = 2L)
  }
  storage.mode(arcs) <- "integer"
  arcs
}

.bnlearn_available <- function() requireNamespace("bnlearn", quietly = TRUE)

#' Constraint-based CPDAG discovery
#'
#' Two-phase structure learning, self-contained in \pkg{mcglmgraph}, that
#' estimates the completed partially directed acyclic graph (CPDAG) of the
#' Markov equivalence class of the directed acyclic graph (DAG) underlying the
#' second-order dependence between the responses.
#' Phase one recovers the skeleton and the separating sets from the
#' between-response covariance via \code{\link{mc_pc_skeleton}}. Phase two
#' orients the edges up to the Markov equivalence class: the native route
#' (\code{orient_method = "pc"}) orients v-structures by the separating-set rule
#' and then applies Meek's rules; the optional route (\code{"hc"} or
#' \code{"tabu"}) delegates the second phase to \pkg{bnlearn}'s score-based
#' search.  A blacklist prevents arcs outside the native skeleton, but the
#' score search may remove native skeleton edges; this route is therefore not
#' a pure orientation of a fixed skeleton.
#'
#' @section Which separating sets the orientation consults:
#' \code{orient_rule} selects between the rule that reads the single set the
#' skeleton recorded (\code{"standard"}) and the two rules that enumerate every
#' set separating the endpoints of each unshielded triple on the final skeleton,
#' \code{"conservative"} and \code{"majority"}. The enumeration is native and
#' uses the same conditional-association test and the same degrees-of-freedom
#' contract as the skeleton, so a difference between two runs is attributable to
#' the rule and not to a change of test. A triple that the rule leaves ambiguous
#' receives no orientation, and Meek's R1, R3 and R4 are withheld where they
#' would settle it by propagation; the count is reported in
#' \code{n_ambiguous}.
#'
#' @section What background knowledge changes:
#' Without it the output is the CPDAG of an equivalence class, the object the
#' second-order separations identify. With it the output is no longer the CPDAG
#' of any class: when no triple is left ambiguous it is the maximally oriented
#' partially directed graph of Meek (1995), the union of the members of the
#' class that satisfy the knowledge. The \code{label} records which of the two
#' it is, and \code{knowledge} carries the arcs that were imposed. Knowledge
#' is refused, rather than silently dropped, when it contradicts itself, when
#' it requires an arc between responses the skeleton does not join, and when
#' it contradicts an edge oriented by a v-structure. A contradiction with an
#' edge that only Meek's rules would orient is not refused beforehand: it
#' surfaces during propagation as a conflict or a cycle and the output then
#' fails the checks below.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param alpha significance level of the second-order separation test.
#' @param orient_method \code{"pc"} (native, default), \code{"hc"} or
#'   \code{"tabu"} (\pkg{bnlearn}). When \pkg{bnlearn} is not installed, the
#'   last two fall back on \code{"pc"} with a warning, and the result records
#'   it in \code{orientation_fallback}.
#' @param orient_rule rule for the unshielded triples of the native route:
#'   \code{"standard"} (default), \code{"conservative"} or \code{"majority"}.
#'   It has no effect when \code{orient_method} is a score-based route.
#' @param required,forbidden optional background knowledge, as two-column
#'   matrices of arcs \eqn{(from, to)} given either by response name or by
#'   index. An arc in \code{required} must be present with that direction, an
#'   arc in \code{forbidden} must not. They are imposed before propagation, and
#'   Meek's fourth rule, idle without them, carries their consequences through
#'   the rest of the graph. Available on the native route only.
#' @param max_ord maximum conditioning order; \code{NULL} uses \eqn{R - 2}.
#' @param df_residual optional residual degrees of freedom \eqn{\nu_e}:
#'   \code{NULL} (default) uses the raw-sample contract based on \eqn{n}; a
#'   number or the result of \code{\link{mc_residual_df}} is used as given;
#'   \code{"auto"} computes it with \code{\link{mc_residual_df}} and refuses a
#'   fit for which it is not exact.
#' @return An object of class \code{"mcg_cpdag"}, a list with the 0/1 matrix
#'   \code{skeleton} of the oriented structure and its \code{edges} (two-column
#'   integer matrix); the constraint-stage \code{pc_skeleton} and \code{sepset},
#'   as returned by \code{\link{mc_pc_skeleton}}; the \code{cpdag} adjacency in
#'   the convention of \code{\link{mc_orient_meek}}, labeled by response; the
#'   \code{vstructures} table (\code{NULL} on a score-based route); the
#'   \code{conflicts} table of edges required in opposite directions; the
#'   named logical vector \code{gates}; the logical \code{valid_cpdag}; the
#'   \code{label} to be used when reporting the output; \code{knowledge} and
#'   \code{n_imposed}, the background knowledge with the checks on it and the
#'   number of arcs it oriented; \code{orient_method}, \code{orient_rule},
#'   \code{orient_method_requested} and \code{orientation_fallback};
#'   \code{sepsets}, the collections enumerated under a rule other than
#'   \code{"standard"} (\code{NULL} otherwise); \code{ambiguous}, the triples
#'   left undecided, and their number \code{n_ambiguous};
#'   \code{meek_iterations} and \code{meek_blocked}, the number of sweeps of
#'   Meek's rules and of orientations withheld by an ambiguous triple;
#'   \code{separation_diagnostics}, the diagnostics of the skeleton search; and
#'   \code{alpha}, \code{df_residual} and \code{response_names}. The
#'   \code{print} method returns \code{x} invisibly.
#' @section What the output may be called:
#' The output is submitted to seven checks, reported in \code{gates}: no pair
#' of colliders or of rules required opposite directions on an edge
#' (\code{no_orientation_conflict}); no separation test needed a ridge to
#' invert its block (\code{unregularized_separation}); the orientation stayed
#' within the learned skeleton (\code{orientation_within_skeleton}) and kept
#' all of it (\code{skeleton_preserved}); the directed part has no cycle
#' (\code{acyclic}); every collider that was written survives Meek propagation
#' (\code{vstructures_preserved}); and a consistent DAG extension exists, by
#' the algorithm of Dor and Tarsi (1992) (\code{consistent_extension}). When
#' all pass, \code{valid_cpdag} is \code{TRUE} and \code{label} is
#' \code{"CPDAG"}, or \code{"maximally oriented graph under background
#' knowledge"} when the knowledge oriented some edge. When any fails,
#' \code{valid_cpdag} is \code{FALSE} and \code{label} becomes
#' \code{"second-order partially directed structure"}: the object still records
#' what the separations supported, but it is not certified as the CPDAG of a
#' Markov equivalence class. The first two checks record facts of the run and
#' not properties of the returned graph, so a refused label does not establish
#' that the output differs from the CPDAG of a class. A triple left ambiguous
#' by the conservative or majority rule stays unoriented, so a certified output
#' can be less oriented than the CPDAG of the class of its extension.
#'
#' @section Assumptions behind the reading:
#' Passing the checks is a statement about internal consistency only. Reading
#' the output as the CPDAG of the second-order structure assumes that the
#' between-response covariance is generated by some DAG and that every zero
#' partial correlation corresponds to a d-separation in it (second-order Markov
#' and faithfulness). Reading it as the Markov equivalence class of the
#' distribution of the responses assumes, in addition, that a zero partial
#' correlation of the residuals is equivalent to conditional independence of
#' the responses given the conditioning responses and the covariates, which
#' holds for jointly Gaussian residuals whose conditional covariance does not
#' depend on the covariates and fails in general otherwise. Recovery of the
#' class is asymptotic, for a fixed number of responses and a level
#' \code{alpha} that decreases with the sample size. A causal reading of the
#' arcs requires further causal sufficiency, absence of selection bias and
#' substantive justification. The procedure checks none of these assumptions.
#' @references Spirtes, P., Glymour, C. and Scheines, R. (2000).
#'   \emph{Causation, Prediction, and Search}, 2nd edition. MIT Press.
#'
#'   Meek, C. (1995). Causal inference and causal explanation with background
#'   knowledge. \emph{Proceedings of the 11th Conference on Uncertainty in
#'   Artificial Intelligence}, 403-410.
#'
#'   Dor, D. and Tarsi, M. (1992). A simple algorithm to construct a consistent
#'   extension of a partially oriented graph. Technical Report R-185, Cognitive
#'   Systems Laboratory, University of California, Los Angeles.
#'
#'   Andersson, S. A., Madigan, D. and Perlman, M. D. (1997). A characterization
#'   of Markov equivalence classes for acyclic digraphs. \emph{The Annals of
#'   Statistics} 25(2), 505-541. \doi{10.1214/aos/1031833662}
#' @seealso \code{\link{mc_pc_skeleton}}, \code{\link{mc_separating_sets}},
#'   \code{\link{mc_orient_meek}}, \code{\link{mc_class_members}},
#'   \code{\link{mc_ges}} for the score-based alternative.
#' @family structure discovery
#' @examples
#' \donttest{
#' set.seed(1); n <- 600
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7*z1 + 0.7*z2 + rnorm(n)
#' d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = rnorm(n))
#' f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#' mc_graph_pipeline(f)$cpdag
#' }
#' @export
mc_graph_pipeline <- function(fit, alpha = 0.05,
                              orient_method = c("pc", "hc", "tabu"),
                              orient_rule = c("standard", "conservative",
                                              "majority"),
                              required = NULL, forbidden = NULL,
                              max_ord = NULL, df_residual = NULL) {
  stopifnot(inherits(fit, "mcglm_graph"))
  orient_method <- match.arg(orient_method)
  orient_rule <- match.arg(orient_rule)
  orient_method_requested <- orient_method
  df_residual <- .resolve_df_residual(df_residual, fit)
  Sigma <- fit$graph$Shat; n <- fit$graph$n
  rn <- fit$graph$response_names
  required <- .knowledge_to_index(required, rn, "required")
  forbidden <- .knowledge_to_index(forbidden, rn, "forbidden")
  if ((!is.null(required) || !is.null(forbidden)) &&
      orient_method_requested != "pc")
    stop("background knowledge applies to the native orientation only; use ",
         "orient_method = \"pc\".", call. = FALSE)
  skel <- mc_pc_skeleton(Sigma, n, alpha = alpha, max_ord = max_ord,
                         df_residual = df_residual)

  if (orient_method %in% c("hc", "tabu") && !.bnlearn_available()) {
    warning("Package 'bnlearn' is not installed; using the native PC orientation.")
    orient_method <- "pc"
  }

  conflicts <- data.frame(from = integer(0), to = integer(0))
  meek_iterations <- 0L
  meek_blocked <- 0L
  sep <- NULL
  knowledge <- NULL
  ambiguous <- matrix(integer(0), nrow = 0L, ncol = 3L)
  if (orient_method == "pc") {
    if (orient_rule != "standard")
      sep <- .sepsets_for_triples(Sigma, n, skel, alpha, max_ord, df_residual)
    res <- .orient_vstructures(skel, rule = orient_rule, sep = sep)
    ambiguous <- res$ambiguous
    amat <- mc_orient_meek(res$amat,
                           ambiguous = if (nrow(ambiguous)) ambiguous else NULL,
                           required = required, forbidden = forbidden)
    meek_conflicts <- attr(amat, "meek_conflicts", exact = TRUE)
    meek_iterations <- attr(amat, "meek_iterations", exact = TRUE)
    meek_blocked <- attr(amat, "meek_blocked", exact = TRUE)
    knowledge <- attr(amat, "knowledge", exact = TRUE)
    vs <- res$vstructures
    conflicts <- unique(rbind(res$conflicts, meek_conflicts))
  } else {
    amat <- .orient_bnlearn(skel$adjacency, fit, orient_method)
    vs <- NULL
  }
  output_skeleton <- (amat == 1L) | (t(amat) == 1L)
  native_skeleton <- skel$adjacency == 1L
  gates <- c(no_orientation_conflict = nrow(conflicts) == 0L,
             unregularized_separation = skel$diagnostics$n_regularized == 0L,
             orientation_within_skeleton = !any(output_skeleton &
                                                  !native_skeleton),
             skeleton_preserved = identical(unname(output_skeleton),
                                             unname(native_skeleton)),
             acyclic = !.has_directed_cycle(amat),
             vstructures_preserved = .vstructures_preserved(amat, vs),
             consistent_extension = .has_consistent_extension(amat))
  valid <- all(gates)
  n_imposed <- if (is.null(knowledge)) 0L else knowledge$n_imposed
  # The diagnostics of the orientation step travel in fields of their own, so
  # the adjacency is handed over as a plain labeled matrix and printing it does
  # not spill them.
  amat <- matrix(as.integer(amat), nrow(amat), ncol(amat))
  dimnames(amat) <- list(rn, rn)
  skeleton <- 1L * output_skeleton
  dimnames(skeleton) <- list(rn, rn)
  pc_skeleton <- skel$adjacency
  dimnames(pc_skeleton) <- list(rn, rn)
  edges <- which(skeleton == 1L & upper.tri(skeleton),
                 arr.ind = TRUE)
  edges <- matrix(as.integer(edges), ncol = 2)
  structure(list(skeleton = skeleton, pc_skeleton = pc_skeleton,
                 sepset = skel$sepset,
                 edges = edges, cpdag = amat, vstructures = vs,
                 conflicts = conflicts, gates = gates, valid_cpdag = valid,
                 label = if (!valid)
                   "second-order partially directed structure"
                   else if (n_imposed > 0L)
                     "maximally oriented graph under background knowledge" else
                       "CPDAG",
                 knowledge = knowledge,
                 n_imposed = n_imposed,
                 orient_method = orient_method,
                 orient_rule = if (orient_method == "pc") orient_rule else NA_character_,
                 orient_method_requested = orient_method_requested,
                 orientation_fallback = orient_method != orient_method_requested,
                 sepsets = sep,
                 ambiguous = ambiguous,
                 n_ambiguous = nrow(ambiguous),
                 meek_iterations = meek_iterations,
                 meek_blocked = meek_blocked,
                 separation_diagnostics = skel$diagnostics,
                 alpha = alpha,
                 df_residual = df_residual,
                 response_names = rn),
            class = "mcg_cpdag")
}

#' @rdname mc_graph_pipeline
#' @param x an object of class \code{"mcg_cpdag"}.
#' @param ... ignored.
#' @export
print.mcg_cpdag <- function(x, ...) {
  cat(sprintf("Structure discovery (orientation = %s%s, alpha = %.3g)\n",
              x$orient_method,
              if (!is.null(x$orient_rule) && !is.na(x$orient_rule) &&
                  x$orient_rule != "standard")
                paste0(", ", x$orient_rule, " rule") else "",
              x$alpha))
  if (!is.null(x$n_imposed) && x$n_imposed > 0L) {
    kn <- x$knowledge
    cat(sprintf("Background knowledge: %d arc%s imposed (%d required, %d forbidden)\n",
                x$n_imposed, if (x$n_imposed == 1L) "" else "s",
                nrow(kn$required), nrow(kn$forbidden)))
    if (!isTRUE(kn$required_kept))
      cat("A required arc did not survive propagation.\n")
  }
  if (!is.null(x$n_ambiguous) && x$n_ambiguous > 0L)
    cat(sprintf("Ambiguous triples: %d | propagations withheld: %d\n",
                x$n_ambiguous,
                if (is.null(x$meek_blocked)) 0L else x$meek_blocked))
  if (isTRUE(x$orientation_fallback))
    cat(sprintf("Requested orientation '%s' was unavailable; fallback used.\n",
                x$orient_method_requested))
  cat(sprintf("Output: %s\n", x$label))
  cat(sprintf("Skeleton edges: %d\n", nrow(x$edges)))
  if (!is.null(x$vstructures) && nrow(x$vstructures) > 0) {
    cat(sprintf("V-structures (colliders): %d",
                sum(x$vstructures$is_collider)))
    if (!is.null(x$vstructures$applied))
      cat(sprintf(", oriented: %d", sum(x$vstructures$applied)))
    cat("\n")
  }
  if (!is.null(x$gates) && !isTRUE(x$valid_cpdag)) {
    failed <- names(x$gates)[!x$gates]
    cat(sprintf("Failed validity gates: %s\n", paste(failed, collapse = ", ")))
    if (!is.null(x$conflicts) && nrow(x$conflicts) > 0)
      cat(sprintf("Edges left undirected by conflicting colliders: %s\n",
                  paste(sprintf("%d-%d", x$conflicts$from, x$conflicts$to),
                        collapse = ", ")))
  }
  if (!is.null(x$separation_diagnostics) &&
      x$separation_diagnostics$n_regularized > 0L)
    cat(sprintf("Regularized separation tests: %d (max ridge %.3g)\n",
                x$separation_diagnostics$n_regularized,
                x$separation_diagnostics$max_ridge))
  cat("Adjacency:\n"); print(x$cpdag)
  invisible(x)
}

# ---- internal: optional bnlearn orientation ---------------------------------
# bnlearn is restricted to arcs in the native PC skeleton by blacklisting every
# arc outside it.  Because no skeleton arc is whitelisted, the score search may
# remove edges and is not a pure orientation procedure.

.orient_bnlearn <- function(skel_adj, fit, method) {
  dat <- as.data.frame(fit$graph$residuals)
  names(dat) <- fit$graph$response_names
  edges <- which(skel_adj == 1L & upper.tri(skel_adj), arr.ind = TRUE)
  edges <- matrix(as.integer(edges), ncol = 2)
  allowed <- mc_skeleton_to_bnlearn(edges, fit$graph$response_names)
  bl <- .complement_arcs(fit$graph$response_names, allowed)
  net <- switch(method,
                hc   = bnlearn::hc(dat, blacklist = bl),
                tabu = bnlearn::tabu(dat, blacklist = bl))
  cp  <- bnlearn::cpdag(net)
  bnlearn::amat(cp)
}

.complement_arcs <- function(nodes, wl) {
  all_arcs <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
  all_arcs <- all_arcs[all_arcs$from != all_arcs$to, ]
  key <- function(d) paste(d$from, d$to, sep = "->")
  all_arcs[!key(all_arcs) %in% key(wl), ]
}
