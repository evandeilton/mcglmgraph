# crossval.R: out-of-sample comparison of candidate supports.
# Part of the mcglmgraph package. License: GPL-3.
#
# An information criterion penalizes the number of parameters on the sample that
# chose them. Splitting the units instead asks the question directly: fit the
# layer on one part and score it on the other. What is scored is the Gaussian
# working criterion of the residual layer, the same scalar the criteria are
# built from, evaluated per held-out unit so that folds of different sizes are
# comparable.

# Candidate supports, as a named list of edge matrices or NULL for complete.
.normalize_supports <- function(supports, fit) {
  g <- fit$graph
  if (is.null(supports)) {
    own <- g$spec$edges
    supports <- list(unstructured = NULL)
    if (!is.null(own) && nrow(as.matrix(own)) <
        choose(g$R, 2L)) supports$fitted <- own
    return(supports)
  }
  if (!is.list(supports) || is.null(names(supports)) || anyNA(names(supports)) ||
      any(names(supports) == ""))
    stop("supports must be a named list of edge matrices, with NULL for the ",
         "complete graph.", call. = FALSE)
  supports
}

#' Out-of-sample comparison of candidate supports
#'
#' Splits the units into folds, fits the graph layer on the units outside each
#' fold and scores it on the units inside it, which compares candidate supports
#' without a penalty term and without scoring a support on the sample that chose
#' it. The assignment of units to folds is random, so call \code{set.seed()}
#' beforehand for a reproducible split.
#'
#' @section What is scored:
#' The Gaussian working criterion of the residual layer, the scalar the
#' information criteria of \code{\link{mc_compute_criterion}} are built from,
#' evaluated on the held-out residuals under the covariance fitted on the rest
#' and divided by the number of held-out units,
#' \deqn{-\tfrac12\bigl[R\log(2\pi)+\log|\hat\Sigma_{\mathrm{train}}|
#'   +\mathrm{tr}(\hat\Sigma_{\mathrm{train}}^{-1}\hat S_{\mathrm{test}})\bigr],}
#' where \eqn{R} is the number of responses,
#' \eqn{\hat\Sigma_{\mathrm{train}}} is the between-response covariance fitted
#' to the sample covariance of the training residuals, and
#' \eqn{\hat S_{\mathrm{test}}} is the covariance of the held-out residuals.
#' Each covariance is centered at the mean of its own part and divided by the
#' number of units in it. Larger is better. The scale is
#' that of a log density per unit, so differences between supports are directly
#' comparable and a difference of a given size means the same at any fold size.
#'
#' @section What is held fixed:
#' The first stage is not refitted: the split is of the stored Pearson
#' residuals, so the comparison is between structures for the residual
#' covariance and not between whole pipelines. That is the question the sparse
#' against unstructured comparison asks, and it keeps every support scored on
#' exactly the same residuals.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param supports a named list of candidate supports, each a two-column matrix
#'   of edges (rows \eqn{(i,j)}, response indices) or \code{NULL} for the
#'   complete graph. The default, \code{NULL}, compares the complete graph
#'   (named \code{"unstructured"}) with the support the fit carries (named
#'   \code{"fitted"}), when that one is sparse.
#' @param folds integer number of folds, at least two and at most the number of
#'   units (default 5).
#' @param workers number of parallel workers (default 1, sequential). The folds
#'   are formed in the parent before any dispatch, so a given seed yields the
#'   same split at any worker count.
#' @return An object of class \code{"mcg_cv"}: a list with \code{per_fold}, a
#'   data frame with one row per support and fold (columns \code{support},
#'   \code{fold}, \code{score}, the held-out score, and \code{npar}, the
#'   number of free parameters of the support, which is the same in every
#'   fold); \code{summary}, a data frame with one row per support ordered by
#'   decreasing score (columns \code{support}, \code{n_folds}, the number of
#'   folds with a usable fit, \code{npar}, \code{score}, the mean over those
#'   folds, \code{se}, the standard deviation of the fold scores divided by
#'   the square root of their number, and \code{delta}, the difference from
#'   the best); \code{failures}, a data frame of the fits that failed
#'   (\code{support}, \code{fold}, \code{error}); and \code{folds}, \code{n}
#'   and \code{supports} (the names of the candidates). \code{print} returns
#'   its argument invisibly.
#' @examples
#' \donttest{
#' set.seed(1); n <- 400
#' K <- diag(4) * 2
#' K[1, 2] <- K[2, 1] <- -0.7; K[3, 4] <- K[4, 3] <- -0.7
#' Y <- matrix(rnorm(n * 4), n) %*% chol(solve(K))
#' d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], y4 = Y[, 4],
#'                 x = rnorm(n))
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
#' sel <- mc_graph_stepwise(f, criterion = "bic")
#' mc_graph_cv(f, supports = list(unstructured = NULL, selected = sel$edges))
#' }
#' @seealso \code{\link{mc_compute_criterion}} for the in-sample criteria,
#'   \code{\link{mc_graph_stepwise}} and \code{\link{mc_select_lambda}} for
#'   supports to compare.
#' @family information criteria
#' @export
mc_graph_cv <- function(fit, supports = NULL, folds = 5L, workers = 1L) {
  stopifnot(inherits(fit, "mcglm_graph"))
  if (length(folds) != 1L || !is.finite(folds) || folds < 2 ||
      folds != as.integer(folds))
    stop("folds must be an integer of at least two.", call. = FALSE)
  folds <- as.integer(folds)
  g <- fit$graph
  n <- g$n
  if (folds > n)
    stop("folds cannot exceed the number of units.", call. = FALSE)
  supports <- .normalize_supports(supports, fit)
  ctrl <- .graph_control(fit)

  E <- as.matrix(g$residuals)
  assign_fold <- sample(rep(seq_len(folds), length.out = n))

  grid <- expand.grid(support = seq_along(supports), fold = seq_len(folds))

  one <- function(k) {
    sname <- names(supports)[grid$support[k]]
    ed <- supports[[grid$support[k]]]
    test <- which(assign_fold == grid$fold[k])
    train <- setdiff(seq_len(n), test)
    err <- NULL
    out <- tryCatch({
      spec <- mc_between_graph(g$type, edges = ed,
                               ordering = if (g$type == "dag")
                                 g$spec$ordering else NULL)
      spec$R <- g$R
      Etr <- E[train, , drop = FALSE]
      Ctr <- scale(Etr, center = TRUE, scale = FALSE)
      Str <- crossprod(Ctr) / length(train)
      b <- .fit_between(Str, length(train), spec, max_iter = ctrl$max_iter,
                        tol = ctrl$tol, score_tol = ctrl$score_tol,
                        rcond_tol = ctrl$rcond_tol, verbose = FALSE,
                        dag_solver = ctrl$dag_solver)
      if (!isTRUE(b$converged))
        stop("the training fit did not converge: ", b$status)
      Ete <- E[test, , drop = FALSE]
      Cte <- scale(Ete, center = TRUE, scale = FALSE)
      Ste <- crossprod(Cte) / length(test)
      Sinv <- solve(b$Sigma_b)
      ld <- as.numeric(determinant(b$Sigma_b, logarithm = TRUE)$modulus)
      score <- -0.5 * (g$R * log(2 * pi) + ld + sum(Sinv * Ste))
      list(score = score, npar = b$npar)
    }, error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(out) || !is.finite(out$score))
      return(list(ok = FALSE, support = sname, fold = grid$fold[k],
                  error = if (is.null(err)) "non-finite score" else err))
    list(ok = TRUE, support = sname, fold = grid$fold[k],
         score = out$score, npar = out$npar)
  }

  res <- .mcg_lapply(seq_len(nrow(grid)), one, workers = workers)

  rows <- list(); failures <- list()
  for (k in seq_len(nrow(grid))) {
    z <- res[[k]]
    msg <- .mcg_error_message(z)
    if (!is.null(msg))
      z <- list(ok = FALSE, support = names(supports)[grid$support[k]],
                fold = grid$fold[k], error = msg)
    if (isTRUE(z$ok))
      rows[[length(rows) + 1L]] <- data.frame(
        support = z$support, fold = z$fold, score = z$score, npar = z$npar,
        stringsAsFactors = FALSE)
    else
      failures[[length(failures) + 1L]] <- data.frame(
        support = z$support, fold = z$fold, error = z$error,
        stringsAsFactors = FALSE)
  }
  if (!length(rows))
    stop("no fold produced a usable fit for any support.", call. = FALSE)
  per_fold <- do.call(rbind, rows)
  fail_df <- if (length(failures)) do.call(rbind, failures) else
    data.frame(support = character(), fold = integer(), error = character(),
               stringsAsFactors = FALSE)

  agg <- lapply(split(per_fold, per_fold$support), function(z) data.frame(
    support = z$support[1L], n_folds = nrow(z),
    npar = z$npar[1L], score = mean(z$score),
    se = if (nrow(z) > 1L) stats::sd(z$score) / sqrt(nrow(z)) else NA_real_,
    stringsAsFactors = FALSE))
  summ <- do.call(rbind, agg)
  summ <- summ[order(-summ$score), , drop = FALSE]
  summ$delta <- summ$score[1L] - summ$score
  rownames(summ) <- NULL

  structure(list(per_fold = per_fold, summary = summ, folds = folds,
                 failures = fail_df, n = n,
                 supports = names(supports)), class = "mcg_cv")
}

#' @rdname mc_graph_cv
#' @param x an object of class \code{"mcg_cv"}.
#' @param ... ignored.
#' @export
print.mcg_cv <- function(x, ...) {
  cat(sprintf("Out-of-sample comparison of %d support%s (%d folds, n = %d)\n",
              length(x$supports), if (length(x$supports) == 1L) "" else "s",
              x$folds, x$n))
  cat("Held-out working criterion per unit, larger is better.\n")
  if (nrow(x$failures))
    cat(sprintf("Fits that failed: %d\n", nrow(x$failures)))
  cat("\n")
  print(format(x$summary, digits = 4))
  invisible(x)
}
