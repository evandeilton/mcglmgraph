# class_members.R: enumerating the members of a Markov equivalence class.
# Part of the mcglmgraph package. License: GPL-3.
#
# A CPDAG stands for a set of DAGs. The fitted covariance, the working criterion
# and the implied concentration graph are properties of that set, but a
# structural coefficient need not be: it is a partial regression coefficient
# given the other parents of its head, and those can differ from one member to
# another. Reporting one value for such a coefficient reports a quantity that
# depends on a choice the data do not make. Enumerating the members is what lets
# the range be reported instead.

# The immoralities of a partially directed graph, as comparable keys.
.vstructure_keys <- function(amat) {
  storage.mode(amat) <- "integer"
  v <- cpp_vstructure_set(amat)
  if (!nrow(v)) return(character(0))
  unique(paste(v[, 1L], v[, 2L], v[, 3L], sep = "->"))
}

#' Members of a Markov equivalence class
#'
#' Enumerates the directed acyclic graphs (DAGs) that a completed partially
#' directed acyclic graph (CPDAG) stands for: the acyclic orientations of its undirected
#' edges that keep every directed edge and introduce no v-structure (immorality)
#' the class does not have. By the characterization of Markov equivalence by
#' skeleton and v-structures (Verma and Pearl, 1990; Andersson, Madigan and
#' Perlman, 1997), these are exactly the DAGs of the equivalence class. The
#' members imply the same separations, so observational data that carry only
#' those constraints do not distinguish among them. For a partially directed
#' graph that is not a CPDAG, such as one oriented further by background
#' knowledge, the same enumeration returns its consistent DAG extensions.
#'
#' @section How they are enumerated:
#' One undirected edge is oriented at a time, both ways, Meek's rules propagate
#' the consequences, and a branch is abandoned as soon as it creates a cycle or
#' an immorality the class does not have. Because the two branches fix opposite
#' directions on the same edge, no member is produced twice. The count grows
#' quickly with the number of undirected edges, so \code{max_members} caps it
#' and the result says whether the cap was reached.
#'
#' @param cpdag a square 0/1 adjacency matrix in the package convention (an
#'   undirected edge has both entries one, an arc \eqn{i \to j} only
#'   \code{[i, j]}), or an object of class \code{"mcg_cpdag"}.
#' @param max_members a single positive integer; the enumeration stops after
#'   this many members and reports the truncation.
#' @return A list of adjacency matrices, one per member, each a DAG in the same
#'   convention and with the dimnames of \code{cpdag}. The list carries the
#'   attribute \code{"truncated"}, \code{TRUE} when \code{max_members} was
#'   reached before the enumeration ended. The list is empty when the graph
#'   admits no consistent extension.
#' @references Verma, T. and Pearl, J. (1990). Equivalence and synthesis of
#'   causal models. \emph{Proceedings of the 6th Conference on Uncertainty in
#'   Artificial Intelligence}, 220-227.
#'
#'   Andersson, S. A., Madigan, D. and Perlman, M. D. (1997). A characterization
#'   of Markov equivalence classes for acyclic digraphs. \emph{The Annals of
#'   Statistics} 25(2), 505-541. \doi{10.1214/aos/1031833662}
#'
#'   Meek, C. (1995). Causal inference and causal explanation with background
#'   knowledge. \emph{Proceedings of the 11th Conference on Uncertainty in
#'   Artificial Intelligence}, 403-410.
#' @family equivalence class
#' @examples
#' # the chain 1 - 2 - 3 has three members, the collider none of them
#' a <- matrix(0L, 3, 3)
#' a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L
#' length(mc_class_members(a))
#' @seealso \code{\link{mc_class_coefficients}} for what the members disagree
#'   about, \code{\link{mc_orient_meek}} for the propagation rules.
#' @export
mc_class_members <- function(cpdag, max_members = 1000L) {
  if (inherits(cpdag, "mcg_cpdag")) cpdag <- cpdag$cpdag
  amat <- as.matrix(cpdag)
  if (!nrow(amat) || nrow(amat) != ncol(amat) || !is.numeric(amat) ||
      any(!amat %in% c(0, 1)) || any(diag(amat) != 0))
    stop("cpdag must be a square 0/1 adjacency with a zero diagonal.",
         call. = FALSE)
  storage.mode(amat) <- "integer"
  if (length(max_members) != 1L || !is.finite(max_members) || max_members < 1)
    stop("max_members must be one positive integer.", call. = FALSE)
  max_members <- as.integer(max_members)

  target <- .vstructure_keys(amat)
  nm <- dimnames(amat)
  out <- list()
  truncated <- FALSE

  recurse <- function(A) {
    if (truncated) return(invisible(NULL))
    und <- which(A == 1L & t(A) == 1L & upper.tri(A), arr.ind = TRUE)
    if (!nrow(und)) {
      if (!.has_directed_cycle(A) && identical(.vstructure_keys(A), target)) {
        if (length(out) >= max_members) { truncated <<- TRUE; return(invisible(NULL)) }
        dimnames(A) <- nm
        out[[length(out) + 1L]] <<- A
      }
      return(invisible(NULL))
    }
    i <- und[1L, 1L]; j <- und[1L, 2L]
    for (dir in list(c(i, j), c(j, i))) {
      B <- A
      B[dir[1L], dir[2L]] <- 1L
      B[dir[2L], dir[1L]] <- 0L
      B <- mc_orient_meek(B)
      attributes(B) <- list(dim = dim(B))
      if (.has_directed_cycle(B)) next
      if (!identical(.vstructure_keys(B), target)) next
      recurse(B)
      if (truncated) return(invisible(NULL))
    }
    invisible(NULL)
  }

  recurse(amat)
  attr(out, "truncated") <- truncated
  out
}

#' What the members of a class disagree about
#'
#' Fits the directed layer on the residual covariance of a fit under every
#' member of a Markov equivalence class and collects, for each arc, the values
#' of its structural coefficient that the class allows. The coefficients are
#' those of the closed-form root, \code{\link{mc_dag_closed_root}}, and the
#' standardized ones are \eqn{\gamma^*_{jk}=\gamma_{jk}(\Sigma_{b,kk}/
#' \Sigma_{b,jj})^{1/2}}, as in \code{\link{mc_dag_coefficients}}, with
#' \eqn{j} the child and \eqn{k} the parent. An arc whose values coincide
#' across members is one whose coefficient the class determines; an arc whose
#' values move is one for which a single number would report a choice the data
#' did not make.
#'
#' @section What the values are:
#' Every member reproduces the same fitted covariance when the class is a
#' CPDAG, and the second moments do not distinguish among them. The direction
#' of a compelled arc is a property of the class, but its coefficient is a
#' partial regression coefficient given the other parents of its head, and
#' those can change from one member to another. The values are coefficients of
#' linear projections of Pearson residuals. Enumerating a quantity over the
#' members of a class follows Maathuis, Kalisch and Buehlmann (2009), who apply
#' it to intervention effects; the coefficients reported here admit that
#' reading only under causal assumptions the procedure does not check: a linear
#' structural model with independent innovations whose graph belongs to the
#' class, causal sufficiency and absence of selection bias.
#'
#' @section Which arcs the class pins down:
#' An arc has the same coefficient in every member whose parent set for its head
#' is the same. That happens for every member when no undirected edge of the
#' CPDAG touches the head of the arc, which is the condition to check before
#' reporting a single value. The returned table flags it in \code{common}.
#'
#' @param fit an object of class \code{"mcglm_graph"}, whose residual covariance
#'   the layer is fitted to.
#' @param cpdag the class, as an adjacency matrix or an \code{"mcg_cpdag"}
#'   object; \code{NULL} runs \code{\link{mc_graph_pipeline}} on the fit.
#' @param max_members cap on the number of members enumerated, as in
#'   \code{\link{mc_class_members}}.
#' @param ... passed to \code{\link{mc_graph_pipeline}} when \code{cpdag} is
#'   \code{NULL}; ignored by the \code{print} method.
#' @return An object of class \code{"mcg_class_coef"}: a list with
#'   \code{arcs}, a data frame with one row per arc present in some member
#'   (\code{from}, \code{to}, \code{parent}, \code{child}, the number of
#'   members that contain it \code{n_members}, the range of its coefficient
#'   \code{gamma_min}, \code{gamma_max} and of its standardized coefficient
#'   \code{gamma_std_min}, \code{gamma_std_max}, and the logical
#'   \code{common}, true when no undirected edge of the class touches the head
#'   of the arc); \code{members}, the list of DAGs; \code{values}, the
#'   per-member table of coefficients (\code{member}, \code{from}, \code{to},
#'   \code{gamma}, \code{gamma_std}); \code{n_members}; \code{truncated},
#'   whether the enumeration was capped; and \code{response_names}. The
#'   \code{print} method returns \code{x} invisibly. An error is raised when the
#'   class has no consistent extension.
#' @references Maathuis, M. H., Kalisch, M. and Buehlmann, P. (2009).
#'   Estimating high-dimensional intervention effects from observational data.
#'   \emph{The Annals of Statistics} 37(6A), 3133-3164.
#'   \doi{10.1214/09-AOS685}
#' @family equivalence class
#' @seealso \code{\link{mc_class_members}}, \code{\link{mc_dag_closed_root}}
#' @examples
#' \donttest{
#' set.seed(1); n <- 400; x <- rnorm(n)
#' z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
#' d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' mc_class_coefficients(f)
#' }
#' @export
mc_class_coefficients <- function(fit, cpdag = NULL, max_members = 1000L, ...) {
  stopifnot(inherits(fit, "mcglm_graph"))
  if (is.null(cpdag)) cpdag <- mc_graph_pipeline(fit, ...)
  cp <- if (inherits(cpdag, "mcg_cpdag")) cpdag$cpdag else as.matrix(cpdag)
  members <- mc_class_members(cp, max_members = max_members)
  if (!length(members))
    stop("the class has no consistent extension, so there is nothing to fit.",
         call. = FALSE)

  g <- fit$graph
  R <- g$R
  rn <- g$response_names
  Se <- g$Shat
  storage.mode(cp) <- "integer"
  # a head free of undirected edges has the same parent set in every member
  free_head <- vapply(seq_len(R), function(j)
    !any(cp[j, ] == 1L & cp[, j] == 1L), logical(1))

  rows <- list()
  for (m in seq_along(members)) {
    A <- members[[m]]
    arcs <- which(A == 1L & t(A) == 0L, arr.ind = TRUE)
    arcs <- matrix(as.integer(arcs), ncol = 2L)
    root <- mc_dag_closed_root(Se, R = R, edges = arcs, directed = TRUE)
    Sd <- diag(root$Sigma_b)
    if (nrow(root$edges))
      rows[[length(rows) + 1L]] <- data.frame(
        member = m,
        from = as.integer(root$edges[, 1L]), to = as.integer(root$edges[, 2L]),
        gamma = root$gamma,
        gamma_std = root$gamma * sqrt(Sd[root$edges[, 1L]] /
                                        Sd[root$edges[, 2L]]),
        stringsAsFactors = FALSE)
  }
  values <- if (length(rows)) do.call(rbind, rows) else
    data.frame(member = integer(), from = integer(), to = integer(),
               gamma = numeric(), gamma_std = numeric())

  key <- paste(values$from, values$to, sep = ":")
  agg <- lapply(split(seq_len(nrow(values)), key), function(ix) {
    z <- values[ix, , drop = FALSE]
    data.frame(
      from = z$from[1L], to = z$to[1L],
      parent = rn[z$from[1L]], child = rn[z$to[1L]],
      n_members = nrow(z),
      gamma_min = min(z$gamma), gamma_max = max(z$gamma),
      gamma_std_min = min(z$gamma_std), gamma_std_max = max(z$gamma_std),
      common = free_head[z$to[1L]],
      stringsAsFactors = FALSE)
  })
  arcs_df <- if (length(agg)) do.call(rbind, agg) else
    data.frame(from = integer(), to = integer(), parent = character(),
               child = character(), n_members = integer(),
               gamma_min = numeric(), gamma_max = numeric(),
               gamma_std_min = numeric(), gamma_std_max = numeric(),
               common = logical())
  rownames(arcs_df) <- NULL
  arcs_df <- arcs_df[order(arcs_df$from, arcs_df$to), , drop = FALSE]

  structure(list(arcs = arcs_df, members = members, values = values,
                 n_members = length(members),
                 truncated = isTRUE(attr(members, "truncated")),
                 response_names = rn), class = "mcg_class_coef")
}

#' @rdname mc_class_coefficients
#' @param x an object of class \code{"mcg_class_coef"}.
#' @export
print.mcg_class_coef <- function(x, ...) {
  cat(sprintf("Coefficients over the %d member%s of the class%s\n",
              x$n_members, if (x$n_members == 1L) "" else "s",
              if (x$truncated) " (enumeration truncated)" else ""))
  cat("An arc marked common has the same value in every member.\n\n")
  show <- x$arcs[, c("parent", "child", "n_members", "gamma_min", "gamma_max",
                     "gamma_std_min", "gamma_std_max", "common")]
  print(format(show, digits = 3))
  invisible(x)
}
