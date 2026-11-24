# parallel.R: optional worker support for the embarrassingly parallel parts.
# Part of the mcglmgraph package. License: GPL-3.
#
# Five routines in this package do genuinely independent repeated work and go
# through .mcg_lapply(): the half-samples of mc_graph_stability(), the penalty
# grid of mc_select_lambda(), the resamples of mc_graph_bootstrap() and of
# mc_cpdag_stability(), and the folds of mc_graph_cv(). All of them take
# workers, opt-in and sequential by default, for one reason: calling code often
# parallelizes at an outer level, with its own PSOCK cluster. Turning on workers
# underneath that would oversubscribe the machine, which is slower than staying
# sequential and is a failure mode that is hard to see from inside a worker.
#
# Reproducibility: every stochastic draw is made in the parent, before any work
# is dispatched. A given seed therefore produces the same subsamples whatever
# the worker count, and `workers = 4` returns exactly what `workers = 1` returns.

# TRUE when the current process is a PSOCK worker: those are started with
# `--args MASTER=<host> PORT=<n> ...`, which nothing else in an ordinary R
# session sets.
.inside_psock_worker <- function() {
  any(grepl("^MASTER=", commandArgs(trailingOnly = FALSE)))
}

# Resolve the requested worker count against the machine and the context.
.resolve_workers <- function(workers) {
  if (is.null(workers)) return(1L)
  if (length(workers) != 1L || !is.finite(workers) || workers < 1 ||
      workers != as.integer(workers))
    stop("workers must be a single positive integer.", call. = FALSE)
  workers <- as.integer(workers)
  if (workers == 1L) return(1L)
  if (.inside_psock_worker()) {
    warning("Running inside a parallel worker; falling back to workers = 1. ",
            "Nesting parallelism oversubscribes the machine. Parallelize at ",
            "one level only.", call. = FALSE)
    return(1L)
  }
  avail <- tryCatch(parallel::detectCores(logical = TRUE),
                    error = function(e) NA_integer_)
  if (is.finite(avail) && workers > avail) {
    warning("workers = ", workers, " exceeds the ", avail,
            " cores detected; using ", avail, ".", call. = FALSE)
    workers <- as.integer(avail)
  }
  workers
}

# lapply with an optional fork-based backend. Falls back to sequential whenever
# forking is unavailable (Windows) so behavior never depends on the platform.
.mcg_lapply <- function(X, FUN, workers = 1L, ...) {
  workers <- .resolve_workers(workers)
  if (workers == 1L || .Platform$OS.type != "unix")
    return(lapply(X, FUN, ...))
  parallel::mclapply(X, FUN, ..., mc.cores = workers,
                     mc.preschedule = FALSE)
}

# mclapply reports a failed element as a "try-error" object rather than raising,
# so callers must convert those back into the error contract they expect.
.mcg_error_message <- function(x) {
  if (!inherits(x, "try-error")) return(NULL)
  cond <- attr(x, "condition")
  if (is.null(cond)) as.character(x) else conditionMessage(cond)
}
