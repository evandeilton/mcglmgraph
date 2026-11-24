# twostage.R: propagation of the first-stage uncertainty into the graph layer.
# Part of the mcglmgraph package. License: GPL-3.
#
# Treating the Pearson residuals as observed keeps the point estimate of the
# graph layer correct but not its standard errors. Stacking the two sets of
# estimating functions gives a block triangular sensitivity whose inversion
# yields the influence function of the second stage; see the section of the same
# name in mc_vcov_graph(). Two routes form it: the shortcut, the reduced form
# evaluated at the fitted covariance, which needs second-stage quantities only
# and applies to a restricted class of first stages; and the full form, which
# rebuilds the first stage unit by unit and evaluates B at the data.

# Variance functions of mcglm whose dispersion parameter multiplies the whole
# variance, so that v_ir = tau_r V_r(mu_ir; p_r). The Poisson-Tweedie function
# is additive, v = mu + tau mu^p, and is multiplicative only at p = 1, where it
# collapses to (1 + tau) mu.
.MULTIPLICATIVE_VARIANCES <- c("constant", "tweedie", "binomialP", "binomialPQ")

# Is the reduced form, and with it the shortcut, available for this fit, and if
# not, why not. Only the first-stage structure is checked here; on a DAG with
# an immorality the shortcut is then first-order only for the parameters of
# the collider and its descendants (see .dag_immoralities).
.two_stage_applicable <- function(fit) {
  g <- fit$graph
  why <- character(0)

  vr <- tryCatch(unlist(fit$variance), error = function(e) NULL)
  if (is.null(vr) || length(vr) != g$R)
    return(list(ok = FALSE, why = "the variance functions of the first stage could not be read"))

  npc <- tryCatch(lengths(fit$matrix_pred), error = function(e) NULL)
  if (is.null(npc) || length(npc) != g$R) {
    why <- c(why, "the covariance components of the first stage could not be read")
  } else if (any(npc != 1L)) {
    bad <- which(npc != 1L)
    why <- c(why, paste0("response", if (length(bad) > 1L) "s " else " ",
                         paste(g$response_names[bad], collapse = ", "),
                         " carr", if (length(bad) > 1L) "y" else "ies",
                         " more than one dispersion parameter"))
  } else if (!.identity_matrix_pred(fit)) {
    # A single dispersion parameter is not enough. The identity the reduced form
    # rests on comes from the first-stage Pearson equations being the score of
    # a working Gaussian model that factorizes over units, which needs the
    # within-response structure to be a multiple of the identity. One
    # covariance component that is not diagonal, a fixed correlation among the
    # rows of a response say, keeps the parameter count at one and breaks the
    # factorization all the same.
    why <- c(why, paste0("the within-response structure is not a multiple of ",
                         "the identity, so the rows of a response are not ",
                         "independent and the first stage does not fix the ",
                         "average of their squared residuals"))
  }

  unknown <- setdiff(vr, c(.MULTIPLICATIVE_VARIANCES, "poisson_tweedie"))
  if (length(unknown))
    why <- c(why, paste0("the variance function ",
                         paste(unique(unknown), collapse = ", "),
                         " is not one the shortcut is established for"))

  pt <- which(vr == "poisson_tweedie")
  if (length(pt)) {
    fixed <- tryCatch(vapply(fit$power_fixed[pt], isTRUE, logical(1)),
                      error = function(e) rep(FALSE, length(pt)))
    pw <- tryCatch(vapply(fit$list_initial$power[pt], function(z) z[1L],
                          numeric(1)),
                   error = function(e) rep(NA_real_, length(pt)))
    ok_pt <- fixed & is.finite(pw) & abs(pw - 1) < 1e-8
    if (!all(ok_pt))
      why <- c(why, paste0("response",
                           if (sum(!ok_pt) > 1L) "s " else " ",
                           paste(g$response_names[pt[!ok_pt]], collapse = ", "),
                           " ", if (sum(!ok_pt) > 1L) "are" else "is",
                           " Poisson-Tweedie with a power that is not fixed at one"))
  }

  # With prior weights that vary, the first stage fixes a weighted average of
  # the squared residuals, not their plain average.
  w <- tryCatch(.first_stage_weights(fit, length(fit$residuals)),
                error = function(e) NULL)
  if (is.null(w)) {
    why <- c(why, "the prior weights of the first stage could not be read")
  } else if (length(w) && diff(range(w)) > 1e-12 * max(abs(w))) {
    why <- c(why, paste0("the prior weights of the first stage are not ",
                         "constant, so it fixes a weighted average of the ",
                         "squared residuals rather than their average"))
  }
  list(ok = length(why) == 0L, why = why)
}

# Colliders of the fitted DAG with two non-adjacent parents (immoralities).
# For the parameters of those nodes and of their descendants the shortcut,
# which reads the fitted covariance where the reduced form reads the sample
# one, is exact only to first order and under correct specification.
.dag_immoralities <- function(fit) {
  g <- fit$graph
  if (!identical(g$type, "dag")) return(integer(0))
  eo <- g$meta$edges_original
  if (is.null(eo)) eo <- .dag_setup(g$spec, g$R)$edges_original
  if (!nrow(eo)) return(integer(0))
  adj <- matrix(FALSE, g$R, g$R)
  adj[eo] <- TRUE
  adj <- adj | t(adj)
  out <- integer(0)
  for (j in unique(eo[, 2L])) {
    pa <- eo[eo[, 2L] == j, 1L]
    if (length(pa) < 2L) next
    pr <- utils::combn(pa, 2L)
    if (any(!adj[cbind(pr[1L, ], pr[2L, ])])) out <- c(out, j)
  }
  sort(as.integer(out))
}

# Route taken by method = "auto": the full form on a DAG with an immorality
# when it is available, the shortcut when it applies, the full form otherwise.
.two_stage_auto_method <- function(fit, cluster = NULL,
                                   sc = .two_stage_applicable(fit),
                                   warn = FALSE) {
  imm <- .dag_immoralities(fit)
  if (length(imm)) {
    fa <- .two_stage_full_applicable(fit, cluster)
    if (fa$ok) return("full")
    if (sc$ok) {
      if (warn)
        warning("the fitted DAG has an immorality at ",
                paste(fit$graph$response_names[imm], collapse = ", "),
                " and the full two-stage form is not available (",
                paste(fa$why, collapse = "; "), "); the shortcut is used, ",
                "and its correction is first-order only for the parameters ",
                "of colliders and their descendants.", call. = FALSE)
      return("shortcut")
    }
  }
  if (sc$ok) "shortcut" else "full"
}

# ---- the full form -----------------------------------------------------------

# Internal routines of mcglm that rebuild its first stage. They are reached
# through the namespace rather than by `:::`, and a missing one is reported by
# name, since a change in mcglm would otherwise surface as an obscure error.
.mcglm_internal <- function(name) {
  f <- tryCatch(utils::getFromNamespace(name, "mcglm"), error = function(e) NULL)
  if (!is.function(f))
    stop("the installed mcglm does not provide the internal routine ", name,
         ", which the full two-stage correction rebuilds the first stage with.",
         call. = FALSE)
  f
}

# The units the contributions are summed over: the clusters when given, the
# rows otherwise. Codes run from 1 to the number of units.
.two_stage_units <- function(fit, cluster) {
  n <- fit$graph$n
  if (is.null(cluster)) return(list(unit = seq_len(n), n_units = n))
  if (length(cluster) != n)
    stop("cluster must have one identifier per observational row.", call. = FALSE)
  if (anyNA(cluster))
    stop("cluster identifiers must not be missing.", call. = FALSE)
  f <- as.factor(cluster)
  if (nlevels(f) < 2L)
    stop("at least two independent clusters are required.", call. = FALSE)
  list(unit = as.integer(f), n_units = nlevels(f))
}

# Does every within-response matrix vanish between distinct units? That is what
# makes the first stage a sum over units, and it is what a cluster has to
# satisfy when the structure is not the identity.
.structure_nested_in <- function(fit, unit) {
  for (r in seq_along(fit$matrix_pred)) for (Z in fit$matrix_pred[[r]]) {
    if (.is_identity_structure(Z)) next
    if (isS4(Z)) {
      Tz <- methods::as(Z, "TsparseMatrix")
      i <- Tz@i + 1L; j <- Tz@j + 1L
      x <- if (methods::.hasSlot(Tz, "x")) Tz@x else rep(1, length(i))
    } else {
      ij <- which(as.matrix(Z) != 0, arr.ind = TRUE)
      i <- ij[, 1L]; j <- ij[, 2L]; x <- rep(1, length(i))
    }
    keep <- x != 0
    if (any(unit[i[keep]] != unit[j[keep]])) return(FALSE)
  }
  TRUE
}

# Is the full form available, and if not, why not. Only the structure is read
# here; whether the rebuilt first stage reproduces the fitted covariance is
# checked when the correction is formed.
.two_stage_full_applicable <- function(fit, cluster = NULL) {
  need <- c("list_initial", "Regression", "Covariance", "Information",
            "list_X", "link", "variance", "covariance", "matrix_pred",
            "power_fixed", "Ntrial", "offset", "sparse", "residuals", "C",
            "linear_pred", "n_obs")
  miss <- need[!need %in% names(fit)]
  if (length(miss))
    return(list(ok = FALSE, why = paste0(
      "the fit does not carry the first-stage components ",
      paste(miss, collapse = ", "))))
  g <- fit$graph
  why <- character(0)
  if (!identical(as.integer(fit$n_obs), as.integer(g$n)) ||
      length(g$residuals) != length(fit$residuals))
    why <- c(why, paste0("the graph layer does not sit on the residuals of ",
                         "this first stage, so the two cannot be stacked"))
  if (!.identity_matrix_pred(fit)) {
    if (is.null(cluster)) {
      why <- c(why, paste0("the within-response structure is not the identity, ",
                           "so the rows are not the units the first stage sums ",
                           "over; pass cluster with one identifier per row, ",
                           "grouping the rows the structure connects"))
    } else {
      u <- tryCatch(.two_stage_units(fit, cluster), error = function(e) NULL)
      if (is.null(u)) {
        why <- c(why, "cluster is not a valid grouping of the rows")
      } else if (!.structure_nested_in(fit, u$unit)) {
        why <- c(why, paste0("the within-response structure connects rows of ",
                             "different clusters; every matrix of matrix_pred ",
                             "must vanish between distinct clusters"))
      }
    }
  }
  list(ok = length(why) == 0L, why = why)
}

# Weights of the first stage, one per stacked row.
.first_stage_weights <- function(fit, N) {
  W <- fit$W
  if (is.null(W)) return(rep(1, N))
  if (inherits(W, "ddiMatrix")) {
    w <- if (identical(attr(W, "diag"), "U")) rep(1, N) else attr(W, "x")
  } else {
    w <- tryCatch(as.numeric(diag(W)), error = function(e) NULL)
  }
  if (length(w) != N || any(!is.finite(w)))
    stop("the weights of the first stage could not be read.", call. = FALSE)
  w
}

.as_dgc <- function(M) {
  if (is.null(M) || inherits(M, "dgCMatrix")) return(M)
  methods::as(methods::as(M, "CsparseMatrix"), "generalMatrix")
}

# The first stage of a fit at a given parameter vector, rebuilt with mcglm's own
# routines: the means and their derivatives, and C with its derivatives.
.first_stage_build <- function(fit, alpha, derivative_beta = TRUE) {
  R <- length(fit$linear_pred)
  inf <- fit$Information
  pb <- length(fit$Regression)
  li <- .mcglm_internal("mc_updateBeta")(fit$list_initial, alpha[seq_len(pb)],
                                         information = inf, n_resp = R)
  li <- .mcglm_internal("mc_updateCov")(li, covariance = alpha[-seq_len(pb)],
                                        list_power_fixed = fit$power_fixed,
                                        information = inf, n_resp = R)
  mu_list <- Map(.mcglm_internal("mc_link_function"), beta = li$regression,
                 offset = fit$offset, X = fit$list_X, link = fit$link)
  Cf <- .mcglm_internal("mc_build_C")(
    list_mu = mu_list, list_Ntrial = fit$Ntrial, rho = li$rho,
    list_tau = li$tau, list_power = li$power, list_Z = fit$matrix_pred,
    list_sparse = fit$sparse, list_variance = fit$variance,
    list_covariance = fit$covariance, list_power_fixed = fit$power_fixed,
    compute_C = TRUE, compute_derivative_beta = derivative_beta,
    compute_derivative_cov = TRUE)
  list(mu_list = mu_list, Cf = Cf)
}

# Did the first stage solve its Pearson equations with mcglm's bias correction?
# A fit that does not record the option used mcglm's default, which is TRUE.
.first_stage_corrected <- function(fit) {
  co <- fit$con$correct
  if (is.null(co)) TRUE else isTRUE(co)
}

# Totals of the first-stage estimating functions at alpha, as the first stage
# solved them: with `correct`, the Pearson functions carry mcglm's bias
# correction. This is the function whose Jacobian the observed sensitivity is.
.first_stage_totals <- function(fit, alpha, y, w,
                                correct = .first_stage_corrected(fit)) {
  b <- .first_stage_build(fit, alpha, derivative_beta = FALSE)
  mu <- unlist(lapply(b$mu_list, function(m) as.numeric(m$mu)))
  U <- as.numeric(cpp_first_stage_totals(
    .as_dgc(b$Cf$inv_C), lapply(b$mu_list, function(m) as.matrix(m$D)),
    lapply(b$Cf$D_C, .as_dgc), y - mu, w))
  if (!correct) return(U)
  # mcglm's own correction term, with the inverse quasi-score sensitivity
  # -(D' C^-1 W D)^-1 evaluated at alpha, as its fitting algorithm forms it
  D <- Matrix::bdiag(lapply(b$mu_list, function(m) m$D))
  Sb <- -as.matrix(Matrix::crossprod(D, b$Cf$inv_C %*%
                                       (Matrix::Diagonal(x = w) %*% D)))
  cr <- .mcglm_internal("mc_correction")(D_C = b$Cf$D_C, inv_J_beta = solve(Sb),
                                         D = D, inv_C = b$Cf$inv_C)
  U + c(rep(0, ncol(D)), as.numeric(cr))
}

# The observed sensitivity of the first stage, d U_1 / d alpha' at the estimate,
# by central differences of the totals, on mcglm's sign convention. A step at
# which the equations cannot be evaluated (a covariance that is not positive
# definite, a singular quasi-score sensitivity) falls back to a one-sided
# difference.
.first_stage_observed <- function(fit) {
  alpha <- c(fit$Regression, fit$Covariance)
  y <- as.numeric(fit$fitted) + as.numeric(fit$residuals)
  w <- .first_stage_weights(fit, length(y))
  correct <- .first_stage_corrected(fit)
  # a missing correction routine is reported by name, not as a failed step
  if (correct) .mcglm_internal("mc_correction")
  U <- function(a) tryCatch(.first_stage_totals(fit, a, y, w, correct),
                            error = function(e) NULL)
  U0 <- NULL
  p <- length(alpha)
  S <- matrix(0, p, p)
  for (m in seq_len(p)) {
    h <- 1e-5 * max(1, abs(alpha[m]))
    up <- alpha; up[m] <- up[m] + h
    dn <- alpha; dn[m] <- dn[m] - h
    Up <- U(up); Dn <- U(dn)
    if (!is.null(Up) && !is.null(Dn)) {
      S[, m] <- (Up - Dn) / (2 * h)
    } else {
      if (is.null(U0)) U0 <- U(alpha)
      if (!is.null(Up)) S[, m] <- (Up - U0) / h
      else if (!is.null(Dn)) S[, m] <- (U0 - Dn) / h
      else stop("the observed sensitivity of the first stage could not be ",
                "formed: the first-stage equations could not be evaluated on ",
                "either side of parameter ", m, ".", call. = FALSE)
    }
  }
  S
}

# The first stage at its estimate, rebuilt with mcglm's own routines and split
# by unit: contributions, expected sensitivity (the cross block included) and
# the derivative of the Pearson residuals.
.first_stage_system <- function(fit, unit, n_units) {
  R <- length(fit$linear_pred)
  inf <- fit$Information
  fb <- .first_stage_build(fit, c(fit$Regression, fit$Covariance))
  mu_list <- fb$mu_list
  Cf <- fb$Cf

  # mcglm returns no derivative in beta for a response whose variance does not
  # depend on the mean, and concatenates the others; they are placed back
  # response by response so that every column lines up with its parameter.
  nb <- as.integer(unlist(inf$n_betas))
  off <- c(0L, cumsum(nb))[seq_len(R)]
  has <- vapply(fit$variance, function(v) !identical(v, "constant"), logical(1))
  idx <- unlist(lapply(which(has), function(r) off[r] + seq_len(nb[r])))
  if (length(idx) != length(Cf$D_C_beta))
    stop("the derivatives of the first-stage covariance in the regression ",
         "parameters do not line up with those parameters.", call. = FALSE)
  DCbeta <- vector("list", sum(nb))
  for (k in seq_along(idx)) DCbeta[[idx[k]]] <- .as_dgc(Cf$D_C_beta[[k]])

  C <- .as_dgc(Cf$C)
  res <- as.numeric(fit$residuals)
  out <- cpp_first_stage_units(
    .as_dgc(Cf$inv_C), C, lapply(mu_list, function(m) as.matrix(m$D)),
    lapply(Cf$D_C, .as_dgc), DCbeta, res,
    .first_stage_weights(fit, length(res)),
    as.integer(unit) - 1L, as.integer(n_units))
  out$r <- as.numeric(out$r)

  # The rebuilt covariance must be the fitted one, or everything below refers
  # to another model.
  gap <- tryCatch(cpp_sp_max_abs_diff(C, .as_dgc(fit$C)),
                  error = function(e) Inf)
  if (!is.finite(gap) || gap > 1e-8 * max(1, max(abs(attr(C, "x")))))
    stop("the first stage rebuilt from the fitted parameters does not ",
         "reproduce the covariance of the fit.", call. = FALSE)
  out
}

# Corrected contributions by the full form, one row per unit:
# c_u = psi_2u - B_n A_n^-1 psi_1u, with A_n = -S the sensitivity of the first
# stage on the accumulated scale.
.two_stage_full_scores <- function(fit, unit, n_units,
                                  sensitivity = c("observed", "expected")) {
  sensitivity <- match.arg(sensitivity)
  g <- fit$graph
  fs <- .first_stage_system(fit, unit, n_units)
  E <- as.matrix(g$residuals)
  if (length(E) != length(fs$r) ||
      max(abs(as.numeric(E) - fs$r)) > 1e-8 * max(1, max(abs(fs$r))))
    stop("the graph layer was not fitted on the Pearson residuals of this ",
         "first stage, so the two stages cannot be stacked.", call. = FALSE)
  Ec <- scale(E, center = TRUE, scale = FALSE)
  Sinv <- solve(g$Sigma_b)
  Bn <- cpp_two_stage_cross(Ec, Sinv, g$Dlist, fs$G)

  # The bias correction of mcglm is one term for the whole sample, not a sum
  # over units, so the unit contributions omit it and sum to minus it at the
  # estimate. Centering them adds it back in equal shares, and a term common
  # to all units does not enter the outer products.
  psi1 <- sweep(fs$psi, 2L, colMeans(fs$psi))
  A <- -(if (sensitivity == "observed") .first_stage_observed(fit) else fs$S)
  M <- tryCatch(t(solve(t(A), t(Bn))), error = function(e) NULL)
  if (is.null(M) || any(!is.finite(M)))
    stop("the sensitivity of the first stage could not be inverted.",
         call. = FALSE)
  Psi2 <- rowsum(.pearson_scores(g$residuals, g$Sigma_b, g$Dlist), unit,
                 reorder = TRUE)
  list(scores = Psi2 - psi1 %*% t(M), B = Bn, A = A)
}

# Corrected sandwich H_n^-1 (sum_u c_u c_u') H_n^-1, by the shortcut or by
# the full form; H_n differs from the model Gram only for a misspecified DAG.
.two_stage_vcov <- function(fit, cluster = NULL, small_sample = FALSE,
                            method = c("auto", "shortcut", "full"),
                            sensitivity = c("observed", "expected")) {
  method <- match.arg(method)
  sensitivity <- match.arg(sensitivity)
  g <- fit$graph
  if (is.null(g$Dlist))
    stop("The fit does not carry the covariance derivatives.", call. = FALSE)
  sc <- .two_stage_applicable(fit)
  if (method == "auto") method <- .two_stage_auto_method(fit, cluster, sc,
                                                         warn = TRUE)
  if (method == "shortcut" && !sc$ok)
    stop("the two-stage shortcut does not apply to this fit: ",
         paste(sc$why, collapse = "; "),
         ". method = \"full\" forms the correction from the first-stage ",
         "contributions instead.", call. = FALSE)
  if (method == "full") {
    fa <- .two_stage_full_applicable(fit, cluster)
    if (!fa$ok)
      stop("the two-stage correction is not available for this fit. ",
           "The shortcut does not apply",
           if (length(sc$why)) paste0(" (", paste(sc$why, collapse = "; "), ")")
           else "",
           ", and the full form does not either: ",
           paste(fa$why, collapse = "; "), ".", call. = FALSE)
  }

  u <- .two_stage_units(fit, cluster)
  extra <- list()
  if (method == "shortcut") {
    Ec <- scale(g$residuals, center = TRUE, scale = FALSE)
    Psi <- .pearson_scores(g$residuals, g$Sigma_b, g$Dlist)
    ts <- cpp_two_stage_scores(Psi, as.matrix(g$K_b), g$Dlist, Ec)
    C <- ts$scores
    if (!is.null(cluster)) C <- rowsum(C, u$unit, reorder = FALSE)
    extra$Q <- ts$Q
    extra$centre <- as.numeric(ts$centre)
  } else {
    fl <- .two_stage_full_scores(fit, u$unit, u$n_units, sensitivity)
    C <- fl$scores
    extra$B <- fl$B
    extra$A <- fl$A
    extra$sensitivity <- sensitivity
    extra$sensitivity_rcond <- tryCatch(rcond(fl$A), error = function(e) NA_real_)
  }
  meat <- crossprod(C)
  if (!is.null(cluster) && isTRUE(small_sample))
    meat <- meat * u$n_units / (u$n_units - 1)

  rct <- g$rcond_tol
  if (is.null(rct) || !is.finite(rct)) rct <- .Machine$double.eps^(2 / 3)
  Jinv <- .safe_solve(.graph_inference_bread(fit), rcond_tol = rct)
  M <- Jinv %*% meat %*% Jinv
  M <- 0.5 * (M + t(M))
  attr(M, "inversion_status") <- attr(Jinv, "inversion_status", exact = TRUE)
  attr(M, "inversion_ridge") <- attr(Jinv, "inversion_ridge", exact = TRUE)
  attr(M, "two_stage_method") <- method
  for (nm in names(extra)) attr(M, nm) <- extra[[nm]]
  M
}

#' Is the two-stage correction available for a fit
#'
#' Reports which of the two routes that propagate the first-stage uncertainty
#' into the graph layer apply to a fit, and if one does not, why. The shortcut
#' is the reduced form evaluated at the fitted covariance (see the section
#' "Propagating the first stage" of \code{\link{mc_vcov_graph}}). The reduced
#' form needs every response to carry a single dispersion parameter that
#' multiplies its variance function, \eqn{v_{ri}=\tau_{r}V_r(\mu_{ri};p_r)},
#' an identity within-response structure and constant prior weights, because
#' together they make the first-stage Pearson equations fix the average of the
#' squared residuals of each response and let the correction be written, to
#' first order, with second-stage quantities alone. It also needs the first
#' stage to estimate the correlation between responses without restriction,
#' which \pkg{mcglm} does and the check therefore does not read. The check
#' below is that structural one: on a
#' DAG with an immorality the shortcut can qualify and still be, for the
#' parameters of the collider and of its descendants, only a first-order
#' approximation to the reduced form, which is why
#' \code{method = "auto"} prefers the full form there. The full form needs
#' only that the first stage be a sum over units, which is what the rows are
#' when the within-response structure is the identity, and what the clusters
#' must be otherwise.
#'
#' @section Which fits qualify for the shortcut:
#' The variance functions \code{"constant"}, \code{"tweedie"},
#' \code{"binomialP"} and \code{"binomialPQ"} are multiplicative in the
#' dispersion parameter and qualify. \code{"poisson_tweedie"} is additive,
#' \eqn{v=\mu+\tau\mu^{p}}, and qualifies only at \eqn{p=1}, where it collapses
#' to \eqn{(1+\tau)\mu}; the power must be fixed there rather than estimated,
#' since an estimated power is a first-stage parameter of its own. A response
#' with more than one covariance component has more than one dispersion
#' parameter and does not qualify. Prior weights must be constant:
#' with weights that vary, the first stage fixes a weighted average of the
#' squared residuals rather than their average.
#'
#' A single dispersion parameter is necessary and not sufficient. The identity
#' comes from the first-stage equations being the score of a working Gaussian
#' model that factorizes over units, which needs the within-response structure
#' to be a multiple of the identity. One covariance component that is not
#' diagonal keeps the parameter count at one and breaks the factorization all
#' the same, so it is refused too.
#'
#' @section Which fits qualify for the full form:
#' Any family, any number of covariance components and any power, estimated or
#' fixed, provided the fit carries its first-stage components and the graph
#' layer sits on its Pearson residuals. With an identity within-response
#' structure the units are the rows, or the clusters when \code{cluster} is
#' given. With any other structure \code{cluster} is required and has to nest
#' it: every matrix of the matrix linear predictor must vanish between rows of
#' distinct clusters, which is checked. A random-effect block per subject nested
#' in a cluster per hunter, say, qualifies with the hunter as the cluster.
#'
#' @section What the check cannot verify:
#' The check reads the structure of the first-stage fit and nothing else: the
#' variance functions, the number of dispersion parameters, whether the
#' within-response structure is the identity or nests in the clusters, whether
#' a Poisson-Tweedie power was fixed at one and whether the prior weights are
#' constant. For the shortcut, that
#' structure makes the first stage fix the average of the squared residuals,
#' which is the algebraic half of it. The statistical half, that the first
#' stage moves the residual covariance along a diagonal rescaling, rests on
#' statements about the data that cannot be read off the fit:
#' \itemize{
#'   \item a correctly specified mean, condition (b) of the section "What the
#'     correction reaches" of \code{\link{mc_vcov_graph}}, which the constant
#'     variance function with a design common to all responses does without;
#'   \item a conditional residual covariance that does not depend on the
#'     covariates, condition (c) there, which the constant variance function
#'     does without, and which a variance function that depends on the mean
#'     needs at least in the weaker form of second moments that change with the
#'     covariates only through a factor common to the responses;
#'   \item independence between the units, that is between the rows once the
#'     within-response structure is the identity, or between the clusters.
#' }
#' The first two concern the shortcut only; the full form, with its observed
#' sensitivity, does not need them. The third concerns every form. A fit can
#' pass the check and still violate them. The spatial example is the plain
#' case: a soil survey fitted with a trend surface in the mean and an identity
#' within-response structure qualifies, while its neighboring sites remain
#' dependent, since coordinates in the mean do not whiten the residuals. There
#' the rows are not independent units, and neither the corrected form nor the
#' sandwich it corrects has the variance it claims.
#'
#' @param fit an object of class \code{"mcglm_graph"}.
#' @param cluster optional vector with one identifier per observational row,
#'   as in \code{\link{mc_vcov_graph}}; default \code{NULL}. It changes only
#'   whether the full form is available for a within-response structure that is
#'   not the identity, and, on a DAG with an immorality, which route
#'   \code{method = "auto"} takes.
#' @return A list with \code{shortcut_ok} (logical) and \code{why_shortcut}
#'   (character vector, empty when the shortcut applies), the verdict and the
#'   reasons for the shortcut; \code{full_ok} and \code{why_full}, the same for
#'   the full form; and \code{method}, a string, \code{"shortcut"} or
#'   \code{"full"}, the route \code{mc_vcov_graph(type = "two_stage")} takes by
#'   default (\code{method = "auto"}), or \code{NA_character_} when neither
#'   applies. The components \code{ok} and \code{why} repeat
#'   \code{shortcut_ok} and \code{why_shortcut}.
#' @seealso \code{\link{mc_vcov_graph}}, whose \code{type = "two_stage"}
#'   consumes this check, and \code{\link{mc_graph_bootstrap}}, which
#'   propagates the first stage by refitting.
#' @family inference on the graph layer
#' @references
#' Murphy, K. M. and Topel, R. H. (1985). Estimation and inference in two-step
#' econometric models. \emph{Journal of Business and Economic Statistics} 3(4),
#' 370--379. \doi{10.1080/07350015.1985.10509471}
#'
#' Newey, W. K. and McFadden, D. (1994). Large sample estimation and hypothesis
#' testing. In Engle, R. F. and McFadden, D. L. (eds), \emph{Handbook of
#' Econometrics}, vol. 4, 2111--2245. Elsevier.
#' \doi{10.1016/S1573-4412(05)80005-4}
#' @examples
#' \donttest{
#' set.seed(1); n <- 200; x <- rnorm(n)
#' d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
#' f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#' mc_two_stage_available(f)
#' }
#' @export
mc_two_stage_available <- function(fit, cluster = NULL) {
  stopifnot(inherits(fit, "mcglm_graph"))
  sc <- .two_stage_applicable(fit)
  fa <- .two_stage_full_applicable(fit, cluster)
  list(ok = sc$ok, why = sc$why,
       shortcut_ok = sc$ok, why_shortcut = sc$why,
       full_ok = fa$ok, why_full = fa$why,
       method = if (sc$ok || fa$ok)
         .two_stage_auto_method(fit, cluster, sc) else NA_character_)
}
