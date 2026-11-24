# Fits for the tests of the route taken by the two-stage correction. Three
# responses with Gaussian-copula dependence, a log-linear mean in one covariate
# and the variance functions of `kind`.
#   "gaussian"   constant variance, identity link
#   "tweedie"    gamma margins, Tweedie variance of power `power` (fixed)
#   "mixed"      first response Gaussian (constant), the others Tweedie
#   "binomialP"  binomial margins with 10 trials, logit link
#   "tweedie0"   Gaussian data under a Tweedie variance of power fixed at 0
route_fit <- function(kind, n = 250, seed = 1, power = 2,
                      between_graph = mc_between_graph("ggm")) {
  set.seed(seed); x <- stats::runif(n)
  Rl <- matrix(c(1, .5, .2, .5, 1, .4, .2, .4, 1), 3)
  U <- stats::pnorm(MASS::mvrnorm(n, rep(0, 3), Rl))
  qg <- function(r, mu) stats::qgamma(U[, r], shape = 2, rate = 2 / mu)
  mu <- cbind(exp(1 + .5 * x), exp(.5 - .3 * x), rep(exp(.8), n))
  lp <- list(y1 ~ x, y2 ~ x, y3 ~ x)
  ini <- function(reg, pw, tau)
    list(regression = reg, power = as.list(pw), tau = as.list(tau),
         rho = c(.5, .2, .4))
  reg <- list(c(1, .5), c(.5, -.3), c(.8, 0))
  if (kind == "tweedie") {
    Y <- sapply(1:3, function(r) qg(r, mu[, r]))
    args <- list(link = rep("log", 3), variance = rep("tweedie", 3),
                 power_fixed = rep(TRUE, 3),
                 control_initial = ini(reg, rep(power, 3), rep(.5, 3)))
  } else if (kind == "count") {
    # overdispersed counts for the Poisson-Tweedie variance with power one
    Y <- sapply(1:3, function(r) stats::qnbinom(U[, r], size = 3, mu = mu[, r]))
    args <- list(link = rep("log", 3), variance = rep("poisson_tweedie", 3),
                 power_fixed = rep(TRUE, 3),
                 control_initial = ini(reg, rep(1, 3), rep(.3, 3)))
  } else if (kind == "mixed") {
    Y <- cbind(stats::qnorm(U[, 1], 1 + x, 1), qg(2, mu[, 2]), qg(3, mu[, 3]))
    args <- list(link = c("identity", "log", "log"),
                 variance = c("constant", "tweedie", "tweedie"),
                 power_fixed = rep(TRUE, 3),
                 control_initial = ini(list(c(1, 1), reg[[2]], reg[[3]]),
                                       c(1, 2, 2), rep(.5, 3)))
  } else if (kind == "binomialP") {
    p <- stats::plogis(cbind(-.5 + x, .3 - x, rep(0, n)))
    # mcglm takes the proportion of successes with the number of trials
    Y <- sapply(1:3, function(r) stats::qbinom(U[, r], 10, p[, r])) / 10
    args <- list(link = rep("logit", 3), variance = rep("binomialP", 3),
                 Ntrial = rep(list(rep(10, n)), 3))
  } else if (kind == "tweedie0") {
    Y <- cbind(1 + x, 2 - x, rep(.5, n)) + MASS::mvrnorm(n, rep(0, 3), Rl)
    args <- list(link = rep("identity", 3), variance = rep("tweedie", 3),
                 power_fixed = rep(TRUE, 3),
                 control_initial = ini(list(c(1, 1), c(2, -1), c(.5, 0)),
                                       rep(0, 3), rep(1, 3)))
  } else {
    Y <- cbind(1 + x, 2 - x, rep(.5, n)) + MASS::mvrnorm(n, rep(0, 3), Rl)
    args <- list()
  }
  d <- data.frame(y1 = Y[, 1], y2 = Y[, 2], y3 = Y[, 3], x = x)
  suppressMessages(do.call(mcglm_graph,
                           c(list(lp, between_graph = between_graph,
                                  data = d), args)))
}
