# Members of a Markov equivalence class, and what they disagree about.

vkeys <- function(A) {
  storage.mode(A) <- "integer"
  v <- mcglmgraph:::cpp_vstructure_set(A)
  if (!nrow(v)) character(0) else
    sort(unique(paste(v[, 1], v[, 2], v[, 3], sep = "->")))
}

has_cycle <- function(A) {
  R <- nrow(A); indeg <- colSums(A); live <- rep(TRUE, R)
  repeat {
    src <- which(live & indeg == 0)
    if (!length(src)) break
    s <- src[1]; live[s] <- FALSE
    indeg[A[s, ] == 1] <- indeg[A[s, ] == 1] - 1
    indeg[s] <- -1
  }
  any(live)
}

test_that("the count of members is the one the class has", {
  # a chain has three orientations without a collider
  chain <- matrix(0L, 3, 3)
  chain[1, 2] <- chain[2, 1] <- chain[2, 3] <- chain[3, 2] <- 1L
  expect_equal(length(mc_class_members(chain)), 3L)

  # a collider is its own class
  coll <- matrix(0L, 3, 3)
  coll[1, 3] <- coll[2, 3] <- 1L
  expect_equal(length(mc_class_members(coll)), 1L)

  # a complete graph on R nodes has one member per topological order
  for (R in 3:4) {
    full <- matrix(1L, R, R); diag(full) <- 0L
    expect_equal(length(mc_class_members(full)), factorial(R))
  }

  # a fully directed CPDAG has exactly itself
  expect_identical(mc_class_members(coll)[[1L]], coll)
})

test_that("every member is an acyclic graph of the same class", {
  set.seed(3)
  for (i in 1:25) {
    R <- sample(3:5, 1)
    D <- matrix(0L, R, R)
    for (a in seq_len(R - 1L)) for (b in (a + 1L):R)
      if (runif(1) < 0.5) D[a, b] <- 1L
    ord <- sample(R); D <- D[ord, ord]
    cp <- mc_true_cpdag(t(D))        # mc_true_cpdag reads arcs from the lower triangle
    mem <- mc_class_members(cp)
    if (!length(mem)) next
    target <- vkeys(cp)
    for (A in mem) {
      expect_false(has_cycle(A))
      expect_true(all((A == 1L) + (t(A) == 1L) <= 1L))    # no undirected edge
      expect_identical(vkeys(A), target)
      expect_identical(unname((A + t(A)) > 0), unname((cp + t(cp)) > 0))
    }
    expect_equal(anyDuplicated(vapply(mem, function(A)
      paste(A, collapse = ""), character(1))), 0L)        # no repetition
  }
})

test_that("the members agree on the fitted covariance and differ on coefficients", {
  # a five-response counterexample:
  # a -> k <- b, k -> j, k -> m, and an edge between j and m
  g <- 0.7; h <- 0.5; c0 <- 0.4
  arcs <- rbind(c(1, 3), c(2, 3), c(3, 4), c(3, 5), c(4, 5))
  gam <- c(0.6, 0.6, g, h, c0)
  T <- diag(5)
  for (e in seq_len(nrow(arcs))) T[arcs[e, 2], arcs[e, 1]] <- -gam[e]
  Ti <- solve(T); Sig <- Ti %*% diag(5) %*% t(Ti)

  cp <- matrix(0L, 5, 5)
  cp[1, 3] <- cp[2, 3] <- 1L
  cp[3, 4] <- cp[3, 5] <- 1L
  cp[4, 5] <- cp[5, 4] <- 1L
  mem <- mc_class_members(cp)
  expect_equal(length(mem), 2L)

  fits <- lapply(mem, function(A) {
    a <- which(A == 1L & t(A) == 0L, arr.ind = TRUE)
    mc_dag_closed_root(Sig, R = 5, edges = matrix(as.integer(a), ncol = 2L),
                       directed = TRUE)
  })
  # every member reproduces the same covariance: the class is what is identified
  expect_equal(fits[[1]]$Sigma_b, fits[[2]]$Sigma_b, tolerance = 1e-10)
  expect_equal(fits[[1]]$Sigma_b, Sig, tolerance = 1e-10)

  val <- function(f, from, to) {
    k <- paste(f$edges[, 1], f$edges[, 2], sep = "->")
    unname(f$gamma[k == paste(from, to, sep = "->")])
  }
  kj <- sort(c(val(fits[[1]], 3, 4), val(fits[[2]], 3, 4)))
  km <- sort(c(val(fits[[1]], 3, 5), val(fits[[2]], 3, 5)))
  # the compelled arc k -> j takes two values, and the closed form gives one
  expect_equal(kj, c((g - c0 * h) / (1 + c0^2), g), tolerance = 1e-9)
  expect_equal(round(kj, 3), c(0.431, 0.700))
  expect_equal(round(km, 3), c(0.500, 0.780))
  # the coefficients into k, whose head no undirected edge touches, do not move
  expect_equal(val(fits[[1]], 1, 3), val(fits[[2]], 1, 3), tolerance = 1e-10)
  expect_equal(val(fits[[1]], 2, 3), val(fits[[2]], 2, 3), tolerance = 1e-10)
})

test_that("the table marks the arcs the class pins down", {
  skip_if_not_installed("MASS")
  set.seed(9); n <- 600; x <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.8 * z1 + 0.8 * z2 + rnorm(n)
  z4 <- 0.6 * z3 + rnorm(n)
  d <- data.frame(y1 = z1, y2 = z2, y3 = z3, y4 = z4, x = x)
  f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x, y4 ~ x), data = d)
  cc <- mc_class_coefficients(f)
  expect_s3_class(cc, "mcg_class_coef")
  expect_gte(cc$n_members, 1L)
  expect_true(all(cc$arcs$gamma_min <= cc$arcs$gamma_max))
  # an arc flagged common has a single value across members
  for (k in which(cc$arcs$common))
    expect_equal(cc$arcs$gamma_min[k], cc$arcs$gamma_max[k], tolerance = 1e-9)
  # and every member fits the same covariance
  Se <- f$graph$Shat
  sig <- lapply(cc$members, function(A) {
    a <- which(A == 1L & t(A) == 0L, arr.ind = TRUE)
    mc_dag_closed_root(Se, R = 4, edges = matrix(as.integer(a), ncol = 2L),
                       directed = TRUE)$Sigma_b
  })
  for (S in sig[-1]) expect_equal(S, sig[[1]], tolerance = 1e-9)
  txt <- capture.output(print(cc))
  expect_true(any(grepl("member", txt)))
})

test_that("the enumeration can be capped", {
  full <- matrix(1L, 5, 5); diag(full) <- 0L
  mem <- mc_class_members(full, max_members = 10L)
  expect_lte(length(mem), 10L)
  expect_true(attr(mem, "truncated"))
  small <- mc_class_members(full, max_members = 1000L)
  expect_equal(length(small), factorial(5))
  expect_false(attr(small, "truncated"))
})

test_that("the contracts of the arguments are enforced", {
  a <- matrix(0L, 3, 3)
  expect_error(mc_class_members(matrix(2L, 3, 3)), "0/1 adjacency")
  expect_error(mc_class_members(matrix(0L, 2, 3)), "square")
  expect_error(mc_class_members(a, max_members = 0), "positive integer")
  expect_error(mc_class_coefficients(structure(list(), class = "lm")),
               "inherits")
})
