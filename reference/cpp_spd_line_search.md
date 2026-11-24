# Backtracking line search preserving positive-definiteness (C++ kernel)

Largest step alpha among 1, 1/2, ..., 2^-max_halvings such that
K_b(eta + alpha\*d) has smallest eigenvalue above tol. When no
admissible step is found, is_spd is false, eta_new is the point tried
with the smallest step 2^-max_halvings (it fails the test), halvings is
max_halvings + 1 and step_size is 2^-(max_halvings + 1), half the step
that produced eta_new.

## Usage

``` r
cpp_spd_line_search(eta, d, bases, max_halvings = 30L, tol = 1e-10)
```

## Arguments

- eta:

  current parameter vector.

- d:

  Newton direction.

- bases:

  numeric cube R x R x q_G of base matrices Z_d, with q_G the number of
  graph parameters.

- max_halvings:

  maximum number of step halvings (default 30).

- tol:

  tolerance on the smallest eigenvalue (default 1e-10).

## Value

list with eta_new (last point tried), step_size, halvings (number of
halvings performed) and is_spd (whether eta_new passed the test).
