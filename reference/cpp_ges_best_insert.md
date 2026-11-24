# Best insertion available on a CPDAG (C++ kernel)

Scans the operators Insert(X, Y, T) of the greedy equivalence search: X
and Y nonadjacent, T a subset of the neighbors of Y not adjacent to X,
valid when NA(Y,X) union T is a clique and every semi-directed path from
Y to X meets it. The operator writes X -\> Y and T' -\> Y for T' in T,
and its gain is the difference of the local scores of Y with and without
X. The operator with the largest gain above 1e-10 is returned.

## Usage

``` r
cpp_ges_best_insert(S, n, amat, penalty, max_set = 12L)
```

## Arguments

- S:

  covariance of the responses.

- n:

  number of units.

- amat:

  the current CPDAG.

- penalty:

  multiplier on the parameter count.

- max_set:

  largest T enumerated, as a guard against the exponential.

## Value

list with the gain, the nodes x and y, the set T, and a flag saying
whether any valid insertion was found.
