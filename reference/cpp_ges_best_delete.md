# Best deletion available on a CPDAG (C++ kernel)

Scans the operators Delete(X, Y, H): X and Y adjacent by X - Y or X -\>
Y, H a subset of NA(Y,X), valid when NA(Y,X) minus H is a clique. The
operator removes the edge and writes Y -\> H' and X -\> H' for H' in H.
The operator with the largest gain above 1e-10 is returned.

## Usage

``` r
cpp_ges_best_delete(S, n, amat, penalty, max_set = 12L)
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

list with the gain, the nodes x and y, the set H, and a flag.
