# Local score of a node given a parent set (C++ kernel)

Returns -(n/2) log(r) - (penalty/2) (\|pa\| + 1), where r = S\[j, j\] -
S\[j, pa\] S\[pa, pa\]^-1 S\[pa, j\] is the residual variance of the
regression of node j on its parents. It is -Inf when the parent block is
numerically singular or r is not positive.

## Usage

``` r
cpp_ges_local_score(S, n, j, pa, penalty)
```

## Arguments

- S:

  covariance of the responses.

- n:

  number of units.

- j:

  node, 1-based.

- pa:

  parents, 1-based, possibly empty.

- penalty:

  multiplier on the parameter count; log(n) gives the BIC.

## Value

the penalized Gaussian working contribution of the node.
