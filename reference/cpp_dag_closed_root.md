# Closed-form root of the directed layer (C++ kernel)

Solves the Pearson estimating equations of the modified Cholesky layer
node by node: gamma_j = S\[P_j, P_j\]^-1 S\[P_j, j\] and delta_j = S\[j,
j\] - S\[j, P_j\] gamma_j. A node with no parent keeps delta_j = S_jj
and contributes no coefficient. Every delta_j is lifted to floor_rel \*
S_jj when it falls below that floor.

## Usage

``` r
cpp_dag_closed_root(
  Shat,
  R,
  edge_idx,
  floor_rel = 1e-10,
  rcond_tol = 3.66685286250104e-11
)
```

## Arguments

- Shat:

  R x R sample covariance of the Pearson residuals, in the positional
  (elimination) order of the nodes.

- R:

  integer number of responses.

- edge_idx:

  integer \|E\| x 2 matrix of edges (i, j), i \< j, 1-based, in
  positional order.

- floor_rel:

  positivity floor of delta_j, relative to S_jj (to max \|S_kk\| when
  S_jj is not positive), so that it follows the scale of each response.

- rcond_tol:

  a parent block S\[P_j, P_j\] whose equilibrated form D^-1 S\[P_j,
  P_j\] D^-1, D = diag(S\[P_j, P_j\])^(1/2), has a reciprocal condition
  number not above this tolerance is near singular: its coefficients are
  not identified and are taken from a solve with a ridge of 1e-8 on the
  unit diagonal of the equilibrated block.

## Value

list with gamma (one per edge, in the order of edge_idx), delta (length
R), the number of nodes whose normal equations needed a ridge, the
number of innovation variances lifted to the positivity floor, and
ill_conditioned, the positional (1-based) nodes whose parent block is
near singular.
