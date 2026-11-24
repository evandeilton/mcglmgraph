# Innovation variances that force a unit diagonal

Returns the innovation variances under which the covariance of the
directed layer has unit diagonal, that is under which it is a
correlation matrix, and says whether the coefficients supplied admit
them.

## Usage

``` r
mc_dag_unit_diagonal(gamma, R, edges = NULL)
```

## Arguments

- gamma:

  finite numeric vector of structural coefficients, exactly one per
  edge.

- R:

  a single positive integer, the number of responses.

- edges:

  optional integer matrix with columns \\(i, j)\\, \\i \< j\\, the
  directed edges \\i \to j\\. If `NULL`, the complete DAG on the natural
  ordering \\1 \< 2 \< \dots \< R\\ is used. Each edge must occur
  exactly once; duplicated rows, including reversed duplicates, are
  rejected because they do not define distinct structural coordinates. A
  row given as \\(j, i)\\ with \\j \> i\\ is read as the same edge \\i
  \to j\\: the direction always follows the natural order of the
  responses. The rows are sorted lexicographically before use, and
  `gamma` is read in that sorted order, not in the order the rows were
  supplied. Sort the matrix yourself, or read the order back from the
  `edges` component of the result, whenever the correspondence matters.

## Value

A list with `delta`, the innovation variances implied by the
coefficients; `admissible`, true when all of them are positive;
`min_delta`; `Sigma_b`, the covariance they produce, which has unit
diagonal when admissible; `r_squared`, the variance of each response
accounted for by its parents; and `to_unit_innovations`, the diagonal
rescaling that sends the covariance to the representative of its orbit
with unit innovations.

## Why this is a constraint and not a reparametrization

The undirected layer can fix the diagonal at one at the cost of nothing.
Here the requirement
\$\$\delta_j=1-\boldsymbol\gamma_j^\top\boldsymbol\Sigma\_{b,P_jP_j}
\boldsymbol\gamma_j,\qquad\text{equivalently}\qquad
\delta_1=1,\quad\delta_j=1-\sum\_{k\<j}L\_{jk}^2\delta_k,\$\$ with
\\\mathbf L=\mathbf T^{-1}\\, leaves \\\|E_D\|\\ free parameters but
confines them to the open set where every \\\delta_j\\ stays positive, a
set defined by \\R\\ nonlinear inequalities. It contains the origin and
is a proper subset of the coefficient space whenever there is an edge.
That is exactly the freedom the layer buys by leaving the diagonal free,
which is why the fitting routines do not impose it: the constraint is
useful for reading a covariance as a correlation and for reproducing a
generator, not for estimating.

## One representative per orbit

Rescaling by \\\mathbf A\\ sends \\\delta_j\\ to \\a_j^2\delta_j\\, so
every orbit \\\\\mathbf A\boldsymbol\Sigma_b\mathbf A\\\\ contains
exactly one matrix with unit innovation variances, reached by
\\a_j=\delta_j^{-1/2}\\. The component `to_unit_innovations` carries
that rescaling.

## See also

[`mc_dag_coefficients`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)
for the scale-free quantities.

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_closed_root()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md)

## Examples

``` r
# a chain: the second and third responses keep less and less of their own
mc_dag_unit_diagonal(c(0.6, 0.6), R = 3, edges = rbind(c(1, 2), c(2, 3)))$delta
#> [1] 1.00 0.64 0.64
# coefficients too large for a correlation matrix are flagged
mc_dag_unit_diagonal(c(0.9, 0.9), R = 3,
                     edges = rbind(c(1, 3), c(2, 3)))$admissible
#> [1] FALSE
```
