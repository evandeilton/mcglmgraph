# Closed-form root of the directed layer

Solves the Pearson estimating equations of the modified Cholesky (LDL)
layer of a directed acyclic graph (DAG) directly, from the covariance of
the Pearson residuals. Because the structural matrix has unit
determinant, the Gaussian working criterion of the residual layer
separates by node, \$\$\ell^p_e=-\frac
n2\sum\_{j=1}^R\left\\\log\delta_j+
\frac{q_j(\boldsymbol\gamma_j)}{\delta_j}\right\\,\qquad q_j=\widehat
S\_{e,jj}-2\boldsymbol\gamma_j^\top\widehat{\mathbf S}\_{e,P_j,j}
+\boldsymbol\gamma_j^\top\widehat{\mathbf
S}\_{e,P_jP_j}\boldsymbol\gamma_j,\$\$ so the system is \\R\\
independent least-squares problems and its root is
\$\$\hat{\boldsymbol\gamma}\_j=\widehat{\mathbf S}\_{e,P_jP_j}^{-1}
\widehat{\mathbf S}\_{e,P_j,j},\qquad \hat\delta_j=\widehat
S\_{e,jj}-\widehat{\mathbf S}\_{e,j,P_j} \widehat{\mathbf
S}\_{e,P_jP_j}^{-1}\widehat{\mathbf S}\_{e,P_j,j}.\$\$ For a fixed
directed support the layer is therefore a sequence of least-squares
regressions of the Pearson residuals of each response on those of its
parents. A response with no parent keeps \\\hat\delta_j=\widehat
S\_{e,jj}\\ and contributes no coefficient. Here \\P_j\\ is the set of
parents of response \\j\\, \\\boldsymbol\gamma_j\\ the vector of its
coefficients \\\gamma\_{jk}\\, \\k \in P_j\\, and \\n\\ the number of
units; the criterion is written up to an additive constant.

## Usage

``` r
mc_dag_closed_root(Shat, R, edges = NULL, ordering = NULL, directed = FALSE)
```

## Arguments

- Shat:

  the \\R \times R\\ covariance of the Pearson residuals, in the order
  of the responses.

- R:

  a single positive integer, the number of responses.

- edges:

  optional two-column integer matrix of response indices. With
  `directed = FALSE` each row is an unordered adjacency and its
  direction follows `ordering` (the earlier response is the parent).
  With `directed = TRUE` each row is an arc (parent, child); each pair
  of responses may appear at most once, and self-loops are refused. If
  `NULL`, the complete DAG on `ordering` is used.

- ordering:

  optional integer permutation of `1:R` giving the order of the
  responses in the factorization, parents before children, so that
  `ordering[1]` has no parent; `NULL` uses `1:R`. With `directed = TRUE`
  it is derived from the arcs when left `NULL`, and checked against them
  when supplied.

- directed:

  logical; `TRUE` reads the columns of `edges` as parent and child
  rather than as an unordered pair.

## Value

A list with `gamma`, `delta`, the implied `Sigma_b` and `K_b` in the
order of the responses, the `edges` actually used as (parent, child)
pairs in response labels with `gamma` aligned to their rows, the
`ordering`, the `parents` of each response, and `status`, which records
how many nodes needed a ridge on their normal equations, how many
innovation variances were lifted to the positivity floor
(\\10^{-10}\widehat S\_{e,jj}\\), and `ill_conditioned`, the responses
whose parent block, rescaled to a unit diagonal, has a reciprocal
condition number at or below `.Machine$double.eps^(2/3)`; their
coefficients are not identified and come from a solve with a ridge of
\\10^{-8}\\ on that unit diagonal. Both tests are invariant to a
rescaling of the responses, and so is the root: rescaling response \\j\\
by \\a_j\\ multiplies \\\gamma\_{jk}\\ by \\a_j/a_k\\ and \\\delta_j\\
by \\a_j^2\\. A nonzero status means the residual covariance is close to
singular on some parent set, and
[`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
then does not certify the fit for inference. `delta` is indexed by
response, like the diagonal of `Sigma_b`.

## Relation to the iterative engine

[`mcglm_graph`](https://evandeilton.github.io/mcglmgraph/reference/mcglm_graph.md)
takes this route by default. Passing
`control = list(dag_solver = "newton")` runs instead the scoring
iteration on \\(\boldsymbol\gamma,\log\boldsymbol\delta)\\ that the
directed and undirected layers share. It solves the same system and
serves as an independent check on the closed form.

## Giving the arcs directly

With `directed = FALSE` the edges are unordered adjacencies and their
direction comes from `ordering`, as elsewhere in the package. That is
awkward when the support to fit is an estimated DAG, whose arcs need not
run along the natural order of the responses. With `directed = TRUE` the
two columns of `edges` are read as parent and child, any acyclic set of
arcs is accepted, and a compatible ordering is derived internally by
Kahn's algorithm. A cyclic set is refused, since it has no modified
Cholesky factorization at all.

## References

Pourahmadi, M. (1999). Joint mean-covariance models with applications to
longitudinal data: unconstrained parameterisation. *Biometrika* 86(3),
677-690.
[doi:10.1093/biomet/86.3.677](https://doi.org/10.1093/biomet/86.3.677)

## See also

[`mc_dag_sensitivity`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_coefficients`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md)

Other directed layer:
[`mc_build_D_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_D_matrix.md),
[`mc_build_T_matrix()`](https://evandeilton.github.io/mcglmgraph/reference/mc_build_T_matrix.md),
[`mc_cholesky_derivatives()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_derivatives.md),
[`mc_cholesky_to_covariance()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_covariance.md),
[`mc_cholesky_to_precision()`](https://evandeilton.github.io/mcglmgraph/reference/mc_cholesky_to_precision.md),
[`mc_dag_coefficients()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
[`mc_dag_sensitivity()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_sensitivity.md),
[`mc_dag_unit_diagonal()`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_unit_diagonal.md)

## Examples

``` r
S <- mc_cholesky_to_covariance(c(0.7, 0.3, 0.5), c(1, 0.64, 0.49),
                               R = 3)$Sigma_b
fit <- mc_dag_closed_root(S, R = 3)
round(fit$gamma, 12)     # recovers the generating coefficients
#> [1] 0.7 0.3 0.5
round(fit$delta, 12)
#> [1] 1.00 0.64 0.49
```
