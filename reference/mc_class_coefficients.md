# What the members of a class disagree about

Fits the directed layer on the residual covariance of a fit under every
member of a Markov equivalence class and collects, for each arc, the
values of its structural coefficient that the class allows. The
coefficients are those of the closed-form root,
[`mc_dag_closed_root`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md),
and the standardized ones are
\\\gamma^\*\_{jk}=\gamma\_{jk}(\Sigma\_{b,kk}/ \Sigma\_{b,jj})^{1/2}\\,
as in
[`mc_dag_coefficients`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_coefficients.md),
with \\j\\ the child and \\k\\ the parent. An arc whose values coincide
across members is one whose coefficient the class determines; an arc
whose values move is one for which a single number would report a choice
the data did not make.

## Usage

``` r
mc_class_coefficients(fit, cpdag = NULL, max_members = 1000L, ...)

# S3 method for class 'mcg_class_coef'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`, whose residual covariance the
  layer is fitted to.

- cpdag:

  the class, as an adjacency matrix or an `"mcg_cpdag"` object; `NULL`
  runs
  [`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  on the fit.

- max_members:

  cap on the number of members enumerated, as in
  [`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md).

- ...:

  passed to
  [`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md)
  when `cpdag` is `NULL`; ignored by the `print` method.

- x:

  an object of class `"mcg_class_coef"`.

## Value

An object of class `"mcg_class_coef"`: a list with `arcs`, a data frame
with one row per arc present in some member (`from`, `to`, `parent`,
`child`, the number of members that contain it `n_members`, the range of
its coefficient `gamma_min`, `gamma_max` and of its standardized
coefficient `gamma_std_min`, `gamma_std_max`, and the logical `common`,
true when no undirected edge of the class touches the head of the arc);
`members`, the list of DAGs; `values`, the per-member table of
coefficients (`member`, `from`, `to`, `gamma`, `gamma_std`);
`n_members`; `truncated`, whether the enumeration was capped; and
`response_names`. The `print` method returns `x` invisibly. An error is
raised when the class has no consistent extension.

## What the values are

Every member reproduces the same fitted covariance when the class is a
CPDAG, and the second moments do not distinguish among them. The
direction of a compelled arc is a property of the class, but its
coefficient is a partial regression coefficient given the other parents
of its head, and those can change from one member to another. The values
are coefficients of linear projections of Pearson residuals. Enumerating
a quantity over the members of a class follows Maathuis, Kalisch and
Buehlmann (2009), who apply it to intervention effects; the coefficients
reported here admit that reading only under causal assumptions the
procedure does not check: a linear structural model with independent
innovations whose graph belongs to the class, causal sufficiency and
absence of selection bias.

## Which arcs the class pins down

An arc has the same coefficient in every member whose parent set for its
head is the same. That happens for every member when no undirected edge
of the CPDAG touches the head of the arc, which is the condition to
check before reporting a single value. The returned table flags it in
`common`.

## References

Maathuis, M. H., Kalisch, M. and Buehlmann, P. (2009). Estimating
high-dimensional intervention effects from observational data. *The
Annals of Statistics* 37(6A), 3133-3164.
[doi:10.1214/09-AOS685](https://doi.org/10.1214/09-AOS685)

## See also

[`mc_class_members`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md),
[`mc_dag_closed_root`](https://evandeilton.github.io/mcglmgraph/reference/mc_dag_closed_root.md)

Other equivalence class:
[`mc_class_members()`](https://evandeilton.github.io/mcglmgraph/reference/mc_class_members.md),
[`mc_sid()`](https://evandeilton.github.io/mcglmgraph/reference/mc_sid.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 400; x <- rnorm(n)
z1 <- rnorm(n); z2 <- rnorm(n); z3 <- 0.7 * z1 + 0.7 * z2 + rnorm(n)
d <- data.frame(y1 = z1, y2 = z2, y3 = z3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
mc_class_coefficients(f)
#> Coefficients over the 1 member of the class
#> An arc marked common has the same value in every member.
#> 
#>   parent child n_members gamma_min gamma_max gamma_std_min gamma_std_max common
#> 1     y1    y3         1     0.559     0.559         0.560         0.560   TRUE
#> 2     y2    y3         1     0.496     0.496         0.497         0.497   TRUE
# }
```
