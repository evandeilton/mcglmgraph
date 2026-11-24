# Residual degrees of freedom of the first stage

Returns \\\nu_e\\, the residual degrees of freedom that the Fisher
factor of the conditional-association test uses in the form
\\\nu_e-\|S\|-2\\, when the first stage fixes it exactly, and says why
it does not otherwise.

## Usage

``` r
mc_residual_df(fit)

# S3 method for class 'mcg_residual_df'
print(x, ...)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"` or `"mcglm"`.

- x:

  an object of class `"mcg_residual_df"`.

- ...:

  ignored.

## Value

An object of class `"mcg_residual_df"`: a list with `df`, the value of
\\\nu_e\\ or `NA`; `exact`, a logical; `n`, the number of observations;
`rank`, the common rank \\h\\ of the designs or `NA`; and `why`, a
character vector with the reasons the correction is not exact, empty
when it is. The `print` method returns `x` invisibly.

## When the correction is exact

The Pearson residuals of a Gaussian response with identity link and
constant variance are the least-squares residuals of its design, scaled.
When every response has that form and the same column space of design,
of rank \\h\\, the residual covariance is the sample covariance after
projecting out a common \\h\\-dimensional space. For Gaussian responses
its partial correlations are then distributed as those of a sample
covariance with \\n-h\\ degrees of freedom, that of \\n-h+1\\
observations with an estimated mean, so \\\nu_e=n-h\\. With the
intercept alone, \\h=1\\ and \\\nu_e-\|S\|-2=n-\|S\|-3\\, the factor
used when no correction is supplied. Outside that case no finite-sample
correction is established and none is offered: the function returns
`exact = FALSE` with the reasons, and the approximation by \\n\\
(`df_residual = NULL` in the separation routines) is the one to use. The
check requires constant variance, the identity link, an identity
within-response structure, unit first-stage weights, and designs of
equal rank spanning the same column space for all responses.

## See also

[`mc_pc_skeleton`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
whose `df_residual` accepts it, and
[`mc_graph_pipeline`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
where `df_residual = "auto"` calls it.

Other structure discovery:
[`mc_find_vstructures()`](https://evandeilton.github.io/mcglmgraph/reference/mc_find_vstructures.md),
[`mc_ges()`](https://evandeilton.github.io/mcglmgraph/reference/mc_ges.md),
[`mc_graph_pipeline()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_pipeline.md),
[`mc_orient_meek()`](https://evandeilton.github.io/mcglmgraph/reference/mc_orient_meek.md),
[`mc_pc_skeleton()`](https://evandeilton.github.io/mcglmgraph/reference/mc_pc_skeleton.md),
[`mc_separating_sets()`](https://evandeilton.github.io/mcglmgraph/reference/mc_separating_sets.md),
[`mc_skeleton_to_bnlearn()`](https://evandeilton.github.io/mcglmgraph/reference/mc_skeleton_to_bnlearn.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 60; x <- rnorm(n); g <- gl(3, 20)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x, g = g)
f <- mcglm_graph(list(y1 ~ x + g, y2 ~ x + g, y3 ~ x + g), data = d)
#> Automatic initial values selected.
mc_residual_df(f)          # 60 - 4
#> Residual degrees of freedom: 56 (n = 60, rank of the design = 4)
# }
```
