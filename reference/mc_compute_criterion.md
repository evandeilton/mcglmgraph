# Information criteria for a graph-structured fit

Computes working AIC and BIC, the quasi-likelihood information criterion
(QIC), a Godambe-adjusted working BIC (GBIC) and the extended BIC
(EBIC), all built on the Gaussian working log-criterion \\\ell^p_e\\ of
[`mc_graph_ploglik`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md).
With \\p\\ the number of free parameters of the between-response
covariance (diagonal and edge parameters) and \\n\\ the number of units
(rows of the residual matrix), \$\$\mathrm{AIC} = -2\ell^p_e + 2p,\qquad
\mathrm{BIC} = -2\ell^p_e + \log(n)\\p,\$\$ \$\$\mathrm{QIC} =
-2\ell^p_e + 2\\d\_{\mathrm{eff}},\qquad \mathrm{GBIC} = -2\ell^p_e +
\log(n)\\d\_{\mathrm{eff}},\$\$ with the effective dimension
\\d\_{\mathrm{eff}}=\mathrm{tr}(H^{-1}V)\\. Here \\H\\ is the curvature
(sensitivity) and \\V\\ the empirical variability of the Pearson
estimating functions, both per unit, so that the trace is invariant to
the scale. Under the Gaussian working model, when \\H = V = J\\, the
effective dimension reduces to the parameter count and QIC and GBIC
reduce to AIC and BIC. `"clbic"` is an exact alias of `"gbic"`; it does
not assert a composite-likelihood factorization. Smaller values indicate
a better working model.

## Usage

``` r
mc_compute_criterion(
  fit,
  which = c("aic", "bic", "qic", "gbic", "clbic", "ebic"),
  ebic_gamma = 0.5,
  cluster = NULL,
  small_sample = FALSE
)
```

## Arguments

- fit:

  an object of class `"mcglm_graph"`.

- which:

  character vector among `"aic"`, `"bic"`, `"qic"`, `"gbic"`, `"clbic"`
  and `"ebic"`. `"clbic"` is an alias of `"gbic"`. When left out, the
  first five are returned; `"ebic"` is computed only when requested.

- ebic_gamma:

  the parameter \\\gamma\ge0\\ of the extended BIC (default 0.5), used
  only when `"ebic"` is requested.

- cluster:

  optional independent-cluster identifier, one per residual row (at
  least two distinct values, no missing values). For QIC and GBIC,
  Pearson contributions are summed within clusters before forming their
  empirical variability. BIC, GBIC and EBIC use the logarithm of the
  number of independent clusters, instead of the number of rows, in
  their penalty. Supplying the identifier does not by itself establish
  independence of the clusters or centering of their scores.

- small_sample:

  logical; multiply the clustered variability by \\G/(G-1)\\, with \\G\\
  the number of clusters; ignored without `cluster`.

## Value

A named numeric vector of the requested criteria, carrying the
attributes `"criterion_status"` (`"godambe"`, `"godambe_pseudoinverse"`,
`"parameter_count_fallback"` or `"not_required"`) and
`"effective_dimension"` (the value of \\d\_{\mathrm{eff}}\\, `NA` when
no Godambe-adjusted criterion was requested). A status other than
`"godambe"` means the effective dimension could not be computed from a
well-conditioned sensitivity (a spectral pseudo-inverse or, last, the
parameter count was used instead) and must be reported as such.

## Details

With `cluster`, the logarithmic multiplier \\\log(n)\\ of BIC, GBIC and
EBIC becomes \\\log(G)\\, with \\G\\ the number of independent clusters,
and the variability \\V\\ is formed from the contributions summed within
clusters. The criterion \\\ell^p_e\\, which remains a sum over the \\n\\
rows, the parameter count \\p\\, AIC and the multiplier 2 of QIC are not
altered.

## The extended BIC

`"ebic"` is the form of Foygel and Drton (2010) for Gaussian graphical
models, \\\mathrm{EBIC}\_\gamma=\mathrm{BIC}+4\gamma\|E\|\log R\\, with
\\\|E\|\\ the number of edge parameters of the fit and \\R\\ the number
of responses. The diagonal parameters are common to every candidate
support, so they enter through the BIC alone and do not alter
differences between supports. \\\gamma=0\\ gives the BIC back. It is
offered as a working criterion: its consistency for the graph is
established for Gaussian likelihoods, and it does not carry over, as
such, to the working criterion of the residual layer.

## Scope of the Godambe penalties

For a sparse DAG, \\H\\ is the observed curvature of the residual
working criterion. It can differ from the model-based Pearson
sensitivity when the fitted covariance is misspecified. For the
linear-precision GGM they agree. The QIC optimism calculation assumes
genuinely known first-stage marginal parameters. With marginals
estimated on the same data, QIC remains a working criterion; the
reported penalty is not an exact two-stage optimism correction. The
empirical variability is \\V=n^{-1}\sum_i\psi_i\psi_i^\top\\, with
\\\psi_i\\ the Pearson contribution of unit \\i\\ (summed within
clusters when `cluster` is given, then divided by \\n\\). It estimates
the variance of the scores only when their means vanish per independent
unit or cluster (or the aggregate contribution of squared means
vanishes). A pooled zero score is insufficient.

## References

Foygel, R. and Drton, M. (2010). Extended Bayesian information criteria
for Gaussian graphical models. *Advances in Neural Information
Processing Systems* 23, 604–612.

Pan, W. (2001). Akaike's information criterion in generalized estimating
equations. *Biometrics* 57(1), 120–125.
[doi:10.1111/j.0006-341X.2001.00120.x](https://doi.org/10.1111/j.0006-341X.2001.00120.x)

Gao, X. and Song, P. X.-K. (2010). Composite likelihood Bayesian
information criteria for model selection in high-dimensional data.
*Journal of the American Statistical Association* 105(492), 1531–1540.
[doi:10.1198/jasa.2010.tm09414](https://doi.org/10.1198/jasa.2010.tm09414)

Varin, C. and Vidoni, P. (2005). A note on composite likelihood
inference and model selection. *Biometrika* 92(3), 519–528.
[doi:10.1093/biomet/92.3.519](https://doi.org/10.1093/biomet/92.3.519)

## See also

[`mc_graph_ploglik`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md),
[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md),
[`mc_graph_cv`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md).

Other information criteria:
[`mc_graph_cv()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_cv.md),
[`mc_graph_ploglik()`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_ploglik.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 150; x <- rnorm(n)
d <- data.frame(y1 = rnorm(n), y2 = rnorm(n), y3 = rnorm(n), x = x)
f <- mcglm_graph(list(y1~x, y2~x, y3~x), data = d)
#> Automatic initial values selected.
mc_compute_criterion(f, c("aic", "bic", "qic", "gbic"))
#>      aic      bic      qic     gbic 
#> 1282.718 1300.782 1283.775 1303.430 
#> attr(,"criterion_status")
#> [1] "godambe"
#> attr(,"effective_dimension")
#> [1] 6.528506
mc_compute_criterion(f, "ebic", ebic_gamma = 0.5)
#>     ebic 
#> 1307.374 
#> attr(,"criterion_status")
#> [1] "not_required"
#> attr(,"effective_dimension")
#> [1] NA
# }
```
