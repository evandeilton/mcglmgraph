# Methods for graph structure selection results

Print, summarize and plot the result of
[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md)
or
[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md).
`print` shows the selector, the criterion, the number and list of
selected edges and, for the L1 route, the selected penalty. `summary`
adds the search trace (stepwise), the head of the penalty path (L1), and
the v-structures and CPDAG when the object carries them. `plot` draws
the network of the selected fit with
[`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md).

## Usage

``` r
# S3 method for class 'mcg_selection'
print(x, ...)

# S3 method for class 'mcg_selection'
summary(object, ...)

# S3 method for class 'mcg_selection'
plot(x, ...)
```

## Arguments

- x, object:

  an object of class `"mcg_selection"`.

- ...:

  for `plot`, further arguments passed to
  [`mc_plot_network`](https://evandeilton.github.io/mcglmgraph/reference/mc_plot_network.md);
  ignored by `print` and `summary`.

## Value

`print` and `summary` return the object invisibly; `plot` returns
invisibly the value of `mc_plot_network` applied to the selected fit (or
the object itself when it carries no fit).

## See also

[`mc_graph_stepwise`](https://evandeilton.github.io/mcglmgraph/reference/mc_graph_stepwise.md),
[`mc_select_lambda`](https://evandeilton.github.io/mcglmgraph/reference/mc_select_lambda.md).

## Examples

``` r
# \donttest{
set.seed(1); n <- 200; x <- rnorm(n)
y1 <- rnorm(n); y2 <- 0.6 * y1 + rnorm(n); y3 <- rnorm(n)
d <- data.frame(y1 = y1, y2 = y2, y3 = y3, x = x)
f <- mcglm_graph(list(y1 ~ x, y2 ~ x, y3 ~ x), data = d)
#> Automatic initial values selected.
sel <- mc_graph_stepwise(f)
print(sel)
#> Graph structure selection (stepwise, criterion = bic)
#> Selected edges: 1
#> Edges:
#>      [,1] [,2]
#> [1,]    1    2
summary(sel)
#> Graph structure selection (stepwise, criterion = bic)
#> Selected edges: 1
#> Edges:
#>      [,1] [,2]
#> [1,]    1    2
#> 
#> Search trace:
#>   step action edge_i edge_j criterion
#> 1    1 remove      1      3  1671.188
#> 2    2 remove      2      3  1666.173
# }
```
