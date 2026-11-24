# Delta-method variance of the between-response partial correlations

For the pair (i, j) the partial correlation is rho_ij = -K_ij /
sqrt(K_ii K_jj). Differentiating it with respect to a covariance
parameter theta_d gives g_d = -Kdot_d(i,j)/sqrt(K_ii K_jj) - (rho_ij/2)
\[ Kdot_d(i,i)/K_ii + Kdot_d(j,j)/K_jj \], where Kdot_d = dK/dtheta_d.
The delta method gives the variance Var(rho_ij) = g' V g, with V the
estimated covariance of the parameter vector.

## Usage

``` r
cpp_pcor_delta(K, Sdot, pairs, V)
```

## Arguments

- K:

  the R x R precision matrix of the fitted layer.

- Sdot:

  list of D matrices dSigma_b/dtheta_d.

- pairs:

  integer matrix q x 2 with the response pairs (i, j), 1-based.

- V:

  estimated D x D covariance matrix of the parameter vector.

## Value

list with pcor (length q), variance (length q) and gradient (the q x D
matrix of gradients g).

## Details

The derivatives of the precision are obtained from the derivatives of
the covariance by Kdot_d = -K Sdot_d K, an identity that holds for any
parametrization. In the GGM layer it returns exactly the canonical base
Z_d; in the DAG layer it supplies the corresponding directions for
(gamma, delta) without changing the formula.

K is assumed positive definite and is not checked. For a pair with K_ii
K_jj \<= 0 the kernel returns rho = 0, a non-finite gradient and a
variance of zero. A quadratic form g' V g that is negative or not finite
is returned as zero.
