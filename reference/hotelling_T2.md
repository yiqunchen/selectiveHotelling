# Two-sample Hotelling T-squared statistic

Computes the classical two-sample Hotelling \\T^2\\ statistic comparing
the means of the rows of `X` assigned to cluster `k1` versus cluster
`k2`, using a pooled (ridge-stabilized) covariance estimate.

## Usage

``` r
hotelling_T2(X, cl, k1, k2, ridge = 1e-08)
```

## Arguments

- X:

  A numeric data matrix with observations in rows.

- cl:

  Integer/character vector of cluster labels, one per row of `X`.

- k1, k2:

  The two cluster labels to compare.

- ridge:

  Small non-negative ridge added to the diagonal of the pooled
  covariance for numerical stability.

## Value

A single numeric \\T^2\\ value, or `NA_real_` if either cluster has
fewer than two members.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(40), 20, 2)
cl <- rep(1:2, each = 10)
hotelling_T2(X, cl, 1, 2)
#> [1] 0.6052868
```
