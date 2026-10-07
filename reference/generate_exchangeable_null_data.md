# Generate exchangeable (equicorrelated) Gaussian null data

Draws an `n` by `p` matrix of multivariate Gaussian observations with a
common mean `mu` and an exchangeable (equicorrelation) covariance:
variance `sig^2` on the diagonal and `sig^2 * rho` off-diagonal. Under
this model all observations share the same mean, so it serves as a null
for cluster-mean comparisons.

## Usage

``` r
generate_exchangeable_null_data(n, p, rho = 0, mu = 0, sig = 1)
```

## Arguments

- n:

  Number of observations (rows).

- p:

  Number of variables (columns).

- rho:

  Equicorrelation parameter (off-diagonal correlation).

- mu:

  Common mean added to every entry.

- sig:

  Common standard deviation.

## Value

A list with elements `dat` (the `n` by `p` data matrix) and `Sig` (the
`p` by `p` covariance matrix used).

## Examples

``` r
set.seed(1)
out <- generate_exchangeable_null_data(30, 2, rho = 0.1)
dim(out$dat)
#> [1] 30  2
```
