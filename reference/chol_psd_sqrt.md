# Symmetric positive-semidefinite matrix square root

Symmetrizes a matrix, eigendecomposes it, thresholds tiny or negative
eigenvalues at `tol`, and returns both the symmetric square root and
inverse square root.

## Usage

``` r
chol_psd_sqrt(M, tol = 1e-10)
```

## Arguments

- M:

  A (numerically) symmetric matrix.

- tol:

  Lower bound applied to eigenvalues before taking square roots.

## Value

A list with elements `sqrt` (a symmetric square root) and `invsqrt` (its
inverse).
