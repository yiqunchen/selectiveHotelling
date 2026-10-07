# Build polar-decomposition objects for a selected cluster pair

Constructs the fixed reduction used by all samplers. For the rows of `X`
in clusters `k1` and `k2`, it forms the projector onto the orthogonal
complement of the all-ones vector, an orthonormal basis `A` of that
complement, the Gram matrix `S = X_sub' P_perp X_sub` with its symmetric
square root, the observed polar factor `Q_obs` on the Stiefel manifold,
and the part of `X` that is held fixed while only `Q` varies.

## Usage

``` r
prepare_pair_objects(X, cl, k1, k2)
```

## Arguments

- X:

  A numeric data matrix with observations in rows.

- cl:

  Cluster-assignment vector for the rows of `X`.

- k1, k2:

  The two cluster labels defining the tested pair.

## Value

A list with elements `idx1`, `idx2`, `m`, `n`, `r`, `A`, `P_perp`,
`S_sqrt`, `X_fixed`, and `Q_obs`.
