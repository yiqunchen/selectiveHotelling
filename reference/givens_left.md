# Givens rotation acting on two rows of Q

Applies a Givens rotation of angle `theta` to rows `i` and `j` of `Q` (a
left action on those two rows), used by the geodesic hit-and-run / slice
samplers.

## Usage

``` r
givens_left(Q, i, j, theta)
```

## Arguments

- Q:

  An `r` by `p` matrix with orthonormal columns.

- i, j:

  Row indices to rotate.

- theta:

  Rotation angle.

## Value

The rotated matrix.
