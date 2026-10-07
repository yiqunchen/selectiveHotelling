# One left-multiplicative Cayley random-walk step on the Stiefel manifold

Proposes a new Stiefel point by left-multiplying `Q` by an orthogonal
matrix generated from a random skew-symmetric matrix via the Cayley
transform, with step size `eps`.

## Usage

``` r
cayley_left_step(Q, eps)
```

## Arguments

- Q:

  An `r` by `p` matrix with orthonormal columns.

- eps:

  Step size controlling proposal scale.

## Value

A new `r` by `p` matrix with orthonormal columns.
