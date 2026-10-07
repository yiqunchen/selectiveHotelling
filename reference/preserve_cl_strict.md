# Check preservation of the final cluster pair

Legacy/final-only feasibility check: requires that the observations in
original clusters `k1` and `k2` form exactly those two clusters (up to
relabeling) in the candidate final assignment `cl1`.

## Usage

``` r
preserve_cl_strict(cl0, cl1, k1, k2)
```

## Arguments

- cl0, cl1:

  Original and candidate final cluster-assignment vectors.

- k1, k2:

  The two original cluster labels being tested.

## Value

Logical scalar.
