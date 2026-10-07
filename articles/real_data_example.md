# Real data: Palmer penguins

``` r

library(selectiveHotelling)
library(palmerpenguins)
```

The Palmer penguins data contain four body measurements for penguins of
three species. We standardize the measurements, cluster the penguins
with $`k`$-means, and test every pair of clusters.

``` r

peng <- na.omit(penguins[, c("species", "bill_length_mm", "bill_depth_mm",
                             "flipper_length_mm", "body_mass_g")])
X <- scale(as.matrix(peng[, -1]))
species <- peng$species
dim(X)
#> [1] 342   4
```

``` r

set.seed(1)
init <- sample.int(nrow(X), 3)
cl_fun <- function(x) kmeans_with_fixed_init(x, init_centers_idx = init)$cluster
cl <- cl_fun(X)
table(cluster = cl, species)
#>        species
#> cluster Adelie Chinstrap Gentoo
#>       1      7        63      0
#>       2      0         0    123
#>       3    144         5      0
```

``` r

naive_pval <- function(X, cl, k1, k2) {
  n1 <- sum(cl == k1); n2 <- sum(cl == k2); p <- ncol(X); r <- n1 + n2
  F_stat <- (r - p - 1) / ((r - 2) * p) * hotelling_T2(X, cl, k1, k2)
  pf(F_stat, p, r - p - 1, lower.tail = FALSE)
}
pairs <- rbind(c(1, 2), c(1, 3), c(2, 3))
out <- t(apply(pairs, 1, function(k) {
  set.seed(10 * k[1] + k[2])
  c(k1 = k[1], k2 = k[2],
    naive = naive_pval(X, cl, k[1], k[2]),
    selective = selective_hotelling_bc(X, cl_fun, k[1], k[2], M = 1000)$pval)
}))
signif(out, 3)
#>      k1 k2     naive selective
#> [1,]  1  2 5.08e-117   0.00091
#> [2,]  1  3  4.27e-68   0.20800
#> [3,]  2  3 1.03e-167   0.00837
```

The clusters match the species: cluster 1 is mostly Chinstrap, cluster 2
is Gentoo, and cluster 3 is mostly Adelie. The naive test rejects every
null hypothesis with a vanishingly small $`p`$-value. The selective test
rejects the two contrasts that involve the Gentoo cluster, the species
that differs most in these measurements, but not the contrast between
the Chinstrap and Adelie clusters.
