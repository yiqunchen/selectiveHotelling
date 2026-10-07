# selectiveHotelling: Selective Hotelling T-squared Testing After Clustering with Unknown Covariance

Tests whether two clusters found by a clustering algorithm have the same
mean, with Hotelling's \\T^2\\ and an unknown covariance matrix
\\\Sigma\\. The p-values condition on the clustering and on the
statistics that remove the unknown means and \\\Sigma\\, so they control
the selective Type I error.

## Main functions

- [`selective_hotelling_test()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_test.md) -
  clusters the data (k-means or hierarchical clustering) and runs the
  main test.

- [`selective_hotelling_bc()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_hotelling_bc.md) -
  the main test for any clustering function `cl_fun`: the Besag and
  Clifford (1989) construction with the random Givens kernel on the
  event that the clustering returns the tested pair. The p-value is
  exactly uniform under the null for any path length and any reversible
  kernel; the kernel and the path length affect only power.

- [`selective_truncated_F()`](https://yiqunchen.github.io/selectiveHotelling/reference/selective_truncated_F.md) -
  the one-dimensional test, which also fixes all of the polar factor
  except one angle. It is computed by the same construction on the angle
  for any clustering function, and in closed form for k-means when the
  conditioning includes the whole Lloyd trajectory.

- [`kmeans_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_fun.md)
  /
  [`kmeans_with_fixed_init()`](https://yiqunchen.github.io/selectiveHotelling/reference/kmeans_with_fixed_init.md) -
  k-means with its trajectory.

- [`hclust_fun()`](https://yiqunchen.github.io/selectiveHotelling/reference/hclust_fun.md) -
  hierarchical clustering (squared-distance default).

- [`hotelling_T2()`](https://yiqunchen.github.io/selectiveHotelling/reference/hotelling_T2.md) -
  the two-sample Hotelling \\T^2\\ statistic.

- [`generate_exchangeable_null_data()`](https://yiqunchen.github.io/selectiveHotelling/reference/generate_exchangeable_null_data.md) -
  simulate null data.

## References

Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance
tests. *Biometrika*, 76(4), 633-642.

Chen, Y. T. and Witten, D. M. (2023). Selective inference for k-means
clustering. *Journal of Machine Learning Research*.

Gao, L. L., Bien, J. and Witten, D. (2024). Selective inference for
hierarchical clustering. *Journal of the American Statistical
Association*.

## See also

Useful links:

- <https://github.com/yiqunchen/selectiveHotelling>

- Report bugs at
  <https://github.com/yiqunchen/selectiveHotelling/issues>

## Author

**Maintainer**: Yiqun Chen <yiqun.t.chen@gmail.com>
