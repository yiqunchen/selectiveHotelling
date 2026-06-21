# =============================================================================
# hotelling_T2: two-sample Hotelling T^2 statistic.
# =============================================================================

test_that("hotelling_T2 equals a hand-computed value on a tiny 2-group example", {
  # Two groups of three points in R^2, chosen so the pooled covariance is
  # full-rank (avoids a singular solve()).
  X1 <- rbind(c(0, 0), c(2, 1), c(1, 3))
  X2 <- rbind(c(5, 4), c(7, 5), c(6, 8))
  X <- rbind(X1, X2)
  cl <- c(1, 1, 1, 2, 2, 2)

  got <- hotelling_T2(X, cl, 1, 2, ridge = 0)

  # Hand computation, mirroring the classical pooled-covariance formula.
  n1 <- 3
  n2 <- 3
  m1 <- colMeans(X1)
  m2 <- colMeans(X2)
  S1 <- stats::cov(X1)
  S2 <- stats::cov(X2)
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  d <- m1 - m2
  hand <- as.numeric((n1 * n2) / (n1 + n2) * t(d) %*% solve(Sp, d))

  expect_equal(got, hand, tolerance = 1e-10)
  expect_equal(got, 39.13514, tolerance = 1e-4)
})

test_that("hotelling_T2 returns NA when a group has fewer than two rows", {
  X1 <- rbind(c(0, 0), c(2, 1), c(1, 3), c(2, 2))
  X2 <- rbind(c(5, 4))  # only one observation
  X <- rbind(X1, X2)
  cl <- c(1, 1, 1, 1, 2)

  res <- hotelling_T2(X, cl, 1, 2)
  expect_true(is.na(res))
  expect_identical(res, NA_real_)

  # Symmetric: the OTHER group being a singleton also yields NA.
  cl2 <- c(2, 1, 1, 1, 1)
  expect_true(is.na(hotelling_T2(X, cl2, 1, 2)))
})

test_that("hotelling_T2 is symmetric in its two group labels", {
  X1 <- rbind(c(0, 0), c(2, 1), c(1, 3))
  X2 <- rbind(c(5, 4), c(7, 5), c(6, 8))
  X <- rbind(X1, X2)
  cl <- c(1, 1, 1, 2, 2, 2)

  expect_equal(hotelling_T2(X, cl, 1, 2, ridge = 0),
               hotelling_T2(X, cl, 2, 1, ridge = 0),
               tolerance = 1e-12)
})

test_that("hotelling_T2 is non-negative", {
  set.seed(11)
  X <- generate_exchangeable_null_data(20, 2, rho = 0.1)$dat
  cl <- rep(1:2, each = 10)
  val <- hotelling_T2(X, cl, 1, 2)
  expect_true(is.finite(val))
  expect_gte(val, 0)
})
