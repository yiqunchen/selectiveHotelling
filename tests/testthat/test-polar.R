# =============================================================================
# Polar decomposition / reconstruction.
#
# prepare_pair_objects() performs the exact polar reduction X_sub = A Q S^{1/2}
# (plus a fixed part). Reconstructing the data from the OBSERVED polar factor
# Q_obs must return the original X on the rows of the selected cluster pair.
# This mirrors the X_from_Q closure used inside every sampler.
# =============================================================================

test_that("prepare_pair_objects reconstructs X exactly from Q_obs on the pair rows", {
  prepare_pair_objects <- sh_internal("prepare_pair_objects")

  d <- make_null_clustering(seed = 42)
  X <- d$X
  cl <- d$cl

  tab <- table(cl)
  usable <- as.integer(names(tab)[tab >= 2])
  expect_gte(length(usable), 2L)
  k1 <- usable[1]
  k2 <- usable[2]

  prep <- prepare_pair_objects(X, cl, k1, k2)

  # Reconstruct using the same algebra the samplers use.
  XQ <- x_from_Q(prep, prep$Q_obs)

  pair_rows <- c(prep$idx1, prep$idx2)
  expect_equal(XQ[pair_rows, , drop = FALSE],
               X[pair_rows, , drop = FALSE],
               tolerance = 1e-7)
  expect_lt(max(abs(XQ[pair_rows, ] - X[pair_rows, ])), 1e-7)

  # Rows NOT in the pair are held fixed and must equal X identically.
  other_rows <- setdiff(seq_len(nrow(X)), pair_rows)
  if (length(other_rows) > 0) {
    expect_equal(XQ[other_rows, , drop = FALSE],
                 X[other_rows, , drop = FALSE],
                 tolerance = 1e-10)
  }
})

test_that("Q_obs lies on the Stiefel manifold (orthonormal columns)", {
  prepare_pair_objects <- sh_internal("prepare_pair_objects")

  d <- make_null_clustering(seed = 7)
  X <- d$X
  cl <- d$cl
  tab <- table(cl)
  usable <- as.integer(names(tab)[tab >= 2])
  k1 <- usable[1]
  k2 <- usable[2]

  prep <- prepare_pair_objects(X, cl, k1, k2)
  G <- t(prep$Q_obs) %*% prep$Q_obs
  expect_equal(G, diag(nrow(G)), tolerance = 1e-8)
})
