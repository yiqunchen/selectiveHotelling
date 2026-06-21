// kmeans_cpp_dev.cpp
// Bit-for-bit drop-in C++ (Rcpp) replacement for
// selectiveHotelling::kmeans_with_fixed_init().
//
// Replicates Lloyd's k-means from a FIXED initialization with full trajectory
// tracking, matching the pure-R reference EXACTLY:
//   * centers initialized to rows x[init_centers_idx, ] (1-based indices)
//   * assignment: cl[i] = which.min_j sum((x[i,] - center[j,])^2)
//     -- R's which.min breaks ties by LOWEST index (first minimum); we use a
//        strict "<" so the first (lowest j) minimum wins.
//   * trajectory[[t]] = assignment vector at iteration t
//   * convergence: if new assignment == previous assignment, STOP and trim
//   * update: center[j,] = colMeans(x[members,]); EMPTY cluster keeps its
//     previous center (not updated)
//   * cap at max_iter (default 50)
//   * labels are 1..k (1-based)
//
// Compile via: Rcpp::sourceCpp("kmeans_cpp_dev.cpp")
//
// FLOATING-POINT NOTE (bit-for-bit): R computes the squared Euclidean distance
// as rowSums((x - center)^2): each elementwise product (x-c)^2 is rounded to a
// double BEFORE summation, and the sum is a plain left-to-right double accumulate.
// A naive C++ "d2 += diff*diff" can be contracted by the compiler into a single
// fused multiply-add (FMA), which does NOT round the product first and therefore
// yields a 1-ULP-different distance. On exact ties that flips which.min and breaks
// bit-for-bit agreement. We disable FP contraction (pragma + -ffp-contract=off)
// and force the squared term through a named double so it is rounded like R.

// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
using namespace Rcpp;

// Disable fused-multiply-add / FP contraction so (diff*diff) rounds to double
// before accumulation, matching R's (x - center)^2 elementwise rounding.
#ifdef __clang__
#pragma STDC FP_CONTRACT OFF
#endif

// [[Rcpp::export]]
List kmeans_fixed_init_cpp(NumericMatrix x,
                           IntegerVector init_centers_idx,
                           int max_iter = 50) {
  const int n = x.nrow();
  const int p = x.ncol();
  const int k = init_centers_idx.size();

  // centers: k x p. Initialize from x[init_centers_idx, ] (1-based -> 0-based).
  NumericMatrix centers(k, p);
  for (int j = 0; j < k; ++j) {
    int ridx = init_centers_idx[j] - 1;  // 1-based to 0-based
    for (int d = 0; d < p; ++d) {
      centers(j, d) = x(ridx, d);
    }
  }

  List trajectory(max_iter);
  List center_history(max_iter);

  // cl: current assignment (1-based labels). cl_prev: previous iteration.
  IntegerVector cl(n);
  IntegerVector cl_prev;       // empty => "NULL" (no previous yet)
  bool have_prev = false;

  int t_run = 0;               // number of iterations actually executed

  for (int t = 0; t < max_iter; ++t) {
    // Record center snapshot at start of this iteration (clone: independent copy)
    center_history[t] = clone(centers);

    // ---- Assignment step ----
    // For each point i, find j minimizing squared Euclidean distance to center j.
    // Ties broken by lowest j (first minimum), matching which.min.
    for (int i = 0; i < n; ++i) {
      double best_dist = R_PosInf;
      int best_j = 0;
      for (int j = 0; j < k; ++j) {
        double d2 = 0.0;
        for (int dcol = 0; dcol < p; ++dcol) {
          double diff = x(i, dcol) - centers(j, dcol);
          // Force rounding of the squared term to double BEFORE accumulating,
          // so the result is bit-identical to R's rowSums((x - center)^2).
          double sq = diff * diff;
          d2 += sq;
        }
        // strict "<": keeps the FIRST (lowest-index) minimum, like which.min
        if (d2 < best_dist) {
          best_dist = d2;
          best_j = j;
        }
      }
      cl[i] = best_j + 1;  // 1-based label
    }

    // Store a copy of the assignment vector for this iteration.
    trajectory[t] = clone(cl);
    t_run = t + 1;

    // ---- Convergence check ----
    // Stop if assignment identical to previous iteration's.
    if (have_prev) {
      bool same = true;
      for (int i = 0; i < n; ++i) {
        if (cl[i] != cl_prev[i]) { same = false; break; }
      }
      if (same) {
        break;  // converged; t_run already = t+1, trim happens below
      }
    }
    cl_prev = clone(cl);
    have_prev = true;

    // ---- Update step ----
    // center[j,] = colMeans(x[members,]); empty cluster keeps previous center.
    // Accumulate sums and counts in one pass.
    NumericMatrix sums(k, p);   // zero-initialized
    std::vector<int> counts(k, 0);
    for (int i = 0; i < n; ++i) {
      int j = cl[i] - 1;
      counts[j] += 1;
      for (int dcol = 0; dcol < p; ++dcol) {
        sums(j, dcol) += x(i, dcol);
      }
    }
    for (int j = 0; j < k; ++j) {
      if (counts[j] > 0) {
        for (int dcol = 0; dcol < p; ++dcol) {
          centers(j, dcol) = sums(j, dcol) / counts[j];
        }
      }
      // else: empty cluster -> leave centers(j, .) unchanged
    }
  }

  // Trim trajectory and center_history to the iterations actually run.
  List traj_out(t_run);
  List ch_out(t_run);
  for (int t = 0; t < t_run; ++t) {
    traj_out[t] = trajectory[t];
    ch_out[t]   = center_history[t];
  }

  return List::create(
    _["cluster"]          = cl,
    _["trajectory"]       = traj_out,
    _["center_history"]   = ch_out,
    _["init_centers_idx"] = init_centers_idx,
    _["n_iter"]           = t_run
  );
}
