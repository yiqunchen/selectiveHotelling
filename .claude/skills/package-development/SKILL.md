---
name: package-development
description: Develop, test, document, and release the selectiveHotelling R package and its pkgdown site. Use when changing code in R/, adding a function or kernel, editing tutorials in vignettes/articles/, running R CMD check, or deploying the website.
---

# Developing selectiveHotelling

Pure R, no compiled code. `devtools` may be unavailable; use `pkgload`, `testthat`, `roxygen2`.

## Loop

```bash
Rscript -e 'pkgload::load_all(".")'                                     # load
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'  # tests
Rscript -e 'roxygen2::roxygenise()'                                     # after editing roxygen comments
R CMD build . && R CMD check --no-manual selectiveHotelling_*.tar.gz    # before pushing
```

Run the tests after every change to `R/`. A change is done when `R CMD check` reports
`Status: OK`.

## Rules that protect validity

The Besag-Clifford p-value is exact for any `M`, any `L`, and any reversible kernel **only if**
these hold. Never break them, and add a test when touching the samplers.

- The kernel, `M`, `L` and the step-size distribution may depend on the data only through
  `P2 X`, `H` and the clustering, never on `Q_obs`. No step-size tuning from trial moves at the
  observed point (that breaks exactness; see Appendix C of the paper). A step size drawn afresh
  from a fixed distribution (`eps = "mixture"`) is fine.
- The selection event is "the clustering returns the tested pair"; feasibility is checked by
  rerunning `cl_fun`.
- `cl_fun` must be deterministic; k-means uses initial centers drawn once
  (`kmeans_with_fixed_init`).
- No grid approximation of the selection event anywhere (author rule).
- The test requires `r >= p + 2`.

## Layout

- `R/besag_clifford.R`: `selective_hotelling_bc()`, the main test (random Givens kernel by default).
- `R/test.R`: `selective_hotelling_test()` wrapper for k-means / hclust, plus print and summary.
- `R/truncated_f.R`, `R/truncated_f_mc.R`: `selective_truncated_F()`, the one-dimensional tests.
- `R/harmonic_arc.R`: closed-form arc solvers. `R/polar.R`: polar decomposition, Givens and Cayley moves.
- `R/clustering.R`, `R/kmeans.R`, `R/utils.R`: helpers.
- `mcp/server.R`: MCP server exposing the tests to agents (update its tools if signatures change).
- `vignettes/articles/`: pkgdown tutorials (not built with the package).

## Website

Built with pkgdown from `_pkgdown.yml`, `pkgdown/extra.css`, and `vignettes/articles/`; the logo
is drawn by `pkgdown/logo.R`. Every exported function must appear in the `reference:` section of
`_pkgdown.yml`, or the build fails.

```bash
Rscript -e 'pkgdown::build_site(preview = FALSE)'      # writes docs/ (gitignored)
```

Deploy by replacing the contents of the `gh-pages` branch with `docs/` (keep a `.nojekyll`
file) and pushing. Check the live site at https://yiqunchen.github.io/selectiveHotelling/.

## Writing style for tutorials and docs

American spelling. No math in section headings (it breaks the table of contents). Color-blind-safe
colors in plots (e.g. `#0072B2`, `#E69F00`, `#009E73`). Plain words over jargon. Numbers quoted in
text must come from the rendered output.
