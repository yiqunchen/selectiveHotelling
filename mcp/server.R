#!/usr/bin/env Rscript
# MCP server for selectiveHotelling (Model Context Protocol, stdio transport).
#
# Lets an AI agent cluster a data set and compute selective p-values with the package.
# Pure R: the only extra dependency is jsonlite. Register it with an MCP client, e.g. in
# Claude Code from the package root:
#   claude mcp add selectiveHotelling -- Rscript mcp/server.R
# or use the .mcp.json file at the package root. Run from a source checkout the server loads
# the package with pkgload::load_all(); otherwise it uses the installed package.
#
# Tools: cluster_data, selective_test, onedim_test, package_help.
# Every tool takes the data either as `data_path` (a CSV file) or as `data` (an array of rows),
# and re-clusters it from `seed`, so the same arguments always give the same clustering.

suppressPackageStartupMessages(library(jsonlite))

# ---- load the package ---------------------------------------------------------------------
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
pkg_root <- dirname(script_dir)
if (file.exists(file.path(pkg_root, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE)) {
  suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE, export_all = FALSE))
} else {
  suppressPackageStartupMessages(library(selectiveHotelling))
}
ns <- asNamespace("selectiveHotelling")
PKG_VERSION <- as.character(utils::packageVersion("selectiveHotelling"))

# ---- data and clustering ------------------------------------------------------------------
read_data <- function(a) {
  if (!is.null(a$data_path)) {
    df <- utils::read.csv(path.expand(a$data_path), check.names = FALSE)
  } else if (!is.null(a$data)) {
    df <- as.data.frame(do.call(rbind, lapply(a$data, function(r) as.numeric(unlist(r)))))
  } else stop("Give the data as `data_path` (a CSV file) or `data` (an array of rows).")
  if (!is.null(a$columns)) {
    cols <- unlist(a$columns)
    miss <- setdiff(cols, names(df))
    if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "))
    df <- df[, cols, drop = FALSE]
  }
  num <- vapply(df, is.numeric, logical(1))
  dropped <- names(df)[!num]
  X <- as.matrix(df[, num, drop = FALSE])
  ok <- stats::complete.cases(X)
  X <- X[ok, , drop = FALSE]
  if (isTRUE(a$scale)) X <- scale(X)
  if (ncol(X) < 1 || nrow(X) < 3) stop("Need at least 3 complete rows and 1 numeric column.")
  list(X = unname(X), features = colnames(df)[num], dropped_columns = dropped,
       rows_removed_na = sum(!ok))
}

make_clustering <- function(X, a) {
  method <- if (is.null(a$clustering)) "kmeans" else a$clustering
  K <- if (is.null(a$K)) 3L else as.integer(a$K)
  seed <- if (is.null(a$seed)) 1L else as.integer(a$seed)
  if (method == "kmeans") {
    set.seed(seed)
    km <- ns$kmeans_with_trajectory(X, k = K)          # draws the initial centers once
    init <- km$init_centers_idx
    cl_fun <- function(x) kmeans_with_fixed_init(x, init_centers_idx = init)$cluster
    list(method = "kmeans", K = K, seed = seed, cl = as.integer(km$cluster), cl_fun = cl_fun,
         init = init)
  } else if (method == "hclust") {
    linkage <- if (is.null(a$linkage)) "average" else a$linkage
    cl_fun <- ns$make_hclust_cl_fun(k = K, method = linkage, square = TRUE)
    list(method = paste0("hclust (", linkage, " linkage)"), K = K, seed = seed,
         cl = as.integer(cl_fun(X)), cl_fun = cl_fun, init = NULL, linkage = linkage)
  } else stop("`clustering` must be \"kmeans\" or \"hclust\".")
}

naive_pval <- function(X, cl, k1, k2) {
  n1 <- sum(cl == k1); n2 <- sum(cl == k2); p <- ncol(X); r <- n1 + n2
  if (r < p + 2) return(NA_real_)
  T2 <- hotelling_T2(X, cl, k1, k2)
  stats::pf((r - p - 1) / ((r - 2) * p) * T2, p, r - p - 1, lower.tail = FALSE)
}

pairs_from <- function(a, K) {
  if (identical(a$pairs, "all") || (is.null(a$k1) && is.null(a$pairs))) return(utils::combn(K, 2))
  if (!is.null(a$pairs)) return(do.call(cbind, lapply(a$pairs, function(pr) as.integer(unlist(pr)))))
  matrix(c(as.integer(a$k1), as.integer(a$k2)), 2)
}

# ---- tools --------------------------------------------------------------------------------
tool_cluster_data <- function(a) {
  d <- read_data(a); X <- d$X; cl <- make_clustering(X, a)
  sizes <- tabulate(cl$cl, nbins = cl$K)
  means <- t(vapply(seq_len(cl$K), function(k) colMeans(X[cl$cl == k, , drop = FALSE]), numeric(ncol(X))))
  pr <- utils::combn(cl$K, 2)
  testable <- apply(pr, 2, function(kk) sum(sizes[kk]) >= ncol(X) + 2)
  out <- list(n = nrow(X), p = ncol(X), features = I(d$features), dropped_columns = I(d$dropped_columns),
              rows_removed_na = d$rows_removed_na, clustering = cl$method, K = cl$K, seed = cl$seed,
              cluster_sizes = sizes, cluster_means = unname(split(round(means, 4), row(means))),
              testable_pairs = lapply(which(testable), function(j) pr[, j]),
              untestable_pairs = lapply(which(!testable), function(j) pr[, j]))
  if (nrow(X) <= 2000) out$labels <- cl$cl
  out$note <- paste("Pass the same data, clustering, K and seed to selective_test to test this clustering.",
                    "A pair is testable when its two clusters hold at least p + 2 observations.")
  out
}

tool_selective_test <- function(a) {
  d <- read_data(a); X <- d$X; cl <- make_clustering(X, a)
  M <- if (is.null(a$M)) 5000L else as.integer(a$M)
  L <- if (is.null(a$L)) 2L else as.integer(a$L)
  pr <- pairs_from(a, cl$K); p <- ncol(X)
  res <- lapply(seq_len(ncol(pr)), function(j) {
    k1 <- pr[1, j]; k2 <- pr[2, j]; n1 <- sum(cl$cl == k1); n2 <- sum(cl$cl == k2)
    row <- list(k1 = k1, k2 = k2, n1 = n1, n2 = n2)
    if (n1 + n2 < p + 2) return(c(row, list(error = sprintf("n1 + n2 = %d < p + 2 = %d: not testable.", n1 + n2, p + 2))))
    set.seed(cl$seed + j)
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(selective_hotelling_bc(X, cl$cl_fun, k1, k2, M = M, L = L, cl_obs = cl$cl),
                  error = function(e) e)
    if (inherits(r, "error")) return(c(row, list(error = conditionMessage(r))))
    c(row, list(T2 = signif(r$T_obs, 6), naive_pval = signif(naive_pval(X, cl$cl, k1, k2), 4),
                selective_pval = signif(r$pval, 4), mc_se = signif(r$se, 3),
                move_fraction = round(r$move_frac, 3), seconds = round(proc.time()[["elapsed"]] - t0, 1)))
  })
  list(test = "p_polar: selective Hotelling test, Besag-Clifford with the random Givens kernel",
       n = nrow(X), p = p, clustering = cl$method, K = cl$K, seed = cl$seed,
       cluster_sizes = tabulate(cl$cl, nbins = cl$K), M = M, L = L, results = res,
       note = paste("The selective p-value is exact for any M and L; mc_se is its Monte Carlo",
                    "standard error. The naive p-value ignores that the clusters were estimated."))
}

tool_onedim_test <- function(a) {
  d <- read_data(a); X <- d$X; cl <- make_clustering(X, a)
  method <- if (is.null(a$method)) "bc" else a$method
  M <- if (is.null(a$M)) 5000L else as.integer(a$M)
  L <- if (is.null(a$L)) 2L else as.integer(a$L)
  pr <- pairs_from(a, cl$K)
  res <- lapply(seq_len(ncol(pr)), function(j) {
    k1 <- pr[1, j]; k2 <- pr[2, j]
    row <- list(k1 = k1, k2 = k2, n1 = sum(cl$cl == k1), n2 = sum(cl$cl == k2))
    set.seed(cl$seed + j)
    r <- tryCatch(
      if (method == "exact") {
        if (cl$method != "kmeans") stop("method = \"exact\" is available for k-means only.")
        selective_truncated_F(X, k1, k2, clustering = "kmeans", method = "exact", K = cl$K,
                              init_centers_idx = cl$init)
      } else {
        selective_truncated_F(X, k1, k2, clustering = "custom", method = "bc", onedim_kernel = "givens",
                              cl_fun = cl$cl_fun, cl_obs = cl$cl, M = M, L = L)
      }, error = function(e) e)
    if (inherits(r, "error")) return(c(row, list(error = conditionMessage(r))))
    c(row, list(pval = signif(r$pval, 4), mc_se = if (is.na(r$se)) NULL else signif(r$se, 3)))
  })
  list(test = if (method == "exact") "p_polar,1d^traj: closed form, also conditions on the k-means trajectory"
              else "p_polar,1d: one-dimensional selective test (Algorithm 2)",
       clustering = cl$method, K = cl$K, seed = cl$seed, results = res)
}

tool_package_help <- function(a) {
  topic <- if (is.null(a$topic)) "selective_hotelling_test" else a$topic
  db <- if (dir.exists(file.path(pkg_root, "man"))) tools::Rd_db(dir = pkg_root) else tools::Rd_db("selectiveHotelling")
  key <- paste0(topic, ".Rd")
  if (!key %in% names(db)) return(list(error = paste("No help topic", topic), topics = sub("\\.Rd$", "", names(db))))
  txt <- utils::capture.output(tools::Rd2txt(db[[key]], options = list(underline_titles = FALSE)))
  list(topic = topic, text = paste(txt, collapse = "\n"))
}

data_props <- list(
  data_path  = list(type = "string", description = "Path to a CSV file with a header row. Non-numeric columns are dropped."),
  data       = list(type = "array", items = list(type = "array", items = list(type = "number")),
                    description = "The data as an array of rows (use instead of data_path)."),
  columns    = list(type = "array", items = list(type = "string"), description = "Columns of the CSV to use (default: all numeric columns)."),
  scale      = list(type = "boolean", description = "Standardize each column first (default false)."),
  clustering = list(type = "string", enum = list("kmeans", "hclust"), description = "Clustering algorithm (default kmeans)."),
  K          = list(type = "integer", description = "Number of clusters (default 3)."),
  linkage    = list(type = "string", description = "hclust linkage (default average)."),
  seed       = list(type = "integer", description = "Seed for the k-means initial centers and the Monte Carlo (default 1).")
)
pair_props <- list(
  k1 = list(type = "integer", description = "First cluster to test."),
  k2 = list(type = "integer", description = "Second cluster to test."),
  pairs = list(description = "Instead of k1/k2: \"all\" or an array of [k1, k2] pairs. Default: all pairs."),
  M = list(type = "integer", description = "Path length (default 5000; 1000 is enough for exploration). Fix it before looking at the data."),
  L = list(type = "integer", description = "Steps between points of the path (default 2).")
)

TOOLS <- list(
  list(name = "cluster_data", fn = tool_cluster_data,
       description = paste("Cluster a data set with k-means or hierarchical clustering and report the cluster sizes,",
                           "means and which pairs can be tested. Use before selective_test to see the clusters."),
       inputSchema = list(type = "object", properties = data_props)),
  list(name = "selective_test", fn = tool_selective_test,
       description = paste("Test whether two estimated clusters have the same mean, with a selective p-value",
                           "(p_polar) that accounts for the clusters having been estimated from the same data and",
                           "needs no estimate of the covariance matrix. Also returns the naive Hotelling p-value.",
                           "Assumes Gaussian rows with a common covariance. Takes seconds to minutes per pair."),
       inputSchema = list(type = "object", properties = c(data_props, pair_props))),
  list(name = "onedim_test", fn = tool_onedim_test,
       description = paste("The one-dimensional selective test: p_polar,1d (method \"bc\", any clustering) or its",
                           "closed form p_polar,1d^traj (method \"exact\", k-means only, instant). Both condition on",
                           "more than selective_test and have less power."),
       inputSchema = list(type = "object", properties = c(data_props, pair_props,
         list(method = list(type = "string", enum = list("bc", "exact"), description = "Default bc."))))),
  list(name = "package_help", fn = tool_package_help,
       description = "Return the help page of a selectiveHotelling function (default selective_hotelling_test).",
       inputSchema = list(type = "object", properties = list(topic = list(type = "string"))))
)

# ---- JSON-RPC over stdio ------------------------------------------------------------------
send <- function(x) { cat(toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, na = "null"), "\n", sep = ""); flush(stdout()) }
reply <- function(id, result) send(list(jsonrpc = "2.0", id = id, result = result))
fail  <- function(id, code, msg) send(list(jsonrpc = "2.0", id = id, error = list(code = code, message = msg)))

handle <- function(req) {
  id <- req$id; m <- req$method
  if (is.null(id)) return(invisible())                       # notifications need no reply
  if (m == "initialize") {
    pv <- req$params$protocolVersion
    reply(id, list(protocolVersion = if (is.null(pv)) "2025-06-18" else pv,
                   capabilities = list(tools = setNames(list(), character(0))),
                   serverInfo = list(name = "selectiveHotelling", version = PKG_VERSION),
                   instructions = paste("Selective inference after clustering with unknown covariance.",
                                        "Fix the pair to test (or test all pairs) and M before looking at results.")))
  } else if (m == "ping") {
    reply(id, setNames(list(), character(0)))
  } else if (m == "tools/list") {
    reply(id, list(tools = lapply(TOOLS, function(t) t[c("name", "description", "inputSchema")])))
  } else if (m == "tools/call") {
    nm <- req$params$name; args <- req$params$arguments
    tool <- Filter(function(t) t$name == nm, TOOLS)
    if (!length(tool)) return(fail(id, -32602, paste("Unknown tool:", nm)))
    out <- NULL
    err <- tryCatch({ utils::capture.output(out <- suppressWarnings(tool[[1]]$fn(if (is.null(args)) list() else args))); NULL },
                    error = function(e) conditionMessage(e))
    if (!is.null(err)) {
      reply(id, list(content = list(list(type = "text", text = paste("Error:", err))), isError = TRUE))
    } else {
      txt <- toJSON(out, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA, na = "null")
      reply(id, list(content = list(list(type = "text", text = as.character(txt))), isError = FALSE))
    }
  } else {
    fail(id, -32601, paste("Method not found:", m))
  }
}

con <- file("stdin", open = "r")
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if (!nzchar(trimws(line))) next
  req <- tryCatch(fromJSON(line, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(req)) { fail(NULL, -32700, "Parse error"); next }
  tryCatch(handle(req), error = function(e) fail(req$id, -32603, conditionMessage(e)))
}
