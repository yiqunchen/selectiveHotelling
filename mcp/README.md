# MCP server

`server.R` is a [Model Context Protocol](https://modelcontextprotocol.io) server that lets an AI
agent (Claude Code, Claude Desktop, or any MCP client) run the selective tests in this package.
It is plain R and needs only `jsonlite`. The package itself is documented at
<https://yiqunchen.github.io/selectiveHotelling/>. Run from a clone of this repository, it loads the
package source with `pkgload`; otherwise it uses the installed package.

## Setup

In Claude Code, start a session in a clone of this repository: `.mcp.json` registers the server.
Elsewhere:

```bash
claude mcp add selectiveHotelling -- Rscript /path/to/selectiveHotelling/mcp/server.R
```

For other clients, the command is `Rscript /path/to/selectiveHotelling/mcp/server.R` (stdio transport).

## Tools

| Tool | What it returns |
|---|---|
| `cluster_data` | cluster sizes and means, and which pairs have enough observations to test |
| `selective_test` | for each pair: the selective p-value, its Monte Carlo standard error, and the naive p-value |
| `onedim_test` | the one-dimensional test: `method = "bc"` (any clustering) or `"exact"` (k-means closed form) |
| `package_help` | the help page of a package function |

The data go in as `data_path` (a CSV file; non-numeric columns are dropped) or `data` (an array
of rows). The clustering is set by `clustering` (`"kmeans"` or `"hclust"`), `K`, `linkage`, and
`seed`; the same arguments always give the same clustering, so `cluster_data` and
`selective_test` see the same clusters. `M` (default 5000) and `L` (default 2) set the
Monte Carlo path. All pairs are tested unless `k1` and `k2` or `pairs` are given.

Example request an agent might make: *"Cluster `penguins.csv` into 3 groups with k-means and
test every pair of clusters."*

## Check that it works

```bash
printf '%s\n' '{"jsonrpc":"2.0","id":1,"method":"tools/list"}' | Rscript mcp/server.R
```
