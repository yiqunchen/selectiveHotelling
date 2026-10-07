# Agent skills

Skills are instructions that an AI agent (Claude Code and other tools that read
`.claude/skills/`) loads when a task matches their description. Each one is a folder with a
`SKILL.md` file.

| Skill | Used when | What it covers |
|---|---|---|
| `selective-test-after-clustering` | someone wants to test differences between estimated clusters | the checks to do first (pair fixed in advance, `n1 + n2 >= p + 2`, deterministic clustering, `M` chosen in advance, roughly Gaussian rows), which of the three p-values to use, and how to report the results |
| `package-development` | someone changes the package code, tutorials, or website | the test and check loop, the rules that keep the p-values exact, the file layout, and how to build and deploy the website |

The first skill uses the [MCP server](../../mcp/README.md) when it is connected, and R code
otherwise. To use the skills in another project, copy their folders into that project's
`.claude/skills/` (or `~/.claude/skills/` for all projects).
