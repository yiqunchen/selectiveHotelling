# aggregate.R -- collect the shards of the high-dimensional study and summarise
# Type I error (null cells) and power (all pairs, and pairs with pair_alt = 1).
#   Rscript dev/highdim/cluster/aggregate.R      (from the package root)
SHARDS <- Sys.getenv("HD_SHARDS", "dev/highdim/cluster/shards")
OUT    <- Sys.getenv("HD_OUT", "dev/highdim/cluster/results")
fs <- list.files(SHARDS, "^shard_\\d+\\.rds$", recursive = TRUE, full.names = TRUE)
fs <- fs[!grepl("_smoke", fs)]
r  <- do.call(rbind, lapply(fs, readRDS))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(r, file.path(OUT, "highdim_raw.rds"))
pc <- grep("^p_", names(r), value = TRUE)
summ <- function(d) {
  rej <- vapply(pc, function(m) mean(d[[m]] <= 0.05, na.rm = TRUE), numeric(1))
  data.frame(n_pairs = nrow(d), n_data = length(unique(d$rep)), pair_alt = mean(d$pair_alt),
             cpu_pair = mean(d$cpu), t(rej), check.names = FALSE)
}
keys <- c("alg", "n", "p", "ratio", "sig", "signal", "s")
size  <- do.call(rbind, lapply(split(r[r$signal == "null", ], r[r$signal == "null", keys], drop = TRUE),
                              function(d) cbind(d[1, keys], summ(d))))
power <- do.call(rbind, lapply(split(r[r$signal != "null", ], r[r$signal != "null", keys], drop = TRUE),
                              function(d) cbind(d[1, keys], summ(d))))
alt   <- r[r$signal != "null" & r$pair_alt == 1, ]
cpow  <- do.call(rbind, lapply(split(alt, alt[, keys], drop = TRUE), function(d) cbind(d[1, keys], summ(d))))
write.csv(size,  file.path(OUT, "highdim_size.csv"),  row.names = FALSE)
write.csv(power, file.path(OUT, "highdim_power.csv"), row.names = FALSE)
write.csv(cpow,  file.path(OUT, "highdim_power_alt.csv"), row.names = FALSE)
cat(sprintf("%d shards, %d pairs. Type I over null cells (median [max]):\n", length(fs), nrow(r)))
print(round(apply(size[, pc], 2, function(z) c(median = median(z), max = max(z))), 3))
cat("Wrote", OUT, "\n")
