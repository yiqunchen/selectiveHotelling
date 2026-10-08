args <- commandArgs(TRUE); f <- if (length(args)) args[1] else "dev/highdim/sim_raw.rds"
r <- readRDS(f); r <- r[r$ok, ]
meth <- c("hotelling", "l2", "cq", "diag", "diag_w", "max", "ridge", "proj", "adapt", "oracle")
rej <- function(d) sapply(meth, function(m) { x <- d[[m]]; if (all(is.na(x))) NA else round(mean(x <= 0.05, na.rm = TRUE), 3) })
cat("== TYPE I (kappa = 0), all reps ==\n")
s <- do.call(rbind, lapply(split(r[r$kappa == 0, ], ~ p + sig, drop = TRUE), function(d)
  data.frame(p = d$p[1], sig = d$sig[1], nrep = nrow(d), med_r = median(d$r), t(rej(d)))))
print(s[order(s$sig, s$p), ], row.names = FALSE)
cat("\n== POWER (unconditional), by signal/sig/p/kappa ==\n")
pw <- r[r$kappa > 0, ]
s <- do.call(rbind, lapply(split(pw, ~ signal + sig + p + kappa, drop = TRUE), function(d)
  data.frame(signal = d$signal[1], sig = d$sig[1], p = d$p[1], kappa = d$kappa[1],
             pure = round(mean(d$pure), 2), t(rej(d)))))
print(s[order(s$signal, s$sig, s$p, s$kappa), ], row.names = FALSE)
