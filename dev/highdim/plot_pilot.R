# Pilot figure: Type I (left) and power vs kappa (right), n = 50, p = 100.
suppressMessages({ library(ggplot2); library(patchwork) })
m  <- c("l2", "cq", "diag", "diag_w", "max", "ridge", "proj", "adapt", "oracle")
a  <- readRDS("dev/highdim/sim_raw.rds"); b <- readRDS("dev/highdim/sim_big_raw.rds")
l  <- readRDS("dev/highdim/sim_line_raw.rds")
r  <- rbind(a, b[b$kappa > 0, ], l[l$kappa > 0, ]); r <- r[r$ok, ]
r$signal <- factor(r$signal, levels = c("line", "sparse", "block", "dense"),
                   labels = c("1-sparse (line)", "5-sparse", "1/4-sparse (block)", "dense"))
long <- do.call(rbind, lapply(m, function(s) data.frame(r[, c("sig", "signal", "kappa", "pure")], stat = s, rej = r[[s]] <= 0.05)))
long$stat <- factor(long$stat, levels = m)
cols <- c(l2 = "#0072B2", cq = "#56B4E9", diag = "#E69F00", diag_w = "#D55E00", max = "#CC79A7",
          ridge = "#009E73", proj = "#F0E442", adapt = "#999999", oracle = "black")

t1 <- aggregate(rej ~ sig + stat, long[long$kappa == 0, ], function(x) c(m = mean(x), n = length(x)))
t1 <- data.frame(t1[1:2], rate = t1$rej[, "m"], n = t1$rej[, "n"])
se <- sqrt(0.05 * 0.95 / median(t1$n))
p1 <- ggplot(t1, aes(sig, rate, colour = stat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05 - 2 * se, ymax = 0.05 + 2 * se, alpha = 0.12) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_point(position = position_dodge(0.7), size = 1.8) +
  scale_colour_manual(values = cols) + coord_cartesian(ylim = c(0, 0.15)) +
  labs(x = "covariance", y = "Type I error at 0.05", colour = NULL) + theme_bw(10) + guides(colour = "none")

pw <- long[long$kappa > 0 & long$sig %in% c("iso", "cs", "hetero"), ]
pa <- aggregate(rej ~ sig + signal + kappa + stat, pw, mean)
pu <- aggregate(pure ~ sig + signal + kappa, pw[pw$stat == "l2", ], mean)
p2 <- ggplot(pa, aes(kappa, rej, colour = stat)) +
  geom_line(data = pu, aes(kappa, pure), inherit.aes = FALSE, linetype = 3, colour = "grey40") +
  geom_line(linewidth = 0.5) + geom_point(size = 1) +
  facet_grid(sig ~ signal) + scale_x_log10(breaks = c(4, 16, 64, 256)) +
  scale_colour_manual(values = cols) + ylim(0, 1) +
  labs(x = expression(kappa ~ "(signal size, log scale)"), y = "power (dotted: P(pair recovered))", colour = NULL) +
  theme_bw(10)

g <- p1 + p2 + plot_layout(widths = c(1, 2), guides = "collect")
ggsave("dev/highdim/pilot_typeI_power.png", g, width = 14, height = 5.5, dpi = 150)
