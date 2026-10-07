# Hex logo: the circle of rotation angles (Figure 3(c) of the paper). The thick arcs are the
# angles at which the clustering still returns the tested pair; the black dot is the observed
# data (angle 0) and the triangle an accepted step. Run from the package root:
#   Rscript pkgdown/logo.R          -> man/figures/logo.png
library(ggplot2)
FEAS <- "#6A51A3"; INFEAS <- "#D4D4D4"; CL1 <- "#3B6FA0"; CL2 <- "#D9A13B"; INK <- "#2D1F4F"

hex <- data.frame(x = 0.95 * cos(pi / 6 + (0:5) * pi / 3), y = 0.95 * sin(pi / 6 + (0:5) * pi / 3))
arc <- function(a, b, n = 200) { t <- seq(a, b, length.out = n); data.frame(x = cos(t), y = sin(t)) }
R <- 0.42; cy <- 0.12                                   # circle radius and centre
sc <- function(d) transform(d, x = R * x, y = R * y + cy)
arcs <- list(c(-0.55, 0.75), c(1.9, 2.7), c(3.6, 4.6))  # feasible angles (radians)
feas <- do.call(rbind, lapply(seq_along(arcs), function(i) cbind(sc(arc(arcs[[i]][1], arcs[[i]][2])), g = i)))
obs <- sc(data.frame(x = 1, y = 0))
acc <- sc(data.frame(x = cos(2.3), y = sin(2.3)))
set.seed(3)                                             # two small tested clusters inside the circle
pts <- rbind(data.frame(x = rnorm(8, -0.13, 0.04), y = rnorm(8, cy + 0.03, 0.04), c = CL1),
             data.frame(x = rnorm(8,  0.13, 0.04), y = rnorm(8, cy - 0.03, 0.04), c = CL2))

p <- ggplot() +
  geom_polygon(data = hex, aes(x, y), fill = "white", colour = FEAS, linewidth = 1.6) +
  geom_path(data = sc(arc(0, 2 * pi, 400)), aes(x, y), colour = INFEAS, linewidth = 1.6) +
  geom_path(data = feas, aes(x, y, group = g), colour = FEAS, linewidth = 3.2, lineend = "round") +
  geom_point(data = pts, aes(x, y), colour = pts$c, size = 1.5) +
  geom_point(data = obs, aes(x, y), colour = "black", size = 3) +
  geom_point(data = acc, aes(x, y), shape = 24, fill = FEAS, colour = "white", size = 3.6, stroke = 0.6) +
  annotate("text", x = 0, y = -0.48, label = "selectiveHotelling", colour = INK,
           size = 3.6, fontface = "bold", family = "sans") +
  coord_equal(xlim = c(-0.9, 0.9), ylim = c(-1.02, 1.02), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "transparent", colour = NA))
ggsave("man/figures/logo.png", p, width = 1.73, height = 2, units = "in", dpi = 300, bg = "transparent")
