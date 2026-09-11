# FoF algorithm illustration: fit separately for Tianjin and Xi'an origins.
# Required inputs / variable definitions
# F0/F1: wide-format matrices; Step 2 expands F1 to long format
# while retaining the full F0 curve for each output time point.
# F0 / F1: complete matrices, one retained token per row, 11 columns for time 0:10.
# F0 = interpolated within-speaker z-scored f0; F1 = within-speaker z-scored F1.
# meta: one row per token, in exactly the same order as F0 and F1.
# durationZ = standardized duration; speaker = ordinary speaker-ID factor.
# DiphVar.ord = diphthong x variety; sex.ord = sex; type.ord = carrier sentence;
# order.ord = reading order; leftSegment.ord = onset.
# All .ord columns are ordered factors with original reference levels and treatment contrasts.

library(mgcv)
library(itsadug)
library(ggplot2)

# Preserve the original reference levels and ordered-factor smooth coding.
for (column in c("DiphVar.ord", "sex.ord", "type.ord", "order.ord", "leftSegment.ord")) {
  stopifnot(is.ordered(meta[[column]]))
  contrasts(meta[[column]]) <- "contr.treatment"
}

# 1. Center f0 at each input time and set trapezoidal integration weights.
n <- nrow(F0)
K <- 11
s_grid <- 0:10
t_grid <- 0:10
F0_centered <- sweep(F0, 2, colMeans(F0), "-")
weights <- c(0.5, rep(1, 9), 0.5) / 10

# 2. Each model row represents one token at one output time.
dat <- meta[rep(1:n, each = K), ]
dat$Y <- as.vector(t(F1))
dat$t_out <- rep(t_grid, times = n)
dat$start <- rep(c(TRUE, rep(FALSE, 10)), times = n)
S <- matrix(rep(s_grid, times = n * K), ncol = K, byrow = TRUE)
Tmat <- matrix(dat$t_out, nrow = n * K, ncol = K)
L <- F0_centered[rep(1:n, each = K), ]
L <- sweep(L, 2, weights, "*")
dat$S <- I(S)
dat$Tmat <- I(Tmat)
dat$L <- I(L)

# 3. Original FoF formula: the tensor term represents beta(s,t).
fof_formula <- Y ~ DiphVar.ord + sex.ord + type.ord + order.ord + leftSegment.ord +
  s(t_out, bs = "cr", k = 11, m = 1) +
  s(durationZ, bs = "cr", k = 11, m = 1) +
  s(t_out, by = DiphVar.ord, bs = "cr", k = 11, m = 1) +
  s(t_out, by = sex.ord, bs = "cr", k = 11, m = 1) +
  s(t_out, by = type.ord, bs = "cr", k = 11, m = 1) +
  s(t_out, by = order.ord, bs = "cr", k = 11, m = 1) +
  s(t_out, by = leftSegment.ord, bs = "cr", k = 11, m = 1) +
  te(S, Tmat, by = L, bs = c("cr", "cr"), k = c(11, 11)) +
  s(t_out, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1)

# 4. First pass without AR(1); itsadug estimates rho; second pass uses AR(1).
fof_noAR <- bam(fof_formula, data = dat, method = "fREML", discrete = TRUE)
rho_fof <- itsadug::start_value_rho(fof_noAR)
fof_model <- update(fof_noAR, rho = rho_fof, AR.start = dat$start)

# 5. Identify the coefficient block for the FoF surface.
for (sm in fof_model$smooth) {
  if (sm$label == "te(S,Tmat):L") surface_smooth <- sm
}
cols <- surface_smooth$first.para:surface_smooth$last.para
b <- coef(fof_model)[cols]
V <- vcov(fof_model)[cols, cols]
surface <- expand.grid(input_time = s_grid, output_time = t_grid)
surface$beta <- NA_real_
surface$se <- NA_real_

# 6. At each (s,t), use a unit input selector and the lpmatrix.
for (i in 1:nrow(surface)) {
  s <- surface$input_time[i]
  t <- surface$output_time[i]
  nd <- dat[1, ]
  nd$t_out <- t
  nd$S <- I(matrix(s_grid, nrow = 1))
  nd$Tmat <- I(matrix(t, nrow = 1, ncol = K))
  selector <- matrix(0, nrow = 1, ncol = K)
  selector[1, s + 1] <- 1   # No quadrature weight: extract beta itself.
  nd$L <- I(selector)
  X <- predict(fof_model, nd, type = "lpmatrix")
  X <- X[, cols, drop = FALSE]
  surface$beta[i] <- as.numeric(X %*% b)
  surface$se[i] <- sqrt(max(0, as.numeric(X %*% V %*% t(X))))
}

# 7. Display the 11 x 11 coefficient surface.
surface_plot <- ggplot(surface, aes(input_time, output_time, fill = beta)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  scale_x_continuous(breaks = 0:10, labels = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = 0:10, labels = seq(0, 1, 0.1)) +
  labs(x = "Input f0 time", y = "Output F1 time", fill = "beta") +
  theme_classic()
