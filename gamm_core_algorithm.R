# GAMM algorithm illustration: fit separately for Tianjin and Xi'an origins.
# Required inputs / variable definitions
# dat: long-format data, one row per token and time point.
# data sorted by token and measurement.no (0:10).
# f1Z / f0_10Z = within-speaker z-scored F1 / smoothed f0; durationZ = standardized duration.
# speaker / character = ordinary factors identifying speaker / lexical item.
# ToneVarDiphInt.ord = tone x variety x diphthong (e.g. "1 SM aj").
# DiphVar.ord = diphthong x variety (e.g. "aj SM"); variety.ord = SM or local variety.
# sex.ord, leftSegment.ord, type.ord, order.ord = sex, onset, carrier, reading order.
# All .ord columns MUST be ordered factors with original reference levels and
# treatment contrasts, e.g. contrasts(dat$sex.ord) <- "contr.treatment".
# start: TRUE at the first retained row of each token, FALSE otherwise. 
# (Based on retained rows included in the model, might be different between models.)
# nd: prepared prediction grid with original factor levels, all realized covariate
# cells, time points, and condition-cell mean durationZ; IDs are placeholders.


library(mgcv)
library(itsadug)

# Keep the preprocessed reference-level ordering; set treatment contrasts explicitly.
for (column in grep("[.]ord$", names(dat), value = TRUE)) {
  stopifnot(is.ordered(dat[[column]]))
  contrasts(dat[[column]]) <- "contr.treatment"
}

# 1. Original baseline formula (the f0 baseline has the same right-hand side).
baseline_formula <- f1Z ~ ToneVarDiphInt.ord + order.ord + sex.ord + leftSegment.ord + type.ord + 

    s(measurement.no, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = ToneVarDiphInt.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = order.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = sex.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = leftSegment.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = type.ord, bs = "cr", k = 11, m = 1) + 

    s(measurement.no, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1) + 
    s(measurement.no, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1, by = ToneVarDiphInt.ord) + 

    s(measurement.no, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1) + 
    s(measurement.no, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1, by = sex.ord) + 
    s(measurement.no, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1, by = order.ord) + 
    s(measurement.no, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1, by = type.ord) + 
    s(measurement.no, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1, by = variety.ord)

baseline_f0_formula <- update(baseline_formula, f0_10Z ~ .)

# 2. Original mechanistic formula. 
mechanistic_formula <- f1Z ~ DiphVar.ord + sex.ord + type.ord + order.ord + leftSegment.ord + 

    s(measurement.no, bs = "cr", k = 11, m = 1) + 
    s(f0_10Z, bs = "cr", k = 11, m = 1) + 
    s(durationZ, bs = "cr", k = 11, m = 1) + 

    ti(measurement.no, f0_10Z, k = c(11, 11)) + 
    ti(measurement.no, durationZ, k = c(11, 11)) + 
    ti(f0_10Z, durationZ, k = c(11, 11)) + 
    ti(measurement.no, durationZ, f0_10Z, k = c(11, 11, 11)) + 

    ti(measurement.no, f0_10Z, by = DiphVar.ord, k = c(11, 11)) + 
    ti(measurement.no, f0_10Z, by = sex.ord, k = c(11, 11)) + 

    s(measurement.no, by = DiphVar.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = type.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = order.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = sex.ord, bs = "cr", k = 11, m = 1) + 
    s(measurement.no, by = leftSegment.ord, bs = "cr", k = 11, m = 1) + 

    s(f0_10Z, by = DiphVar.ord, bs = "cr", k = 11, m = 1) + 
    s(f0_10Z, by = sex.ord, bs = "cr", k = 11, m = 1) + 

    s(measurement.no, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1) + 
    s(measurement.no, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1, by = DiphVar.ord) + 

    s(measurement.no, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1) + 

    s(f0_10Z, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1) + 
    s(f0_10Z, character, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1) + 

    s(durationZ, speaker, bs = "fs", xt = list(bs = "tp"), k = 11, m = 1)

# 3. First-pass fits without AR(1), then estimate rho from their residuals.
f1_noAR <- bam(baseline_formula, data = dat, method = "fREML", discrete = TRUE)
f0_noAR <- bam(baseline_f0_formula, data = dat, method = "fREML", discrete = TRUE)
mech_noAR <- bam(mechanistic_formula, data = dat, method = "fREML", discrete = TRUE)
rho_baseline_f1 <- itsadug::start_value_rho(f1_noAR)
rho_baseline_f0 <- itsadug::start_value_rho(f0_noAR)
rho_mechanistic <- itsadug::start_value_rho(mech_noAR)

# 4. Second-pass AR(1) fits, with exactly the same formulas. 
# For start, TRUE should be at the first retained row of each token. 
# Since baseline_f0 and mechanistic_f1 have f0 values in the models, 
# the start should be based on retained observations (with f0 values) of each token.
baseline_f1 <- update(f1_noAR, rho = rho_baseline_f1, AR.start = dat$start)
baseline_f0 <- update(f0_noAR, rho = rho_baseline_f0, AR.start = dat$start)
mechanistic_f1 <- update(mech_noAR, rho = rho_mechanistic, AR.start = dat$start)

# 5. Population-level direct and cascade predictions (exclude random smooths).
random_terms_f1 <- character()
for (sm in baseline_f1$smooth) {
  if (grepl("speaker|character", sm$label)) random_terms_f1 <- c(random_terms_f1, sm$label)
}
random_terms_f0 <- character()
for (sm in baseline_f0$smooth) {
  if (grepl("speaker|character", sm$label)) random_terms_f0 <- c(random_terms_f0, sm$label)
}
random_terms_mech <- character()
for (sm in mechanistic_f1$smooth) {
  if (grepl("speaker|character", sm$label)) random_terms_mech <- c(random_terms_mech, sm$label)
}
direct_f1 <- predict(baseline_f1, nd, exclude = random_terms_f1)
predicted_f0 <- predict(baseline_f0, nd, exclude = random_terms_f0)
nd_cascade <- nd
nd_cascade$f0_10Z <- as.numeric(predicted_f0)
cascade_f1 <- predict(mechanistic_f1, nd_cascade, exclude = random_terms_mech)

# 6. Reconstruction error: first calculate RMSE within each full condition cell.
curves <- nd
curves$squared_error <- (as.numeric(direct_f1) - as.numeric(cascade_f1))^2
cell_mse <- aggregate(squared_error ~ variety.ord + sex.ord + diphthong + tone +
                        leftSegment.ord + type.ord + order.ord,
                      data = curves, FUN = mean)
cell_mse$RMSE <- sqrt(cell_mse$squared_error)

# Equal mean of cell RMSEs, then equal mean across sex x diphthong for each tone.
condition_rmse <- aggregate(RMSE ~ variety.ord + sex.ord + diphthong + tone,
                            data = cell_mse, FUN = mean)
tone_rmse <- aggregate(RMSE ~ variety.ord + tone, data = condition_rmse, FUN = mean)
