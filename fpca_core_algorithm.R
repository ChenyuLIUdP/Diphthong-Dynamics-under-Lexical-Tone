# fPCA algorithm illustration: common bases across both origins and both varieties.
# Required inputs / variable definitions
# F0 / F1: wide-format matrices, complete normalized curves, one token per row, 11 time columns (0:10).
# F0 uses interpolated f0_10Z; F1 uses f1Z. Rows align exactly with meta.
# meta: token_id, speaker_id, origin, variety_status, tone, sex, type, order,
# leftSegment, diphthong, duration_z, ToneOriginVar, cell_label.
# type = carrier sentence; order = reading order; leftSegment = onset.
# duration_z = standardized duration; ToneOriginVar = tone x origin x variety_status;
# cell_label identifies the four origin x variety_status cells.
# Categorical columns retain the original ordinary factors and reference levels.
# This presents the primary pooled /ai/ + /au/ algorithm.

library(fda)
library(lmerTest)
library(emmeans)

# 1. Separate fPCAs, using the original basis, smoothing and 99% retention rule.
basis <- create.bspline.basis(rangeval = c(0, 10), nbasis = 11, norder = 4)
fd_par <- fdPar(basis, Lfdobj = int2Lfd(2), lambda = 1e-6)
curves <- list(F0 = F0, F1 = F1)
fpca <- list()
for (signal in c("F0", "F1")) {
  Y <- curves[[signal]]
  smoothed <- smooth.basis(0:10, t(Y), fd_par)
  nharm_max <- max(3, min(10, nrow(Y) - 1))
  pc <- pca.fd(smoothed$fd, nharm = nharm_max, centerfns = TRUE)
  n_keep <- which(cumsum(pc$values) / sum(pc$values) >= 0.99)[1]
  if (is.na(n_keep)) n_keep <- length(pc$values)
  n_keep <- min(max(3, n_keep), length(pc$values))
  fpca[[signal]] <- list(
    mu = as.numeric(eval.fd(0:10, pc$meanfd)),
    efunctions = eval.fd(0:10, pc$harmonics[1:n_keep]),
    scores = pc$scores[, 1:n_keep, drop = FALSE],
    evalues_full = pc$values
  )
}

# 2. Attach the first three scores used by the score and bridge regressions.
model_df <- meta
for (j in 1:3) {
  model_df[[paste0("F0_PC", j)]] <- fpca$F0$scores[, j]
  model_df[[paste0("F1_PC", j)]] <- fpca$F1$scores[, j]
}

# 3. Primary mixed models (same predictors for F1-PC1 and F1-PC2).
score_pc1 <- lmer(F1_PC1 ~ F0_PC1 + F0_PC2 + F0_PC3 + duration_z +
  tone * origin * variety_status + sex + type + order + leftSegment +
  diphthong + (1 | speaker_id), data = model_df, REML = FALSE)
score_pc2 <- update(score_pc1, F1_PC2 ~ .)

# 4. Full-versus-reduced comparison for the complete condition block.
# The reduced models omit tone, origin, variety_status and all of their
# interactions. The likelihood-ratio tests therefore assess the joint
# 15-parameter condition block rather than a tone-only effect.
reduced_pc1 <- lmer(F1_PC1 ~ F0_PC1 + F0_PC2 + F0_PC3 + duration_z +
  sex + type + order + leftSegment + diphthong + (1 | speaker_id),
  data = model_df, REML = FALSE)
reduced_pc2 <- update(reduced_pc1, F1_PC2 ~ .)

compare_condition_block <- function(reduced_model, full_model, outcome) {
  reduced_logLik <- logLik(reduced_model)
  full_logLik <- logLik(full_model)
  reduced_npar <- attr(reduced_logLik, "df")
  full_npar <- attr(full_logLik, "df")
  chisq <- 2 * (as.numeric(full_logLik) - as.numeric(reduced_logLik))
  df <- full_npar - reduced_npar
  data.frame(
    outcome = outcome,
    reduced_npar = reduced_npar,
    full_npar = full_npar,
    reduced_AIC = AIC(reduced_model),
    full_AIC = AIC(full_model),
    chisq = chisq,
    df = df,
    p_value = pchisq(chisq, df = df, lower.tail = FALSE)
  )
}

condition_block_lrt <- rbind(
  compare_condition_block(reduced_pc1, score_pc1, "F1_PC1"),
  compare_condition_block(reduced_pc2, score_pc2, "F1_PC2")
)

# 5. Duration contribution: companion OLS R-squared, as in the source script.
duration_without <- lm(F1_PC1 ~ F0_PC1 + F0_PC2 + F0_PC3 +
  tone * origin * variety_status + sex + type + order + leftSegment +
  diphthong, data = model_df)
duration_with <- update(duration_without, . ~ . + duration_z)
duration_without_pc2 <- update(duration_without, F1_PC2 ~ .)
duration_with_pc2 <- update(duration_with, F1_PC2 ~ .)
duration_delta_R2 <- c(
  PC1 = summary(duration_with)$r.squared - summary(duration_without)$r.squared,
  PC2 = summary(duration_with_pc2)$r.squared - summary(duration_without_pc2)$r.squared
)

# 6. Four adjusted score models used for the PC-shape visualization and
# relative tonal configuration.
adjusted_scores <- unique(model_df[c("ToneOriginVar", "cell_label", "tone")])
adjusted_score_models <- list()
for (outcome in c("F0_PC1", "F0_PC2", "F1_PC1", "F1_PC2")) {
  score_formula <- reformulate(c("ToneOriginVar", "sex", "type", "order",
    "leftSegment", "diphthong", "duration_z", "(1 | speaker_id)"), response = outcome)
  fit <- lmer(score_formula, data = model_df, REML = FALSE)
  adjusted_score_models[[outcome]] <- fit
  emm <- as.data.frame(emmeans(fit, ~ ToneOriginVar))
  adjusted_scores[[outcome]] <- emm$emmean[match(adjusted_scores$ToneOriginVar,
                                                  emm$ToneOriginVar)]
}

# 7. Score-weighted PC1 and PC2 contributions used for the adjusted time-shape curves.
adjusted_pc_shapes <- list()
for (i in seq_len(nrow(adjusted_scores))) {
  a <- adjusted_scores[i, ]
  for (signal in c("F0", "F1")) {
    for (pc_number in 1:2) {
      score_name <- paste0(signal, "_PC", pc_number)
      adjusted_pc_shapes[[length(adjusted_pc_shapes) + 1]] <- data.frame(
        cell_label = a$cell_label,
        tone = a$tone,
        signal = signal,
        PC = paste0("PC", pc_number),
        measurement.no = 0:10,
        contribution = a[[score_name]] * fpca[[signal]]$efunctions[, pc_number]
      )
    }
  }
}
adjusted_pc_shapes <- do.call(rbind, adjusted_pc_shapes)

# 8. Reconstruct adjusted trajectories with PC1-PC2.
reconstructed <- list()
for (i in seq_len(nrow(adjusted_scores))) {
  a <- adjusted_scores[i, ]
  reconstructed[[i]] <- data.frame(cell_label = a$cell_label, tone = a$tone,
    measurement.no = 0:10,
    F0 = fpca$F0$mu + as.vector(fpca$F0$efunctions[, 1:2] %*% c(a$F0_PC1, a$F0_PC2)),
    F1 = fpca$F1$mu + as.vector(fpca$F1$efunctions[, 1:2] %*% c(a$F1_PC1, a$F1_PC2)))
}
relative_configuration <- do.call(rbind, reconstructed)

# 9. Center across four tones at each time, then correlate ordering and spacing.
relative_configuration$F0 <- with(relative_configuration,
  F0 - ave(F0, cell_label, measurement.no, FUN = mean))
relative_configuration$F1 <- with(relative_configuration,
  F1 - ave(F1, cell_label, measurement.no, FUN = mean))
relative_similarity <- unique(relative_configuration[c("cell_label", "measurement.no")])
for (i in seq_len(nrow(relative_similarity))) {
  dat <- subset(relative_configuration,
    cell_label == relative_similarity$cell_label[i] &
    measurement.no == relative_similarity$measurement.no[i])
  relative_similarity$rank_cor[i] <- cor(dat$F0, dat$F1, method = "spearman")
  relative_similarity$level_cor[i] <- cor(dat$F0, dat$F1, method = "pearson")
}
relative_summary <- aggregate(cbind(rank_cor, level_cor) ~ cell_label,
                              data = relative_similarity, FUN = mean)

# 10. Origin-specific bridges: three regressions in the SAME common PC bases.
bridge_by_origin <- list()
bridge_summary <- list()
for (origin_value in levels(model_df$origin)) {
  dat <- subset(model_df, origin == origin_value)
  bridge_pc1 <- lm(F1_PC1 ~ F0_PC1 + F0_PC2 + F0_PC3, data = dat)
  bridge_pc2 <- update(bridge_pc1, F1_PC2 ~ .)
  bridge_pc3 <- update(bridge_pc1, F1_PC3 ~ .)
  dat$F1_PC1_bridged <- predict(bridge_pc1, dat)
  dat$F1_PC2_bridged <- predict(bridge_pc2, dat)
  dat$F1_PC3_bridged <- predict(bridge_pc3, dat)
  bridge_by_origin[[origin_value]] <- dat
  # Discrepancy is computed between condition-mean scores, as in the source.
  means <- aggregate(cbind(F1_PC1, F1_PC2, F1_PC1_bridged, F1_PC2_bridged) ~
    ToneOriginVar, data = dat, FUN = mean)
  distance <- sqrt((means$F1_PC1 - means$F1_PC1_bridged)^2 +
                   (means$F1_PC2 - means$F1_PC2_bridged)^2)
  bridge_summary[[origin_value]] <- c(
    R2_PC1 = summary(bridge_pc1)$r.squared,
    R2_PC2 = summary(bridge_pc2)$r.squared,
    R2_PC3 = summary(bridge_pc3)$r.squared,
    mean_distance_PC12 = mean(distance)
  )
}

# 11. Separate POOLED bridge for the adjusted direct-versus-bridged comparison.
bridge_emmeans <- list()
bridge_contrasts <- list()
for (outcome in c("F1_PC1", "F1_PC2", "F1_PC3")) {
  bridge_formula <- reformulate(c("F0_PC1", "F0_PC2", "F0_PC3"), response = outcome)
  pooled_bridge <- lm(bridge_formula, data = model_df)
  direct <- model_df
  direct$score <- model_df[[outcome]]
  direct$method <- "direct"
  bridged <- model_df
  bridged$score <- predict(pooled_bridge, model_df)
  bridged$method <- "bridged"
  paired <- rbind(direct, bridged)
  paired$method <- factor(paired$method, levels = c("direct", "bridged"))
  # Original paired model: REML, no sex term, bobyqa optimizer.
  paired_fit <- lmer(score ~ method * ToneOriginVar + leftSegment + order +
    type + diphthong + (1 | speaker_id) + (1 | token_id), data = paired,
    REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  bridge_emmeans[[outcome]] <- emmeans(paired_fit, ~ method | ToneOriginVar)
  bridge_contrasts[[outcome]] <- contrast(bridge_emmeans[[outcome]], method = "pairwise")
}
