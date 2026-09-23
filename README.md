# Core analysis algorithms and supplementary documents

This [GitHub repository](https://github.com/ChenyuLIUdP/Diphthong-Dynamics-under-Lexical-Tone) accompanies *Diphthong Dynamics under Lexical Tone: Cross-Dialectal Evidence for Category-Specific f0–F1 Coupling*. It contains two supplementary documents and three R scripts illustrating the core analyses.

## Repository contents

| File | Contents |
| --- | --- |
| [SuppPub1.html](https://ChenyuLIUdP.github.io/Diphthong-Dynamics-under-Lexical-Tone/SuppPub1.html) | Supplementary methods, data-screening summaries, GAMM, fPCA and FoF results, and robustness analyses. |
| [SuppPub2.html](https://ChenyuLIUdP.github.io/Diphthong-Dynamics-under-Lexical-Tone/SuppPub2.html) | Example data frames, detailed GAMM summaries, diagnostics and model comparisons. |
| [gamm_core_algorithm.R](gamm_core_algorithm.R) | Baseline and mechanistic GAMMs, cascade predictions and reconstruction RMSE summaries. |
| [fpca_core_algorithm.R](fpca_core_algorithm.R) | Common-basis fPCA, score models, bridge comparisons and relative tonal configuration. |
| [fof_core_algorithm.R](fof_core_algorithm.R) | Function-on-function regression and coefficient-surface extraction. |

## Inputs and analysis scope

### GAMM

The analysis is performed separately for each speaker origin. All three models use one long-format table, `dat`, with one row per token and time point, sorted by token and time (0–10). A separate prediction grid, `nd`, contains the realized covariate combinations and condition-cell mean duration.

`ToneVarDiphInt.ord` represents tone × variety × diphthong, and `DiphVar.ord` represents diphthong × variety. Variables ending in `.ord` are ordered factors with treatment contrasts; speaker and lexical-item identifiers are ordinary factors.

The script fits F1 and *f*0 baseline models and a mechanistic F1 model. For each model, `itsadug::start_value_rho()` estimates rho from a first fit without AR(1), followed by a second fit with AR(1). Population-level direct and cascade F1 predictions exclude random smooths. RMSE between the direct and cascade F1 predictions is calculated within each full condition cell, averaged equally across realized onset, carrier-sentence and reading-order combinations, and then averaged equally across sexes and diphthongs for each tone and variety.

### fPCA

Inputs are complete wide-format matrices, `F0` and `F1`, with one token per row and 11 time columns, plus an aligned token-level metadata table, `meta`. Both origins, both varieties and both diphthongs are pooled. Separate fPCAs establish common bases for *f*0 and F1, retaining at least three PCs under the 99% variance rule. Categorical predictors use ordinary factors.

The script includes score mixed models, comparative analysis between full and reduced models, auxiliary duration analyses, origin-specific bridges for fit and condition-mean discrepancy summaries, and a separate pooled bridge for paired direct-versus-bridged comparisons. Covariate-adjusted PC1 and PC2 marginal means are multiplied by their corresponding eigenfunctions to obtain score-weighted component curves over time. PC1–PC2 reconstructions of *f*0 and F1 are then centered across tones at each time point within each origin × variety-status cell to compare relative tonal configurations.

### FoF

The analysis is performed separately for each speaker origin. Inputs are complete wide-format matrices, `F0` and `F1`, with one token per row and 11 time columns, plus aligned metadata, `meta`. The script expands F1 into long format while retaining the full *f*0 input curve for each output time point. Ordered factors use treatment contrasts.

The matrix-valued tensor term approximates the functional linear effect ∫ xᵢ(s) β(s,t) ds, where s indexes input *f*0 time and t indexes output F1 time. Here, xᵢ(s) is the *f*0 input curve centered at each input time point by subtracting the across-token mean within the speaker-origin subset. The matrix inputs are:

- `S`: input-time matrix; each row contains the 11 *f*0 time points (0–10).
- `Tmat`: output-time matrix; each row repeats its F1 output time across all 11 columns.
- `L`: the token's centered *f*0 curve multiplied by trapezoidal integration weights, repeated for each output time. The model term `te(S, Tmat, by = L)` sums the weighted contributions across input times.

The script includes input centering, trapezoidal integration weights, a first fit without AR(1), rho estimation with `itsadug::start_value_rho()`, and a second fit with AR(1). The coefficient surface and its standard errors are extracted using the model's `lpmatrix`.

## Required R packages

| Script | Packages |
| --- | --- |
| GAMM | `mgcv`, `itsadug` |
| fPCA | `fda`, `lmerTest`, `emmeans` |
| FoF | `mgcv`, `itsadug`, `ggplot2` |
