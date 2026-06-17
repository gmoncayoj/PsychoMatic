# PsychoMatic complete example workflow
# Copy, paste, and run this script after installing/loading PsychoMatic.
# The workflow uses bundled datasets and temporary files for exports.

library(PsychoMatic)

# -------------------------------------------------------------------------
# 0. Load bundled example datasets
# -------------------------------------------------------------------------
# psychomatic_ordinal: five-category ordinal items.
# psychomatic_continuous: continuous item responses.
# psychomatic_multigroup: continuous items plus a group variable.

data(psychomatic_ordinal)
data(psychomatic_continuous)
data(psychomatic_multigroup)

ordinal_items <- psychomatic_ordinal
continuous_items <- psychomatic_continuous
multigroup_items <- psychomatic_multigroup

# -------------------------------------------------------------------------
# 1. desc_auto(): automated item descriptives
# -------------------------------------------------------------------------
# Produces item-level descriptive statistics and response percentages.

desc_results <- desc_auto(
  ordinal_items,
  digits = 2,
  language = "eng",
  report = FALSE
)
print(desc_results)

# -------------------------------------------------------------------------
# 2. screen_items(): item-quality screening
# -------------------------------------------------------------------------
# Screens missingness, zero variance, floor/ceiling effects, item-rest
# correlations, and probable reverse-keyed items.

screen_results <- screen_items(ordinal_items)
print(screen_results)

# -------------------------------------------------------------------------
# 3. reverse_items(): reverse-score selected items
# -------------------------------------------------------------------------
# Use this when theory or item wording indicates reverse-keyed items.

reversed_items <- reverse_items(
  ordinal_items,
  items = "ord2",
  min_score = 1,
  max_score = 5,
  overwrite = FALSE
)
head(reversed_items)

# -------------------------------------------------------------------------
# 4. score_scale(): compute scale scores
# -------------------------------------------------------------------------
# Computes mean or sum scores with a minimum valid-response rule.

score_results <- score_scale(
  ordinal_items,
  keys = list(total = names(ordinal_items)),
  method = "mean",
  min_valid = 0.80,
  impute = "person_mean"
)
head(score_results$scores)
score_results$reliability

# -------------------------------------------------------------------------
# 5. cormat(): inter-item correlation matrix
# -------------------------------------------------------------------------
# Pearson is fast and always available for numeric items. Polychoric and
# tetrachoric options are also supported when appropriate.

pearson_cor <- cormat(
  continuous_items,
  type = "pearson",
  report = FALSE
)
print(round(pearson_cor, 3))

# -------------------------------------------------------------------------
# 6. efa_auto(): automated exploratory factor analysis
# -------------------------------------------------------------------------
# Runs estimator selection, parallel analysis, rotation, item refinement,
# reliability, and structured summaries.

efa_result <- efa_auto(
  ordinal_items,
  rotation = "oblique",
  max_iter = 3,
  verbose = FALSE,
  language = "eng"
)
summary(efa_result)

# -------------------------------------------------------------------------
# 7. export_efa(): export EFA results
# -------------------------------------------------------------------------
# Export to a temporary Excel file to avoid writing uncontrolled files to the
# working directory.

efa_xlsx <- tempfile("psychomatic-efa-", fileext = ".xlsx")
export_efa(
  efa_result,
  format = "excel",
  file_name = efa_xlsx
)
efa_xlsx

# -------------------------------------------------------------------------
# 8. compare_models_auto(): compare CFA model candidates with ML
# -------------------------------------------------------------------------
# ML workflows report fit indices plus AIC/BIC when available. The
# recommendation is based on fit indices, not only on information criteria.

cfa_models <- list(
  one_factor = "general =~ item1 + item2 + item3 + item4 + item5 + item6",
  two_factor = "
    factor1 =~ item1 + item2 + item3
    factor2 =~ item4 + item5 + item6
  "
)

model_comparison_ml <- compare_models_auto(
  continuous_items,
  models = cfa_models,
  estimator = "ML",
  language = "eng"
)
model_comparison_ml$fit
model_comparison_ml$recommendation

# -------------------------------------------------------------------------
# 9. compare_models_auto(): compare ordinal CFA models with WLSMV
# -------------------------------------------------------------------------
# WLSMV does not use AIC/BIC for model comparison. PsychoMatic compares scaled
# or robust fit indices when available: higher CFI/TLI and lower RMSEA/SRMR.

ordinal_cfa_models <- list(
  one_factor = "general =~ ord1 + ord2 + ord3 + ord4 + ord5 + ord6",
  two_factor = "
    factor1 =~ ord1 + ord2 + ord3
    factor2 =~ ord4 + ord5 + ord6
  "
)

model_comparison_wlsmv <- compare_models_auto(
  ordinal_items,
  models = ordinal_cfa_models,
  estimator = "WLSMV",
  language = "eng"
)
model_comparison_wlsmv$fit
model_comparison_wlsmv$recommendation

# -------------------------------------------------------------------------
# 10. cfa_auto(): automated confirmatory factor analysis
# -------------------------------------------------------------------------
# Adds estimator handling, fit interpretation, reliability, warnings, and
# structured summary tables.

cfa_model <- "
  factor1 =~ item1 + item2 + item3
  factor2 =~ item4 + item5 + item6
"

cfa_result <- cfa_auto(
  continuous_items,
  model = cfa_model,
  estimator = "ML",
  language = "eng"
)
summary(cfa_result)

# -------------------------------------------------------------------------
# 11. export_cfa(): export CFA results
# -------------------------------------------------------------------------
# Export to a temporary Excel file.

cfa_xlsx <- tempfile("psychomatic-cfa-", fileext = ".xlsx")
export_cfa(
  cfa_result,
  format = "excel",
  file_name = cfa_xlsx
)
cfa_xlsx

# -------------------------------------------------------------------------
# 12. factorial_invariance_auto(): multi-group measurement invariance
# -------------------------------------------------------------------------
# Fits configural, metric, scalar, and strict invariance models.

invariance_model <- "
  factor1 =~ mg1 + mg2 + mg3
  factor2 =~ mg4 + mg5 + mg6
"

invariance_result <- factorial_invariance_auto(
  multigroup_items,
  group = "group",
  model = invariance_model,
  estimator = "ML",
  language = "eng",
  report = FALSE
)
summary(invariance_result)

# -------------------------------------------------------------------------
# 13. inv_align_auto(): alignment invariance
# -------------------------------------------------------------------------
# Matrix input is useful when loading and intercept estimates are already
# available from prior models.

lambda <- matrix(
  c(0.70, 0.80, 0.75,
    0.72, 0.78, 0.74),
  nrow = 2,
  byrow = TRUE
)
nu <- matrix(
  c(0.00, 0.00, 0.00,
    0.10, 0.05, 0.08),
  nrow = 2,
  byrow = TRUE
)

alignment_result <- inv_align_auto(
  lambda = lambda,
  nu = nu,
  language = "eng"
)
summary(alignment_result)

# -------------------------------------------------------------------------
# 14. irt_auto(): basic IRT models with interpretation
# -------------------------------------------------------------------------
# Fits Rasch, 2PL, or graded response models through mirt and returns item-level
# interpretation for discrimination, difficulty/thresholds, and review flags.

irt_result <- irt_auto(
  ordinal_items[, 1:4],
  model = "grm",
  language = "eng"
)
irt_result$item_parameters
irt_result$model_summary
irt_result$interpretation
irt_result$references

# -------------------------------------------------------------------------
# 15. dif_auto(): differential item functioning screening with interpretation
# -------------------------------------------------------------------------
# Performs lightweight logistic DIF screening by group and returns item-level
# interpretation plus an overall summary.

dif_results <- dif_auto(
  multigroup_items[paste0("mg", 1:6)],
  group = multigroup_items$group,
  alpha = 0.01
)
dif_results$results
dif_results$summary
dif_results$interpretation
dif_results$references

# -------------------------------------------------------------------------
# 16. psychomatic_report(): integrated report
# -------------------------------------------------------------------------
# Creates an Excel report from available PsychoMatic result objects.

full_report <- tempfile("psychomatic-full-report-", fileext = ".xlsx")
psychomatic_report(
  descriptives = desc_results,
  screening = screen_results,
  scores = score_results,
  efa = efa_result,
  cfa = cfa_result,
  cfa_comparison_ml = model_comparison_ml,
  cfa_comparison_wlsmv = model_comparison_wlsmv,
  invariance = invariance_result,
  alignment = alignment_result,
  irt = irt_result,
  dif = dif_results,
  format = "excel",
  file_name = full_report,
  language = "eng"
)
full_report
