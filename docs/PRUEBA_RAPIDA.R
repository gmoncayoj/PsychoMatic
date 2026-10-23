# PsychoMatic: reproducible checks of the simplified workflow.
# Run each section separately in RStudio. No print()/summary() calls are needed.
library(PsychoMatic)
data <- lavaan::HolzingerSwineford1939
items <- data[paste0("x", 1:9)]
model <- 'visual =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed =~ x7 + x8 + x9'

desc_auto(items, digits = 4)
cormat(items, type = "pearson", digits = 4)
efa_auto(items, seed = 916, refine = "suggest")
cfa_auto(data, model)  # automatic estimator selection
cfa_auto(data, model, estimator = "MLR", fit_variant = "scaled")
cfa_auto(data, model, estimator = "MLR", fit_variant = "robust")
reliability_auto(items[1:3])
factorial_invariance_auto(data, "school", model)

# CFA-based reliability: save the object only for reuse.
# CFA and reliability print their results when run.
fit <- cfa_auto(data, model, estimator = "MLR")
reliability_auto(items, model = fit$fit)

# Export directly from each function; existing files are protected.
output_dir <- tempfile("PsychoMatic-tests-")
dir.create(output_dir)
efa_auto(items, seed = 916, refine = "suggest", report = TRUE,
         format = c("excel", "word"), file_name = file.path(output_dir, "EFA"))
cfa_auto(data, model, report = TRUE, format = "word",
         file_name = file.path(output_dir, "CFA.docx"))
reliability_auto(items[1:3], report = TRUE,
                 file_name = file.path(output_dir, "Reliability.xlsx"))
factorial_invariance_auto(data, "school", model, report = TRUE,
                         file_name = file.path(output_dir, "Invariance.xlsx"))

# Optional: use only when relevant to the research question.
dimensionality_sensitivity(items, n_iter = 100, seed = 916)
dif_auto(items, group = data$school)

# Partial invariance requires a deliberate release and justification.
# Teaching example, not a recommendation to modify this model:
invariance_partial(data, "school", model, free = "x2 ~ 1",
                   rationale = "Teaching example: assess sensitivity of the x2 intercept.",
                   estimator = "ML")

# For silent programmatic workflows:
result <- cfa_auto(data, model, verbose = FALSE)

# If an ordinal category is empty within a group, inspect its frequencies:
# table(my_data$R5, my_data$group, useNA = "ifany")
# Do not automatically remove or merge categories to force a fit.
