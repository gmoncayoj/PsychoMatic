# PsychoMatic 1.0.0: complete measurement-invariance sequence

grouped_items <- lavaan::simulateData(
  "general =~ 0.7*mg1 + 0.7*mg2 + 0.7*mg3 + 0.7*mg4 + 0.7*mg5 + 0.7*mg6",
  sample.nobs = 400, seed = 2402
)
grouped_items$group <- factor(rep(c("A", "B"), each = 200))

model <- "general =~ mg1 + mg2 + mg3 + mg4 + mg5 + mg6"
complete_invariance <- PsychoMatic::factorial_invariance_auto(
  grouped_items,
  group = "group",
  model = model,
  estimator = "ML",
  report = FALSE,
  language = "eng"
)
