# PsychoMatic 1.0.0: concise continuous CFA

continuous_items <- lavaan::simulateData(
  "general =~ 0.7*item1 + 0.7*item2 + 0.7*item3 + 0.7*item4 + 0.7*item5 + 0.7*item6",
  sample.nobs = 200, seed = 2401
)

model <- "general =~ item1 + item2 + item3 + item4 + item5 + item6"
cfa <- PsychoMatic::cfa_auto(
  continuous_items,
  model = model,
  estimator = "ML",
  language = "eng"
)
