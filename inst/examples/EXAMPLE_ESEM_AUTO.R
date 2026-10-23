# Ready-to-use templates: replace MY_DATA with your data frame.
# Run the alternative matching your data and research question.

# Geomin with ordinal items in columns 2:13; column 1 is not an item.
result <- PsychoMatic::esem_auto(
  data = MY_DATA,
  n_factors = 3,
  items = 2:13,
  rotation = "geomin",
  ordered = TRUE,
  language = "esp"
)
result$manuscript_text

# Continuous indicators, when MY_DATA contains only item columns.
result <- PsychoMatic::esem_auto(
  data = MY_DATA,
  n_factors = 3,
  items = names(MY_DATA),
  rotation = "geomin",
  ordered = FALSE,
  language = "esp"
)

# Target: replace these item names/memberships with your theoretical structure.
factors <- list(
  Factor1 = c("R1", "R2", "R3", "R4"),
  Factor2 = c("R5", "R6", "R7", "R8"),
  Factor3 = c("R9", "R10", "R11", "R12")
)
result <- PsychoMatic::esem_auto(
  data = MY_DATA,
  factors = factors,
  rotation = "target",
  ordered = TRUE,
  language = "esp"
)

# Optional exports:
# PsychoMatic::psychomatic_report(result, format = c("word", "excel"), language = "esp")
