# PsychoMatic 1.0.0: graded-response IRT and conservative DIF screening

set.seed(2402)
latent <- matrix(stats::rnorm(600 * 3), ncol = 3)
group <- factor(rep(c("A", "B"), each = 300))
ordinal_groups <- as.data.frame(sapply(seq_len(9), function(j) {
  response <- 0.8 * latent[, ceiling(j / 3)] + stats::rnorm(600, sd = 0.6)
  values <- integer(600)
  for (g in levels(group)) {
    thresholds <- c(-Inf, -1, -0.3, 0.3, 1, Inf)
    if (j == 3 && g == "B") thresholds[4] <- 0.8
    values[group == g] <- cut(response[group == g], thresholds, labels = FALSE)
  }
  values
}))
names(ordinal_groups) <- paste0("i", 1:9)
ordinal_groups[] <- lapply(ordinal_groups, ordered, levels = 1:5)
ordinal_groups$group <- group
grouped_items <- lavaan::simulateData(
  "general =~ 0.7*mg1 + 0.7*mg2 + 0.7*mg3 + 0.7*mg4 + 0.7*mg5 + 0.7*mg6",
  sample.nobs = 400, seed = 2402
)
grouped_items$group <- factor(rep(c("A", "B"), each = 200))

item_names <- paste0("i", 1:9)
numeric_items <- as.data.frame(lapply(
  ordinal_groups[item_names],
  as.integer
))

if (!requireNamespace("mirt", quietly = TRUE)) {
  message("Skipped IRT: optional package 'mirt' is not installed.")
} else {
  irt_items <- paste0("i", 1:3)
  irt_rows <- stats::complete.cases(numeric_items[irt_items])
  irt <- PsychoMatic::irt_auto(
    numeric_items[irt_rows, irt_items],
    model = "grm",
    language = "eng"
  )
  irt$model_summary
  utils::head(irt$item_parameters)
}

# A separate binary teaching fixture is derived from continuous simulated data.
# Do not dichotomize ordinal study items merely to make this screening applicable.
binary_items <- as.data.frame(lapply(
  grouped_items[paste0("mg", 1:6)],
  function(item) as.integer(item > stats::median(item))
))
dif <- PsychoMatic::dif_auto(
  binary_items,
  group = grouped_items$group,
  items = names(binary_items),
  polytomous = "not_evaluable",
  language = "eng"
)
dif$summary
utils::head(dif$interpretation)
