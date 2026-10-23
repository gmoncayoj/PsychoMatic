# PsychoMatic 1.0.0: researcher-directed partial invariance

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

model <- paste(
  "factor1 =~ i1 + i2 + i3",
  "factor2 =~ i4 + i5 + i6",
  "factor3 =~ i7 + i8 + i9",
  sep = "\n"
)

if (!requireNamespace("semTools", quietly = TRUE)) {
  message("Skipped partial invariance: optional package 'semTools' is not installed.")
} else {
  partial <- PsychoMatic::invariance_partial(
    ordinal_groups,
    group = "group",
    model = model,
    free = "i3 | t3",
    rationale = "The simulation preregisters a group-B shift in the third i3 threshold.",
    stage = "thresholds",
    workflow = "ordinal",
    measurement = "ordinal",
    estimator = "WLSMV",
    ordered = paste0("i", 1:9),
    min_anchor_items = 2,
    report = FALSE,
    language = "eng"
  )
}
