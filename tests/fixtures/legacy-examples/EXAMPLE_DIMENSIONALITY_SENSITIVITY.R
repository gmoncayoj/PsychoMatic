# PsychoMatic 0.4.0: sensitivity of factor-retention decisions

set.seed(2401)
latent <- matrix(stats::rnorm(300 * 3), ncol = 3)
ordinal_items <- as.data.frame(sapply(seq_len(12), function(j) {
  as.integer(cut(0.8 * latent[, ceiling(j / 4)] + stats::rnorm(300, sd = 0.6),
                 breaks = c(-Inf, -1, -0.3, 0.3, 1, Inf), labels = FALSE))
}))
names(ordinal_items) <- sprintf("item%02d", 1:12)
attr(ordinal_items, "reverse_keys") <- c("item03", "item07", "item11")
ordinal_items[attr(ordinal_items, "reverse_keys")] <-
  lapply(ordinal_items[attr(ordinal_items, "reverse_keys")], function(x) 6L - x)
attr(ordinal_items, "factor_structure") <- list(
  factor1 = names(ordinal_items)[1:4], factor2 = names(ordinal_items)[5:8],
  factor3 = names(ordinal_items)[9:12]
)

items <- PsychoMatic::reverse_items(
  ordinal_items,
  items = attr(ordinal_items, "reverse_keys"),
  min_score = 1,
  max_score = 5,
  overwrite = TRUE
)
specification <- PsychoMatic::psychomatic_spec(
  items,
  level = "ordinal",
  missing = "available",
  language = "eng"
)
sensitivity <- PsychoMatic::dimensionality_sensitivity(
  specification,
  correlations = "auto",
  quantiles = c(0.90, 0.95),
  n_iter = 100,
  seed = 2401,
  language = "eng"
)

print(sensitivity)
summary(sensitivity)$configurations

plot_path <- tempfile("psychomatic-dimensionality-", fileext = ".pdf")
grDevices::pdf(plot_path)
plot(sensitivity)
grDevices::dev.off()
