# PsychoMatic 1.0.0: specification, item checks, and scoring

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

item_names <- names(ordinal_items)
item_ranges <- stats::setNames(rep(list(c(1, 5)), length(item_names)), item_names)
specification <- PsychoMatic::psychomatic_spec(
  ordinal_items,
  items = item_names,
  level = "ordinal",
  reverse = attr(ordinal_items, "reverse_keys"),
  ranges = item_ranges,
  missing = "available",
  language = "eng"
)
scores <- PsychoMatic::score_scale(
  specification,
  keys = attr(ordinal_items, "factor_structure"),
  reverse = attr(ordinal_items, "reverse_keys"),
  min_score = 1,
  max_score = 5,
  min_valid = 0.80
)

summary(specification)$items
utils::head(scores$scores)
scores$reliability
