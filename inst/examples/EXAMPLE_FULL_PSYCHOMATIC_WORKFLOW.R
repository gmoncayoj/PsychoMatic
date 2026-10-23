# PsychoMatic 1.0.0: compact end-to-end psychometric workflow

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

reverse_keys <- attr(ordinal_items, "reverse_keys")
domains <- attr(ordinal_items, "factor_structure")
items <- PsychoMatic::reverse_items(
  ordinal_items,
  items = reverse_keys,
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
scores <- PsychoMatic::score_scale(
  specification,
  keys = domains,
  min_valid = 0.80
)
reliability <- PsychoMatic::reliability_auto(
  items,
  items = domains$factor1,
  level = "ordinal",
  language = "eng"
)
efa <- PsychoMatic::efa_auto(
  specification,
  refine = "suggest",
  parallel_iter = 100,
  seed = 2401,
  language = "eng"
)
utils::head(scores$scores)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("Skipped report export: optional package 'openxlsx' is not installed.")
} else {
  report_path <- tempfile("psychomatic-workflow-", fileext = ".xlsx")
  PsychoMatic::psychomatic_report(
    list(reliability = reliability, efa = efa),
    format = "excel",
    file_name = report_path,
    language = "eng"
  )
  message("Temporary report: ", report_path)
}
