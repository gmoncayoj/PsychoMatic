# Historical identified ordinal sequences are still required by the partial
# invariance machinery. Exercise those private primitives separately from the
# two public workflows, whose current behavior is covered by the split tests.
legacy_ordinal_identification <- function(
  data, group, model, ordered = NULL, estimator = "auto", report = FALSE,
  missing = "auto"
) {
  .ifa_auto(data, group, model, language = "eng", estimator = estimator,
    ordered = ordered, report = report, criteria = "auto", missing = missing,
    min_group_n = 20L, overwrite = FALSE, digits = 3, file_name = NULL,
    format = "excel", verbose = FALSE, group_expression = substitute(group),
    group_environment = parent.frame(), original_call = match.call(),
    workflow = "legacy")
}
