test_that("the public interface excludes developer utilities and datasets", {
  exports <- getNamespaceExports("PsychoMatic")
  expect_false(any(c("check_items", "dimensionality_sensitivity", "audit_report",
                    "likert_demo", "invariance_demo", "psychomatic_continuous",
                    "psychomatic_ordinal", "psychomatic_multigroup") %in% exports))
  for (name in c("check_items", "dimensionality_sensitivity", "audit_report")) {
    expect_true(is.function(get(name, asNamespace("PsychoMatic"))))
  }
  catalog <- utils::data(package = "PsychoMatic")$results
  expect_length(intersect(catalog[, "Item"], c(
    "likert_demo", "invariance_demo", "psychomatic_continuous",
    "psychomatic_ordinal", "psychomatic_multigroup"
  )), 0L)
})
