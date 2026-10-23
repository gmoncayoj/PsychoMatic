.audit_has_unsafe_value <- function(x, depth = 0L) {
  if (depth > 12L) return(TRUE)
  if (is.environment(x) || is.function(x) || isS4(x) ||
      typeof(x) == "externalptr" || inherits(x, "connection") ||
      inherits(x, c("lavaan", "SingleGroupClass", "mirt", "sirt", "psych", "fa"))) {
    return(TRUE)
  }
  if (!is.list(x)) return(FALSE)
  any(vapply(x, .audit_has_unsafe_value, logical(1L), depth = depth + 1L))
}

.audit_regression_result <- function() {
  fit <- data.frame(metric = c("CFI", "RMSEA"), value = c(0.95, 0.04))
  anchors <- data.frame(item = c("i1", "i2"), retained = c(TRUE, TRUE))
  out <- PsychoMatic:::.pm_result(
    class = "psychomatic_invariance_partial",
    call = quote(invariance_partial(data, group, model, free, rationale)),
    requested_settings = list(stage = "scalar"),
    effective_settings = list(stage = "scalar", language = "eng"),
    data_summary = list(n = 300L, groups = c(A = 150L, B = 150L)),
    results = list(
      fit_table = fit,
      anchors = list(items = anchors),
      comparison = list(fit = fit)
    ),
    warnings = NULL,
    errors = NULL,
    diagnostics = list(fit_variants = fit),
    decision_log = data.frame(code = "partial", decision = "retained"),
    conclusion = list(status = "interpretable", reason = NA_character_)
  )
  out$tables <- list(fit = fit, anchors = anchors)
  out
}

test_that("audit paths distinguish reserved literal names from positions", {
  value <- structure(list(1, 2, 3), names = c("[[2]]", "", "[[2]]"))
  flattened <- PsychoMatic:::.pm_audit_flatten(value)

  expect_identical(nrow(flattened), 3L)
  expect_identical(length(unique(flattened$path)), 3L)
  expect_true("[[2]]" %in% flattened$path)
  expect_true(any(grepl("\\[\\[2\\]\\]", flattened$path, fixed = TRUE)))
  expect_equal(sort(vapply(flattened$value_raw, identity, numeric(1L))), 1:3)
})

test_that("diagnostic data frames retain values as safe canonical tables", {
  object <- .audit_regression_result()
  expected <- object$diagnostics$fit_variants
  audit <- audit_report(object)

  diagnostic_row <- audit$diagnostics[
    audit$diagnostics$path == "fit_variants",
    ,
    drop = FALSE
  ]
  expect_identical(nrow(diagnostic_row), 1L)
  expect_false(diagnostic_row$truncated[[1L]])
  expect_identical(diagnostic_row$value_raw[[1L]], expected)
  diagnostic_tables <- grep(
    "diagnostics_fit_variants$",
    names(audit$results),
    value = TRUE
  )
  expect_length(diagnostic_tables, 1L)
  expect_identical(audit$results[[diagnostic_tables]], expected)
})

test_that("audit recursively removes unsafe live values from every section", {
  object <- .audit_regression_result()
  unsafe <- new.env(parent = emptyenv())
  object$tables$unsafe <- data.frame(id = 1L)
  object$tables$unsafe$payload <- I(list(unsafe))
  object$results$unsafe <- object$tables$unsafe
  object$diagnostics$unsafe <- object$tables$unsafe
  object$decision_log$payload <- I(list(unsafe))
  object$references <- data.frame(reference = "x")
  object$references$payload <- I(list(unsafe))
  object$warnings <- data.frame(
    code = "unsafe",
    message = "probe",
    stringsAsFactors = FALSE
  )
  object$warnings$condition <- I(list(unsafe))

  audit <- audit_report(object)

  expect_false(.audit_has_unsafe_value(audit))
  unsafe_tables <- audit$results[grepl("unsafe", names(audit$results))]
  expect_true(length(unsafe_tables) >= 1L)
  expect_true(all(vapply(unsafe_tables, function(table) {
    all(vapply(table$payload, function(value) {
      is.character(value) && length(value) == 1L && grepl("omitted", value)
    }, logical(1L)))
  }, logical(1L))))
})

test_that("audit uses canonical object tables and retains unique nested evidence", {
  object <- .audit_regression_result()
  before <- serialize(object, NULL)
  audit <- audit_report(object)
  labels <- names(audit$results)

  expect_identical(sum(labels == "result_001_fit"), 1L)
  expect_identical(sum(labels == "result_001_anchors"), 1L)
  expect_length(grep("diagnostics_fit_variants$", labels, value = TRUE), 1L)
  expect_true(any(grepl("comparison_fit$", labels)))
  expect_match(audit$objects$call_text, "invariance_partial", fixed = TRUE)
  expect_false(grepl("audit_report", audit$objects$call_text, fixed = TRUE))
  expect_identical(serialize(object, NULL), before)
})
