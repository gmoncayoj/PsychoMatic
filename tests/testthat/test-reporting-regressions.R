reporting_regression_result <- function() {
  fit <- data.frame(metric = c("CFI", "RMSEA"), value = c(0.956789, 0.041234))
  loadings <- data.frame(item = c("i1", "i2"), loading = c(0.712345, 0.823456))
  out <- PsychoMatic:::.pm_result(
    class = "cfa_auto",
    call = quote(cfa_auto(items, model = "f =~ i1 + i2")),
    requested_settings = list(estimator = "auto"),
    effective_settings = list(estimator = "MLR", language = "eng"),
    data_summary = list(n = 250L, n_items = 2L),
    results = list(fit = fit, loadings = loadings),
    warnings = NULL,
    errors = NULL,
    diagnostics = list(fit_variants = data.frame(
      metric = fit$metric,
      base = fit$value,
      scaled = fit$value + c(0.001, -0.001)
    )),
    decision_log = PsychoMatic:::.pm_decision(
      "estimator", "auto", "MLR", "non-normality", "information"
    ),
    conclusion = list(status = "interpretable", reason = NA_character_)
  )
  out$tables <- list(fit = fit, loadings = loadings)
  out$fit <- fit
  out$loadings <- loadings
  out
}

test_that("report collectors reject every unsupported positional member", {
  objects <- list(
    good = data.frame(value = 1),
    bad = 42
  )

  expect_error(
    PsychoMatic:::.pm_collect_report_tables(objects),
    "index 2|object 2|unsupported",
    ignore.case = TRUE
  )

  path <- tempfile("psychomatic-invalid-collection-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  expect_error(
    psychomatic_report(objects, format = "excel", file_name = path),
    "index 2|object 2|unsupported",
    ignore.case = TRUE
  )
  expect_false(file.exists(path))
})

test_that("common result adapters select canonical tables without duplicates", {
  object <- reporting_regression_result()
  before <- serialize(object, NULL)
  tables <- PsychoMatic:::.pm_object_tables(object)
  fit <- object$tables$fit
  loadings <- object$tables$loadings

  expect_identical(
    sum(vapply(tables, identical, logical(1), y = fit)),
    1L
  )
  expect_identical(
    sum(vapply(tables, identical, logical(1), y = loadings)),
    1L
  )
  expect_identical(serialize(object, NULL), before)
})

test_that("Excel digits are presentation-only numeric formats", {
  skip_if_not_installed("openxlsx")
  raw <- pi
  table <- data.frame(value = raw, count = 3L)
  path <- tempfile("psychomatic-digits-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)

  PsychoMatic:::.pm_write_excel_report(
    list(values = table),
    path,
    digits = 2L
  )

  reopened <- openxlsx::read.xlsx(path)
  expect_true(is.numeric(reopened$value))
  expect_equal(reopened$value[[1L]], raw, tolerance = 1e-14)
  styles <- openxlsx::getStyles(openxlsx::loadWorkbook(path))
  formats <- vapply(styles, function(style) {
    value <- style$numFmt$formatCode
    if (is.null(value)) "" else as.character(value)
  }, character(1L))
  expect_true("0.00" %in% formats)
  expect_true("0" %in% formats)
})

test_that("export_cfa rejects duplicate formats before writing", {
  result <- reporting_regression_result()
  base <- tempfile("psychomatic-cfa-duplicate-")
  paths <- paste0(base, c(".xlsx", ".docx"))
  on.exit(unlink(paths), add = TRUE)

  expect_error(
    export_cfa(
      result,
      format = c("excel", "excel"),
      file_name = base
    ),
    "duplicate",
    ignore.case = TRUE
  )
  expect_false(any(file.exists(paths)))
})

test_that("export_cfa uses the common mapped writer and returns paths", {
  skip_if_not_installed("openxlsx")
  result <- reporting_regression_result()
  before <- serialize(result, NULL)
  path <- tempfile("psychomatic-cfa-common-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)

  visible <- withVisible(export_cfa(
    result,
    format = "excel",
    file_name = path,
    digits = 2L
  ))

  expect_false(visible$visible)
  expect_identical(as.character(visible$value), normalizePath(path, winslash = "/"))
  expect_s3_class(attr(visible$value, "sheet_map", exact = TRUE), "data.frame")
  expect_identical(
    attr(visible$value, "sheet_map", exact = TRUE)$effective_sheet,
    openxlsx::getSheetNames(path)
  )
  expect_identical(serialize(result, NULL), before)
})
