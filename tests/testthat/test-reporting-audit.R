reporting_probe <- function() {
  out <- PsychoMatic:::.pm_result(
    class = "psychomatic_probe",
    call = quote(probe(x)),
    requested_settings = list(method = "auto", digits = 2L),
    effective_settings = list(method = "spearman", dimensions = 2L),
    data_summary = list(
      n = 10L,
      p = 2L,
      variables = c("i1", "i2"),
      missing = "pairwise"
    ),
    results = list(value = 0.4123456789),
    warnings = data.frame(
      code = "fallback",
      message = "fallback used",
      stringsAsFactors = FALSE
    ),
    errors = data.frame(
      code = "polychoric_error",
      message = "polychoric failed",
      stringsAsFactors = FALSE
    ),
    diagnostics = list(
      original_error = "polychoric failed",
      nested = list(reason = "sparse cells")
    ),
    decision_log = PsychoMatic:::.pm_decision(
      "fallback", "polychoric", "spearman", "sparse cells", "warning"
    ),
    conclusion = list(
      status = "not_evaluable",
      reason = "The preferred correlation was not evaluable."
    )
  )
  out$references <- data.frame(
    reference = "Example reference",
    identifier = "doi:10.0000/example",
    stringsAsFactors = FALSE
  )
  out$versions <- list(PsychoMatic = "0.4.0", psych = "2.5.6")
  out$report_path <- "C:/tmp/probe.xlsx"
  out
}

test_that("audit_report has the exact public contract and preserves evidence", {
  expect_identical(names(formals(audit_report)), c("result", "language"))
  expect_null(eval(formals(audit_report)$language))

  probe <- reporting_probe()
  before <- serialize(probe, NULL)
  visible <- withVisible(audit_report(probe, language = "eng"))
  out <- visible$value

  expect_true(visible$visible)
  expect_identical(class(out)[1:2], c("psychomatic_audit", "psychomatic_result"))
  expect_identical(out$effective_settings$language, "eng")
  required <- c(
    "objects", "settings", "data", "warnings", "errors", "decisions",
    "diagnostics", "conclusions", "references", "versions", "artifacts",
    "results"
  )
  expect_true(all(required %in% names(out)))
  expect_true(all(vapply(
    out[setdiff(required, "results")],
    is.data.frame,
    logical(1)
  )))
  method <- out$settings[out$settings$path == "method", , drop = FALSE]
  expect_identical(method$requested, "auto")
  expect_identical(method$effective, "spearman")
  expect_true("fallback used" %in% out$warnings$message)
  expect_true("polychoric failed" %in% out$errors$message)
  expect_true("polychoric failed" %in% out$diagnostics$value_display)
  expect_identical(out$conclusions$status, "not_evaluable")
  expect_true("doi:10.0000/example" %in% out$references$identifier)
  expect_true("psych" %in% out$versions$package)
  expect_true("C:/tmp/probe.xlsx" %in% out$artifacts$path)
  expect_identical(serialize(probe, NULL), before)

  captured <- capture.output(returned <- print(out))
  expect_identical(returned, out)
  expect_true(length(captured) > 0L)
  summarized <- summary(out)
  expect_true(all(required %in% names(summarized)))
  expect_identical(summarized$conclusions, out$conclusions)
})

test_that("audit collections are positional and preserve duplicate labels", {
  probe <- reporting_probe()
  collection <- structure(
    list(probe, probe, probe),
    names = c("Analysis", "analysis", "")
  )
  before <- serialize(collection, NULL)
  out <- audit_report(collection)

  expect_identical(out$objects$object_id, sprintf("result_%03d", 1:3))
  expect_identical(out$objects$display_label[1:2], c("Analysis", "analysis"))
  expect_false(any(out$objects$label_generated[1:2]))
  expect_true(out$objects$label_generated[3L])
  expect_identical(out$objects$display_label[3L], "result_003")
  expect_identical(out$settings$object_id[match(1:3, out$settings$object_index)],
                   sprintf("result_%03d", 1:3))
  expect_identical(serialize(collection, NULL), before)

  single <- audit_report(probe)
  expect_identical(nrow(single$objects), 1L)
  expect_error(audit_report(data.frame(x = 1)), "unsupported|data.frame")
  expect_error(audit_report(list(bad = 1)), "unsupported|object 1|index 1")
})

test_that("audit extraction is bounded and empty sections stay typed", {
  probe <- reporting_probe()
  probe$warnings <- PsychoMatic:::.pm_empty_condition_table()
  probe$errors <- PsychoMatic:::.pm_empty_condition_table()
  probe$diagnostics$large <- as.list(seq_len(800L))
  probe$diagnostics$unsafe_environment <- new.env(parent = emptyenv())

  out <- audit_report(probe)
  expect_s3_class(out$warnings, "data.frame")
  expect_s3_class(out$errors, "data.frame")
  expect_identical(nrow(out$warnings), 0L)
  expect_identical(nrow(out$errors), 0L)
  expect_true(any(out$diagnostics$truncated))
  expect_false(any(vapply(
    out$diagnostics$value_raw,
    function(value) is.environment(value) || inherits(value, "externalptr"),
    logical(1)
  )))
})

test_that("nested audit paths are unambiguous", {
  probe <- reporting_probe()
  probe$diagnostics <- list("a.b" = 1, a = list(b = 2))
  out <- audit_report(probe)

  expect_true(all(c("a\\.b", "a.b") %in% out$diagnostics$path))
  literal <- out$diagnostics$value_raw[[match("a\\.b", out$diagnostics$path)]]
  nested <- out$diagnostics$value_raw[[match("a.b", out$diagnostics$path)]]
  expect_identical(literal, 1)
  expect_identical(nested, 2)
})

test_that("audit adapters cover legacy item screening and scoring objects", {
  items <- data.frame(
    i1 = c(1, 2, 3, 4, 5, 4),
    i2 = c(2, 2, 3, 4, 5, 3),
    i3 = c(1, 2, 2, 3, 4, 5)
  )
  screened <- screen_items(items)
  scored <- score_scale(items, items = names(items))
  screened_before <- serialize(screened, NULL)
  scored_before <- serialize(scored, NULL)

  screened_audit <- audit_report(screened, language = "eng")
  scored_audit <- audit_report(scored, language = "eng")

  expect_identical(screened_audit$objects$primary_class, "psychomatic_item_screen")
  expect_identical(scored_audit$objects$primary_class, "psychomatic_scores")
  expect_true(any(grepl("item_screen", names(screened_audit$results))))
  expect_true(any(grepl("scores", names(scored_audit$results))))
  expect_true(any(grepl("reliability", names(scored_audit$results))))
  expect_identical(serialize(screened, NULL), screened_before)
  expect_identical(serialize(scored, NULL), scored_before)
})

test_that("audit preserves scalar results and omits heavyweight live values", {
  probe <- reporting_probe()
  probe$effective_settings$language <- "eng"
  probe$conclusion <- NA_character_
  probe$results <- list(value = pi, fit = new.env(parent = emptyenv()))
  probe$diagnostics$callback <- function() NULL
  probe$warnings$unsafe_detail <- I(list(new.env(parent = emptyenv())))

  out <- audit_report(probe)
  expect_identical(out$effective_settings$language, "eng")
  expect_identical(out$conclusions$status, "evaluable")
  values <- out$results[[grep("values$", names(out$results), value = TRUE)[1L]]]
  scalar <- values[values$path == "value", , drop = FALSE]
  expect_identical(scalar$value_raw[[1L]], pi)
  expect_true(any(values$path == "fit" & values$truncated))
  expect_true(any(out$diagnostics$path == "callback" & out$diagnostics$truncated))
  raw_values <- c(out$diagnostics$value_raw, values$value_raw)
  expect_false(any(vapply(
    raw_values,
    function(value) is.environment(value) || is.function(value) || isS4(value),
    logical(1)
  )))
  warning_list_columns <- out$warnings[vapply(out$warnings, is.list, logical(1))]
  expect_false(any(vapply(
    unlist(warning_list_columns, recursive = FALSE),
    is.environment,
    logical(1)
  )))
})

test_that("audit language resolution is deterministic for one or many objects", {
  english <- reporting_probe()
  english$effective_settings$language <- "eng"
  spanish <- reporting_probe()
  spanish$effective_settings$language <- "esp"

  expect_identical(audit_report(english)$effective_settings$language, "eng")
  old <- options(PsychoMatic.language = "esp")
  on.exit(options(old), add = TRUE)
  expect_identical(
    audit_report(list(english = english, spanish = spanish))$effective_settings$language,
    "esp"
  )
  expect_identical(
    audit_report(list(english, spanish), language = "eng")$effective_settings$language,
    "eng"
  )
})

test_that("report collection is positional, lossless, and mapped", {
  first <- data.frame(value = 1.123456789, label = "first")
  second <- data.frame(value = 9.876543211, label = "second")
  objects <- structure(list(first, second), names = c("Analysis", "analysis"))
  before <- serialize(objects, NULL)

  tables <- PsychoMatic:::.pm_collect_report_tables(objects, digits = 2)
  expect_identical(length(tables), 2L)
  expect_identical(length(unique(tolower(names(tables)))), 2L)
  expect_equal(tables[[1L]]$value, first$value, tolerance = 0)
  expect_equal(tables[[2L]]$value, second$value, tolerance = 0)
  map <- attr(tables, "table_map", exact = TRUE)
  expect_s3_class(map, "data.frame")
  expect_identical(map$object_index, 1:2)
  expect_identical(map$display_label, c("Analysis", "analysis"))
  expect_identical(map$object_id, c("result_001", "result_002"))
  expect_identical(map$requested_sheet, c("Analysis_table", "analysis_table"))
  expect_identical(serialize(objects, NULL), before)

  probe_tables <- PsychoMatic:::.pm_collect_report_tables(
    list(probe = reporting_probe())
  )
  expect_true(all(vapply(
    c("settings", "warnings", "errors", "decisions", "diagnostics", "conclusions"),
    function(section) any(grepl(section, names(probe_tables), fixed = TRUE)),
    logical(1)
  )))

  skip_if_not_installed("openxlsx")
  path <- tempfile("psychomatic-list-report-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  psychomatic_report(objects, format = "excel", file_name = path, detail = "audit")
  expect_true(file.exists(path))
  expect_identical(length(openxlsx::getSheetNames(path)), 2L)
})

test_that("alignment item evidence survives audit and report collection", {
  skip_if_not_installed("sirt")
  lambda <- matrix(
    c(0.70, 0.80, 0.75, 0.72, 0.78, 0.74),
    nrow = 2,
    byrow = TRUE
  )
  nu <- matrix(
    c(0.00, 0.00, 0.00, 0.10, 0.05, 0.08),
    nrow = 2,
    byrow = TRUE
  )
  result <- inv_align_auto(lambda = lambda, nu = nu, language = "eng")
  audit <- audit_report(result)
  tables <- PsychoMatic:::.pm_collect_report_tables(list(alignment = result))

  expected <- c(
    "item_results_loadings_summary",
    "item_results_loadings_parm_est",
    "item_results_loadings_noninvariant_matrix",
    "item_results_intercepts_summary",
    "item_results_intercepts_parm_est"
  )
  expect_true(all(vapply(
    expected,
    function(label) any(grepl(label, names(audit$results), fixed = TRUE)),
    logical(1)
  )))
  expect_true(all(vapply(
    expected,
    function(label) any(grepl(label, names(tables), fixed = TRUE)),
    logical(1)
  )))
  parm_table <- audit$results[[grep(
    "item_results_loadings_parm_est$", names(audit$results), value = TRUE
  )[1L]]]
  expect_true(any(vapply(parm_table, is.numeric, logical(1))))
})

test_that("Excel sheet names are safe, unique, mapped, and numerically lossless", {
  skip_if_not_installed("openxlsx")
  labels <- c(
    "'A/very*long?invalid:[worksheet] name that exceeds thirty one characters'",
    "'a/VERY*LONG?INVALID:[worksheet] name that exceeds thirty one characters'",
    "''''",
    "[]:*?/\\"
  )
  tables <- structure(
    lapply(seq_along(labels), function(i) data.frame(value = pi * i)),
    names = labels
  )
  path <- tempfile("psychomatic-sheet-map-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  written <- PsychoMatic:::.pm_write_excel_report(tables, path)
  sheet_map <- attr(written, "sheet_map", exact = TRUE)
  sheets <- openxlsx::getSheetNames(path)

  expect_identical(unname(sheet_map$effective_sheet), sheets)
  expect_true(all(nchar(sheets) <= 31L))
  expect_false(any(grepl("[\\[\\]*?/\\\\:]", sheets)))
  expect_false(any(startsWith(sheets, "'") | endsWith(sheets, "'")))
  expect_identical(length(unique(tolower(sheets))), length(sheets))
  expect_identical(sheet_map$requested_sheet, labels)
  reopened <- openxlsx::read.xlsx(path, sheet = sheets[1L])
  # OOXML stores numeric cells as decimal text; reopening can differ by a few
  # machine epsilons while remaining a numeric, full-precision cell.
  expect_equal(reopened$value[[1L]], pi, tolerance = 1e-14)
})

test_that("Excel keeps finite numeric cells and labels missing or nonfinite values", {
  skip_if_not_installed("openxlsx")
  table <- data.frame(value = c(pi, NA_real_, NaN, Inf, -Inf))
  path <- tempfile("psychomatic-special-numbers-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)

  PsychoMatic:::.pm_write_excel_report(list(values = table), path)
  reopened <- openxlsx::read.xlsx(path, sheet = 1L, skipEmptyRows = FALSE)
  expect_true(is.numeric(reopened$value))
  expect_equal(reopened$value[[1L]], pi, tolerance = 1e-14)
  marker <- reopened$value__missing_or_nonfinite
  expect_identical(marker[2:5], c("<NA>", "<NaN>", "<Inf>", "<-Inf>"))
})

test_that("psychomatic_report validates inputs and preflights all outputs atomically", {
  expect_identical(
    names(formals(psychomatic_report)),
    c("...", "format", "file_name", "language", "digits", "overwrite", "detail")
  )
  expect_identical(eval(formals(psychomatic_report)$overwrite), FALSE)
  table <- data.frame(value = pi)
  expect_error(psychomatic_report(), "At least one")
  expect_error(psychomatic_report(x = table, digits = -1), "digits")
  expect_error(psychomatic_report(x = table, digits = 1.5), "integer")
  expect_error(psychomatic_report(x = table, overwrite = NA), "overwrite")
  expect_error(psychomatic_report(x = table, format = "pdf"), "arg")
  expect_error(psychomatic_report(x = table, language = "fr"), "arg|esp|eng")
  expect_error(
    psychomatic_report(x = table, format = c("excel", "excel")),
    "duplicate"
  )

  base <- tempfile("psychomatic-atomic-")
  xlsx <- paste0(base, ".xlsx")
  docx <- paste0(base, ".docx")
  writeLines("occupied", docx)
  on.exit(unlink(c(xlsx, docx)), add = TRUE)
  expect_error(
    psychomatic_report(
      x = table,
      format = c("excel", "word"),
      file_name = base,
      overwrite = FALSE
    ),
    "exists"
  )
  expect_false(file.exists(xlsx))

  excel <- tempfile("psychomatic-report-", fileext = ".xlsx")
  on.exit(unlink(excel), add = TRUE)
  visible <- withVisible(psychomatic_report(
    x = table,
    format = "excel",
    file_name = excel
  ))
  expect_false(visible$visible)
  expect_report_test_path(as.character(visible$value), excel)
  expect_true(file.exists(excel))
  expect_error(psychomatic_report(x = table, file_name = excel), "exists")
})

test_that("integrated Excel reports retain audit evidence and source objects", {
  skip_if_not_installed("openxlsx")
  probe <- reporting_probe()
  before <- serialize(probe, NULL)
  path <- tempfile("psychomatic-evidence-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)

  written <- psychomatic_report(
    probe = probe,
    detail = "audit",
    format = "excel",
    file_name = path,
    language = "eng"
  )
  sheets <- openxlsx::getSheetNames(path)
  expect_true(any(grepl("settings", sheets, ignore.case = TRUE)))
  expect_true(any(grepl("warnings", sheets, ignore.case = TRUE)))
  expect_true(any(grepl("errors", sheets, ignore.case = TRUE)))
  expect_true(any(grepl("conclusions", sheets, ignore.case = TRUE)))
  map <- attr(written, "sheet_map", exact = TRUE)
  expect_identical(map$effective_sheet, sheets)
  expect_identical(map$requested_sheet, attr(
    PsychoMatic:::.pm_collect_report_tables(list(probe = probe)),
    "table_map",
    exact = TRUE
  )$requested_sheet)
  expect_identical(serialize(probe, NULL), before)
})

test_that("integrated Word reports contain evidence headings", {
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  probe <- reporting_probe()
  path <- tempfile("psychomatic-evidence-", fileext = ".docx")
  on.exit(unlink(path), add = TRUE)

  psychomatic_report(
    probe = probe,
    detail = "audit",
    format = "word",
    file_name = path,
    language = "eng",
    digits = 4L
  )
  text <- officer::docx_summary(officer::read_docx(path))$text
  expect_true(any(grepl("PsychoMatic Report", text, fixed = TRUE)))
  expect_true(any(grepl("warnings", text, ignore.case = TRUE)))
  expect_true(any(grepl("conclusions", text, ignore.case = TRUE)))
})

test_that("CFA export adds a trailing fail-closed overwrite argument", {
  expect_identical(
    names(formals(export_cfa)),
    c("result", "format", "file_name", "digits", "overwrite")
  )
  expect_identical(eval(formals(export_cfa)$overwrite), FALSE)

  skip_if_not_installed("openxlsx")
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  result <- structure(list(), class = "cfa_auto")
  base <- tempfile("psychomatic-cfa-atomic-")
  xlsx <- paste0(base, ".xlsx")
  docx <- paste0(base, ".docx")
  writeLines("occupied", docx)
  on.exit(unlink(c(xlsx, docx)), add = TRUE)
  expect_error(
    export_cfa(
      result,
      format = c("excel", "word"),
      file_name = base,
      overwrite = FALSE
    ),
    "exists"
  )
  expect_false(file.exists(xlsx))
})
