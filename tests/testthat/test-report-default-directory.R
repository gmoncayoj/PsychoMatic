with_report_project <- function(code) {
  directory <- tempfile("proyecto reportes-\u00e1-")
  dir.create(directory)
  previous <- getwd()
  on.exit({ setwd(previous); unlink(directory, recursive = TRUE) }, add = TRUE)
  setwd(directory)
  eval.parent(substitute(code))
}

test_that("desc_auto saves a readable default workbook in the current project", {
  skip_if_not_installed("openxlsx")
  with_report_project({
    items <- lavaan::HolzingerSwineford1939[paste0("x", 1:3)]
    result <- desc_auto(items, report = TRUE, digits = 2, language = "esp", verbose = FALSE)
    path <- unname(attr(result, "report_paths")[["excel"]])
    expect_identical(dirname(path), normalizePath(getwd(), winslash = "/"))
    expect_true(file.exists(path))
    expect_true(length(openxlsx::getSheetNames(path)) > 0L)
    expect_match(basename(path), "^psychomatic_report_[0-9]{8}_[0-9]{6}_.*\\.xlsx$")
    before <- readBin(path, "raw", file.info(path)$size)
    second <- desc_auto(items, report = TRUE, verbose = FALSE)
    expect_false(identical(path, unname(attr(second, "report_paths")[["excel"]])))
    expect_identical(readBin(path, "raw", file.info(path)$size), before)
  })
})

test_that("integrated reader and audit reports use paired default project paths", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  with_report_project({
    for (detail in c("report", "audit")) {
      paths <- psychomatic_report(data.frame(value = pi),
        format = c("excel", "word"), detail = detail)
      expect_identical(names(paths), c("excel", "word"))
      expect_true(all(file.exists(paths)))
      expect_true(all(dirname(paths) == normalizePath(getwd(), winslash = "/")))
      expect_identical(tools::file_path_sans_ext(basename(paths[["excel"]])),
        tools::file_path_sans_ext(basename(paths[["word"]])))
      expect_true(nrow(openxlsx::read.xlsx(paths[["excel"]])) > 0L)
      expect_true(nrow(officer::docx_summary(officer::read_docx(paths[["word"]]))) > 0L)
    }
    descriptive <- desc_auto(lavaan::HolzingerSwineford1939[paste0("x", 1:3)],
      report = TRUE, format = c("excel", "word"), verbose = FALSE)
    paths <- attr(descriptive, "report_paths")
    expect_true(all(file.exists(paths)))
    expect_true(all(dirname(paths) == normalizePath(getwd(), winslash = "/")))
    expect_identical(tools::file_path_sans_ext(basename(paths[["excel"]])),
      tools::file_path_sans_ext(basename(paths[["word"]])))
  })
})

test_that("EFA and CFA exports do not override the project destination with temporary paths", {
  skip_if_not_installed("openxlsx")
  with_report_project({
    items <- lavaan::HolzingerSwineford1939[paste0("x", 1:9)]
    efa <- efa_auto(items, refine = "suggest", parallel_iter = 5, seed = 206, report = TRUE, verbose = FALSE)
    cfa <- cfa_auto(items,
      model = "visual =~ x1+x2+x3\ntextual =~ x4+x5+x6\nspeed =~ x7+x8+x9",
      estimator = "ML", report = TRUE, verbose = FALSE)
    for (result in list(efa, cfa)) {
      paths <- attr(result, "report_paths")
      expect_length(paths, 1L)
      expect_true(file.exists(paths))
      expect_identical(dirname(unname(paths[[1L]])), normalizePath(getwd(), winslash = "/"))
    }
    before <- list.files(getwd(), full.names = TRUE)
    expect_message(export_efa(efa, format = "excel"), "EFA report exported")
    efa_paths <- setdiff(list.files(getwd(), full.names = TRUE), before)
    expect_length(efa_paths, 1L)
    if (length(efa_paths)) expect_true(length(openxlsx::getSheetNames(efa_paths[1])) > 0L)
    cfa_path <- export_cfa(cfa, format = "excel")
    expect_identical(dirname(as.character(cfa_path)), normalizePath(getwd(), winslash = "/"))
    expect_true(file.exists(cfa_path))
  })
})

test_that("shared and legacy report options resolve the default destination at call time", {
  skip_if_not_installed("openxlsx")
  with_report_project({
    first <- PsychoMatic:::.pm_output_options(report = TRUE)
    legacy <- PsychoMatic:::.ifa_preflight_report(TRUE, overwrite = FALSE)
    expected <- normalizePath(getwd(), winslash = "/")
    expect_identical(dirname(unname(first$paths[[1L]])), expected)
    expect_identical(dirname(legacy), expected)
    dir.create("otro proyecto")
    setwd("otro proyecto")
    next_output <- PsychoMatic:::.pm_output_options(report = TRUE)
    expect_identical(dirname(unname(next_output$paths[[1L]])), normalizePath(getwd(), winslash = "/"))
    expect_false(identical(dirname(first$paths[[1L]]), dirname(next_output$paths[[1L]])))
  })
})

test_that("explicit report paths remain authoritative and existing files remain protected", {
  skip_if_not_installed("openxlsx")
  with_report_project({
    dir.create("explicit destination")
    explicit <- normalizePath("explicit destination", winslash = "/")
    target <- file.path(explicit, "descriptives.xlsx")
    items <- lavaan::HolzingerSwineford1939[paste0("x", 1:3)]
    first <- desc_auto(items, report = TRUE, file_name = target, verbose = FALSE)
    expect_identical(unname(attr(first, "report_paths")[["excel"]]), target)
    expect_error(desc_auto(items, report = TRUE, file_name = target, verbose = FALSE), "exists")
    expect_no_error(desc_auto(items, report = TRUE, file_name = target, overwrite = TRUE, verbose = FALSE))
    relative <- psychomatic_report(data.frame(value = pi), file_name = "explicit destination/relative")
    expect_true(file.exists(relative))
    expect_identical(normalizePath(as.character(relative), winslash = "/"),
      file.path(explicit, "relative.xlsx"))
    disabled <- desc_auto(items, report = FALSE, verbose = FALSE)
    expect_null(attr(disabled, "report_paths"))
  })
})
