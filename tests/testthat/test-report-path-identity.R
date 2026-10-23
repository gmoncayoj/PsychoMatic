test_that("report path comparisons distinguish identity from textual aliases", {
  directory <- tempfile("report path identity ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  target <- file.path(directory, "report.xlsx")
  writeLines("original report", target)
  expect_identical(canonical_report_test_path(target, must_exist = TRUE),
    normalizePath(target, winslash = "/", mustWork = TRUE))
  expect_error(canonical_report_test_path(file.path(directory, "missing.xlsx"),
    must_exist = TRUE), "exist")
  expect_error(canonical_report_test_path(file.path(directory, "absent", "report.xlsx")),
    "exist")
  expect_false(identical(canonical_report_test_path(target),
    canonical_report_test_path(file.path(directory, "another.xlsx"))))
  dir.create(file.path(directory, "other"))
  expect_false(identical(canonical_report_test_path(target),
    canonical_report_test_path(file.path(directory, "other", "report.xlsx"))))
  if (.Platform$OS.type == "windows") {
    alias <- utils::shortPathName(directory)
  } else {
    alias <- paste0(directory, "-link")
    expect_true(file.symlink(directory, alias))
    on.exit(unlink(alias), add = TRUE)
  }
  expect_identical(canonical_report_test_path(file.path(alias, "future.xlsx")),
    canonical_report_test_path(file.path(directory, "future.xlsx")))
  expect_identical(canonical_report_test_path(file.path(alias, "report.xlsx"), must_exist = TRUE),
    canonical_report_test_path(target, must_exist = TRUE))
  expect_identical(readLines(target), "original report")
})
