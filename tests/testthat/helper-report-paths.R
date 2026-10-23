# Normalize filesystem identities in assertions, not user-facing report paths.
canonical_report_test_path <- function(path, must_exist = FALSE) {
  unname(vapply(as.character(path), function(value) {
    if (isTRUE(must_exist) && !file.exists(value)) {
      stop("Expected report path does not exist: ", value, call. = FALSE)
    }
    if (file.exists(value)) return(normalizePath(value, winslash = "/", mustWork = TRUE))
    parent <- dirname(value)
    if (!dir.exists(parent)) stop("Parent directory does not exist: ", parent, call. = FALSE)
    file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(value))
  }, character(1L)))
}

expect_report_test_path <- function(actual, expected, must_exist = TRUE) {
  testthat::expect_identical(canonical_report_test_path(actual, must_exist),
    canonical_report_test_path(expected, must_exist))
}
