# Expensive native failure probes and exhaustive export/permutation regressions
# remain executable, but require an explicit opt-in (set in the CI workflow).
.pm_extended_tests <- function() {
  identical(Sys.getenv("PSYCHOMATIC_EXTENDED_TESTS", unset = ""), "true")
}

.pm_skip_extended <- function() {
  if (!.pm_extended_tests()) {
    testthat::skip("Extended regression: set PSYCHOMATIC_EXTENDED_TESTS=true")
  }
  invisible(NULL)
}
