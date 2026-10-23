load_data_fixture <- function(name, envir = parent.frame()) {
  load(testthat::test_path("..", "fixtures", paste0(name, ".rda")), envir = envir)
  invisible(NULL)
}
