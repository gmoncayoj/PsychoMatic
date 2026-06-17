#' Export EFA results to Excel or Word
#'
#' Exports the structured tables from an `efa_auto()` result. The function uses
#' the package-level report writer so EFA exports stay consistent with the
#' integrated `psychomatic_report()` workflow.
#'
#' @param result Object of class `"efa_auto"` returned by `efa_auto()`.
#' @param format Output format: `"excel"`, `"word"`, or both.
#' @param file_name Optional output path or base path. If `NULL`, a temporary
#'   file path is generated.
#'
#' @return The input `result`, invisibly.
#'
#' @examples
#' if (identical(Sys.getenv("PSYCHOMATIC_RUN_SLOW_EXAMPLES"), "true") &&
#'     requireNamespace("openxlsx", quietly = TRUE)) {
#'   set.seed(1)
#'   factor_score <- rnorm(80)
#'   items <- data.frame(
#'     item1 = factor_score + rnorm(80, sd = 0.4),
#'     item2 = factor_score + rnorm(80, sd = 0.4),
#'     item3 = factor_score + rnorm(80, sd = 0.4),
#'     item4 = factor_score + rnorm(80, sd = 0.4)
#'   )
#'   result <- efa_auto(items, max_iter = 1, verbose = FALSE, language = "eng")
#'   export_efa(result, format = "excel", file_name = tempfile("efa-", fileext = ".xlsx"))
#' }
#'
#' @export
export_efa <- function(result, format = c("excel", "word"), file_name = NULL) {
  if (!inherits(result, "efa_auto")) {
    stop("`result` must be an object of class `efa_auto`.", call. = FALSE)
  }

  format <- match.arg(format, choices = c("excel", "word"), several.ok = TRUE)

  output_path <- file_name
  if (is.null(output_path)) {
    output_path <- file.path(
      tempdir(),
      paste0("efa_results_", base::format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }

  language <- if (!is.null(result$language)) result$language else "eng"

  paths <- psychomatic_report(
    efa = result,
    format = format,
    file_name = output_path,
    language = language
  )

  for (path in paths) {
    message("EFA report exported: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  invisible(result)
}

# Internal legacy alias for `export_efa()`.
#'
#' @param resultado Deprecated. Use `result`.
#' @param formato Deprecated. Use `format`.
#' @param archivo Deprecated. Use `file_name`.
#'
#' @return The input `resultado`, invisibly.
#' @noRd
exportar_efa <- function(resultado, formato = c("excel", "word"), archivo = NULL) {
  warning("`exportar_efa()` is deprecated; use `export_efa()` instead.", call. = FALSE)
  export_efa(result = resultado, format = formato, file_name = archivo)
}
