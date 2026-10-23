#' Export EFA results to Excel or Word
#'
#' Exports the structured tables from an `efa_auto()` result. The function uses
#' the package-level report writer so EFA exports stay consistent with the
#' integrated `psychomatic_report()` workflow.
#'
#' @param result Object of class `"efa_auto"` returned by `efa_auto()`.
#' @param format Output format: `"excel"`, `"word"`, or both.
#' @param file_name Optional output path or base path. If `NULL`, a unique
#'   timestamped name is generated in the working directory (`getwd()`).
#' @param overwrite If `TRUE`, existing output files may be replaced. The
#'   default is fail-closed.
#'
#' @return The input `result`, invisibly.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' result <- efa_auto(MY_DATA[, 3:14])
#' export_efa(result, format = "excel", file_name = "EFA.xlsx")
#' }
#'
#' @export
export_efa <- function(result,
                       format = c("excel", "word"),
                       file_name = NULL,
                       overwrite = FALSE) {
  if (!inherits(result, "efa_auto")) {
    stop("`result` must be an object of class `efa_auto`.", call. = FALSE)
  }

  format <- match.arg(format, choices = c("excel", "word"), several.ok = TRUE)

  language <- if (!is.null(result$language)) result$language else "eng"

  paths <- psychomatic_report(
    efa = result,
    format = format,
    file_name = file_name,
    language = language,
    overwrite = overwrite
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
#' @param sobrescribir Deprecated. Use `overwrite`.
#'
#' @return The input `resultado`, invisibly.
#' @noRd
exportar_efa <- function(resultado,
                         formato = c("excel", "word"),
                         archivo = NULL,
                         sobrescribir = FALSE) {
  warning("`exportar_efa()` is deprecated; use `export_efa()` instead.", call. = FALSE)
  export_efa(
    result = resultado,
    format = formato,
    file_name = archivo,
    overwrite = sobrescribir
  )
}
