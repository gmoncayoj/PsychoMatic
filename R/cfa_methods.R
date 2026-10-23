# METODO PRINT
# =============================================================================

#' Print a `cfa_auto` object
#'
#' @param x Object returned by `cfa_auto()`.
#' @param digits Number of decimal places shown in the printed report.
#' @param ... Additional unused arguments.
#'
#' @rdname print.cfa_auto
#' @method print cfa_auto
#' @keywords internal
#' @export
print.cfa_auto <- function(x, digits = NULL, ...) {
  .pm_console_report(x, digits)
}

.cfa_respecifications_to_df <- function(x) {
  respec <- if (!is.null(x$respecification_suggestions)) {
    x$respecification_suggestions
  } else {
    x$respecifications
  }

  if (is.null(respec) || length(respec) == 0) {
    return(data.frame())
  }

  do.call(
    rbind,
    lapply(respec, function(item) {
      data.frame(
        Type = if (!is.null(item$type)) item$type else item$tipo,
        Syntax = if (!is.null(item$syntax)) item$syntax else item$sintaxis,
        MI = item$mi,
        EPC = item$epc,
        SameFactor = if (!is.null(item$same_factor)) item$same_factor else item$mismo_factor,
        Note = if (!is.null(item$note)) item$note else item$nota,
        stringsAsFactors = FALSE
      )
    })
  )
}

.cfa_export_excel <- function(result,
                              file_name,
                              digits = 3L,
                              overwrite = FALSE) {
  tables <- .pm_object_tables(result)
  if (!length(tables)) {
    stop("No structured CFA tables were available for Excel export.", call. = FALSE)
  }
  sheet_labels <- c(
    fit = "Fit Indices",
    loadings = "Loadings",
    reliability = "Reliability",
    decisions = "Decisions",
    higher_order_loadings = "Higher Order",
    modification_indices = "Mod Indices",
    respecification_suggestions = "Respec",
    factor_correlation_table = "Factor Corr"
  )
  matched <- match(names(tables), names(sheet_labels))
  use_label <- !is.na(matched)
  names(tables)[use_label] <- unname(sheet_labels[matched[use_label]])
  .pm_write_excel_report(
    tables,
    file_name,
    overwrite = overwrite,
    digits = digits
  )
}

.cfa_export_word <- function(result, file_name, digits = 3, overwrite = FALSE) {
  tables <- .pm_object_tables(result)
  if (!length(tables)) {
    stop("No structured CFA tables were available for Word export.", call. = FALSE)
  }
  language <- if (!is.null(result$language)) result$language else "eng"
  .pm_write_word_report(
    tables,
    file_name,
    language = .pm_normalize_language(language),
    digits = digits,
    overwrite = overwrite
  )
}

#' Export CFA results to Excel or Word
#'
#' @param result Object of class `"cfa_auto"` returned by `cfa_auto()`.
#' @param format `"excel"` or `"word"` (you can provide both).
#' @param file_name Optional base file name without extension. If `NULL`, a
#'   unique timestamped name is generated in the working directory (`getwd()`).
#' @param digits Number of decimal places used to format report values. Raw
#'   numeric values remain unchanged in Excel files.
#' @param overwrite If `TRUE`, existing output files may be replaced.
#'
#' @return Normalized output path, invisibly. When both formats are requested,
#'   returns a named character vector; Excel output includes a `sheet_map`
#'   attribute.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' model <- "F1 =~ R1 + R2 + R3 + R4"
#' result <- cfa_auto(MY_DATA, model = model)
#' export_cfa(result, format = "word", file_name = "CFA.docx")
#' }
#'
#' @export
export_cfa <- function(result,
                       format = c("excel", "word"),
                       file_name = NULL,
                       digits = 3,
                       overwrite = FALSE) {
  if (!inherits(result, "cfa_auto")) {
    stop("'result' must be an object of class 'cfa_auto'.", call. = FALSE)
  }

  format <- match.arg(format, c("excel", "word"), several.ok = TRUE)
  if (anyDuplicated(format)) {
    stop("`format` must not contain duplicate output formats.", call. = FALSE)
  }
  digits <- .pm_validate_scalar_number(digits, "digits", lower = 0, integer = TRUE)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(file_name)) {
    file_name <- .pm_report_path(NULL, format[[1L]])
  }
  .pm_check_report_dependencies(format)

  paths <- vapply(format, function(fmt) .pm_report_path(file_name, fmt), character(1))
  paths <- vapply(seq_along(paths), function(i) {
    extension <- if (format[[i]] == "excel") ".xlsx" else ".docx"
    .pm_safe_output_path(paths[[i]], extension = extension, overwrite = overwrite)
  }, character(1))

  outputs <- lapply(seq_along(format), function(i) {
    fmt <- format[[i]]
    if (fmt == "excel") {
      .cfa_export_excel(
        result,
        paths[[i]],
        digits = digits,
        overwrite = overwrite
      )
    } else {
      .cfa_export_word(result, paths[[i]], digits = digits, overwrite = overwrite)
    }
  })

  if (length(outputs) == 1L) {
    return(invisible(outputs[[1L]]))
  }
  output_paths <- vapply(outputs, as.character, character(1L))
  names(output_paths) <- format
  sheet_maps <- lapply(outputs, attr, which = "sheet_map", exact = TRUE)
  sheet_maps <- Filter(Negate(is.null), sheet_maps)
  if (length(sheet_maps)) {
    attr(output_paths, "sheet_map") <- do.call(rbind, sheet_maps)
  }
  invisible(output_paths)
}

# Internal legacy alias for `export_cfa()`
#'
#' @param resultado Deprecated. Use `result`.
#' @param formato Deprecated. Use `format`.
#' @param archivo Deprecated. Use `file_name`.
#' @param digits Number of decimal places used when exporting the Word report.
#'
#' @return Normalized output path(s), invisibly.
#' @noRd
exportar_cfa <- function(resultado,
                         formato = c("excel", "word"),
                         archivo = NULL,
                         digits = 3,
                         overwrite = FALSE) {
  warning("`exportar_cfa()` is deprecated; use `export_cfa()` instead.", call. = FALSE)
  export_cfa(
    result = resultado,
    format = formato,
    file_name = archivo,
    digits = digits,
    overwrite = overwrite
  )
}
