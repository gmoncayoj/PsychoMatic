#' Create an integrated PsychoMatic report
#'
#' Collects tables from PsychoMatic result objects and writes an Excel workbook
#' and/or Word document. The function is intentionally conservative: it only
#' writes to the user-supplied path, or to `tempdir()` when `file_name = NULL`.
#'
#' @param ... Named PsychoMatic result objects or data frames. Supported objects
#'   include outputs from `desc_auto()`, `screen_items()`, `score_scale()`,
#'   `efa_auto()`, `cfa_auto()`, `factorial_invariance_auto()`,
#'   `inv_align_auto()`, `irt_auto()`, `dif_auto()`, and `compare_models_auto()`.
#' @param format Output format: `"excel"`, `"word"`, or both.
#' @param file_name Output path or base path. If `NULL`, a temporary file is
#'   created.
#' @param language Report language: `"eng"` or `"esp"`.
#' @param digits Number of decimal places for numeric tables.
#'
#' @return Output path, invisibly. When both formats are requested, returns a
#'   named character vector.
#'
#' @examples
#' if (requireNamespace("openxlsx", quietly = TRUE)) {
#'   items <- data.frame(item1 = 1:5, item2 = c(2, 2, 3, 4, 5))
#'   path <- tempfile("psychomatic-report-", fileext = ".xlsx")
#'   psychomatic_report(screening = screen_items(items), format = "excel", file_name = path)
#' }
#'
#' @export
psychomatic_report <- function(...,
                               format = "excel",
                               file_name = NULL,
                               language = c("eng", "esp"),
                               digits = 3) {
  objects <- list(...)
  if (length(objects) == 0L) {
    stop("At least one result object or table must be supplied.", call. = FALSE)
  }
  object_names <- names(objects)
  if (is.null(object_names) || any(!nzchar(object_names))) {
    names(objects) <- paste0("result_", seq_along(objects))
  }

  format <- match.arg(format, choices = c("excel", "word"), several.ok = TRUE)
  language <- match.arg(language)
  tables <- .pm_collect_report_tables(objects, digits = digits)
  if (length(tables) == 0L) {
    stop("No tabular content was found in the supplied objects.", call. = FALSE)
  }

  paths <- vapply(format, function(fmt) {
    path <- .pm_report_path(file_name, fmt)
    if (fmt == "excel") {
      .pm_write_excel_report(tables, path)
    } else {
      .pm_write_word_report(tables, path, language = language)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }, character(1))

  if (length(paths) == 1L) {
    invisible(unname(paths))
  } else {
    invisible(paths)
  }
}

.pm_collect_report_tables <- function(objects, digits = 3) {
  tables <- list()
  for (nm in names(objects)) {
    obj <- objects[[nm]]
    obj_tables <- .pm_object_tables(obj)
    for (tbl_nm in names(obj_tables)) {
      tbl <- obj_tables[[tbl_nm]]
      if (is.data.frame(tbl) && nrow(tbl) > 0L) {
        numeric_cols <- vapply(tbl, is.numeric, logical(1))
        tbl[numeric_cols] <- lapply(tbl[numeric_cols], function(x) round(x, digits))
        tables[[.pm_unique_name(tables, paste(nm, tbl_nm, sep = "_"))]] <- tbl
      }
    }
  }
  tables
}

.pm_object_tables <- function(obj) {
  if (is.data.frame(obj)) {
    return(list(table = obj))
  }
  if (inherits(obj, "psychomatic_scores")) {
    return(list(scores = obj$scores, reliability = obj$reliability))
  }
  if (inherits(obj, "psychomatic_irt")) {
    return(list(
      model_summary = obj$model_summary,
      item_parameters = obj$item_parameters,
      interpretation = obj$interpretation,
      references = obj$references
    ))
  }
  if (inherits(obj, "psychomatic_dif")) {
    return(list(
      summary = obj$summary,
      results = obj$results,
      interpretation = obj$interpretation,
      references = obj$references
    ))
  }
  if (inherits(obj, "psychomatic_model_comparison")) {
    return(list(fit = obj$fit))
  }
  if (inherits(obj, c("efa_auto", "cfa_auto", "factorial_invariance_auto"))) {
    smry <- summary(obj)
    return(smry[vapply(smry, is.data.frame, logical(1))])
  }
  if (inherits(obj, "inv_align_auto_result")) {
    return(list(summary = summary(obj)))
  }
  if (is.list(obj)) {
    return(obj[vapply(obj, is.data.frame, logical(1))])
  }
  list()
}

.pm_report_path <- function(file_name, format) {
  extension <- if (format == "excel") ".xlsx" else ".docx"
  if (is.null(file_name)) {
    return(file.path(
      tempdir(),
      paste0("psychomatic_report_", base::format(Sys.time(), "%Y%m%d_%H%M%S"), extension)
    ))
  }
  if (grepl(paste0("\\", extension, "$"), file_name, ignore.case = TRUE)) {
    return(file_name)
  }
  base_path <- sub("\\.(xlsx|docx)$", "", file_name, ignore.case = TRUE)
  paste0(base_path, extension)
}

.pm_write_excel_report <- function(tables, path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package `openxlsx` is required for Excel reports.", call. = FALSE)
  }
  wb <- openxlsx::createWorkbook()
  sheet_names <- character(0)
  for (nm in names(tables)) {
    sheet <- .pm_excel_sheet_name(nm, sheet_names)
    sheet_names <- c(sheet_names, sheet)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, tables[[nm]])
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

.pm_write_word_report <- function(tables, path, language = "eng") {
  if (!requireNamespace("officer", quietly = TRUE) ||
      !requireNamespace("flextable", quietly = TRUE)) {
    stop("Packages `officer` and `flextable` are required for Word reports.",
         call. = FALSE)
  }
  title <- if (language == "esp") "Reporte PsychoMatic" else "PsychoMatic Report"
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  for (nm in names(tables)) {
    doc <- officer::body_add_par(doc, gsub("_", " ", nm), style = "heading 2")
    doc <- flextable::body_add_flextable(doc, flextable::autofit(flextable::flextable(tables[[nm]])))
  }
  print(doc, target = path)
  path
}

.pm_excel_sheet_name <- function(name, existing) {
  base <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", name)
  base <- substr(base, 1L, 31L)
  candidate <- base
  idx <- 1L
  while (candidate %in% existing) {
    suffix <- paste0("_", idx)
    candidate <- paste0(substr(base, 1L, 31L - nchar(suffix)), suffix)
    idx <- idx + 1L
  }
  candidate
}

.pm_unique_name <- function(x, name) {
  candidate <- name
  idx <- 1L
  while (candidate %in% names(x)) {
    idx <- idx + 1L
    candidate <- paste0(name, "_", idx)
  }
  candidate
}
