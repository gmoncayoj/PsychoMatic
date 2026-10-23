#' Create an integrated PsychoMatic report
#'
#' Collects tables from PsychoMatic result objects and writes an Excel workbook
#' and/or Word document. Reports use the user-supplied path or the current
#' working directory (`getwd()`) when `file_name = NULL`.
#'
#' @param ... Named PsychoMatic result objects or data frames. Supported objects
#'   include outputs from `desc_auto()`, `screen_items()`, `score_scale()`,
#'   `reliability_auto()`, `efa_auto()`, `cfa_auto()`, `esem_auto()`,
#'   `factorial_invariance_auto()`, `ordinal_invariance_auto()`,
#'   `inv_align_auto()`, `irt_auto()`, `dif_auto()`, and `compare_models_auto()`.
#' @param format Output format: `"excel"`, `"word"`, or both.
#' @param file_name Output path or base path. If `NULL`, a unique timestamped
#'   name is generated in `getwd()`. Explicit paths remain authoritative.
#' @param language Report language: `"eng"` or `"esp"`.
#' @param digits Number of decimal places for numeric tables.
#' @param overwrite If `TRUE`, existing report files may be replaced.
#' @param detail `"report"` (default) organizes reader-facing results by analysis.
#'   `"audit"` retains the detailed technical tables and execution metadata.
#'
#' @return Output path, invisibly. When both formats are requested, returns a
#'   named character vector. Excel output carries a `sheet_map` attribute that
#'   records requested and effective worksheet names.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' result <- reliability_auto(MY_DATA[, 3:6])
#' psychomatic_report(reliability = result, format = "word", file_name = "Report.docx")
#' }
#'
#' @export
psychomatic_report <- function(...,
                               format = "excel",
                               file_name = NULL,
                               language = c("eng", "esp"),
                               digits = 3,
                               overwrite = FALSE,
                               detail = c("report", "audit")) {
  objects <- list(...)
  if (length(objects) == 0L) {
    stop("At least one result object or table must be supplied.", call. = FALSE)
  }
  if (length(objects) == 1L && .pm_is_report_collection(objects[[1L]])) {
    objects <- objects[[1L]]
  }
  format <- match.arg(format, choices = c("excel", "word"), several.ok = TRUE)
  if (anyDuplicated(format)) {
    stop("`format` must not contain duplicate output formats.", call. = FALSE)
  }
  language <- .pm_normalize_language(match.arg(language))
  detail <- match.arg(detail)
  digits <- .pm_validate_scalar_number(digits, "digits", lower = 0, integer = TRUE)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(file_name)) file_name <- .pm_report_path(NULL, format[[1L]])
  paths <- vapply(
    format,
    function(fmt) .pm_report_path(file_name, fmt),
    character(1)
  )
  paths <- vapply(seq_along(paths), function(i) {
    extension <- if (format[[i]] == "excel") ".xlsx" else ".docx"
    .pm_safe_output_path(paths[[i]], extension = extension, overwrite = overwrite)
  }, character(1))
  names(paths) <- format
  .pm_check_report_dependencies(format)

  if (identical(detail, "report") && !any(vapply(objects, inherits, logical(1), "psychomatic_audit"))) {
    reports <- lapply(objects, .pm_analysis_sections, language = language, destination = "report")
    labels <- names(objects)
    if (!is.null(labels)) for (i in seq_along(reports)) {
      if (!is.na(labels[i]) && nzchar(labels[i])) reports[[i]]$title <- labels[i]
    }
    outputs <- lapply(seq_along(format), function(i) {
      if (format[i] == "excel") .pm_write_analysis_excel(reports, paths[i], digits, overwrite)
      else .pm_write_analysis_word(reports, paths[i], digits, overwrite)
    })
    if (length(outputs) == 1L) return(invisible(outputs[[1L]]))
    names(outputs) <- format
    return(invisible(vapply(outputs, as.character, character(1))))
  }

  tables <- .pm_collect_report_tables(objects, digits = digits)
  if (length(tables) == 0L) {
    stop("No tabular content was found in the supplied objects.", call. = FALSE)
  }

  outputs <- lapply(seq_along(format), function(i) {
    fmt <- format[[i]]
    path <- paths[[i]]
    if (fmt == "excel") {
      .pm_write_excel_report(
        tables,
        path,
        overwrite = overwrite,
        digits = digits
      )
    } else {
      .pm_write_word_report(
        tables,
        path,
        language = language,
        digits = digits,
        overwrite = overwrite
      )
    }
  })

  if (length(outputs) == 1L) {
    return(invisible(outputs[[1L]]))
  }

  paths <- vapply(outputs, as.character, character(1))
  names(paths) <- format
  sheet_maps <- lapply(outputs, attr, which = "sheet_map", exact = TRUE)
  sheet_maps <- Filter(Negate(is.null), sheet_maps)
  if (length(sheet_maps)) attr(paths, "sheet_map") <- do.call(rbind, sheet_maps)
  invisible(paths)
}

.pm_collect_report_tables <- function(objects, digits = 3) {
  .pm_validate_scalar_number(digits, "digits", lower = 0, integer = TRUE)
  if (!is.list(objects)) stop("`objects` must be a list.", call. = FALSE)
  tables <- list()
  mapping <- list()
  labels <- names(objects)
  if (is.null(labels)) labels <- rep("", length(objects))
  labels[is.na(labels)] <- ""

  for (i in seq_along(objects)) {
    object_id <- sprintf("result_%03d", i)
    display_label <- labels[[i]]
    if (!nzchar(trimws(display_label))) display_label <- object_id
    obj <- objects[[i]]
    if (!.pm_is_supported_report_object(obj)) {
      stop(
        sprintf(
          "Unsupported report object at index %d (class: %s).",
          i,
          paste(class(obj), collapse = "/")
        ),
        call. = FALSE
      )
    }
    obj_tables <- .pm_object_tables(obj)
    if (!length(obj_tables)) {
      stop(
        sprintf("No tabular content was found for report object at index %d.", i),
        call. = FALSE
      )
    }
    sections <- names(obj_tables)
    if (is.null(sections)) sections <- rep("table", length(obj_tables))
    sections[is.na(sections) | !nzchar(sections)] <- paste0(
      "table_", which(is.na(sections) | !nzchar(sections))
    )
    for (j in seq_along(obj_tables)) {
      tbl <- obj_tables[[j]]
      if (!is.data.frame(tbl) || nrow(tbl) == 0L) next
      requested_name <- paste(display_label, sections[[j]], sep = "_")
      internal_name <- .pm_unique_name(tables, requested_name, ignore_case = TRUE)
      tables[[internal_name]] <- tbl
      mapping[[length(mapping) + 1L]] <- data.frame(
        object_index = i,
        object_id = object_id,
        display_label = display_label,
        section = sections[[j]],
        requested_sheet = requested_name,
        effective_sheet = NA_character_,
        table_key = internal_name,
        stringsAsFactors = FALSE
      )
    }
  }
  table_map <- if (length(mapping)) do.call(rbind, mapping) else data.frame(
    object_index = integer(), object_id = character(), display_label = character(),
    section = character(), requested_sheet = character(), effective_sheet = character(),
    table_key = character(), stringsAsFactors = FALSE
  )
  rownames(table_map) <- NULL
  attr(tables, "table_map") <- table_map
  tables
}

.pm_object_tables <- function(obj) {
  if (is.data.frame(obj) && !inherits(obj, "psychomatic_item_screen")) {
    return(list(table = obj))
  }

  if (!.pm_is_supported_audit_result(obj)) {
    stop(
      sprintf("Unsupported report object (class: %s).", paste(class(obj), collapse = "/")),
      call. = FALSE
    )
  }

  tables <- list()
  append_tables <- function(values, prefix = NULL) {
    if (!is.list(values)) return(invisible(NULL))
    keep <- vapply(values, is.data.frame, logical(1))
    values <- values[keep]
    if (!length(values)) return(invisible(NULL))
    labels <- names(values)
    if (is.null(labels)) labels <- rep("table", length(values))
    labels[is.na(labels) | !nzchar(labels)] <- paste0(
      "table_", which(is.na(labels) | !nzchar(labels))
    )
    if (!is.null(prefix)) labels <- paste(prefix, labels, sep = "_")
    for (i in seq_along(values)) {
      key <- .pm_unique_name(tables, labels[[i]], ignore_case = TRUE)
      tables[[key]] <<- values[[i]]
    }
    invisible(NULL)
  }

  audited <- if (inherits(obj, "psychomatic_audit")) obj else audit_report(obj)
  audit_sections <- audited[c(
    "objects", "settings", "data", "warnings", "errors", "decisions",
    "diagnostics", "conclusions", "references", "versions", "artifacts"
  )]
  append_tables(audit_sections, prefix = "audit")
  result_tables <- audited$results
  if (!inherits(obj, "psychomatic_audit") && length(result_tables)) {
    names(result_tables) <- sub("^result_001_", "", names(result_tables))
  }
  append_tables(result_tables)

  tables
}

.pm_report_path <- function(file_name, format) {
  extension <- if (format == "excel") ".xlsx" else ".docx"
  if (is.null(file_name)) {
    return(tempfile(
      paste0("psychomatic_report_", base::format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
      tmpdir = getwd(), fileext = extension
    ))
  }
  if (grepl(paste0("\\", extension, "$"), file_name, ignore.case = TRUE)) {
    return(file_name)
  }
  base_path <- sub("\\.(xlsx|docx)$", "", file_name, ignore.case = TRUE)
  paste0(base_path, extension)
}

.pm_write_excel_report <- function(tables,
                                   path,
                                   overwrite = FALSE,
                                   digits = 3L) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package `openxlsx` is required for Excel reports.", call. = FALSE)
  }
  path <- .pm_safe_output_path(path, extension = ".xlsx", overwrite = overwrite)
  digits <- .pm_validate_scalar_number(digits, "digits", lower = 0, integer = TRUE)
  wb <- openxlsx::createWorkbook()
  sheet_names <- character(0)
  formatted_references <- character()
  input_map <- attr(tables, "table_map", exact = TRUE)
  requested_labels <- names(tables)
  if (is.null(requested_labels)) requested_labels <- rep("", length(tables))
  missing_labels <- is.na(requested_labels) | !nzchar(requested_labels)
  requested_labels[missing_labels] <- paste0("result_", which(missing_labels))
  if (is.data.frame(input_map) && nrow(input_map) == length(tables)) {
    sheet_map <- input_map
    requested_labels <- as.character(sheet_map$requested_sheet)
    missing_requested <- is.na(requested_labels) | !nzchar(requested_labels)
    requested_labels[missing_requested] <- names(tables)[missing_requested]
    sheet_map$effective_sheet <- NA_character_
  } else {
    sheet_map <- data.frame(
      object_index = seq_along(tables),
      display_label = requested_labels,
      requested_sheet = requested_labels,
      effective_sheet = rep(NA_character_, length(tables)),
      stringsAsFactors = FALSE
    )
  }
  for (i in seq_along(tables)) {
    sheet <- .pm_excel_sheet_name(requested_labels[i], sheet_names)
    sheet_names <- c(sheet_names, sheet)
    sheet_map$effective_sheet[i] <- sheet
    openxlsx::addWorksheet(wb, sheet)
    export_table <- .pm_prepare_table_for_export(
      tables[[i]],
      preserve_numeric = TRUE
    )
    openxlsx::writeData(
      wb,
      sheet,
      export_table
    )
    for (column in .pm_reference_columns(export_table)) formatted_references <- c(formatted_references,as.character(export_table[[column]]))
    numeric_columns <- which(vapply(export_table, is.numeric, logical(1L)))
    if (length(numeric_columns) && nrow(export_table)) {
      integer_columns <- numeric_columns[vapply(
        export_table[numeric_columns],
        is.integer,
        logical(1L)
      )]
      decimal_columns <- setdiff(numeric_columns, integer_columns)
      if (length(integer_columns)) {
        openxlsx::addStyle(
          wb,
          sheet,
          openxlsx::createStyle(numFmt = "0"),
          rows = seq_len(nrow(export_table)) + 1L,
          cols = integer_columns,
          gridExpand = TRUE,
          stack = TRUE
        )
      }
      if (length(decimal_columns)) {
        decimal_format <- if (digits == 0L) {
          "0"
        } else {
          paste0("0.", paste(rep("0", digits), collapse = ""))
        }
        openxlsx::addStyle(
          wb,
          sheet,
          openxlsx::createStyle(numFmt = decimal_format),
          rows = seq_len(nrow(export_table)) + 1L,
          cols = decimal_columns,
          gridExpand = TRUE,
          stack = TRUE
        )
      }
    }
  }
  .pm_reference_excel(wb,formatted_references[!is.na(formatted_references)])
  openxlsx::saveWorkbook(wb, path, overwrite = overwrite)
  attr(path, "sheet_map") <- sheet_map
  invisible(path)
}

.pm_write_word_report <- function(tables,
                                  path,
                                  language = "eng",
                                  digits = 3L,
                                  overwrite = FALSE) {
  if (!requireNamespace("officer", quietly = TRUE) ||
      !requireNamespace("flextable", quietly = TRUE)) {
    stop("Packages `officer` and `flextable` are required for Word reports.",
         call. = FALSE)
  }
  path <- .pm_safe_output_path(path, extension = ".docx", overwrite = overwrite)
  digits <- .pm_validate_scalar_number(digits, "digits", lower = 0, integer = TRUE)
  title <- if (language == "esp") "Reporte PsychoMatic" else "PsychoMatic Report"
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  for (nm in names(tables)) {
    doc <- officer::body_add_par(doc, gsub("_", " ", nm), style = "heading 2")
    display_table <- .pm_prepare_table_for_export(tables[[nm]], digits = digits)
    ft <- .pm_reference_flextable(flextable::flextable(display_table),display_table)
    doc <- flextable::body_add_flextable(doc, flextable::autofit(ft))
  }
  print(doc, target = path)
  path
}

.pm_excel_sheet_name <- function(name, existing) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("Excel sheet names must be character scalars.", call. = FALSE)
  }
  base <- name
  for (invalid in c("[", "]", "*", "?", "/", "\\", ":")) {
    base <- gsub(invalid, "_", base, fixed = TRUE)
  }
  base <- gsub("[[:cntrl:]]", "_", base)
  base <- sub("^'+", "", base)
  base <- sub("'+$", "", base)
  base <- trimws(base)
  if (!nzchar(base)) base <- "Sheet"
  base <- substr(base, 1L, 31L)
  candidate <- base
  idx <- 2L
  existing_folded <- tolower(existing)
  while (tolower(candidate) %in% existing_folded) {
    suffix <- paste0("_", idx)
    candidate <- paste0(substr(base, 1L, 31L - nchar(suffix)), suffix)
    idx <- idx + 1L
  }
  candidate
}

.pm_safe_output_path <- function(path, extension = NULL, overwrite = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`path` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(extension)) {
    if (!is.character(extension) || length(extension) != 1L ||
        is.na(extension) || !nzchar(extension)) {
      stop("`extension` must be NULL or a non-empty character scalar.", call. = FALSE)
    }
    if (!startsWith(extension, ".")) extension <- paste0(".", extension)
    if (grepl("[\\/\\:*?\"<>|]", extension)) {
      stop("`extension` contains invalid path characters.", call. = FALSE)
    }
    if (!grepl(paste0("\\Q", extension, "\\E$"), path, ignore.case = TRUE, perl = TRUE)) {
      path <- paste0(path, extension)
    }
  }

  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parent <- dirname(normalized)
  if (!dir.exists(parent)) {
    stop("Parent directory does not exist: ", parent, call. = FALSE)
  }
  if (file.exists(normalized) && !isTRUE(overwrite)) {
    stop("Output path already exists: ", normalized, call. = FALSE)
  }
  normalized
}

.pm_unique_name <- function(x, name, ignore_case = FALSE) {
  candidate <- name
  idx <- 1L
  existing <- names(x)
  if (is.null(existing)) existing <- character()
  is_used <- function(value) {
    if (isTRUE(ignore_case)) tolower(value) %in% tolower(existing) else value %in% existing
  }
  while (is_used(candidate)) {
    idx <- idx + 1L
    candidate <- paste0(name, "_", idx)
  }
  candidate
}

.pm_prepare_table_for_export <- function(table,
                                         digits = NULL,
                                         preserve_numeric = FALSE) {
  if (!is.data.frame(table)) stop("Report tables must be data frames.", call. = FALSE)
  out <- table
  original_names <- names(out)
  for (name in original_names) {
    column <- out[[name]]
    if (is.list(column)) {
      out[[name]] <- vapply(column, .pm_audit_display, character(1))
    } else if (is.numeric(column)) {
      special <- is.na(column) | is.nan(column) | is.infinite(column)
      if (any(special) && isTRUE(preserve_numeric)) {
        status <- rep("", length(column))
        status[is.na(column)] <- "<NA>"
        status[is.nan(column)] <- "<NaN>"
        status[is.infinite(column) & column > 0] <- "<Inf>"
        status[is.infinite(column) & column < 0] <- "<-Inf>"
        column[is.nan(column) | is.infinite(column)] <- NA_real_
        out[[name]] <- column
        status_name <- .pm_unique_column_name(
          names(out), paste0(name, "__missing_or_nonfinite")
        )
        out[[status_name]] <- status
      } else if (any(special)) {
        out[[name]] <- vapply(column, function(value) {
          if (is.nan(value)) return("NaN")
          if (is.na(value)) return("NA")
          if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
          as.character(value)
        }, character(1))
      } else if (!is.null(digits)) {
        out[[name]] <- round(column, digits)
      }
    }
  }
  out
}

.pm_unique_column_name <- function(existing, requested) {
  candidate <- requested
  index <- 1L
  while (tolower(candidate) %in% tolower(existing)) {
    index <- index + 1L
    candidate <- paste0(requested, "_", index)
  }
  candidate
}

.pm_check_report_dependencies <- function(format) {
  if ("excel" %in% format && !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package `openxlsx` is required for Excel reports.", call. = FALSE)
  }
  if ("word" %in% format &&
      (!requireNamespace("officer", quietly = TRUE) ||
       !requireNamespace("flextable", quietly = TRUE))) {
    stop("Packages `officer` and `flextable` are required for Word reports.", call. = FALSE)
  }
  invisible(TRUE)
}

.pm_is_report_collection <- function(x) {
  if (!is.list(x) || is.data.frame(x) || .pm_is_common_result(x) ||
      inherits(x, c("psychomatic_scores", "psychomatic_audit")) ||
      length(x) == 0L) {
    return(FALSE)
  }
  TRUE
}

.pm_is_supported_report_object <- function(x) {
  (is.data.frame(x) && !inherits(x, "psychomatic_item_screen")) ||
    .pm_is_supported_audit_result(x)
}
