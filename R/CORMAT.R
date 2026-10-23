#' Inter-item correlation matrix
#'
#' Computes an inter-item correlation matrix using Pearson, polychoric, or
#' tetrachoric correlations. Optionally exports the matrix to an Excel file.
#'
#' @param data Data frame or matrix containing the items.
#' @param type Defaults to `"auto"`: dichotomous items use tetrachoric;
#'   ordinal items with 2-5 categories use polychoric. With 6+ categories,
#'   ordinal items use polychoric when absolute moment skewness exceeds 1,
#'   absolute excess kurtosis exceeds 2, or endpoint mass is at least 0.50;
#'   otherwise Pearson is used. Numeric integer items with at most 10 observed
#'   categories are inferred ordinal; ordered factors declare this explicitly.
#'   These shape flags are operational rules, not universal literature cutoffs.
#'   Mixed low-category/continuous inputs require explicit selection. Manual
#'   options are `"pearson"`, `"poly"`, and `"tetra"`.
#' @param ordinal_levels Optional uniquely named list giving the ordered levels
#'   for unordered categorical items used with `type = "poly"` or `"tetra"`.
#' @inheritParams efa_auto
#' @param report If `TRUE`, save a report in the requested format.
#' @param file_name Optional Excel file name. If `NULL`, a timestamped name is
#'   generated in the working directory (`getwd()`).
#' @param digits Number of decimal places used in the Excel report.
#' @param overwrite If `TRUE`, an existing report file may be replaced.
#'
#' @return A correlation matrix with pairwise sample sizes, requested and
#'   effective correlation types, and validation diagnostics as attributes.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' cormat(MY_DATA[, 3:14])
#' }
#'
#' @export
cormat <- function(data,
                   type = c("auto", "pearson", "poly", "tetra"),
                   report = FALSE,
                   file_name = NULL,
                   digits = 3,
                   ordinal_levels = NULL,
                   overwrite = FALSE, format = "excel", verbose = TRUE,
                   language = "eng") {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  type <- match.arg(type)
  requested_type <- type
  language <- .pm_normalize_language(language)

  if (inherits(data, "psychomatic_spec")) {
    specification <- data
    data <- specification$data[specification$item_index]
    names(data) <- specification$items
  } else {
    data <- .pm_as_data_frame(data)
  }

  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("'report' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  digits <- .pm_validate_scalar_number(
    digits,
    name = "digits",
    lower = 0,
    integer = TRUE
  )

  if (ncol(data) < 2L) {
    stop("'data' must contain at least two items.", call. = FALSE)
  }

  selection <- .cormat_select(data, type, ordinal_levels, language)
  type <- selection$type
  prepared <- selection$data
  validation <- check_items(prepared, language = "eng")
  .cormat_validate_observed_variation(prepared, validation$issues)
  pairwise_n <- validation$diagnostics$pairwise_n

  correlation_matrix <- switch(
    type,
    pearson = stats::cor(prepared, use = "pairwise.complete.obs"),
    poly = .cormat_extract_rho(polychoric(prepared, smooth = FALSE, progress = FALSE)),
    tetra = .cormat_extract_rho(tetrachoric(prepared, smooth = FALSE))
  )

  correlation_matrix <- as.matrix(correlation_matrix)
  rownames(correlation_matrix) <- colnames(prepared)
  colnames(correlation_matrix) <- colnames(prepared)
  attr(correlation_matrix, "pairwise_n") <- pairwise_n
  attr(correlation_matrix, "requested_type") <- requested_type
  attr(correlation_matrix, "effective_type") <- type
  attr(correlation_matrix, "diagnostics") <- validation$issues
  attr(correlation_matrix, "selection") <- selection$profile
  attr(correlation_matrix, "rationale") <- selection$rationale
  attr(correlation_matrix, "references") <- .cormat_references(type)

  .pm_finish_analysis(correlation_matrix, output, language)
}

.cormat_prepare_data <- function(data, type, ordinal_levels = NULL) {
  data <- .pm_as_data_frame(data)
  if (is.null(names(data)) || any(!nzchar(names(data)))) {
    missing_names <- is.na(names(data)) | !nzchar(names(data))
    names(data)[missing_names] <- paste0("Item_", which(missing_names))
  }

  if (type == "pearson") {
    numeric_cols <- vapply(data, function(x) is.numeric(x) || is.logical(x), logical(1))
    if (!all(numeric_cols)) {
      stop("Pearson correlations require all columns in 'data' to be numeric or logical.", call. = FALSE)
    }
    return(as.data.frame(lapply(data, as.numeric), check.names = FALSE))
  }

  ordinal <- .pm_apply_ordinal_levels(
    data,
    seq_along(data),
    names(data),
    ordinal_levels
  )
  data <- ordinal$data
  prepared <- as.data.frame(
    lapply(data, .cormat_as_ordered_numeric),
    check.names = FALSE,
    optional = TRUE
  )
  names(prepared) <- names(data)

  if (type == "tetra") {
    binary_cols <- vapply(prepared, function(x) length(unique(stats::na.omit(x))) <= 2L, logical(1))
    if (!all(binary_cols)) {
      stop("Tetrachoric correlations require every item to have at most two observed categories.", call. = FALSE)
    }
  }

  prepared
}

.cormat_as_ordered_numeric <- function(x) {
  if (is.ordered(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x) && !inherits(x, c("Date", "POSIXt"))) {
    return(as.numeric(x))
  }

  if (identical(.pm_measurement_level(x), "ambiguous")) {
    stop(
      paste0(
        "Unordered factors, character, and logical items require explicit ",
        "`ordinal_levels` mapping or ordered factors."
      ),
      call. = FALSE
    )
  }

  stop("Ordinal items must be numeric codes or ordered factors.", call. = FALSE)
}

.cormat_extract_rho <- function(result) {
  if (is.matrix(result)) {
    return(result)
  }

  if (is.list(result) && !is.null(result$rho)) {
    return(result$rho)
  }

  stop("Could not extract a correlation matrix from the selected method.", call. = FALSE)
}

.cormat_validate_observed_variation <- function(data, diagnostics = NULL) {
  if (is.null(diagnostics)) {
    diagnostics <- check_items(data, language = "eng")$issues
  }
  if (any(diagnostics$code %in% c("all_missing", "zero_variance"))) {
    stop(
      "All items must have at least two observed values.",
      call. = FALSE
    )
  }
  if (any(diagnostics$code == "non_finite")) {
    stop("All items must contain only finite values or missing values.", call. = FALSE)
  }

  invisible(data)
}

.cormat_export_excel <- function(correlation_matrix,
                                 type,
                                 file_name = NULL,
                                 digits = 3,
                                 overwrite = FALSE) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required to export Excel reports.", call. = FALSE)
  }

  createWorkbook <- openxlsx::createWorkbook
  addWorksheet <- openxlsx::addWorksheet
  writeData <- openxlsx::writeData
  freezePane <- openxlsx::freezePane
  setColWidths <- openxlsx::setColWidths
  saveWorkbook <- openxlsx::saveWorkbook

  if (is.null(file_name)) {
    file_name <- paste0("Correlation_Matrix_", type, "_", base::format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  file_name <- .pm_safe_output_path(file_name, extension = ".xlsx", overwrite = overwrite)

  report_matrix <- round(correlation_matrix, digits = digits)
  report_df <- data.frame(Item = rownames(report_matrix), report_matrix, check.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "Correlation Matrix")
  writeData(wb, "Correlation Matrix", report_df)
  freezePane(wb, "Correlation Matrix", firstRow = TRUE, firstCol = TRUE)
  setColWidths(wb, "Correlation Matrix", cols = seq_len(ncol(report_df)), widths = "auto")
  saveWorkbook(wb, file_name, overwrite = overwrite)

  invisible(file_name)
}
