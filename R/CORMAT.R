#' Inter-item correlation matrix
#'
#' Computes an inter-item correlation matrix using Pearson, polychoric, or
#' tetrachoric correlations. Optionally exports the matrix to an Excel file.
#'
#' @param data Data frame or matrix containing the items.
#' @param type Type of correlation to compute. Options are `"pearson"`, `"poly"`
#'   for polychoric correlations, and `"tetra"` for tetrachoric correlations.
#' @param report If `TRUE`, exports the matrix to an Excel file.
#' @param file_name Optional Excel file name. If `NULL`, a timestamped name is
#'   generated automatically.
#' @param digits Number of decimal places used in the Excel report.
#'
#' @return A correlation matrix.
#'
#' @examples
#' # cormat(my_data, type = "pearson")
#' # cormat(my_data, type = "poly", report = TRUE)
#'
#' @export
cormat <- function(data,
                   type = c("pearson", "poly", "tetra"),
                   report = FALSE,
                   file_name = NULL,
                   digits = 3) {
  type <- match.arg(type)

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix.", call. = FALSE)
  }

  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("'report' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("'digits' must be a non-negative number.", call. = FALSE)
  }

  data <- as.data.frame(data)

  if (ncol(data) < 2L) {
    stop("'data' must contain at least two items.", call. = FALSE)
  }

  prepared <- .cormat_prepare_data(data, type)

  correlation_matrix <- switch(
    type,
    pearson = stats::cor(prepared, use = "pairwise.complete.obs"),
    poly = .cormat_extract_rho(polychoric(prepared)),
    tetra = .cormat_extract_rho(tetrachoric(prepared))
  )

  correlation_matrix <- as.matrix(correlation_matrix)
  rownames(correlation_matrix) <- colnames(prepared)
  colnames(correlation_matrix) <- colnames(prepared)

  if (isTRUE(report)) {
    .cormat_export_excel(correlation_matrix, type = type, file_name = file_name, digits = digits)
  }

  correlation_matrix
}

.cormat_prepare_data <- function(data, type) {
  if (is.null(names(data)) || any(!nzchar(names(data)))) {
    names(data) <- paste0("Item_", seq_len(ncol(data)))
  }

  if (type == "pearson") {
    numeric_cols <- vapply(data, function(x) is.numeric(x) || is.logical(x), logical(1))
    if (!all(numeric_cols)) {
      stop("Pearson correlations require all columns in 'data' to be numeric or logical.", call. = FALSE)
    }
    return(as.data.frame(lapply(data, as.numeric), check.names = FALSE))
  }

  prepared <- as.data.frame(lapply(data, .cormat_as_ordered_numeric), check.names = FALSE)
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
  if (is.logical(x)) {
    return(as.numeric(x))
  }

  if (is.factor(x) || is.character(x)) {
    return(as.numeric(factor(x, ordered = TRUE)))
  }

  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }

  stop("Items must be numeric, logical, factor, ordered factor, or character.", call. = FALSE)
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

.cormat_export_excel <- function(correlation_matrix, type, file_name = NULL, digits = 3) {
  if (is.null(file_name)) {
    file_name <- paste0("Correlation_Matrix_", type, "_", base::format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  if (!grepl("\\.xlsx$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".xlsx")
  }

  report_matrix <- round(correlation_matrix, digits = digits)
  report_df <- data.frame(Item = rownames(report_matrix), report_matrix, check.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "Correlation Matrix")
  writeData(wb, "Correlation Matrix", report_df)
  freezePane(wb, "Correlation Matrix", firstRow = TRUE, firstCol = TRUE)
  setColWidths(wb, "Correlation Matrix", cols = seq_len(ncol(report_df)), widths = "auto")
  saveWorkbook(wb, file_name, overwrite = TRUE)

  invisible(normalizePath(file_name))
}
