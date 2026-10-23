##########################################################

#' Automated descriptives for psychometric items
#'
#' Computes compact descriptive statistics for numeric and ordered item columns.
#' Ordinal response percentages are retained as a separate long table. Results
#' can optionally be exported to Excel.
#'
#' @param data Data frame, matrix, or `psychomatic_spec` containing item
#'   responses.
#' @param digits Number of decimal places used in the report. Returned numeric
#'   values are not rounded.
#' @param language Report language: `"esp"` or `"eng"`.
#' @inheritParams efa_auto
#' @param report If `TRUE`, save a report in the requested format.
#' @param file_name Optional Excel output path. If `NULL`, a timestamped name is
#'   generated in the working directory (`getwd()`) when `report = TRUE`.
#' @param overwrite If `TRUE`, an existing report file may be replaced.
#'
#' @return A compact data frame with one row per item. A long data frame with
#'   columns `item`, `category`, `n`, and `percent` is available from the
#'   `response_percentages` attribute.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' desc_auto(MY_DATA[, 3:14])
#' }
#'
#' @export
desc_auto <- function(data,
                      digits = 3,
                      language = "eng",
                      report = FALSE,
                      file_name = NULL,
                      overwrite = FALSE, format = "excel", verbose = TRUE) {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  language <- .pm_normalize_language(language)
  digits <- .pm_validate_scalar_number(
    digits,
    name = "digits",
    lower = 0,
    integer = TRUE
  )
  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("`report` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  lbl <- if (language == "esp") {
    list(
      item = "Item", n = "N", missing = "Perdidos",
      media = "Media", de = "DE", min = "M\u00ednimo", max = "M\u00e1ximo",
      asimetria = "Asimetr\u00eda", curtosis = "Curtosis",
      perdidos = "Casos perdidos/celdas en blanco",
      hoja_desc = "Descriptivos", hoja_pct = "Porcentajes", hoja_notas = "Notas",
      col_nota = "Nota", guardado = "Reporte guardado en: ",
      no_num = "No se encontraron columnas num\u00e9ricas en la base de datos."
    )
  } else {
    list(
      item = "Item", n = "N", missing = "Missing",
      media = "Mean", de = "SD", min = "Min", max = "Max",
      asimetria = "Skewness", curtosis = "Kurtosis",
      perdidos = "Missing values/blank cells",
      hoja_desc = "Descriptives", hoja_pct = "Response percentages", hoja_notas = "Notes",
      col_nota = "Note", guardado = "Report saved at: ",
      no_num = "No numeric columns found in the dataset."
    )
  }

  if (inherits(data, "psychomatic_spec")) {
    specification <- data
    item_data <- specification$data[specification$item_index]
    names(item_data) <- specification$items
    measurement_level <- unname(specification$measurement_level)
  } else {
    input <- .pm_as_data_frame(data)
    eligible <- vapply(
      input,
      function(x) is.numeric(x) || is.ordered(x),
      logical(1)
    )
    item_data <- input[eligible]
    measurement_level <- vapply(item_data, .pm_measurement_level, character(1))
  }

  if (ncol(item_data) == 0L) stop(lbl$no_num, call. = FALSE)

  stats_rows <- vector("list", length(item_data))
  percentage_rows <- list()
  percentage_index <- 0L
  for (i in seq_along(item_data)) {
    x <- item_data[[i]]
    available <- .pm_available_values(x)
    numeric_x <- if (is.ordered(x)) as.integer(x) else as.numeric(x)
    observed <- numeric_x[available]
    item_n <- length(observed)

    stats_rows[[i]] <- data.frame(
      item = names(item_data)[i],
      n = as.integer(item_n),
      missing = as.integer(length(x) - item_n),
      mean = if (item_n) mean(observed) else NA_real_,
      sd = if (item_n > 1L) stats::sd(observed) else NA_real_,
      min = if (item_n) min(observed) else NA_real_,
      max = if (item_n) max(observed) else NA_real_,
      skewness = .pm_skewness(observed),
      kurtosis = .pm_kurtosis(observed),
      stringsAsFactors = FALSE
    )

    if (identical(measurement_level[[i]], "ordinal")) {
      categories <- if (is.ordered(x)) {
        levels(x)
      } else {
        as.character(sort(unique(x[available])))
      }
      counts <- table(
        factor(as.character(x), levels = categories),
        useNA = "no"
      )
      percentage_index <- percentage_index + 1L
      percentage_rows[[percentage_index]] <- data.frame(
        item = rep(names(item_data)[i], length(categories)),
        category = categories,
        n = as.integer(counts),
        percent = if (sum(counts)) as.numeric(counts) / sum(counts) * 100 else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  internal <- do.call(rbind, stats_rows)
  rownames(internal) <- NULL
  response_percentages <- if (length(percentage_rows)) {
    out <- do.call(rbind, percentage_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      item = character(),
      category = character(),
      n = integer(),
      percent = numeric(),
      stringsAsFactors = FALSE
    )
  }

  result <- internal
  names(result) <- c(
    lbl$item, lbl$n, lbl$missing, lbl$media, lbl$de,
    lbl$min, lbl$max, lbl$asimetria, lbl$curtosis
  )
  attr(result, "response_percentages") <- response_percentages

  .pm_finish_analysis(result, output, language)
}
