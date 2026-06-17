##########################################################

#' Automated descriptives for numeric items
#'
#' Computes descriptive statistics and response percentages for numeric item
#' columns in a data frame. Results can optionally be exported to Excel.
#'
#' @param data Data frame containing the item responses.
#' @param digits Number of decimal places used in the report.
#' @param language Report language: `"esp"` or `"eng"`.
#' @param report If `TRUE`, exports an Excel file.
#'
#' @return Data frame with descriptive statistics and category percentages.
#'
#' @examples
#' items <- data.frame(
#'   item1 = c(1, 2, 3, 4, 5),
#'   item2 = c(2, 2, 3, 4, 5)
#' )
#' desc_auto(items, language = "eng")
#'
#' @export
desc_auto <- function(data, digits = 2, language = "esp", report = FALSE) {

  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("Package `psych` is required. Install it with install.packages('psych').",
         call. = FALSE)
  }

  language <- match.arg(language, choices = c("esp", "eng"))

  # Labels by output language.
  lbl <- if (language == "esp") {
    list(
      item = "Item", media = "Media", de = "DE",
      asimetria = "Asimetr\u00eda", curtosis = "Curtosis",
      pct_prefix = "% ", perdidos = "Casos perdidos/celdas en blanco",
      hoja_desc = "Descriptivos", hoja_notas = "Notas",
      col_nota = "Nota", guardado = "Reporte guardado en: ",
      no_num = "No se encontraron columnas num\u00e9ricas en la base de datos."
    )
  } else {
    list(
      item = "Item", media = "Mean", de = "SD",
      asimetria = "Skewness", curtosis = "Kurtosis",
      pct_prefix = "% ", perdidos = "Missing values/blank cells",
      hoja_desc = "Descriptives", hoja_notas = "Notes",
      col_nota = "Note", guardado = "Report saved at: ",
      no_num = "No numeric columns found in the dataset."
    )
  }

  # Select numeric columns only.
  data_num <- data[, sapply(data, is.numeric), drop = FALSE]

  if (ncol(data_num) == 0) stop(lbl$no_num)

  # Descriptive statistics from psych::describe().
  desc <- psych::describe(data_num)

  stats <- data.frame(
    V1 = rownames(desc),
    V2 = round(desc$mean, digits),
    V3 = round(desc$sd, digits),
    V4 = round(desc$skew, digits),
    V5 = round(desc$kurtosis, digits),
    row.names = NULL
  )
  colnames(stats) <- c(lbl$item, lbl$media, lbl$de, lbl$asimetria, lbl$curtosis)

  # Response percentages by category.
  all_values <- sort(unique(unlist(data_num)))

  pct <- sapply(data_num, function(x) {
    tbl <- table(factor(x, levels = all_values))
    # Use at least two decimals for percentages.
    round(prop.table(tbl) * 100, max(digits, 2))
  })

  # Keep matrix shape when there is a single item.
  if (!is.matrix(pct)) pct <- matrix(pct, nrow = length(all_values))

  pct_df <- as.data.frame(t(pct))
  colnames(pct_df) <- paste0(lbl$pct_prefix, all_values)
  pct_df[[lbl$item]] <- names(data_num)
  rownames(pct_df) <- NULL

  # Combine statistics and percentages.
  result <- merge(stats, pct_df, by = lbl$item)
  result <- result[match(names(data_num), result[[lbl$item]]), ]
  rownames(result) <- NULL

  # Count missing values / blank cells.
  n_perdidos <- sum(is.na(data_num))
  message(lbl$perdidos, " = ", n_perdidos)

  # Export to Excel if report = TRUE.
  if (report) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package `writexl` is required for export. Install it with install.packages('writexl').",
           call. = FALSE)
    }
    nota <- setNames(
      data.frame(paste0(lbl$perdidos, " = ", n_perdidos)),
      lbl$col_nota
    )
    file_prefix <- if (language == "esp") "Descriptivos_" else "Descriptives_"
    archivo <- paste0(file_prefix, format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    writexl::write_xlsx(setNames(list(result, nota), c(lbl$hoja_desc, lbl$hoja_notas)), archivo)
    message(lbl$guardado, archivo)
  }

  return(result)
}
