##########################################################

#' Descriptivos automatizados para items numericos
#'
#' Calcula estadisticos descriptivos y porcentajes de respuesta para las
#' columnas numericas de una base de datos. Opcionalmente exporta el resultado
#' a Excel.
#'
#' @param data Data frame con los datos.
#' @param digits Numero de decimales para el reporte.
#' @param language Idioma del reporte: `"esp"` o `"eng"`.
#' @param report Si es `TRUE`, exporta un archivo Excel.
#'
#' @return Data frame con los descriptivos y porcentajes por categoria.
#'
#' @examples
#' # desc_auto(mi_datos)
#'
#' @export
desc_auto <- function(data, digits = 2, language = "esp", report = FALSE) {

  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("El paquete 'psych' es necesario. Inst\u00e1lelo con install.packages('psych').")
  }

  language <- match.arg(language, choices = c("esp", "eng"))

  # Etiquetas seg\u00fan idioma
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

  # Seleccionar solo columnas num\u00e9ricas
  data_num <- data[, sapply(data, is.numeric), drop = FALSE]

  if (ncol(data_num) == 0) stop(lbl$no_num)

  # Estad\u00edsticos descriptivos con psych::describe
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

  # Porcentaje de respuestas por categor\u00eda
  all_values <- sort(unique(unlist(data_num)))

  pct <- sapply(data_num, function(x) {
    tbl <- table(factor(x, levels = all_values))
    # M\u00ednimo 2 decimales en porcentajes
    round(prop.table(tbl) * 100, max(digits, 2))
  })

  # Garantizar formato matricial (caso de un solo \u00edtem)
  if (!is.matrix(pct)) pct <- matrix(pct, nrow = length(all_values))

  pct_df <- as.data.frame(t(pct))
  colnames(pct_df) <- paste0(lbl$pct_prefix, all_values)
  pct_df[[lbl$item]] <- names(data_num)
  rownames(pct_df) <- NULL

  # Combinar estad\u00edsticos y porcentajes
  result <- merge(stats, pct_df, by = lbl$item)
  result <- result[match(names(data_num), result[[lbl$item]]), ]
  rownames(result) <- NULL

  # Detectar casos perdidos / celdas en blanco
  n_perdidos <- sum(is.na(data_num))
  message(lbl$perdidos, " = ", n_perdidos)

  # Exportar a Excel si report = TRUE
  if (report) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("El paquete 'writexl' es necesario para exportar. Inst\u00e1lelo con install.packages('writexl').")
    }
    nota <- setNames(
      data.frame(paste0(lbl$perdidos, " = ", n_perdidos)),
      lbl$col_nota
    )
    archivo <- paste0("Descriptivos_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    writexl::write_xlsx(setNames(list(result, nota), c(lbl$hoja_desc, lbl$hoja_notas)), archivo)
    message(lbl$guardado, archivo)
  }

  return(result)
}
