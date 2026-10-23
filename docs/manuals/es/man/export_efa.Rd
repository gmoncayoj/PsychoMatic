% Traducción editorial al español.
% Fuente: man/export_efa.Rd
\name{export_efa}
\alias{export_efa}
\title{Exportación de resultados de AFE a Excel o Word}
\usage{
export_efa(
  result,
  format = c("excel", "word"),
  file_name = NULL,
  overwrite = FALSE
)
}
\arguments{
\item{result}{Objeto de clase \code{"efa_auto"} devuelto por \code{efa_auto()}.}

\item{format}{Formato de salida: \code{"excel"}, \code{"word"} o ambos.}

\item{file_name}{Ruta opcional de salida o ruta base. Si es \code{NULL}, se genera un nombre único con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{overwrite}{Si es \code{TRUE}, se permite reemplazar los archivos de salida existentes. De forma predeterminada, se impide el reemplazo.}
}
\value{
El objeto de entrada \code{result}, de forma invisible.
}
\description{
Exporta las tablas estructuradas de un resultado de \code{efa_auto()}. La función utiliza el generador de informes del paquete para que las exportaciones de AFE se mantengan coherentes con el flujo integrado de \code{psychomatic_report()}.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
result <- efa_auto(MY_DATA[, 3:14])
export_efa(result, format = "excel", file_name = "EFA.xlsx")
}

}
