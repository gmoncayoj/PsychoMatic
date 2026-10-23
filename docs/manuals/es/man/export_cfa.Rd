% Traducción editorial al español.
% Fuente: man/export_cfa.Rd
\name{export_cfa}
\alias{export_cfa}
\title{Exportación de resultados de AFC a Excel o Word}
\usage{
export_cfa(
  result,
  format = c("excel", "word"),
  file_name = NULL,
  digits = 3,
  overwrite = FALSE
)
}
\arguments{
\item{result}{Objeto de clase \code{"cfa_auto"} devuelto por \code{cfa_auto()}.}

\item{format}{\code{"excel"} o \code{"word"} (pueden proporcionarse ambos).}

\item{file_name}{Nombre base opcional del archivo, sin extensión. Si es \code{NULL}, se genera un nombre único con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{digits}{Número de decimales utilizados para dar formato a los valores del informe. Los valores numéricos originales se conservan sin cambios en los archivos Excel.}

\item{overwrite}{Si es \code{TRUE}, se permite reemplazar los archivos de salida existentes.}
}
\value{
Ruta de salida normalizada, de forma invisible. Cuando se solicitan ambos formatos, devuelve un vector de caracteres con nombres; la salida Excel incluye un atributo \code{sheet_map}.
}
\description{
Exporta resultados de AFC a Excel o Word.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
model <- "F1 =~ R1 + R2 + R3 + R4"
result <- cfa_auto(MY_DATA, model = model)
export_cfa(result, format = "word", file_name = "CFA.docx")
}

}
