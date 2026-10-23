% Traducción editorial al español.
% Fuente: man/desc_auto.Rd
\name{desc_auto}
\alias{desc_auto}
\title{Estadísticos descriptivos automatizados para ítems psicométricos}
\usage{
desc_auto(
  data,
  digits = 3,
  language = "eng",
  report = FALSE,
  file_name = NULL,
  overwrite = FALSE,
  format = "excel",
  verbose = TRUE
)
}
\arguments{
\item{data}{Data frame, matriz u objeto \code{psychomatic_spec} que contiene respuestas a ítems.}

\item{digits}{Número de decimales utilizados en el informe. Los valores numéricos devueltos no se redondean.}

\item{language}{Idioma del informe: \code{"esp"} o \code{"eng"}.}

\item{report}{Si es \code{TRUE}, guarda un informe en el formato solicitado.}

\item{file_name}{Ruta opcional de salida Excel. Si es \code{NULL}, se genera un nombre con marca temporal en el directorio de trabajo (\code{getwd()}) cuando \code{report = TRUE}.}

\item{overwrite}{Si es \code{TRUE}, se permite reemplazar un archivo de informe existente.}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{verbose}{Si es \code{TRUE}, muestra los resultados directamente, incluso cuando se asignan a un objeto.}
}
\value{
Un data frame compacto con una fila por ítem. En el atributo \code{response_percentages} se encuentra disponible un data frame en formato largo con las columnas \code{item}, \code{category}, \code{n} y \code{percent}.
}
\description{
Calcula estadísticos descriptivos compactos para columnas de ítems numéricos y ordenados. Los porcentajes de respuesta ordinal se conservan en una tabla separada en formato largo. Los resultados pueden exportarse opcionalmente a Excel.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
desc_auto(MY_DATA[, 3:14])
}

}
