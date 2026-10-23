% Traducción editorial al español.
% Fuente: man/irt_auto.Rd
\name{irt_auto}
\alias{irt_auto}
\title{Ajuste automatizado de modelos básicos de TRI}
\usage{
irt_auto(
  data,
  model = c("auto", "rasch", "2pl", "grm"),
  factors = 1,
  itemtype = NULL,
  technical = list(),
  language = c("eng", "esp"),
  digits = 3,
  report = FALSE,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Data frame o matriz que contiene respuestas a ítems.}

\item{model}{Modelo de TRI: \code{"auto"}, \code{"rasch"}, \code{"2pl"} o \code{"grm"}.}

\item{factors}{Número de dimensiones latentes que se pasa a \code{mirt::mirt()}.}

\item{itemtype}{Tipo de ítem opcional de \code{mirt}. Si es \code{NULL}, se selecciona a partir de \code{model}.}

\item{technical}{Argumentos técnicos opcionales que se pasan a \code{mirt::mirt()}.}

\item{language}{Idioma de la interpretación: \code{"eng"} o \code{"esp"}.}

\item{digits}{Número de decimales utilizados en las tablas de interpretación.}

\item{report}{Si es \code{TRUE}, también guarda un informe dirigido al lector.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, muestra los resultados directamente, incluso cuando se asignan a un objeto.}
}
\value{
Una lista de clase \code{psychomatic_irt} con el modelo ajustado, los parámetros de los ítems, la interpretación por ítem, un resumen del modelo y referencias de respaldo en formato APA 7.
}
\description{
Ajusta un modelo unidimensional de Rasch, 2PL o de respuesta graduada mediante \code{mirt}. La dependencia es opcional y solo se carga cuando se llama a esta función.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
irt_auto(MY_DATA[, 3:14])
}

}
