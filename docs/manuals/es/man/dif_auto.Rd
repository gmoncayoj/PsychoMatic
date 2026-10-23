% Traducción editorial al español.
% Fuente: man/dif_auto.Rd
\name{dif_auto}
\alias{dif_auto}
\title{Evaluación preliminar del funcionamiento diferencial de los ítems por grupo}
\usage{
dif_auto(
  data,
  group,
  items = NULL,
  alpha = 0.01,
  language = c("eng", "esp"),
  digits = 3,
  polytomous = c("not_evaluable", "median"),
  report = FALSE,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Data frame o matriz que contiene respuestas a ítems. Si \code{group} es el nombre de una columna de \code{data}, esa columna se utiliza como variable de agrupación.}

\item{group}{Vector de agrupación o nombre de la columna de agrupación en \code{data}.}

\item{items}{Nombres opcionales de los ítems. De forma predeterminada, se utilizan todas las columnas de ítems numéricos excepto la columna de agrupación.}

\item{alpha}{Nivel de significación utilizado para señalar DIF.}

\item{language}{Idioma de la interpretación: \code{"eng"} o \code{"esp"}.}

\item{digits}{Número de decimales utilizados en las tablas de interpretación.}

\item{polytomous}{Tratamiento de los ítems focales politómicos. La opción predeterminada, \code{"not_evaluable"}, no los dicotomiza. Utilice \code{"median"} únicamente para activar una evaluación preliminar de compatibilidad, que se identifica como tal y se acompaña de una advertencia metodológica.}

\item{report}{Si es \code{TRUE}, también guarda un informe dirigido al lector.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, muestra los resultados directamente, incluso cuando se asignan a un objeto.}
}
\value{
Una lista de clase \code{psychomatic_dif} con los resultados de la evaluación preliminar, la interpretación por ítem (\code{dif_type}, decision, severity y recommendation), un resumen general de DIF y referencias de respaldo en formato APA 7.
}
\description{
Realiza una evaluación preliminar sencilla del DIF mediante regresión logística. De forma predeterminada, los ítems politómicos se informan como no evaluables porque la dicotomización puede ocultar DIF específico de las categorías. La dicotomización por la mediana solo está disponible como opción explícita de compatibilidad y debe complementarse con un análisis específico de DIF mediante TRI o modelos ordinales.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
dif_auto(MY_DATA[, 3:14], group = MY_DATA$Country)
}

}
