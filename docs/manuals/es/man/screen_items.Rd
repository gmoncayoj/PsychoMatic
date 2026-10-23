% Traducción editorial al español.
% Fuente: man/screen_items.Rd
\name{screen_items}
\alias{screen_items}
\title{Evaluación de la calidad de los datos por ítem}
\usage{
screen_items(
  data,
  digits = 3,
  floor_cutoff = 0.8,
  ceiling_cutoff = 0.8,
  min_item_rest = -0.1,
  language = "eng",
  report = FALSE,
  format = "excel",
  file_name = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Data frame o matriz que contiene respuestas a ítems.}

\item{digits}{Número de decimales utilizados en la tabla devuelta.}

\item{floor_cutoff}{Proporción en la puntuación mínima utilizada para señalar efectos suelo.}

\item{ceiling_cutoff}{Proporción en la puntuación máxima utilizada para señalar efectos techo.}

\item{min_item_rest}{Correlación ítem-resto por debajo de la cual un ítem se señala como probablemente de clave inversa o inconsistente con los demás ítems.}

\item{language}{Idioma de salida: \code{"eng"} (predeterminado) o \code{"esp"}.}

\item{report}{Genera un informe o permite proporcionar una ruta de salida.}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{file_name}{Nombre del archivo de salida, sin extensión cuando se utilizan ambos formatos.}

\item{overwrite}{Permite reemplazar un informe existente (valor predeterminado: \code{FALSE}).}

\item{verbose}{Muestra los resultados directamente en la consola (valor predeterminado: \code{TRUE}).}
}
\value{
Un data frame de clase \code{psychomatic_item_screen}.
}
\description{
Calcula una tabla compacta de evaluación preliminar de ítems psicométricos numéricos: datos ausentes, varianza, forma de la distribución, proporciones de suelo y techo, correlaciones ítem-resto e ítems probablemente de clave inversa.
}
\details{
Los porcentajes mostrados se expresan en una escala de 0--100. El data frame devuelto conserva sus nombres originales de columnas en inglés y las proporciones de 0--1 por compatibilidad.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
screen_items(MY_DATA[, 3:14])
}

}
