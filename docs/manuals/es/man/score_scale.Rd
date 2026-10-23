% Traducción editorial al español.
% Fuente: man/score_scale.Rd
\name{score_scale}
\alias{score_scale}
\title{Cálculo de puntuaciones de una o más escalas psicométricas}
\usage{
score_scale(
  data,
  keys = NULL,
  items = NULL,
  method = c("mean", "sum"),
  min_valid = 0.8,
  impute = c("none", "person_mean"),
  reverse = NULL,
  min_score = NULL,
  max_score = NULL
)
}
\arguments{
\item{data}{Data frame o matriz que contiene respuestas a ítems.}

\item{keys}{Lista con nombres en la que cada elemento contiene los nombres de los ítems de una escala. Si es \code{NULL}, se calculan las puntuaciones de \code{items} como una única escala denominada \code{total}.}

\item{items}{Nombres opcionales de los ítems utilizados cuando \code{keys = NULL}.}

\item{method}{Método de puntuación: \code{"mean"} o \code{"sum"}.}

\item{min_valid}{Número o proporción mínima de respuestas válidas a los ítems que se requiere para calcular una puntuación de escala.}

\item{impute}{Tratamiento de datos ausentes: \code{"none"} o \code{"person_mean"}.}

\item{reverse}{Vector de caracteres opcional con los ítems que se invertirán antes de calcular las puntuaciones.}

\item{min_score, max_score}{Rango opcional de puntuaciones utilizado cuando se proporciona \code{reverse}.}
}
\value{
Una lista de clase \code{psychomatic_scores} con tablas de puntuaciones y fiabilidad, y un \code{scale_map} que vincula claves únicas con etiquetas de visualización.
}
\description{
Calcula puntuaciones de escala mediante suma o media, con inversión opcional de puntuaciones, imputación opcional mediante la media de la persona y una regla de respuestas válidas mínimas.
}
\details{
Las filas incompletas admisibles se prorratean: las puntuaciones medias utilizan la media de los ítems disponibles y las puntuaciones sumadas multiplican esa media por la longitud de la escala. Los valores ausentes pueden permanecer como \code{NA}, pero las respuestas no finitas se rechazan antes de calcular las puntuaciones de cualquier escala. Las etiquetas de visualización duplicadas reciben claves únicas estables que se registran en \code{scale_map}.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
score_scale(MY_DATA, items = c("R1", "R2", "R3", "R4"))
}

}
