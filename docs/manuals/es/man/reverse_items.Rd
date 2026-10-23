% Traducción editorial al español.
% Fuente: man/reverse_items.Rd
\name{reverse_items}
\alias{reverse_items}
\title{Inversión de las puntuaciones de ítems seleccionados}
\usage{
reverse_items(
  data,
  items = NULL,
  min_score = NULL,
  max_score = NULL,
  suffix = "_r",
  overwrite = FALSE,
  detect = FALSE
)
}
\arguments{
\item{data}{Data frame o matriz que contiene respuestas a ítems.}

\item{items}{Vector de caracteres con los nombres de los ítems que se invertirán. Es obligatorio salvo que \code{detect = TRUE}.}

\item{min_score, max_score}{Rango teórico opcional de puntuaciones. Si se omite, se utiliza el rango observado entre los ítems seleccionados.}

\item{suffix}{Sufijo utilizado para las columnas invertidas cuando \code{overwrite = FALSE}.}

\item{overwrite}{Si es \code{TRUE}, reemplaza las columnas originales. Si es \code{FALSE}, añade columnas invertidas.}

\item{detect}{Si es \code{TRUE}, invierte los ítems señalados por \code{screen_items()} como probablemente de clave inversa.}
}
\value{
Un data frame con los atributos \code{reversed_items} y \code{reverse_range}.
}
\description{
Invierte las puntuaciones de los ítems mediante un rango de puntuaciones definido por el usuario o el rango observado de los ítems seleccionados. Los ítems pueden seleccionarse explícitamente o detectarse a partir de correlaciones ítem-resto negativas.
}
\details{
Las selecciones repetidas se aplican una sola vez. Las respuestas no ausentes deben ser finitas y encontrarse dentro del rango efectivo de puntuación. Al añadir columnas invertidas, un nombre de destino existente se trata como una colisión y nunca se sobrescribe de forma silenciosa.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
reverse_items(MY_DATA, items = c("R2", "R4"), min_score = 1, max_score = 5)
}

}
