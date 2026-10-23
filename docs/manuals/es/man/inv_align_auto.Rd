% Traducción editorial al español; fuente canónica: man/inv_align_auto.Rd.
\name{inv_align_auto}
\alias{inv_align_auto}
\title{Invariancia mediante alineación automatizada}
\usage{
inv_align_auto(
  lambda = NULL,
  nu = NULL,
  data = NULL,
  group = NULL,
  sampling_weights = NULL,
  wgt = NULL,
  language = c("eng", "esp"),
  model = c("2PM", "1PM"),
  align.scale = c(lambda = 0.4, nu = 0.2),
  align.pow = c(lambda = 0.25, nu = 0.25),
  parm_tol = NULL,
  noninvariance_cutoff = 25,
  center = FALSE,
  optimizer = "optim",
  fixed = NULL,
  meth = 1,
  eps = 0.001,
  vcov = NULL,
  digits = 3,
  config_args = list(),
  alignment_args = list(),
  report = FALSE,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{lambda, nu}{Matrices de parámetros de cargas e interceptos.}

\item{data}{Datos a nivel de ítem utilizados cuando no se proporcionan \code{lambda} y \code{nu}.}

\item{group}{Vector de agrupación asociado a \code{data}.}

\item{sampling_weights}{Pesos muestrales opcionales, uno por fila.}

\item{wgt}{Pesos opcionales para el procedimiento de alineación.}

\item{language}{Idioma del informe: \code{"esp"} o \code{"eng"}.}

\item{model}{Tipo de modelo para la etapa configural de alineación.}

\item{align.scale, align.pow, parm_tol}{Parámetros de ajuste de la alineación.}

\item{noninvariance_cutoff}{Porcentaje utilizado para señalar la falta de invariancia.}

\item{center, optimizer, fixed, meth, eps, vcov}{Parámetros avanzados del algoritmo.}

\item{digits}{Número de lugares decimales utilizados en el informe.}

\item{config_args}{Argumentos adicionales para la etapa configural.}

\item{alignment_args}{Argumentos adicionales para la etapa de alineación.}

\item{report}{Si es \code{TRUE}, también guarda un informe destinado al lector.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca de tiempo en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, imprime los resultados directamente, incluso cuando se asignan a un objeto.}
}
\value{
Objeto de la clase \code{inv_align_auto_result}.
}
\description{
Ejecuta la optimización de alineación a partir de matrices de cargas e interceptos (\code{lambda} y
\code{nu}) o de datos a nivel de ítem con un vector de agrupación.
}
\examples{
# Reemplace MY_DATA por su marco de datos y adapte los nombres o las columnas de los ítems.
\dontrun{
inv_align_auto(data = MY_DATA[, 28:39], group = MY_DATA$Country,
               language = "eng")
}

}
