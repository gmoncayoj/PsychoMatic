% Traducción editorial al español.
% Fuente: man/compare_models_auto.Rd
\name{compare_models_auto}
\alias{compare_models_auto}
\title{Comparación de modelos alternativos de AFC}
\usage{
compare_models_auto(
  data,
  models,
  ordered = NULL,
  estimator = NULL,
  language = c("eng", "esp"),
  digits = 3,
  report = FALSE,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE,
  fit_variant = c("scaled", "robust", "base")
)
}
\arguments{
\item{data}{Data frame que contiene las variables observadas.}

\item{models}{Lista con nombres o vector de caracteres con nombres que contiene cadenas de modelos lavaan.}

\item{ordered}{Variables ordenadas opcionales que se pasan a \code{lavaan::cfa()}. Si es \code{NULL} y se utiliza un estimador categórico como WLSMV, los ítems observados referenciados en \code{models} se tratan como ordenados.}

\item{estimator}{Estimador que se pasa a \code{lavaan::cfa()}. \code{NULL} (predeterminado) o \code{"auto"} utiliza los mismos criterios de número de categorías, tamaño muestral y Mardia que \code{cfa_auto()}, aplicados a los indicadores observados del conjunto de modelos candidatos.}

\item{language}{Idioma de la recomendación: \code{"eng"} o \code{"esp"}.}

\item{digits}{Decimales en la consola y los informes; los valores almacenados conservan toda su precisión.}

\item{report}{Si es \code{TRUE}, también guarda un informe dirigido al lector.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, muestra los resultados directamente, incluso cuando se asignan a un objeto.}

\item{fit_variant}{Variante de los índices de ajuste para el informe y la comparación: \code{"scaled"} (predeterminado), \code{"robust"} o \code{"base"}. SRMR se mantiene en su versión base. Los estimadores sin corrección utilizan índices base con la configuración predeterminada scaled. En los demás casos, los índices solicitados no disponibles permanecen ausentes, sin sustitución.}
}
\value{
Una lista de clase \code{psychomatic_model_comparison}.
}
\description{
Ajusta un conjunto de modelos de AFC de lavaan con nombres y devuelve una tabla comparativa compacta. La recomendación es únicamente descriptiva y favorece los modelos convergentes con mejores índices de ajuste comparativo. Los criterios de información se presentan cuando están disponibles, pero no son obligatorios ni se utilizan para estimadores categóricos del tipo WLSMV.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
models <- list(one_factor = "F =~ R1 + R2 + R3 + R4 + R5 + R6",
               two_factors = "F1 =~ R1 + R2 + R3\nF2 =~ R4 + R5 + R6")
compare_models_auto(MY_DATA, models = models)
compare_models_auto(MY_DATA, models = models, fit_variant = "robust")
}

}
