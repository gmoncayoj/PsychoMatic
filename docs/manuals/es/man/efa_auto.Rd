% Traducción editorial al español; fuente canónica: man/efa_auto.Rd.
\name{efa_auto}
\alias{efa_auto}
\title{Análisis factorial exploratorio automatizado}
\usage{
efa_auto(
  data,
  rotation = c("oblique", "orthogonal"),
  min_loading = 0.3,
  min_communality = 0.3,
  min_cross_loading_diff = 0.15,
  max_iter = 20,
  verbose = TRUE,
  export_format = NULL,
  file_name = NULL,
  language = "eng",
  seed = NULL,
  refine = c("remove", "suggest"),
  parallel_iter = 500L,
  parallel_quantile = 0.95,
  repair = c("none", "nearPD"),
  digits = 3,
  report = FALSE,
  format = "excel",
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{data}{Marco de datos, matriz o \code{psychomatic_spec} que contiene las respuestas a los ítems.}

\item{rotation}{Familia de rotación: \code{"oblique"} u \code{"orthogonal"}.}

\item{min_loading}{Carga factorial absoluta mínima aceptable.}

\item{min_communality}{Comunalidad mínima aceptable.}

\item{min_cross_loading_diff}{Diferencia mínima entre las dos cargas más grandes utilizada para detectar cargas cruzadas.}

\item{max_iter}{Número máximo de pasos de eliminación de ítems.}

\item{verbose}{Si es \code{TRUE}, imprime los resultados directamente, incluso cuando se asignan a un objeto.}

\item{export_format}{Formato de exportación opcional que se pasa a \code{export_efa()}.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca de tiempo en el directorio de trabajo (\code{getwd()}).}

\item{language}{Idioma de salida: \code{"esp"} o \code{"eng"}.}

\item{seed}{Semilla entera no negativa opcional para el análisis paralelo.}

\item{refine}{\code{"remove"} conserva el flujo histórico de eliminación automatizada; \code{"suggest"} conserva todos los ítems e informa las señales de alerta detectadas.}

\item{parallel_iter}{Número de simulaciones utilizadas por el análisis paralelo.}

\item{parallel_quantile}{Cuantil de referencia utilizado por el análisis paralelo.}

\item{repair}{Reparación de la matriz de correlaciones. \code{"none"} nunca modifica la matriz; \code{"nearPD"} activa explícitamente la reparación mediante la matriz definida positiva más cercana.}

\item{digits}{Lugares decimales en la consola y los informes; los valores almacenados conservan toda su precisión.}

\item{report}{Si es \code{TRUE}, también guarda un informe destinado al lector.}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{...}{Alias de argumentos en español declarados obsoletos.}
}
\value{
Un objeto de las clases \code{efa_auto} y \code{psychomatic_result}. Los campos
heredados siguen disponibles durante el ciclo de compatibilidad de la versión 1.x.
}
\description{
Ejecuta un flujo compacto de análisis factorial exploratorio: comprobación de ítems,
selección de correlaciones, diagnósticos de Mardia, análisis paralelo, ajuste factorial,
evaluación de ítems y resúmenes de fiabilidad. La evaluación puede sugerir cambios o
eliminar un ítem por vez. Cada modificación va seguida de un reajuste completo
que reclasifica el conjunto actual de ítems y determina su método de correlación.
Omega utiliza el mismo modelo unifactorial independiente por dimensión que
\code{\link[=reliability_auto]{reliability_auto()}}, no las cargas de la solución exploratoria.
}
\examples{
# Reemplace MY_DATA por su marco de datos y adapte los nombres o las columnas de los ítems.
\dontrun{
efa_auto(MY_DATA[, 3:14])
}

}
