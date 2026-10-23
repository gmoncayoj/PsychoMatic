% Traducción editorial al español.
% Fuente: man/reliability_auto.Rd
\name{reliability_auto}
\alias{reliability_auto}
\title{Análisis automatizado de la fiabilidad}
\usage{
reliability_auto(
  data,
  items = NULL,
  level = "auto",
  model = NULL,
  report = FALSE,
  language = "eng",
  digits = 3,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Un data frame, una matriz o un objeto \code{\link[=psychomatic_spec]{psychomatic_spec()}}.}

\item{items}{Nombres opcionales de los ítems o índices positivos. Debe ser \code{NULL} cuando la entrada es un objeto \code{psychomatic_spec}.}

\item{level}{Nivel de medición que se pasa a \code{\link[=psychomatic_spec]{psychomatic_spec()}}: \code{"auto"}, \code{"continuous"} u \code{"ordinal"}.}

\item{model}{Modelo de AFC opcional ajustado con \code{lavaan}.}

\item{report}{Si es \code{TRUE}, envía el resultado estructurado a \code{\link[=psychomatic_report]{psychomatic_report()}}. No se escribe ningún archivo cuando es \code{FALSE}.}

\item{language}{Idioma de salida: \code{"eng"} de forma predeterminada; utilice \code{"esp"} para español. Un \code{NULL} explícito hereda el idioma de una especificación o de la opción del paquete.}

\item{digits}{Decimales en la consola y los informes; los valores almacenados conservan toda su precisión.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, muestra los resultados directamente, incluso cuando se asignan a un objeto.}
}
\value{
Un objeto de clases \code{psychomatic_reliability} y \code{psychomatic_result}.
}
\description{
Calcula coeficientes de fiabilidad manteniendo explícitos el estimando estadístico, la matriz de respuestas, la decisión de inversión de claves y el estado de evaluabilidad. Para las respuestas continuas se calcula el alfa de las puntuaciones observadas; para las respuestas ordinales también se calcula alfa a partir de una matriz policórica con nombres. Un modelo ajustado con \code{lavaan} añade la fiabilidad compuesta y la varianza media extraída mediante las API actuales de \code{semTools}, cuando ese paquete sugerido está instalado. Si no se proporciona un modelo, omega ajusta un único factor congenérico (ML para ítems continuos; WLSMV para ítems ordinales). El omega ordinal se refiere a la puntuación de respuesta observada utilizando la corrección de Green y Yang; el alfa ordinal se refiere a respuestas latentes. Este supuesto requiere una evaluación sustantiva. Alfa utiliza covarianzas por pares. Omega utiliza casos completos y se vuelve a ajustar tras la eliminación de cada ítem; las correlaciones ítem-resto utilizan casos completos en todos los ítems. El alfa tras eliminar un ítem utiliza la submatriz correspondiente de la matriz de la escala completa.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
reliability_auto(MY_DATA[, 3:14])
}

}
