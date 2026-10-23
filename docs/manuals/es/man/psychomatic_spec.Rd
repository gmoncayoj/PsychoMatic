% Traducción editorial al español.
% Fuente: man/psychomatic_spec.Rd
\name{psychomatic_spec}
\alias{psychomatic_spec}
\title{Definición de una especificación de datos psicométricos}
\usage{
psychomatic_spec(
  data,
  items = NULL,
  metadata = NULL,
  level = "auto",
  ordinal_levels = NULL,
  reverse = NULL,
  ranges = NULL,
  group = NULL,
  missing = "available",
  language = "eng"
)
}
\arguments{
\item{data}{Un data frame o una matriz que contiene respuestas a ítems y, opcionalmente, columnas de metadatos o agrupación.}

\item{items}{Nombres de las columnas de ítems o índices enteros positivos. Si es \code{NULL}, se consideran todas las columnas admisibles excepto \code{metadata} y \code{group}; los posibles identificadores provocan un error en lugar de seleccionarse de forma silenciosa.}

\item{metadata}{Nombres opcionales de las columnas de metadatos o índices enteros positivos.}

\item{level}{Nivel de medición: \code{"auto"}, \code{"ordinal"} o \code{"continuous"}. Puede ser un escalar, estar alineado con \code{items} o ser un vector con nombres que especifique excepciones por ítem.}

\item{ordinal_levels}{Lista opcional con nombres únicos que especifica los niveles ordenados de los ítems seleccionados. Los niveles declarados tienen prioridad sobre \code{level} y deben cubrir cada valor observado no ausente.}

\item{reverse}{Nombres opcionales de los ítems o índices enteros positivos que identifican ítems con clave inversa. Registra las claves; no transforma las respuestas.}

\item{ranges}{Lista opcional con nombres únicos de vectores numéricos \code{c(min, max)}, o una matriz o data frame de dos columnas con nombres de ítems en las filas.}

\item{group}{Nombre opcional de una única columna de agrupación o índice entero positivo.}

\item{missing}{Escalar de caracteres que registra la política de datos ausentes. El valor predeterminado es \code{"available"}.}

\item{language}{Idioma de salida: \code{"eng"} de forma predeterminada; \code{"esp"} para español. Un \code{NULL} explícito utiliza \code{getOption("PsychoMatic.language", "eng")}.}
}
\value{
Un objeto de clases \code{psychomatic_spec} y \code{psychomatic_result}. Además de los campos comunes del contrato de resultados, los campos directos exponen \code{data}, \code{raw_data} antes de la correspondencia de valores, \code{items}, \code{measurement_level}, \code{ordinal_levels}, \code{reverse}, \code{ranges}, \code{group}, \code{metadata}, \code{missing} y \code{language}.
}
\description{
Selecciona columnas de ítems sin tratar de forma silenciosa los identificadores o las categorías no ordenadas como mediciones. Registra el orden de los ítems, los niveles de medición, las codificaciones ordinales, las claves inversas, los rangos, la agrupación, los metadatos, la política de datos ausentes, el idioma y las decisiones utilizadas para construir la especificación.
}
\details{
Con \code{level = "auto"}, los factores ordenados son ordinales. Las columnas numéricas ordinarias con entre dos y siete valores enteros finitos distintos son ordinales; las demás columnas numéricas que presentan variación son continuas. Los factores no ordenados, los vectores de caracteres y los vectores lógicos son ambiguos y nunca se ordenan alfabéticamente. Los \code{ordinal_levels} explícitos tienen prioridad sobre un \code{level} explícito por ítem, que a su vez tiene prioridad sobre la clasificación automática.

Cuando \code{items = NULL}, una variable candidata se trata como un posible identificador si su nombre normalizado contiene \code{id}, \code{case}, \code{participant} o \code{respondent}, si al menos el 95 por ciento de sus valores numéricos finitos son únicos, o si esos valores forman una secuencia entera estrictamente monótona. La función se detiene entonces y solicita una declaración explícita de \code{items} o \code{metadata}.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
psychomatic_spec(MY_DATA, items = c("R1", "R2", "R3", "R4"))
}

}
