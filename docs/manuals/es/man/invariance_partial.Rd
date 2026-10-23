% Traducción editorial al español; fuente canónica: man/invariance_partial.Rd.
\name{invariance_partial}
\alias{invariance_partial}
\title{Ajuste de invariancia de medida parcial dirigido por el investigador}
\usage{
invariance_partial(
  data,
  group,
  model,
  free,
  rationale,
  stage = "scalar",
  measurement = "auto",
  estimator = "auto",
  ordered = NULL,
  criteria = "auto",
  min_anchor_items = 2L,
  report = FALSE,
  language = "eng",
  digits = 3,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE,
  fit_variant = c("scaled", "robust"),
  workflow = c("factorial", "ordinal"),
  missing = "auto",
  min_group_n = 20L
)
}
\arguments{
\item{data}{Marco de datos, matriz o \code{psychomatic_spec}.}

\item{group}{Nombre, símbolo o índice de la columna de agrupación, o un vector con un valor por fila.}

\item{model}{Un modelo de medida de primer orden de lavaan.}

\item{free}{Expresiones canónicas de los parámetros que se liberarán, proporcionadas como vector de caracteres o como lista con nombres de familias.}

\item{rationale}{Justificación sustantiva obligatoria, ya sea un único valor global o un vector de caracteres cuyos nombres coincidan exactamente con las expresiones.}

\item{stage}{Etapa objetivo: \code{"metric"}, \code{"scalar"} o \code{"strict"} en el flujo factorial predeterminado. El flujo ordinal explícito admite además \code{"thresholds"} y \code{"loadings"}; los datos ordinales binarios requieren etapas conjuntas \code{"scalar"} o \code{"strict"}.}

\item{measurement}{Declaración (\code{"auto"}, \code{"continuous"}, \code{"ordinal"} o \code{"mixed"}) que se contrasta con los tipos de ítems preparados.}

\item{estimator}{Estimador que se pasa al motor de invariancia completa.}

\item{ordered}{Declaración de ítems ordenados aceptada por el motor de invariancia completa.}

\item{criteria}{Criterios de cambio en el ajuste aceptados por \code{factorial_invariance_auto()}.}

\item{min_anchor_items}{Número mínimo de ítems de anclaje con todas las restricciones por factor de primer orden; debe ser un entero de al menos dos.}

\item{report}{\code{FALSE}, \code{TRUE} o una ruta de salida de Excel. Los archivos existentes se protegen, salvo que \code{overwrite = TRUE} autorice explícitamente su reemplazo.}

\item{language}{Idioma de presentación: \code{"eng"} por defecto; \code{"esp"} para español. Un \code{NULL} explícito hereda una especificación o la opción del paquete.}

\item{digits}{Lugares decimales en la consola y los informes; los valores almacenados conservan toda su precisión.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca de tiempo en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, imprime los resultados directamente, incluso cuando se asignan a un objeto.}

\item{fit_variant}{\code{"scaled"} (predeterminada) o \code{"robust"} para CFI y RMSEA, sus cambios y las decisiones. Los valores robustos no disponibles no se sustituyen. El chi-cuadrado, los grados de libertad y el valor p corregidos siguen siendo escalados; SRMR sigue siendo de base.}

\item{workflow}{\code{"factorial"} (predeterminado) sigue la configuración nativa de lavaan y la secuencia de restricciones de cuatro etapas de \code{factorial_invariance_auto()}, con liberaciones \code{group.partial} especificadas por el investigador. WLSMV/ULSMV utilizan por defecto indicadores ordenados, identificación mediante un indicador de referencia y parametrización delta. \code{"ordinal"} conserva el procedimiento identificado theta/std.lv/Wu-Estabrook de semTools, incluidas las liberaciones de umbrales. Una liberación sin efecto o una comparación sin grados de libertad adicionales no pueden sustentar la invariancia.}

\item{missing}{Opción de datos faltantes aceptada por \code{factorial_invariance_auto()}.}

\item{min_group_n}{Tamaño efectivo mínimo por grupo; 20 por defecto.}
}
\value{
De forma invisible, un objeto de las clases
\code{psychomatic_invariance_partial} y \code{psychomatic_result}.
}
\description{
Ajusta una etapa objetivo de invariancia completa y una etapa parcial equivalente utilizando
liberaciones explícitas de parámetros justificadas por el investigador. La función automatiza
la extensa sintaxis de lavaan/semTools y evita deliberadamente las búsquedas automatizadas
de índices de modificación o la liberación iterativa.
}
\examples{
# Reemplace MY_DATA por su marco de datos y adapte los nombres o las columnas de los ítems.
\dontrun{
model <- "F1 =~ R1 + R2 + R3 + R4"
# Indicadores continuos; justifique el intercepto liberado antes de la prueba.
invariance_partial(MY_DATA, group = "Country", model = model,
                   free = "R2 ~ 1", estimator = "ML",
                   rationale = "A preregistered item-wording hypothesis")
}

}
