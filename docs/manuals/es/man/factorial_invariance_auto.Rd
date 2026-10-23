% Traducción editorial al español; fuente canónica: man/factorial_invariance_auto.Rd.
\name{factorial_invariance_auto}
\alias{factorial_invariance_auto}
\title{Invariancia de medida completa automatizada}
\usage{
factorial_invariance_auto(
  data,
  group,
  model,
  language = c("eng", "esp"),
  estimator = "auto",
  ordered = NULL,
  report = FALSE,
  criteria = "auto",
  missing = "auto",
  min_group_n = 20L,
  overwrite = FALSE,
  digits = 3,
  file_name = NULL,
  format = "excel",
  verbose = TRUE,
  fit_variant = c("scaled", "robust")
)
}
\arguments{
\item{data}{Marco de datos, matriz o \code{psychomatic_spec} que contiene los indicadores y la variable de agrupación (salvo que esta se proporcione como vector).}

\item{group}{Nombre o símbolo de la columna de agrupación, un índice positivo de columna o un vector con un valor por cada fila original.}

\item{model}{Una cadena con un modelo de AFC de lavaan.}

\item{language}{Idioma de salida: español (\code{"esp"}) o inglés (\code{"eng"}).}

\item{estimator}{\code{"auto"}, \code{"ML"}, \code{"MLR"}, \code{"ULS"}, \code{"WLSMV"} o \code{"ULSMV"}. WLSMV/ULSMV explícitos utilizan por defecto \code{ordered = TRUE} y eliminación por lista. La selección automática utiliza WLSMV para respuestas detectadas como ordinales; en los demás casos, ML cuando todos los grupos presentan evidencia favorable de Mardia y MLR en caso contrario. ML/MLR/ULS analizan las puntuaciones numéricas como continuas.}

\item{ordered}{\code{NULL}, \code{TRUE}/\code{FALSE} o nombres/índices de indicadores. Con WLSMV/ULSMV explícitos, \code{NULL} declara ordenados todos los indicadores del modelo. \code{FALSE} solicita explícitamente la estimación a partir de puntuaciones numéricas. Las respuestas de caracteres y de factores no ordenados requieren niveles declarados. Para una secuencia identificada de invariancia de umbrales y cargas, utilice \code{\link[=ordinal_invariance_auto]{ordinal_invariance_auto()}}.}

\item{report}{\code{FALSE}, \code{TRUE} para un informe en \code{getwd()} o una ruta explícita de Excel.}

\item{criteria}{\code{"auto"} o un vector numérico no negativo, con nombres para todos sus elementos, que incluya \code{cfi}, \code{rmsea}, \code{srmr_loading} y \code{srmr_other}; las tres entradas opcionales \verb{configural_*} permiten personalizar el criterio inicial de ajuste configural.}

\item{missing}{\code{"auto"}, \code{"fiml"} o \code{"listwise"}. Los análisis categóricos requieren eliminación por lista.}

\item{min_group_n}{Tamaño analítico mínimo operativo por grupo. El valor predeterminado de 20 no garantiza potencia, identificación ni estabilidad de la estimación.}

\item{overwrite}{Indica si puede reemplazarse una ruta explícita de informe ya existente.}

\item{digits}{Lugares decimales en la consola y los informes; los valores almacenados conservan toda su precisión.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca de tiempo en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{verbose}{Si es \code{TRUE}, imprime los resultados directamente, incluso cuando se asignan a un objeto.}

\item{fit_variant}{\code{"scaled"} (predeterminada) o \code{"robust"} para CFI y RMSEA, sus cambios, las decisiones y los informes. Los estimadores corregidos siempre utilizan chi-cuadrado, grados de libertad y valor p escalados; SRMR no se corrige. ML/ULS sin corrección utilizan por defecto índices de base. Los índices robustos solicitados que no estén disponibles permanecen como faltantes; nunca se sustituyen silenciosamente por índices escalados o de base.}
}
\value{
Un objeto de las clases \code{factorial_invariance_auto} y
\code{psychomatic_result}. Los valores originales no se redondean. Se conservan los alias directos heredados;
\verb{$warnings} es ahora un marco de datos estructurado y \verb{$warnings_list}
(también \verb{$notes}) es la representación de compatibilidad como vector de caracteres.
}
\description{
Ajusta modelos de AFC multigrupo configurales, métricos, escalares y estrictos mediante
llamadas nativas a CFA de lavaan y su identificación predeterminada mediante un indicador de referencia. WLSMV/ULSMV
utilizan por defecto respuestas ordenadas; solo se ajustan las cuatro etapas solicitadas. La
función automatiza la selección del estimador y la construcción de sintaxis, y mantiene
la identificación, las variantes de ajuste, la admisibilidad y la evidencia de comparación de modelos
disponibles para su inspección. Los puntos de corte del cambio en el ajuste son heurísticas descriptivas,
no una prueba de invariancia.
}
\examples{
# Reemplace MY_DATA por su marco de datos y adapte los nombres o las columnas de los ítems.
\dontrun{
model <- "F1 =~ R1 + R2 + R3 + R4"
factorial_invariance_auto(MY_DATA, group = "Country", model = model)
}

}
