% Traducción editorial al español; fuente canónica: man/cfa_auto.Rd.
\name{cfa_auto}
\alias{cfa_auto}
\title{AFC automatizado con informes psicométricos}
\usage{
cfa_auto(
  data,
  model,
  ordered = NULL,
  estimator = NULL,
  std.lv = FALSE,
  mi_threshold = 10,
  n_mi = 10,
  alpha_norm = 0.05,
  language = "eng",
  bifactor_orthogonal = TRUE,
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
\item{data}{Marco de datos que contiene las variables observadas.}

\item{model}{Sintaxis del modelo escrita para \code{lavaan}.}

\item{ordered}{Vector opcional que indica qué ítems deben tratarse como ordinales.}

\item{estimator}{Estimador opcional. Si es \code{NULL}, la función selecciona uno automáticamente.}

\item{std.lv}{Indica si deben estandarizarse las variables latentes.}

\item{mi_threshold}{Umbral mínimo para informar índices de modificación.}

\item{n_mi}{Número máximo de índices de modificación que se mostrarán.}

\item{alpha_norm}{Nivel alfa utilizado para evaluar la normalidad multivariada.}

\item{language}{Idioma del informe: \code{"esp"} o \code{"eng"}.}

\item{bifactor_orthogonal}{Escalar lógico utilizado únicamente en los modelos bifactor detectados. El valor predeterminado conserva el ajuste bifactor ortogonal histórico; establezca \code{FALSE} únicamente cuando se pretendan correlaciones entre los factores latentes por razones sustantivas.}

\item{digits}{Lugares decimales en la consola y los informes; los valores almacenados conservan toda su precisión.}

\item{report}{Si es \code{TRUE}, también guarda un informe destinado al lector.}

\item{file_name}{Ruta base opcional para los archivos exportados. Cuando es \code{NULL}, los informes utilizan un nombre único con marca de tiempo en el directorio de trabajo (\code{getwd()}).}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Si es \code{TRUE}, permite reemplazar un informe existente.}

\item{verbose}{Si es \code{TRUE}, imprime los resultados directamente, incluso cuando se asignan a un objeto.}

\item{fit_variant}{Variante preferida de los índices de ajuste: \code{"scaled"} (predeterminada), \code{"robust"} o \code{"base"}. Cuando los índices solicitados no están disponibles, se recurre primero a los escalados y después a los de base. Las estimaciones de base, escaladas y robustas siguen disponibles en el resultado.}
}
\value{
Un objeto de las clases \code{cfa_auto} y \code{psychomatic_result}. El objeto
\code{lavaan} ajustado está disponible en \code{fit}; las tablas ordenadas y reutilizables están disponibles en
\code{tables}, \code{fit_indices}, \code{loadings_table} y \code{reliability_table}.
\code{warnings} es un marco de datos estructurado con metadatos de las condiciones;
\code{warnings_list} es la representación transitoria como vector de caracteres que se conserva
para el código escrito para PsychoMatic 0.3.x.
}
\description{
Realiza un análisis factorial confirmatorio con selección heurística del estimador,
diagnósticos de normalidad multivariada, índices de ajuste, estimaciones de fiabilidad y
sugerencias de reespecificación del modelo.
\code{lavaan} realiza el ajuste del modelo. Cuando está instalado, \code{MVN} proporciona el
diagnóstico de Mardia sin transformar utilizado para la selección automática del estimador y \code{semTools}
proporciona la fiabilidad compuesta y la varianza media extraída (AVE). Su ausencia se registra
explícitamente y no impide ajustar el AFC; cuando no puede evaluarse la normalidad,
la selección automática utiliza conservadoramente MLR para los datos continuos.
Los índices de ajuste conservan sus valores de base, escalados y robustos, y la variante efectiva
se identifica por separado para cada índice.
}
\section{Compatibilidad}{

Desde PsychoMatic 1.0.0, utilice \code{result$warnings} para inspeccionar programáticamente
las fuentes y los mensajes de las advertencias. El código existente que espera un vector
de caracteres puede migrar a \code{result$warnings_list}; este alias de compatibilidad
se conservará durante un ciclo de obsolescencia. Por ejemplo,
\code{warning_messages <- result$warnings_list} reemplaza el código que antes
trataba \code{result$warnings} como datos de caracteres.
}

\examples{
# Reemplace MY_DATA por su marco de datos y adapte los nombres o las columnas de los ítems.
\dontrun{
model <- "F1 =~ R1 + R2 + R3 + R4"
cfa_auto(MY_DATA, model = model)
}

}
