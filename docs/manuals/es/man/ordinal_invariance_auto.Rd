% Traducción editorial al español; fuente canónica: man/ordinal_invariance_auto.Rd.
\name{ordinal_invariance_auto}
\alias{ordinal_invariance_auto}
\title{Invariancia de medida ordinal automatizada}
\usage{
ordinal_invariance_auto(
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

\item{estimator}{\code{"auto"} (WLSMV), \code{"WLSMV"} o \code{"ULSMV"}. ML, MLR y otros estimadores se rechazan antes del ajuste.}

\item{ordered}{\code{NULL} detecta los indicadores ordinales; \code{TRUE} declara ordenados todos los indicadores del modelo; los nombres/índices declaran una selección explícita. Todos los indicadores del modelo deben ser ordenados. Se rechaza \code{FALSE}.}

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
Un resultado con las clases \code{ordinal_invariance_auto},
\code{factorial_invariance_auto} y \code{psychomatic_result}, que incluye los modelos ajustados originales,
los diagnósticos de identificación, los cambios en el ajuste e informes bilingües.
}
\description{
Ajusta modelos configurales, con igualdad de umbrales y con igualdad de umbrales y cargas,
utilizando la parametrización theta y la identificación de Wu y Estabrook (2016) mediante
\code{semTools::measEq.syntax()}, como ilustran Svetina et al. (2020).
Todos los indicadores deben ser ordinales. Los indicadores binarios o de tres categorías pueden
impedir que una comparación separada de umbrales sea contrastable; dichas etapas se marcan explícitamente
como no identificadas de manera independiente, en lugar de aceptarse como invariantes.
}
\details{
Los cambios son diferencias con signo (modelo actual menos el modelo evaluable anterior).
Las directrices de Chen (2007) son extrapolaciones descriptivas a la estimación ordinal,
no puntos de corte validados específicamente para datos ordinales: una disminución del CFI superior a .010,
junto con un aumento del RMSEA superior a .015 o un aumento del SRMR superior a .030
(al añadir cargas) / .010 (umbrales), señala falta de apoyo. Primero se comprueban el ajuste
configural y la admisibilidad. Pueden proporcionarse criterios personalizados.
}
\examples{
\dontrun{
ordinal_invariance_auto(MY_DATA, group = "Sexo",
  model = "F =~ R1 + R2 + R3 + R4", ordered = TRUE,
  estimator = "WLSMV", report = TRUE, format = c("word", "excel"))
}
}
\references{
Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor analysis
models of different levels of invariance for ordered categorical outcomes.
\emph{Psychometrika, 81}(4), 1014-1045. \doi{10.1007/s11336-016-9506-0}

Svetina, D., Rutkowski, L., & Rutkowski, D. (2020). Multiple-group invariance
with categorical outcomes using updated guidelines: An illustration using
Mplus and the lavaan/semTools packages. \emph{Structural Equation Modeling:
A Multidisciplinary Journal, 27}(1), 111-130.
\doi{10.1080/10705511.2019.1602776}

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
measurement invariance. \emph{Structural Equation Modeling: A Multidisciplinary
Journal, 14}(3), 464-504. \doi{10.1080/10705510701301834}
}
