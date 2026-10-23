% Traducción editorial al español; fuente canónica: man/esem_auto.Rd.
\name{esem_auto}
\alias{esem_auto}
\title{Modelo de medida de ecuaciones estructurales exploratorio con rotación transparente}
\usage{
esem_auto(
  data,
  factors = NULL,
  n_factors = NULL,
  items = NULL,
  rotation = "target",
  ordered = NULL,
  estimator = NULL,
  missing = "listwise",
  rotation_args = list(),
  seed = NULL,
  compare_cfa = FALSE,
  fit_variant = c("scaled", "robust", "base"),
  language = "eng",
  digits = 3,
  report = FALSE,
  file_name = NULL,
  format = "excel",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Marco de datos, matriz o \code{psychomatic_spec}.}

\item{factors}{Lista con nombres de los conjuntos de ítems primarios. Cada ítem seleccionado aparece una sola vez. Target requiere esta lista. Con geomin, solo proporciona la selección, el número de factores y una referencia de AFC opcional, no un objetivo de rotación.}

\item{n_factors}{Entero >= 2, obligatorio para geomin sin \code{factors}.}

\item{items}{Nombres explícitos de columnas o índices positivos de columnas para datos originales con geomin sin \code{factors}; por ejemplo, \code{items = 2:13} o \code{items = names(MY_DATA)[2:13]}. Proporcione selectores de columnas, no un subconjunto de datos: \code{MY_DATA[, 2:13]} es un marco de datos y no constituye un argumento \code{items} válido. Este argumento es obligatorio para los datos originales con geomin, incluso si todas las columnas son ítems; en ese caso, utilice \code{items = names(MY_DATA)}. NULL con \code{factors} o una especificación.}

\item{rotation}{\code{"target"} (predeterminada) o \code{"geomin"}.}

\item{ordered}{NULL hereda o detecta heurísticamente los niveles de medida; TRUE trata todos los ítems seleccionados como ordinales, FALSE solicita un tratamiento numérico continuo, o un vector de caracteres identifica los ítems ordinales. Los factores nominales requieren una correspondencia ordinal explícita; los factores ordenados no pueden convertirse silenciosamente en puntuaciones continuas.}

\item{estimator}{NULL selecciona WLSMV con ítems ordinales y MLR en los demás casos. Opciones explícitas: WLSMV/ULSMV para ítems ordinales o mixtos, ML/MLR para ítems continuos.}

\item{missing}{\code{"listwise"} (predeterminado), \code{"fiml"} para datos continuos con ML/MLR, o \code{"pairwise"} para datos ordinales con WLSMV/ULSMV. No se realiza imputación automática.}

\item{rotation_args}{Opciones con nombres: \code{rstarts} (30), \code{geomin_epsilon} (0.001, solo geomin), \code{gpa_tol} (1e-5), \code{max_iter} (10000).}

\item{seed}{NULL (predeterminado) deja la gestión del generador de números aleatorios (RNG) a lavaan, sin imponer una semilla desde la función envolvente. Un entero no negativo explícito hace reproducibles los inicios de rotación y restaura el estado previo del RNG después del ajuste.}

\item{compare_cfa}{Indica si se ajusta el AFC de estructura simple especificado por \code{factors} con los mismos datos y opciones. La comparación es descriptiva, sin selección automática de modelos, prueba de razón de verosimilitudes (LRT) ni alineación factorial.}

\item{fit_variant}{\code{"scaled"}, \code{"robust"} o \code{"base"}. Los índices solicitados que no están disponibles permanecen en NA con una explicación. El chi-cuadrado, los grados de libertad y el valor p corregidos utilizan variantes escaladas, salvo que se seleccione la variante de base; SRMR siempre utiliza la de base.}

\item{language}{Idioma de salida: \code{"eng"} o \code{"esp"}.}

\item{digits}{Número de decimales que se mostrarán.}

\item{report}{Indica si se escribe un informe, o especifica una ruta explícita para el informe.}

\item{file_name}{Ruta o base del informe; NULL utiliza un nombre único en \code{getwd()}.}

\item{format}{\code{"excel"}, \code{"word"} o ambos.}

\item{overwrite}{Indica si pueden reemplazarse los archivos de informe existentes.}

\item{verbose}{Indica si se imprimen los resultados y el borrador del manuscrito.}
}
\value{
Un objeto \code{esem_auto}/\code{psychomatic_result}. \code{fit} contiene el ajuste nativo de
lavaan con alias seguros para los ítems; \code{results$item_map} establece la correspondencia con los nombres originales.
\code{results} contiene cargas e IC del 95\%, diagnósticos de ítems, correlaciones factoriales,
variantes de ajuste, residuos de covarianza sin transformar, matriz objetivo, correspondencias y sintaxis.
Las tablas de cargas distinguen los intervalos nativos del 95\% en la escala original (\code{raw_ci_lower}, \code{raw_ci_upper})
de los estandarizados (\code{ci_lower}, \code{ci_upper}).
\code{effective_settings$preparation} registra el orden de las categorías, la puntuación inversa declarada
y los límites efectivos, incluidos los inferidos. La codificación del informe corresponde
al orden de entrada de las categorías antes de la puntuación inversa, que se registra por separado.
\code{manuscript_text} y \code{results$manuscript} contienen un borrador basado en los datos que aparece
al final, bajo Interpretación, en la salida de consola, Word y Excel destinada al lector.
La preparación, las correspondencias de ítems y factores y las versiones se conservan en el objeto y en las exportaciones,
pero se omiten en la consola. El texto no constituye una validación independiente
ni una interpretación sustantiva de factores sin etiquetar. Los ajustes fallidos o
inadmisibles generan limitaciones en lugar de conclusiones favorables.
}
\description{
Ajusta un bloque de medida exploratorio directamente en \code{lavaan}, conservando todas
las cargas cruzadas dentro del bloque y la incertidumbre de los parámetros rotados. Sin
regresiones estructurales, esto es AFE en la terminología de \code{lavaan}; ESEM también se
utiliza para estos modelos de medida en la literatura psicométrica.
}
\section{Elección de una rotación}{

\strong{Target} resulta útil cuando la teoría o la investigación previa identifica la dimensión
primaria de cada ítem, permitiendo a la vez cargas secundarias. \code{factors} define
el patrón objetivo. Las cargas secundarias se orientan hacia cero, pero no se fijan
en cero. Una entrada de 1 en la matriz objetivo deja la carga primaria sin objetivo; no
fija esa carga en 1. Un objetivo derivado de los datos de la misma muestra sigue siendo
exploratorio y requiere validación independiente para sustentar una mayor generalización.

\strong{Geomin} resulta útil cuando la organización de los ítems por factores es incierta. Busca
un patrón relativamente simple sin prescribir pertenencias sustantivas.
\code{n_factors} sigue siendo una decisión del investigador. Los factores conservan identificadores
neutros porque su orden y sus signos no establecen su interpretación sustantiva.
El épsilon inicial es 0.001, no un óptimo universal; especifique alternativas
mediante \code{rotation_args} e informe el valor elegido. Ambas opciones son oblicuas.

Las rotaciones del mismo bloque sin restricciones conservan el ajuste global del modelo; CFI y
RMSEA no seleccionan la rotación preferible. Considere la coherencia teórica,
los patrones de cargas y la estabilidad. Treinta inicios aleatorios no establecen un óptimo
global ni la estabilidad muestral. Las cargas estandarizadas informadas son coeficientes
de patrón, no la matriz de estructura oblicua.
}

\examples{
\dontrun{
# Reemplace MY_DATA por su marco de datos.
# Ítems ordinales en las columnas 2:13; la columna 1 no es un ítem.
result <- esem_auto(data = MY_DATA,
                    n_factors = 3,
                    items = 2:13,
                    rotation = "geomin",
                    ordered = TRUE,
                    language = "esp")
result$manuscript_text

# Indicadores continuos: aquí MY_DATA debe contener solo columnas de ítems.
result <- esem_auto(data = MY_DATA,
                    n_factors = 3,
                    items = names(MY_DATA),
                    rotation = "geomin",
                    ordered = FALSE,
                    language = "esp")

# Target: reemplace los nombres de los ítems y sus pertenencias primarias según su teoría.
factors <- list(Factor1 = c("R1", "R2", "R3", "R4"),
                Factor2 = c("R5", "R6", "R7", "R8"),
                Factor3 = c("R9", "R10", "R11", "R12"))
result <- esem_auto(data = MY_DATA,
                    factors = factors,
                    rotation = "target",
                    ordered = TRUE,
                    language = "esp")
}
}
\references{
Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation modeling.
Structural Equation Modeling, 16(3), 397–438. doi:10.1080/10705510903008204.

Marsh, H., & Alamer, A. (2024). When and how to use set-exploratory structural
equation modelling to test structural models: A tutorial using the R package
lavaan. British Journal of Mathematical and Statistical Psychology, 77(3),
459–476. doi:10.1111/bmsp.12336.
}
