% Traducción editorial al español.
% Fuente: man/network_auto.Rd
\name{network_auto}
\alias{network_auto}
\title{Análisis exploratorio automatizado de redes psicométricas}
\usage{
network_auto(
  data,
  items = NULL,
  group = NULL,
  reverse = NULL,
  ranges = NULL,
  uva = NULL,
  uva_cutoff = 0.25,
  corr = c("spearman", "pearson", "polychoric"),
  gamma = 0.5,
  missing = c("listwise", "pairwise"),
  wording = TRUE,
  hierarchy = TRUE,
  bootstrap = TRUE,
  invariance = TRUE,
  n_boot = 1000L,
  n_perm = 5000L,
  consensus_iter = 1000L,
  n_cores = 1L,
  seed = NULL,
  language = "eng",
  report = FALSE,
  format = "excel",
  file_name = NULL,
  overwrite = FALSE,
  digits = 3L,
  verbose = TRUE,
  model = c("glasso", "TMFG", "BGGM"),
  algorithm = c("louvain", "leiden", "walktrap"),
  bootstrap_type = c("resampling", "parametric")
)
}
\arguments{
\item{data}{Datos de ítems numéricos u ordenados, o un objeto \code{psychomatic_spec}.}

\item{items}{Columnas de ítems; de forma predeterminada, se utilizan las columnas admisibles, excluyendo \code{group}.}

\item{group}{Columna de agrupación (p. ej., sexo), o la declarada en la especificación.}

\item{reverse}{Columnas de ítems conocidos con clave inversa. Nunca se infieren a partir de los signos.}

\item{ranges}{Rangos teóricos de respuesta de los ítems inversos.}

\item{uva}{NULL (predeterminado) omite por completo UVA. \code{"diagnostic"} informa los pares redundantes, wTO y propuestas justificadas de eliminación sin eliminar ítems. \code{"removed"} autoriza la selección nativa de representantes y la eliminación, con registro de criterios, puntuaciones y razones. Ningún modo solicita entradas interactivas.}

\item{uva_cutoff}{Umbral de solapamiento topológico ponderado; el valor predeterminado es .25.}

\item{corr}{\code{"spearman"} (predeterminado), \code{"pearson"} o \code{"polychoric"} para ítems ordinales (incluidos los binarios), con un máximo de 11 categorías por ítem en EGAnet. Las correlaciones se vuelven a estimar en las remuestras.}

\item{gamma}{Penalización EBIC; el valor predeterminado es .5.}

\item{missing}{\code{"listwise"} (predeterminado) o \code{"pairwise"} para las correlaciones.}

\item{wording}{Ejecuta riEGA sobre las respuestas originales, sin recodificar.}

\item{hierarchy}{Ejecuta hierEGA y evalúa modelos candidatos de AFC identificables.}

\item{bootstrap}{Ejecuta bootEGA y los análisis de precisión y de eliminación de casos de bootnet.}

\item{invariance}{Ejecuta NCT y pruebas de permutación de las cargas EGA si existe un grupo.}

\item{n_boot}{Replicaciones por análisis bootstrap; el valor predeterminado es 1000.}

\item{n_perm}{Permutaciones por comparación; el valor predeterminado es 5000.}

\item{consensus_iter}{Iteraciones de consenso de Louvain; el valor predeterminado es 1000.}

\item{n_cores}{Procesos de trabajo para bootEGA, bootnet y las permutaciones de cargas; el valor predeterminado es 1. NCT permanece secuencial en la versión compatible.}

\item{seed}{Semilla opcional de valor entero. NULL (predeterminado) obtiene una semilla de análisis del generador de números aleatorios de la sesión y la registra en \code{effective_settings$seed}. Una semilla explícita permite reproducir el análisis y restaura el estado del generador de números aleatorios del entorno que realiza la llamada.}

\item{language}{\code{"eng"} (predeterminado) o \code{"esp"}.}

\item{report}{Exporta un informe; valor lógico o ruta de salida.}

\item{format}{Formato o formatos del informe: \code{"excel"} y/o \code{"word"}.}

\item{file_name}{Nombre base o ruta de salida del informe.}

\item{overwrite}{Permite sobrescribir archivos de informe existentes; el valor predeterminado es FALSE.}

\item{digits}{Precisión de visualización; los resultados numéricos conservan toda su precisión.}

\item{verbose}{Muestra el progreso de las etapas, el tiempo transcurrido y el resumen del análisis.}

\item{model}{Estimador de red: \code{"glasso"} (predeterminado), \code{"TMFG"} o \code{"BGGM"}. BGGM requiere una versión de EGAnet que lo admita; EGAnet 2.4.1 lo deshabilita. Ningún estimador no disponible se sustituye de forma silenciosa. Gamma solo se aplica a GLASSO.}

\item{algorithm}{Detección de comunidades: \code{"louvain"} (predeterminado), \code{"leiden"} o \code{"walktrap"}. Se utiliza en la EGA empírica y por grupos, y en ambos órdenes de hierEGA.}

\item{bootstrap_type}{Bootstrap de bootEGA: \code{"resampling"} (predeterminado) o \code{"parametric"}. También se aplica a los bootstrap configural, de riEGA y de hierEGA. El paramétrico genera datos continuos con distribución normal multivariada; no conserva las categorías de respuesta ordinal. No se selecciona automáticamente una opción por su superioridad. Los procedimientos de precisión de aristas y de eliminación de casos de bootnet se mantienen separados.}
}
\value{
Un objeto \code{psychomatic_network} y \code{psychomatic_result}. Las etapas opcionales que fallan se mantienen explícitas en \code{results$stages} y \code{errors}; los resultados parciales nunca se describen como un análisis completo. Las respuestas originales permanecen disponibles en \code{diagnostics$input}, pero se excluyen de los informes dirigidos al lector.
}
\description{
Opcionalmente ejecuta el análisis de variables únicas antes de EGA (con Spearman, EBICglasso y Louvain de consenso de forma predeterminada), cargas de red revisadas, TEFI, diagnósticos de redacción y jerarquía, bootstrap estructural y comparaciones entre grupos. Los paquetes opcionales solo son necesarios para las etapas solicitadas. Consulte la viñeta del flujo de análisis de redes para conocer los supuestos.
}
\details{
La salida predeterminada de consola y los informes destacan las comunidades, la consistencia estructural y las comparaciones entre grupos. La ausencia de significación no establece equivalencia. Las tablas diagnósticas completas permanecen en el objeto y pueden mostrarse con \code{print(fit, detail="audit")} o exportarse con \code{psychomatic_report(fit, detail="audit")}. Las cursivas de las referencias APA se comparten entre los análisis en Word, Excel y las consolas interactivas compatibles; los registros de texto plano conservan el texto sin formato. \code{options(PsychoMatic.ansi=FALSE)} desactiva el estilo de consola; \code{TRUE} lo fuerza en terminales compatibles con cursivas ANSI.
}
\examples{
\dontrun{
fit <- network_auto(MY_DATA, items = 1:12, group = "sex")
plot(fit, type = "network")
psychomatic_report(fit, format = c("word", "excel"), file_name = "network")
}
}
