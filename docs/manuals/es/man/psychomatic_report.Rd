% Traducción editorial al español.
% Fuente: man/psychomatic_report.Rd
\name{psychomatic_report}
\alias{psychomatic_report}
\title{Creación de un informe integrado de PsychoMatic}
\usage{
psychomatic_report(
  ...,
  format = "excel",
  file_name = NULL,
  language = c("eng", "esp"),
  digits = 3,
  overwrite = FALSE,
  detail = c("report", "audit")
)
}
\arguments{
\item{...}{Objetos de resultados de PsychoMatic con nombres o data frames. Los objetos compatibles incluyen las salidas de \code{desc_auto()}, \code{screen_items()}, \code{score_scale()}, \code{reliability_auto()}, \code{efa_auto()}, \code{cfa_auto()}, \code{esem_auto()}, \code{factorial_invariance_auto()}, \code{ordinal_invariance_auto()}, \code{inv_align_auto()}, \code{irt_auto()}, \code{dif_auto()} y \code{compare_models_auto()}.}

\item{format}{Formato de salida: \code{"excel"}, \code{"word"} o ambos.}

\item{file_name}{Ruta de salida o ruta base. Si es \code{NULL}, se genera un nombre único con marca temporal en \code{getwd()}. Las rutas explícitas prevalecen.}

\item{language}{Idioma del informe: \code{"eng"} o \code{"esp"}.}

\item{digits}{Número de decimales para las tablas numéricas.}

\item{overwrite}{Si es \code{TRUE}, se permite reemplazar los archivos de informe existentes.}

\item{detail}{\code{"report"} (predeterminado) organiza por análisis los resultados dirigidos al lector. \code{"audit"} conserva las tablas técnicas detalladas y los metadatos de ejecución.}
}
\value{
Ruta de salida, de forma invisible. Cuando se solicitan ambos formatos, devuelve un vector de caracteres con nombres. La salida Excel incluye un atributo \code{sheet_map} que registra los nombres solicitados y efectivos de las hojas de cálculo.
}
\description{
Recopila tablas de objetos de resultados de PsychoMatic y escribe un libro Excel y/o un documento Word. Los informes utilizan la ruta proporcionada por el usuario o el directorio de trabajo actual (\code{getwd()}) cuando \code{file_name = NULL}.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
result <- reliability_auto(MY_DATA[, 3:6])
psychomatic_report(reliability = result, format = "word", file_name = "Report.docx")
}

}
