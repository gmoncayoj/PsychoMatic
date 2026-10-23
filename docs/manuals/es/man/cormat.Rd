% Traducción editorial al español.
% Fuente: man/cormat.Rd
\name{cormat}
\alias{cormat}
\title{Matriz de correlaciones entre ítems}
\usage{
cormat(
  data,
  type = c("auto", "pearson", "poly", "tetra"),
  report = FALSE,
  file_name = NULL,
  digits = 3,
  ordinal_levels = NULL,
  overwrite = FALSE,
  format = "excel",
  verbose = TRUE,
  language = "eng"
)
}
\arguments{
\item{data}{Data frame o matriz que contiene los ítems.}

\item{type}{El valor predeterminado es \code{"auto"}: para ítems dicotómicos se utilizan correlaciones tetracóricas; para ítems ordinales con 2 a 5 categorías se utilizan correlaciones policóricas. Con 6 o más categorías, se utilizan correlaciones policóricas cuando la asimetría de momentos en valor absoluto supera 1, la curtosis excedente en valor absoluto supera 2 o la proporción de respuestas en los extremos es al menos 0.50; en caso contrario, se utiliza Pearson. Los ítems numéricos enteros con un máximo de 10 categorías observadas se clasifican como ordinales; los factores ordenados lo declaran explícitamente. Estos indicadores de forma son reglas operativas, no puntos de corte universales de la literatura. Las entradas que combinan variables con pocas categorías y variables continuas requieren una selección explícita. Las opciones manuales son \code{"pearson"}, \code{"poly"} y \code{"tetra"}.}

\item{report}{Si es \code{TRUE}, guarda un informe en el formato solicitado.}

\item{file_name}{Nombre opcional del archivo Excel. Si es \code{NULL}, se genera un nombre con marca temporal en el directorio de trabajo (\code{getwd()}).}

\item{digits}{Número de decimales utilizados en el informe Excel.}

\item{ordinal_levels}{Lista opcional con nombres únicos que especifica los niveles ordenados para los ítems categóricos no ordenados utilizados con \code{type = "poly"} o \code{"tetra"}.}

\item{overwrite}{Si es \code{TRUE}, se permite reemplazar un archivo de informe existente.}

\item{format}{Formato del informe: \code{"excel"}, \code{"word"} o ambos.}

\item{verbose}{Si es \code{TRUE}, muestra los resultados directamente, incluso cuando se asignan a un objeto.}

\item{language}{Idioma de salida: \code{"esp"} o \code{"eng"}.}
}
\value{
Una matriz de correlaciones con los tamaños muestrales por pares, los tipos de correlación solicitados y efectivos, y los diagnósticos de validación como atributos.
}
\description{
Calcula una matriz de correlaciones entre ítems mediante correlaciones de Pearson, policóricas o tetracóricas. Opcionalmente, exporta la matriz a un archivo Excel.
}
\examples{
# Sustituya MY_DATA por su data frame y adapte los nombres o las columnas de los ítems.
\dontrun{
cormat(MY_DATA[, 3:14])
}

}
