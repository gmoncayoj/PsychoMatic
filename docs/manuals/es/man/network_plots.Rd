% Traducción editorial al español.
% Fuente: man/network_plots.Rd
\name{network_plots}
\alias{network_plots}
\title{Exportación de figuras para publicación a partir de resultados de redes almacenados}
\usage{
network_plots(
  x,
  directory = NULL,
  types = NULL,
  formats = c("pdf", "png"),
  language = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  overwrite = FALSE,
  label_cex = 1.2
)
}
\arguments{
\item{x}{Un objeto \code{psychomatic_network}.}

\item{directory}{Directorio de salida existente (se crea si no existe).}

\item{types}{Tipos de figuras; NULL selecciona las figuras pertinentes para publicación y omite los diagnósticos de redacción no disponibles y las figuras descriptivas de TEFI y centralidad. Proporcione explícitamente \code{types=network_plots(x)} para incluir todos los tipos de diagnóstico. Cuando hay resultados bootstrap, \code{"ega_stability"} combina la red EGA y el gráfico de replicación de ítems. \code{"network"} e \code{"item_stability"} siguen disponibles de forma explícita como figuras independientes. Las redes de dos grupos señalan las pruebas completas de cargas con BH.}

\item{formats}{Formatos de archivo: uno o varios de \code{"pdf"}, \code{"svg"} y \code{"png"}.}

\item{language}{\code{"eng"} o \code{"esp"}.}

\item{width, height}{Dimensiones de la figura en pulgadas.}

\item{dpi}{Resolución PNG; el valor predeterminado es 300.}

\item{overwrite}{Permite sobrescribir archivos existentes.}

\item{label_cex}{Multiplicador del tamaño de las etiquetas de los nodos; el valor predeterminado es 1.2. Auméntelo para figuras que se reducirán a una anchura pequeña en la publicación.}
}
\value{
De forma invisible, las rutas de las figuras exportadas; sin \code{directory}, los tipos disponibles.
}
\description{
Exporta figuras para publicación a partir de resultados de redes almacenados.
}
