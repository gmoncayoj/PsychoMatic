% Traducción editorial; fuente canónica: man/plot.psychomatic_network.Rd
\name{plot.psychomatic_network}
\alias{plot.psychomatic_network}
\title{Graficar resultados de redes psicométricas ya calculados}
\usage{
\method{plot}{psychomatic_network}(x, type = "network", language = NULL, label_cex = 1.2, ...)
}
\arguments{
\item{x}{Un resultado \code{psychomatic_network}.}

\item{type}{Tipo de figura. Consulte las figuras disponibles en \code{network_plots()}.}

\item{language}{Idioma de la figura; por defecto utiliza el idioma del resultado.}

\item{label_cex}{Multiplicador del tamaño de las etiquetas de los nodos; por defecto, 1.2.
Las etiquetas no se reducen automáticamente para ajustarse a los nodos; utilice identificadores breves.}

\item{...}{Reservado para futuros argumentos gráficos.}
}
\value{
De forma invisible, un objeto qgraph para \code{"network"} y
\code{"loadings"}, o los resultados subyacentes para las demás figuras. No se vuelve a calcular el análisis.
}
\description{
Graficar resultados de redes psicométricas ya calculados
}
