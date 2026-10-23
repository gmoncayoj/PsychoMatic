% Traducción editorial; fuente canónica: man/PsychoMatic-package.Rd
\docType{package}
\name{PsychoMatic-package}
\alias{PsychoMatic}
\alias{PsychoMatic-package}
\title{PsychoMatic: flujos de trabajo psicométricos automatizados}
\description{
\code{PsychoMatic} ofrece flujos reproducibles con sintaxis reducida para
declarar y revisar datos psicométricos, calcular puntuaciones y fiabilidad,
realizar análisis factorial exploratorio y confirmatorio, ajustar modelos
exploratorios de medición con ecuaciones estructurales y redes psicométricas,
evaluar invarianza completa y parcial dirigida por el investigador, aplicar
alineación, explorar TRI/DIF, auditar y elaborar informes.
La automatización reduce la sintaxis de los procedimientos, pero no sustituye
la teoría de la medición, las decisiones del diseño ni la revisión sustantiva.
}
\details{
Utilice \code{\link[=psychomatic_spec]{psychomatic_spec()}} cuando sea necesario
declarar explícitamente las funciones de los ítems o sus codificaciones.
Entre las funciones principales están
\code{\link[=desc_auto]{desc_auto()}}, \code{\link[=score_scale]{score_scale()}},
\code{\link[=reliability_auto]{reliability_auto()}}, \code{\link[=cormat]{cormat()}},
\code{\link[=efa_auto]{efa_auto()}}, \code{\link[=cfa_auto]{cfa_auto()}},
\code{\link[=esem_auto]{esem_auto()}}, \code{\link[=network_auto]{network_auto()}},
\code{\link[=network_plots]{network_plots()}},
\code{\link[=factorial_invariance_auto]{factorial_invariance_auto()}},
\code{\link[=ordinal_invariance_auto]{ordinal_invariance_auto()}},
\code{\link[=invariance_partial]{invariance_partial()}},
\code{\link[=inv_align_auto]{inv_align_auto()}}, \code{\link[=irt_auto]{irt_auto()}},
\code{\link[=dif_auto]{dif_auto()}} y
\code{\link[=compare_models_auto]{compare_models_auto()}}. Para exportar de forma
explícita, utilice \code{\link[=psychomatic_report]{psychomatic_report()}}.

La mayoría de los objetos de análisis de 1.0.0 heredan de
\code{psychomatic_result} y separan configuración solicitada y efectiva,
resultados brutos, condiciones estructuradas, diagnósticos, decisiones y
conclusiones. Los campos numéricos conservan toda la precisión; la presentación
y el formato de los informes pueden respetar \code{digits}.

Revise siempre las codificaciones ordinales, el tratamiento de datos faltantes,
las opciones de estimación e identificación y los estados \code{not_evaluable}.
Los puntos de corte del ajuste, el análisis paralelo, los umbrales de cambio
de ajuste de Chen, las alertas DIF y los resúmenes de alineación son ayudas
para decidir; no constituyen pruebas formales de validez o invarianza.
}
\seealso{
Enlaces útiles:
\itemize{
  \item \url{https://github.com/gmoncayoj/PsychoMatic}
  \item Comunique errores en \url{https://github.com/gmoncayoj/PsychoMatic/issues}
}
}
\author{
\strong{Mantenedor}: Jose Gamarra-Moncayo \email{gamarramoncayoj@gmail.com} (\href{https://orcid.org/0000-0002-0781-3616}{ORCID})

Autores:
\itemize{
  \item Jose Gamarra-Moncayo \email{gamarramoncayoj@gmail.com} (\href{https://orcid.org/0000-0002-0781-3616}{ORCID})
}
}
\keyword{internal}
