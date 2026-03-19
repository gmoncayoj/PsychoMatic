#' PsychoMatic: herramientas psicometricas automatizadas
#'
#' `PsychoMatic` reune funciones para descriptivos, analisis factorial
#' exploratorio, analisis factorial confirmatorio e invarianza de medicion en
#' flujos de trabajo pensados para R y RStudio.
#'
#' Funciones principales: `desc_auto()`, `efa_auto()`, `exportar_efa()`,
#' `cfa_auto()`, `inv_align_auto()` y `factorial_invariance_auto()`.
#'
#' @keywords internal
#' @import flextable
#' @import officer
#' @import openxlsx
#' @import stats
#' @import utils
#' @importFrom lavaan fitMeasures lavaanify lavInspect modificationIndices standardizedSolution
#' @importFrom psych KMO alpha cortest.bartlett fa fa.parallel omega polychoric
"_PACKAGE"
