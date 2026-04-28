#' PsychoMatic: herramientas psicometricas automatizadas
#'
#' `PsychoMatic` reune funciones para descriptivos, analisis factorial
#' exploratorio, analisis factorial confirmatorio e invarianza de medicion en
#' flujos de trabajo pensados para R y RStudio.
#'
#' Funciones principales: `desc_auto()`, `cormat()`, `efa_auto()`,
#' `export_efa()`, `cfa_auto()`, `export_cfa()`, `inv_align_auto()` y
#' `factorial_invariance_auto()`.
#'
#' @keywords internal
#' @import flextable
#' @import officer
#' @import openxlsx
#' @import stats
#' @import utils
#' @importFrom GPArotation oblimin
#' @importFrom lavaan fitMeasures lavaanify lavInspect modificationIndices standardizedSolution
#' @importFrom Matrix nearPD
#' @importFrom psych KMO alpha cortest.bartlett fa fa.parallel omega polychoric tetrachoric
"_PACKAGE"
