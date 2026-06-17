#' PsychoMatic: automated psychometric workflows
#'
#' `PsychoMatic` provides functions for item descriptives, exploratory factor
#' analysis, confirmatory factor analysis, scoring, measurement invariance,
#' alignment invariance, IRT/DIF screening, and psychometric reporting in
#' reproducible R workflows.
#'
#' Main functions include `desc_auto()`, `screen_items()`, `reverse_items()`,
#' `score_scale()`, `cormat()`, `efa_auto()`, `cfa_auto()`,
#' `factorial_invariance_auto()`, `inv_align_auto()`, `irt_auto()`,
#' `dif_auto()`, `compare_models_auto()`, and `psychomatic_report()`.
#'
#' @keywords internal
#' @import stats
#' @import utils
#' @importFrom GPArotation oblimin
#' @importFrom lavaan fitMeasures lavaanify lavInspect modificationIndices standardizedSolution
#' @importFrom Matrix nearPD
#' @importFrom psych KMO alpha cortest.bartlett fa fa.parallel omega polychoric tetrachoric
"_PACKAGE"
