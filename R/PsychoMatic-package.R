#' PsychoMatic: automated psychometric workflows
#'
#' `PsychoMatic` provides low-syntax, reproducible workflows for declaring and
#' checking psychometric data, scoring, reliability, exploratory and
#' confirmatory factor analysis, exploratory structural equation measurement
#' models, psychometric networks, complete and researcher-directed partial
#' measurement invariance, alignment, IRT/DIF screening, auditing, and reports.
#' Automation reduces procedural syntax but does not replace measurement theory,
#' design decisions, or substantive review.
#'
#' Use [psychomatic_spec()] when item roles or
#' encodings should be explicit. Main analysis functions include
#' [desc_auto()], [score_scale()], [reliability_auto()], [cormat()],
#' [efa_auto()], [cfa_auto()], [esem_auto()], [network_auto()], [network_plots()],
#' [factorial_invariance_auto()], [ordinal_invariance_auto()],
#' [invariance_partial()], [inv_align_auto()],
#' [irt_auto()], [dif_auto()], and [compare_models_auto()]. Use
#' [psychomatic_report()] for explicit export.
#'
#' Most 1.0.0 analysis objects inherit from `psychomatic_result` and separate
#' requested/effective settings, raw results, structured conditions,
#' diagnostics, decisions, and conclusions. Numerical analysis fields retain
#' full precision; display and report formatting may honor `digits`.
#'
#' Ordinal encodings, missing-data policies, estimator/identification settings,
#' and `not_evaluable` states should always be inspected. Fit cutoffs, parallel
#' analysis, Chen fit-change thresholds, DIF flags, and alignment summaries are
#' decision aids rather than formal proof of validity or invariance.
#'
#' @keywords internal
#' @import stats
#' @import utils
#' @importFrom GPArotation oblimin
#' @importFrom lavaan fitMeasures lavaanify lavInspect modificationIndices standardizedSolution
#' @importFrom Matrix nearPD
#' @importFrom psych KMO alpha cortest.bartlett fa fa.parallel omega polychoric tetrachoric
"_PACKAGE"
