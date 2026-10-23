#' Automated ordinal measurement invariance
#'
#' Fits configural, equal-threshold, and equal-threshold-and-loading models using
#' theta parameterization and Wu and Estabrook (2016) identification via
#' `semTools::measEq.syntax()`, as illustrated by Svetina et al. (2020).
#' All indicators must be ordinal. Binary or three-category indicators can make
#' a separate threshold comparison untestable; such stages are explicitly marked
#' as not independently identified rather than accepted as invariant.
#'
#' Changes are signed differences (current minus preceding evaluable model).
#' Chen (2007) guidelines are descriptive extrapolations to ordinal estimation,
#' not ordinal-specific validated cutoffs: a CFI decline exceeding .010 together
#' with an RMSEA increase exceeding .015 or SRMR increase exceeding .030
#' (adding loadings) / .010 (thresholds) signals lack of support. Configural fit
#' and admissibility are checked first. Custom criteria can be supplied.
#'
#' @inheritParams factorial_invariance_auto
#' @param estimator `"auto"` (WLSMV), `"WLSMV"`, or `"ULSMV"`. ML, MLR,
#'   and other estimators are rejected before fitting.
#' @param ordered `NULL` detects ordinal indicators; `TRUE` declares all model
#'   indicators ordered; names/indices declare an explicit selection. Every model
#'   indicator must be ordered. `FALSE` is rejected.
#' @return A result with classes `ordinal_invariance_auto`,
#'   `factorial_invariance_auto`, and `psychomatic_result`, including raw fitted
#'   models, identification diagnostics, fit changes, and bilingual reports.
#' @references
#' Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor analysis
#' models of different levels of invariance for ordered categorical outcomes.
#' *Psychometrika, 81*(4), 1014-1045. \doi{10.1007/s11336-016-9506-0}
#'
#' Svetina, D., Rutkowski, L., & Rutkowski, D. (2020). Multiple-group invariance
#' with categorical outcomes using updated guidelines: An illustration using
#' Mplus and the lavaan/semTools packages. *Structural Equation Modeling:
#' A Multidisciplinary Journal, 27*(1), 111-130.
#' \doi{10.1080/10705511.2019.1602776}
#'
#' Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
#' measurement invariance. *Structural Equation Modeling: A Multidisciplinary
#' Journal, 14*(3), 464-504. \doi{10.1080/10705510701301834}
#' @examples
#' \dontrun{
#' ordinal_invariance_auto(MY_DATA, group = "Sexo",
#'   model = "F =~ R1 + R2 + R3 + R4", ordered = TRUE,
#'   estimator = "WLSMV", report = TRUE, format = c("word", "excel"))
#' }
#' @export
ordinal_invariance_auto <- function(
  data, group, model, language = c("eng", "esp"), estimator = "auto",
  ordered = NULL, report = FALSE, criteria = "auto", missing = "auto",
  min_group_n = 20L, overwrite = FALSE, digits = 3, file_name = NULL,
  format = "excel", verbose = TRUE, fit_variant = c("scaled", "robust")
) {
  if (!is.character(estimator) || length(estimator) != 1L || is.na(estimator) ||
      !toupper(trimws(estimator)) %in% c("AUTO", "WLSMV", "ULSMV")) {
    stop("ordinal_invariance_auto only accepts estimator = 'auto', 'WLSMV', or 'ULSMV'; ML/MLR and other estimators are not ordinal estimators for this workflow.", call. = FALSE)
  }
  if (identical(ordered, FALSE)) {
    stop("ordinal_invariance_auto requires ordered indicators; ordered = FALSE is incompatible.", call. = FALSE)
  }
  .ifa_auto(data, group, model, language, estimator, ordered, report, criteria,
    missing, min_group_n, overwrite, digits, file_name, format, verbose,
    group_expression = substitute(group), group_environment = parent.frame(),
    original_call = match.call(), workflow = "ordinal", fit_variant = fit_variant)
}
