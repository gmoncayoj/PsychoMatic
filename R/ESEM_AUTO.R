#' Exploratory structural equation measurement model with transparent rotation
#'
#' Fits one exploratory measurement block directly in `lavaan`, retaining all
#' within-block cross-loadings and rotated parameter uncertainty. With no
#' structural regressions this is EFA in the terminology of `lavaan`; ESEM is
#' also used for these measurement models in the psychometric literature.
#'
#' @section Choosing a rotation:
#' **Target** is useful when theory or previous research identifies the primary
#' dimension for each item, while allowing secondary loadings. `factors` defines
#' the target pattern. Secondary loadings are oriented towards zero, not fixed
#' to zero. A target entry of 1 leaves the primary loading untargeted; it does
#' not fix that loading to 1. A data-derived target from the same sample remains
#' exploratory and requires independent validation for stronger generalization.
#'
#' **Geomin** is useful when the item-factor organization is uncertain. It seeks
#' a relatively simple pattern without prescribing substantive memberships.
#' `n_factors` remains an investigator decision. Factors retain neutral IDs
#' because order and signs do not establish their substantive interpretation.
#' The initial epsilon is 0.001, not a universal optimum; specify alternatives
#' through `rotation_args` and report the chosen value. Both options are oblique.
#'
#' Rotations of the same unrestricted block preserve global model fit; CFI and
#' RMSEA do not select the preferable rotation. Consider theoretical coherence,
#' loading patterns and stability. Thirty random starts do not establish a global
#' optimum or sampling stability. Reported standardized loadings are pattern
#' coefficients, not the oblique structure matrix.
#'
#' @param data Data frame, matrix, or `psychomatic_spec`.
#' @param factors Named list of primary item sets. Every selected item occurs
#'   once. Target requires this list. With geomin it only supplies selection,
#'   factor count and an optional CFA reference, not a rotation target.
#' @param n_factors Integer >= 2, required for geomin without `factors`.
#' @param items Explicit column names or positive column indices for raw geomin
#'   data without `factors`, for example `items = 2:13` or
#'   `items = names(MY_DATA)[2:13]`. Pass column selectors, not a data subset:
#'   `MY_DATA[, 2:13]` is a data frame and is not a valid `items` argument.
#'   This argument is required for raw geomin data, even if all columns are items;
#'   in that case use `items = names(MY_DATA)`. NULL with `factors` or a specification.
#' @param rotation `"target"` (default) or `"geomin"`.
#' @param ordered NULL inherits/detects measurement levels heuristically; TRUE
#'   treats all selected items as ordinal, FALSE requests numeric continuous
#'   treatment, or a character vector identifies ordinal items. Nominal factors
#'   require explicit ordinal mapping; ordered factors cannot silently become
#'   continuous scores.
#' @param estimator NULL selects WLSMV with ordinal items, otherwise MLR.
#'   Explicit choices: WLSMV/ULSMV for ordinal or mixed items, ML/MLR for continuous.
#' @param missing `"listwise"` (default), `"fiml"` for continuous ML/MLR,
#'   or `"pairwise"` for ordinal WLSMV/ULSMV. No automatic imputation.
#' @param rotation_args Named options: `rstarts` (30), `geomin_epsilon` (0.001,
#'   geomin only), `gpa_tol` (1e-5), `max_iter` (10000).
#' @param seed NULL (default) leaves RNG management to lavaan without imposing
#'   a wrapper seed. An explicit nonnegative integer makes rotation starts
#'   reproducible and restores the incoming RNG state after fitting.
#' @param compare_cfa Whether to fit the simple-structure CFA specified by
#'   `factors` on the same data/options. Comparison is descriptive, without
#'   automatic model selection, LRT or factor alignment.
#' @param fit_variant `"scaled"`, `"robust"` or `"base"`. Missing requested
#'   indices remain NA with a reason. Corrected chi-square/df/p use scaled
#'   variants unless base is selected; SRMR is always base.
#' @param language Output language: `"eng"` or `"esp"`.
#' @param digits Number of display decimals.
#' @param report Whether to write a report, or an explicit report path.
#' @param file_name Report path/base; NULL uses a unique name in `getwd()`.
#' @param format `"excel"`, `"word"` or both.
#' @param overwrite Whether existing report files may be replaced.
#' @param verbose Whether to print results and the manuscript draft.
#' @return An `esem_auto`/`psychomatic_result` object. `fit` contains the native
#'   lavaan fit with safe item aliases; `results$item_map` maps original names.
#'   `results` contains loadings/IC95, item diagnostics, factor correlations,
#'   fit variants, raw covariance residuals, target, maps and syntax.
#'   Loading tables distinguish original-scale (`raw_ci_lower`, `raw_ci_upper`)
#'   from standardized (`ci_lower`, `ci_upper`) native 95% intervals.
#'   `effective_settings$preparation` records category order, declared reverse
#'   scoring and the effective bounds, including inferred bounds. Report encoding
#'   is the input category order before reverse scoring, which is recorded separately.
#'   `manuscript_text` and `results$manuscript` contain a data-based draft placed
#'   last under Interpretation in reader-facing console/Word/Excel output.
#'   Preparation, item/factor maps and versions remain in the object and exports
#'   but are omitted from the console. The text is not an independent
#'   validation or substantive interpretation of unlabeled factors. Failed or
#'   inadmissible fits produce limitations rather than favorable conclusions.
#' @references
#' Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation modeling.
#' Structural Equation Modeling, 16(3), 397–438. doi:10.1080/10705510903008204.
#'
#' Marsh, H., & Alamer, A. (2024). When and how to use set-exploratory structural
#' equation modelling to test structural models: A tutorial using the R package
#' lavaan. British Journal of Mathematical and Statistical Psychology, 77(3),
#' 459–476. doi:10.1111/bmsp.12336.
#' @examples
#' \dontrun{
#' # Replace MY_DATA with your data frame.
#' # Ordinal items in columns 2:13; column 1 is not an item.
#' result <- esem_auto(data = MY_DATA,
#'                     n_factors = 3,
#'                     items = 2:13,
#'                     rotation = "geomin",
#'                     ordered = TRUE,
#'                     language = "esp")
#' result$manuscript_text
#'
#' # Continuous indicators: MY_DATA must contain only item columns here.
#' result <- esem_auto(data = MY_DATA,
#'                     n_factors = 3,
#'                     items = names(MY_DATA),
#'                     rotation = "geomin",
#'                     ordered = FALSE,
#'                     language = "esp")
#'
#' # Target: replace the item names and primary memberships with your theory.
#' factors <- list(Factor1 = c("R1", "R2", "R3", "R4"),
#'                 Factor2 = c("R5", "R6", "R7", "R8"),
#'                 Factor3 = c("R9", "R10", "R11", "R12"))
#' result <- esem_auto(data = MY_DATA,
#'                     factors = factors,
#'                     rotation = "target",
#'                     ordered = TRUE,
#'                     language = "esp")
#' }
#' @export
esem_auto <- function(data, factors = NULL, n_factors = NULL, items = NULL,
                       rotation = "target", ordered = NULL, estimator = NULL,
                       missing = "listwise", rotation_args = list(), seed = NULL,
                       compare_cfa = FALSE, fit_variant = c("scaled", "robust", "base"),
                       language = "eng", digits = 3, report = FALSE, file_name = NULL,
                       format = "excel", overwrite = FALSE, verbose = TRUE) {
  fit_variant <- match.arg(fit_variant)
  language <- .pm_normalize_language(language)
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  if (!is.null(seed)) seed <- as.integer(.pm_validate_scalar_number(seed, "seed", 0,
                                                .Machine$integer.max, integer = TRUE))
  if (!is.logical(compare_cfa) || length(compare_cfa) != 1L || is.na(compare_cfa)) {
    stop("`compare_cfa` must be TRUE or FALSE.", call. = FALSE)
  }
  if (compare_cfa && is.null(factors)) stop("`compare_cfa` requires `factors`.", call. = FALSE)
  prepared <- .pm_esem_prepare(data, factors, n_factors, items, rotation, ordered,
                               estimator, missing, rotation_args, language)
  fitted <- .pm_esem_fit(prepared, seed)
  out <- .pm_esem_result(prepared, fitted, match.call(), fit_variant, language)
  out$requested_settings <- list(factors = factors, n_factors = n_factors, items = items,
    rotation = rotation, ordered = ordered, estimator = estimator, missing = missing,
    rotation_args = rotation_args, seed = seed, compare_cfa = compare_cfa,
    fit_variant = fit_variant, language = language)
  out$effective_settings["seed"] <- list(seed)
  out$effective_settings$compare_cfa <- compare_cfa
  if (compare_cfa && identical(out$conclusion$status, "evaluable")) {
    cmp <- .pm_esem_compare_cfa(prepared, out$fit, seed, fit_variant)
    out$results$cfa_comparison <- cmp
    out$warnings <- rbind(out$warnings, .pm_conditions_to_table(cmp$warnings, "warnings"))
  }
  out$manuscript_text <- .pm_esem_manuscript(out, language, digits)
  out$results$manuscript <- data.frame(Text = out$manuscript_text, stringsAsFactors = FALSE)
  .pm_finish_analysis(out, output, language)
}
