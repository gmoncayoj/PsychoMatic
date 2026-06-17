.pm_null_default <- function(x, default) {
  if (is.null(x)) default else x
}

.pm_empty_table <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)[0, , drop = FALSE]
}

.efa_fit_table <- function(x) {
  rmsea <- if (!is.null(x$efa$RMSEA)) x$efa$RMSEA[1] else NA_real_
  bic <- if (!is.null(x$efa$BIC)) x$efa$BIC else NA_real_

  data.frame(
    metric = c("chisq", "df", "pvalue", "tli", "rmsea", "bic", "kmo"),
    value = c(
      x$chi_cuadrado,
      x$grados_libertad,
      x$p_valor,
      .pm_null_default(x$efa$TLI, NA_real_),
      rmsea,
      bic,
      .pm_null_default(x$kmo_final$MSA, NA_real_)
    ),
    stringsAsFactors = FALSE
  )
}

.efa_reliability_table <- function(x) {
  if (is.null(x$fiabilidad) || length(x$fiabilidad) == 0L) {
    return(.pm_empty_table(
      factor = character(),
      n_items = integer(),
      alpha = numeric(),
      omega = numeric(),
      spearman_brown = numeric(),
      items = character()
    ))
  }

  do.call(rbind, lapply(names(x$fiabilidad), function(name) {
    rel <- x$fiabilidad[[name]]
    data.frame(
      factor = name,
      n_items = .pm_null_default(rel$n_items, length(.pm_null_default(rel$items, character()))),
      alpha = .pm_null_default(rel$alpha, NA_real_),
      omega = .pm_null_default(rel$omega, NA_real_),
      spearman_brown = .pm_null_default(rel$spearman_brown, NA_real_),
      items = paste(.pm_null_default(rel$items, character()), collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))
}

.efa_decisions_table <- function(x) {
  data.frame(
    step = c("estimator", "rotation", "solution"),
    decision = c(
      .pm_null_default(x$estimador$nombre, NA_character_),
      .pm_null_default(x$rotation, NA_character_),
      if (isTRUE(x$solucion_optima)) "optimal" else "review_suggestions"
    ),
    rationale = c(
      .pm_null_default(x$estimador$justificacion, NA_character_),
      "Selected by the `rotation` argument.",
      if (isTRUE(x$solucion_optima)) {
        "No additional respecification suggestions were generated."
      } else {
        paste(x$sugerencias, collapse = " | ")
      }
    ),
    stringsAsFactors = FALSE
  )
}

.efa_references <- function() {
  c(
    "Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30(2), 179-185. https://doi.org/10.1007/BF02289447",
    "Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. Biometrika, 57(3), 519-530. https://doi.org/10.1093/biomet/57.3.519",
    "Revelle, W. (2024). psych: Procedures for psychological, psychometric, and personality research. R package. https://CRAN.R-project.org/package=psych"
  )
}

.efa_add_tidy_fields <- function(x, call = NULL) {
  x$call <- call
  x$fit <- x$efa
  x$loadings <- x$cargas
  x$fit_table <- .efa_fit_table(x)
  x$reliability_table <- .efa_reliability_table(x)
  x$decisions <- .efa_decisions_table(x)
  x$warnings <- if (isTRUE(x$solucion_optima)) character(0) else x$sugerencias
  x$references <- .efa_references()
  x$tables <- list(
    fit = x$fit_table,
    loadings = x$loadings,
    reliability = x$reliability_table,
    decisions = x$decisions
  )
  x
}

.cfa_reliability_table <- function(x) {
  rel <- x$reliability
  if (is.null(rel) || identical(rel$type, "error")) {
    return(.pm_empty_table(
      factor = character(),
      metric = character(),
      value = numeric()
    ))
  }

  if (identical(rel$type, "standard") && !is.null(rel$values)) {
    values <- rel$values
    rows <- expand.grid(
      metric = rownames(values),
      factor = colnames(values),
      stringsAsFactors = FALSE
    )
    rows$value <- as.numeric(values[cbind(rows$metric, rows$factor)])
    return(rows[, c("factor", "metric", "value")])
  }

  if (identical(rel$type, "bifactor")) {
    return(data.frame(
      factor = c(rel$general_factor, names(rel$omega_hs), names(rel$h_specific)),
      metric = c("omega_h", rep("omega_hs", length(rel$omega_hs)), rep("h", length(rel$h_specific))),
      value = c(rel$omega_h, as.numeric(rel$omega_hs), as.numeric(rel$h_specific)),
      stringsAsFactors = FALSE
    ))
  }

  .pm_empty_table(factor = character(), metric = character(), value = numeric())
}

.cfa_decisions_table <- function(x) {
  data.frame(
    step = c("estimator", "model_type", "convergence"),
    decision = c(
      x$estimator,
      x$model_type_code,
      if (isTRUE(x$converged)) "converged" else "not_converged"
    ),
    rationale = c(
      x$estimator_reason,
      x$model_type,
      if (isTRUE(x$converged)) "lavaan reported convergence." else "lavaan did not report convergence."
    ),
    stringsAsFactors = FALSE
  )
}

.cfa_add_tidy_fields <- function(x, call = NULL) {
  x$call <- call
  x$reliability_table <- .cfa_reliability_table(x)
  x$decisions <- .cfa_decisions_table(x)
  x$references <- .build_references(x)
  x$tables <- list(
    fit = x$fit_indices,
    loadings = x$loadings,
    reliability = x$reliability_table,
    decisions = x$decisions
  )
  x
}

.invariance_decisions_table <- function(x) {
  fit_table <- x$fit_table
  data.frame(
    level = rownames(fit_table),
    decision = fit_table$decision,
    status = fit_table$status,
    stringsAsFactors = FALSE
  )
}

#' Summarize an automated EFA result
#'
#' @param object Object returned by `efa_auto()`.
#' @param ... Additional unused arguments.
#'
#' @return A list of tidy tables.
#'
#' @method summary efa_auto
#' @export
summary.efa_auto <- function(object, ...) {
  list(
    fit = object$fit_table,
    loadings = object$loadings,
    reliability = object$reliability_table,
    decisions = object$decisions,
    warnings = object$warnings,
    references = object$references
  )
}

#' Summarize an automated CFA result
#'
#' @param object Object returned by `cfa_auto()`.
#' @param ... Additional unused arguments.
#'
#' @return A list of tidy tables.
#'
#' @method summary cfa_auto
#' @export
summary.cfa_auto <- function(object, ...) {
  list(
    fit = object$fit_indices,
    loadings = object$loadings,
    reliability = object$reliability_table,
    decisions = object$decisions,
    warnings = object$warnings,
    references = object$references
  )
}

#' Summarize a factorial invariance result
#'
#' @param object Object returned by `factorial_invariance_auto()`.
#' @param ... Additional unused arguments.
#'
#' @return A list with fit, decisions, and conclusion.
#'
#' @method summary factorial_invariance_auto
#' @export
summary.factorial_invariance_auto <- function(object, ...) {
  list(
    fit = object$display_table,
    decisions = .invariance_decisions_table(object),
    conclusion = object$conclusion,
    warnings = object$notes
  )
}
