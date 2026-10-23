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
      if (identical(x$status, "not_evaluable")) {
        "not_evaluable"
      } else if (isTRUE(x$screening_passed)) {
        "screening_passed"
      } else {
        "review_suggestions"
      }
    ),
    rationale = c(
      .pm_null_default(x$estimador$justificacion, NA_character_),
      "Selected by the `rotation` argument.",
      if (identical(x$status, "not_evaluable")) {
        x$reason
      } else if (isTRUE(x$screening_passed)) {
        "The automated screening criteria were met."
      } else {
        paste(x$sugerencias, collapse = " | ")
      }
    ),
    stringsAsFactors = FALSE
  )
}

.pm_conditions_for_display <- function(x) {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0L) {
    return(data.frame(
      code = character(), message = character(), class = character(),
      call = character(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    code = as.character(x$code),
    message = as.character(x$message),
    class = vapply(x$class, function(value) paste(value, collapse = "/"), character(1)),
    call = vapply(
      x$call,
      function(value) if (is.null(value)) "" else paste(deparse(value), collapse = " "),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

.efa_metadata_table <- function(x) {
  data.frame(
    field = c(
      "status", "refine", "missing", "correlation_requested",
      "correlation_effective", "correlation_repaired", "parallel_iter",
      "parallel_quantile", "parallel_seed", "screening_passed"
    ),
    value = as.character(c(
      x$status,
      x$effective_settings$refine,
      x$effective_settings$missing,
      if (is.null(x$correlation)) NA_character_ else x$correlation$requested_method,
      if (is.null(x$correlation)) NA_character_ else x$correlation$effective_method,
      if (is.null(x$correlation)) FALSE else x$correlation$repaired,
      x$effective_settings$parallel_analysis$n_iter,
      x$effective_settings$parallel_analysis$quantile,
      .pm_null_default(x$effective_settings$parallel_analysis$seed, NA_integer_),
      x$screening_passed
    )),
    stringsAsFactors = FALSE
  )
}

.efa_references <- function() {
  c(
    "Lloret-Segura, S., Ferreres-Traver, A., Hernandez-Baeza, A., & Tomas-Marco, I. (2014). El analisis factorial exploratorio de los items: Una guia practica, revisada y actualizada. Anales de Psicologia, 30(3), 1151-1169. https://doi.org/10.6018/analesps.30.3.199361",
    "Holgado-Tello, F. P., Chacon-Moscoso, S., Barbero-Garcia, I., & Vila-Abad, E. (2010). Polychoric versus Pearson correlations in exploratory and confirmatory factor analysis of ordinal variables. Quality & Quantity, 44, 153-166. https://doi.org/10.1007/s11135-008-9190-y",
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
  x$warnings_display <- .pm_conditions_for_display(x$warnings)
  x$metadata_table <- .efa_metadata_table(x)
  x$conclusion_table <- data.frame(
    status = x$status,
    screening_passed = x$screening_passed,
    conclusion = x$conclusion,
    stringsAsFactors = FALSE
  )
  x$references <- .efa_references()
  x$tables <- list(
    fit = x$fit_table,
    loadings = x$loadings,
    reliability = x$reliability_table,
    decisions = x$decisions,
    metadata = x$metadata_table,
    conclusion = x$conclusion_table,
    warnings = x$warnings_display
  )
  x
}

.cfa_reliability_table <- function(x) {
  rel <- x$reliability
  empty <- function() {
    .pm_empty_table(
      group = character(),
      level = character(),
      factor = character(),
      metric = character(),
      value = numeric(),
      status = character(),
      reason = character()
    )
  }
  if (is.null(rel) || identical(rel$type, "error")) {
    return(empty())
  }

  if (identical(rel$type, "standard") && is.data.frame(rel$rows)) {
    rows <- data.frame(
      group = rel$rows$group,
      level = rel$rows$level,
      factor = rel$rows$factor,
      metric = rel$rows$coefficient,
      value = rel$rows$value,
      status = rel$rows$status,
      reason = rel$rows$reason,
      stringsAsFactors = FALSE
    )
    rownames(rows) <- NULL
    return(rows)
  }

  if (identical(rel$type, "standard") && !is.null(rel$values)) {
    values <- rel$values
    rows <- expand.grid(
      metric = rownames(values),
      factor = colnames(values),
      stringsAsFactors = FALSE
    )
    rows$value <- as.numeric(values[cbind(rows$metric, rows$factor)])
    rows$group <- "group1"
    rows$level <- "L1"
    valid <- is.finite(rows$value) & rows$value >= 0 & rows$value <= 1
    rows$status <- ifelse(valid, "ok", "not_evaluable")
    rows$reason <- ifelse(
      valid,
      NA_character_,
      ifelse(is.finite(rows$value), "outside_mathematical_domain", "non_finite_coefficient")
    )
    rows$value[!valid] <- NA_real_
    return(rows[, c("group", "level", "factor", "metric", "value", "status", "reason")])
  }

  if (identical(rel$type, "second_order") && is.data.frame(rel$rows)) {
    required <- c("group", "level", "factor", "metric", "value", "status", "reason")
    if (all(required %in% names(rel$rows))) {
      rows <- rel$rows[, required, drop = FALSE]
      rownames(rows) <- NULL
      return(rows)
    }
  }

  if (identical(rel$type, "bifactor")) {
    groups <- rel$groups
    if (is.null(groups) && !is.null(rel$omega_total)) {
      groups <- list(group1 = rel)
    }
    rows <- lapply(names(groups), function(group) {
      value <- groups[[group]]
      estimates <- c(
        omega_total = value$omega_total,
        omega_h = value$omega_h,
        ecv = value$ecv,
        puc = value$puc,
        h_general = value$h_general,
        stats::setNames(as.numeric(value$omega_hs), paste0("omega_hs:", names(value$omega_hs))),
        stats::setNames(as.numeric(value$h_specific), paste0("h:", names(value$h_specific)))
      )
      factors <- c(
        rel$general_factor, rel$general_factor, rel$general_factor, rel$general_factor,
        rel$general_factor, names(value$omega_hs), names(value$h_specific)
      )
      valid <- is.finite(estimates) & estimates >= 0 & estimates <= 1
      data.frame(
        group = group,
        level = "L1",
        factor = factors,
        metric = sub(":.*$", "", names(estimates)),
        value = ifelse(valid, as.numeric(estimates), NA_real_),
        status = ifelse(valid, "ok", "not_evaluable"),
        reason = ifelse(
          valid,
          NA_character_,
          ifelse(
            names(estimates) == "puc" && !is.null(value$puc_diagnostic$code),
            value$puc_diagnostic$code,
            ifelse(is.finite(estimates), "outside_mathematical_domain", "non_finite_coefficient")
          )
        ),
        stringsAsFactors = FALSE
      )
    })
    if (length(rows)) return(do.call(rbind, rows))
  }

  empty()
}

.cfa_loadings_table <- function(x) {
  lower <- as.data.frame(x$loadings, stringsAsFactors = FALSE)
  if (nrow(lower)) lower$level <- "L1"
  higher <- x$higher_order_loadings
  if (is.null(higher)) higher <- x$loadings_ho
  if (!is.null(higher) && nrow(higher)) {
    higher <- as.data.frame(higher, stringsAsFactors = FALSE)
    higher$level <- "L2"
  } else {
    higher <- lower[0, , drop = FALSE]
  }
  columns <- union(names(lower), names(higher))
  add_missing <- function(value) {
    missing <- setdiff(columns, names(value))
    for (name in missing) value[[name]] <- NA
    value[, columns, drop = FALSE]
  }
  lower <- add_missing(lower)
  higher <- add_missing(higher)
  rows <- rbind(lower, higher)
  rownames(rows) <- NULL
  rows[, c("level", setdiff(names(rows), "level")), drop = FALSE]
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
  x$loadings_table <- .cfa_loadings_table(x)
  x$reliability_table <- .cfa_reliability_table(x)
  x$decisions <- .cfa_decisions_table(x)
  x$references <- .build_references(x)
  x$tables <- list(
    fit = x$fit_indices,
    loadings = x$loadings_table,
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
#' @keywords internal
#' @export
summary.efa_auto <- function(object, ...) {
  list(
    fit = object$fit_table,
    loadings = object$loadings,
    reliability = object$reliability_table,
    decisions = object$decisions,
    warnings = object$warnings_display,
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
#' @keywords internal
#' @export
summary.cfa_auto <- function(object, ...) {
  list(
    fit = object$fit_indices,
    loadings = object$loadings_table,
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
#' @return A list preserving raw fit evidence, decisions, diagnostics, criteria,
#'   structured conditions, conclusion, and references.
#'
#' @method summary factorial_invariance_auto
#' @keywords internal
#' @export
summary.factorial_invariance_auto <- function(object, ...) {
  list(
    fit = object$fit_table,
    decisions = .invariance_decisions_table(object),
    diagnostics = object$diagnostics,
    criteria = object$criteria,
    conditions = list(warnings = object$warnings, errors = object$errors),
    conclusion = object$conclusion,
    warnings = object$warnings,
    references = object$references
  )
}

#' Summarize a researcher-directed partial invariance result
#'
#' @param object Object returned by `invariance_partial()`.
#' @param ... Additional unused arguments.
#'
#' @return A structured list preserving raw fit evidence, releases, anchors,
#'   comparisons, conditions, and the cautious conclusion.
#'
#' @method summary psychomatic_invariance_partial
#' @keywords internal
#' @export
summary.psychomatic_invariance_partial <- function(object, ...) {
  list(
    fit = object$fit_table,
    freed_constraints = object$freed_constraints,
    rationale = object$rationale_table,
    anchors = object$anchors,
    anchor_counts = object$anchor_counts,
    comparison = object$comparison,
    criteria = object$criteria,
    conditions = list(warnings = object$warnings, errors = object$errors),
    conclusion_status = object$conclusion_status,
    conclusion = object$conclusion_text,
    references = object$references
  )
}
