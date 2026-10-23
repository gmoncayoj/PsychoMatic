.pm_spearman_brown <- function(r) {
  if (!is.numeric(r) || length(r) != 1L || is.na(r) || !is.finite(r)) {
    return(NA_real_)
  }
  denominator <- 1 + r
  if (abs(denominator) <= sqrt(.Machine$double.eps)) return(NA_real_)
  2 * r / denominator
}

.pm_reliability_row <- function(coefficient,
                                estimate,
                                estimand,
                                 matrix,
                                 keying,
                                 factor = NA_character_,
                                 group = NA_character_,
                                 domain = c(-Inf, Inf),
                                 reason = NA_character_) {
  raw_value <- if (is.numeric(estimate) && length(estimate) == 1L) {
    as.numeric(estimate)
  } else {
    NA_real_
  }
  status <- "ok"
  effective_reason <- reason

  if (!is.na(effective_reason)) {
    status <- "not_evaluable"
  } else if (!is.finite(raw_value)) {
    status <- "not_evaluable"
    effective_reason <- "non_finite_coefficient"
  } else if (length(domain) != 2L || anyNA(domain) ||
             raw_value < domain[1L] || raw_value > domain[2L]) {
    status <- "not_evaluable"
    effective_reason <- "outside_mathematical_domain"
  }

  data.frame(
    group = as.character(group)[1L],
    factor = as.character(factor)[1L],
    coefficient = as.character(coefficient)[1L],
    estimate = if (identical(status, "ok")) raw_value else NA_real_,
    estimand = as.character(estimand)[1L],
    matrix = as.character(matrix)[1L],
    keying = as.character(keying)[1L],
    status = status,
    reason = if (identical(status, "ok")) NA_character_ else effective_reason,
    diagnostic_value = I(list(raw_value)),
    stringsAsFactors = FALSE
  )
}

.pm_alpha_from_matrix <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) < 2L ||
      nrow(x) != ncol(x) || any(!is.finite(x))) {
    return(list(value = NA_real_, reason = "invalid_covariance_matrix"))
  }
  denominator <- sum(x)
  if (!is.finite(denominator) || abs(denominator) <= sqrt(.Machine$double.eps)) {
    return(list(value = NA_real_, reason = "zero_total_score_variance"))
  }
  p <- ncol(x)
  list(
    value = p / (p - 1) * (1 - sum(diag(x)) / denominator),
    reason = NA_character_
  )
}

.pm_ordinal_numeric <- function(x) {
  if (is.ordered(x)) return(as.numeric(x))
  if (is.numeric(x)) return(as.numeric(x))
  stop("Ordinal items must be ordered factors or numeric values.", call. = FALSE)
}

.pm_prepare_reliability_input <- function(data, items, level, language) {
  if (inherits(data, "psychomatic_spec")) {
    if (!is.null(items)) {
      stop("`items` must be NULL when `data` is a psychomatic_spec.", call. = FALSE)
    }
    spec <- data
    language_effective <- .pm_normalize_language(language, spec$language)
  } else {
    raw <- .pm_as_data_frame(data)
    selection <- if (is.null(items)) seq_along(raw) else items
    spec <- psychomatic_spec(
      raw,
      items = selection,
      level = level,
      language = language
    )
    language_effective <- spec$language
  }

  item_data <- spec$data[spec$item_index]
  names(item_data) <- spec$items
  levels <- unname(spec$measurement_level)
  numeric_data <- as.data.frame(
    lapply(seq_along(item_data), function(i) {
      if (identical(levels[i], "ordinal")) {
        .pm_ordinal_numeric(item_data[[i]])
      } else {
        as.numeric(item_data[[i]])
      }
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(numeric_data) <- spec$items

  reverse <- spec$reverse
  if (length(reverse)) {
    for (item in reverse) {
      values <- numeric_data[[item]]
      limits <- spec$ranges[[item]]
      if (is.null(limits)) {
        finite <- values[is.finite(values)]
        if (!length(finite)) next
        limits <- range(finite)
      }
      numeric_data[[item]] <- limits[1L] + limits[2L] - values
    }
  }

  effective_level <- if (all(levels == "ordinal")) {
    "ordinal"
  } else if (all(levels == "continuous")) {
    "continuous"
  } else {
    "mixed"
  }

  list(
    spec = spec,
    data = numeric_data,
    level = effective_level,
    language = language_effective,
    keying = if (length(reverse)) paste0("applied:", paste(reverse, collapse = ",")) else "none"
  )
}

.pm_semtools_available <- function() {
  requireNamespace("semTools", quietly = TRUE)
}

.pm_semtools_coefficients <- function(fit, keying) {
  if (!inherits(fit, "lavaan")) {
    stop("`model` must be a fitted lavaan model.", call. = FALSE)
  }
  if (!.pm_semtools_available()) {
    return(list(
      rows = .pm_reliability_row(
        "composite_reliability", NA_real_, "model_based", "model_implied",
        keying, reason = "semtools_not_available", domain = c(0, 1)
      ),
      diagnostics = list(status = "not_available")
    ))
  }

  captured_cr <- .pm_capture_conditions(semTools::compRelSEM(fit))
  captured_ave <- .pm_capture_conditions(semTools::AVE(fit))
  if (!is.null(captured_cr$error) || !is.null(captured_ave$error)) {
    condition <- if (!is.null(captured_cr$error)) captured_cr$error else captured_ave$error
    return(list(
      rows = .pm_reliability_row(
        "composite_reliability", NA_real_, "model_based", "model_implied",
        keying, reason = "semtools_error", domain = c(0, 1)
      ),
      diagnostics = list(status = "not_evaluable", condition = condition)
    ))
  }

  group_labels <- tryCatch(
    lavaan::lavInspect(fit, "group.label"),
    error = function(e) character()
  )
  normalized <- rbind(
    .pm_semtools_long(
      captured_cr$value,
      coefficient = "composite_reliability",
      group_labels = group_labels,
      level = "L1"
    ),
    .pm_semtools_long(
      captured_ave$value,
      coefficient = "ave",
      group_labels = group_labels,
      level = "L1"
    )
  )
  rows <- lapply(seq_len(nrow(normalized)), function(i) {
    coefficient <- normalized$coefficient[[i]]
    .pm_reliability_row(
      coefficient,
      normalized$value[[i]],
      if (identical(coefficient, "ave")) {
        "average_variance_extracted"
      } else {
        "model_based"
      },
      "model_implied",
      keying,
      factor = normalized$factor[[i]],
      group = normalized$group[[i]],
      domain = c(0, 1)
    )
  })
  if (!length(rows)) {
    rows <- list(.pm_reliability_row(
      "composite_reliability", NA_real_, "model_based", "model_implied",
      keying, reason = "semtools_output_unrecognized", domain = c(0, 1)
    ))
  }

  list(
    rows = do.call(rbind, rows),
    diagnostics = list(
      status = "ok",
      compRelSEM = captured_cr$value,
      AVE = captured_ave$value,
      warnings = c(captured_cr$warnings, captured_ave$warnings),
      messages = c(captured_cr$messages, captured_ave$messages)
    )
  )
}

#' Automated reliability analysis
#'
#' Computes reliability coefficients while keeping the statistical estimand,
#' response matrix, reverse-keying decision, and evaluability status explicit.
#' Continuous responses receive observed-score alpha; ordinal responses also
#' receive alpha from a named polychoric matrix. A fitted `lavaan` model adds
#' composite reliability and average variance extracted through current
#' `semTools` APIs when that suggested package is installed.
#' Without a supplied model, omega fits a single congeneric factor (ML for
#' continuous items; WLSMV for ordinal items). Ordinal omega refers to the
#' observed response score using the Green and Yang correction; ordinal alpha
#' refers to latent responses. This assumption requires substantive assessment.
#' Alpha uses pairwise covariances. Omega uses complete cases, refitted after
#' each item deletion; item-rest correlations use complete cases on all items.
#' Deletion alpha uses the corresponding submatrix of the full-scale matrix.
#'
#' @inheritParams efa_auto
#' @param data A data frame, matrix, or [psychomatic_spec()] object.
#' @param items Optional item names or positive indices. Must be `NULL` for a
#'   `psychomatic_spec` input.
#' @param level Measurement level passed to [psychomatic_spec()]: `"auto"`,
#'   `"continuous"`, or `"ordinal"`.
#' @param model Optional fitted `lavaan` CFA model.
#' @param report If `TRUE`, send the structured result to
#'   [psychomatic_report()]. No file is written when `FALSE`.
#' @param language Output language, `"eng"` by default; use `"esp"` for Spanish.
#'   Explicit `NULL` inherits from a specification or the
#'   package option.
#'
#' @return An object of classes `psychomatic_reliability` and
#'   `psychomatic_result`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' reliability_auto(MY_DATA[, 3:14])
#' }
#'
#' @export
reliability_auto <- function(data,
                             items = NULL,
                             level = "auto",
                             model = NULL,
                             report = FALSE,
                             language = "eng",
                             digits = 3, file_name = NULL, format = "excel",
                             overwrite = FALSE, verbose = TRUE) {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  call <- match.call()
  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("`report` must be TRUE or FALSE.", call. = FALSE)
  }
  prepared <- .pm_prepare_reliability_input(data, items, level, language)
  x <- prepared$data
  p <- ncol(x)
  keying <- prepared$keying
  rows <- list()
  diagnostics <- list(correlation_matrix = NULL, covariance_matrix = NULL)
  captured_conditions <- list()

  nonfinite <- any(vapply(x, function(value) any(!is.finite(value) & !is.na(value)), logical(1)))
  if (nonfinite) {
    rows[[1L]] <- .pm_reliability_row(
      "alpha_observed", NA_real_, "observed_score", "raw_covariance", keying,
      reason = "non_finite_input", domain = c(-Inf, 1)
    )
  } else if (p < 2L) {
    rows[[1L]] <- .pm_reliability_row(
      "alpha_observed", NA_real_, "observed_score", "raw_covariance", keying,
      reason = "insufficient_items", domain = c(-Inf, 1)
    )
  } else if (identical(prepared$level, "mixed")) {
    rows[[1L]] <- .pm_reliability_row(
      "alpha_observed", NA_real_, "observed_score", "raw_covariance", keying,
      reason = "mixed_measurement_levels", domain = c(-Inf, 1)
    )
  } else {
    covariance <- stats::cov(x, use = "pairwise.complete.obs")
    dimnames(covariance) <- list(names(x), names(x))
    diagnostics$covariance_matrix <- covariance
    alpha_raw <- .pm_alpha_from_matrix(covariance)
    rows[[length(rows) + 1L]] <- .pm_reliability_row(
      "alpha_observed", alpha_raw$value, "observed_score", "raw_covariance",
      keying, domain = c(-Inf, 1), reason = alpha_raw$reason
    )

    if (identical(prepared$level, "ordinal")) {
      poly <- NULL
      invisible(utils::capture.output({
        poly <- .pm_capture_conditions(psych::polychoric(x))
      }))
      captured_conditions$polychoric <- poly
      if (is.null(poly$error) && !is.null(poly$value$rho)) {
        correlation <- as.matrix(poly$value$rho)
        dimnames(correlation) <- list(names(x), names(x))
        diagnostics$correlation_matrix <- correlation
        alpha_ordinal <- .pm_alpha_from_matrix(correlation)
        rows[[length(rows) + 1L]] <- .pm_reliability_row(
          "alpha_ordinal", alpha_ordinal$value, "latent_response",
          "polychoric", keying, domain = c(-Inf, 1), reason = alpha_ordinal$reason
        )
      } else {
        rows[[length(rows) + 1L]] <- .pm_reliability_row(
          "alpha_ordinal", NA_real_, "latent_response", "polychoric", keying,
          domain = c(-Inf, 1), reason = "polychoric_not_evaluable"
        )
      }
    } else {
      correlation <- stats::cor(x, use = "pairwise.complete.obs")
      dimnames(correlation) <- list(names(x), names(x))
      diagnostics$correlation_matrix <- correlation
    }

    if (p == 2L) {
      correlation_available <- is.matrix(diagnostics$correlation_matrix) &&
        nrow(diagnostics$correlation_matrix) >= 2L &&
        ncol(diagnostics$correlation_matrix) >= 2L
      r <- if (correlation_available) {
        diagnostics$correlation_matrix[1L, 2L]
      } else {
        NA_real_
      }
      sb <- .pm_spearman_brown(r)
      sb_reason <- if (!correlation_available && identical(prepared$level, "ordinal")) {
        "polychoric_not_evaluable"
      } else if (!is.finite(r)) {
        "non_finite_correlation"
      } else if (abs(1 + r) <= sqrt(.Machine$double.eps)) {
        "zero_denominator"
      } else {
        NA_character_
      }
      rows[[length(rows) + 1L]] <- .pm_reliability_row(
        "spearman_brown", sb, "two_item_composite",
        if (identical(prepared$level, "ordinal")) "polychoric" else "pearson",
        keying, domain = c(-Inf, 1), reason = sb_reason
      )
    }
  }

  semtools <- NULL
  omega <- NULL
  if (!is.null(model)) {
    semtools <- .pm_semtools_coefficients(model, keying)
    rows[[length(rows) + 1L]] <- semtools$rows
  } else {
    omega <- .pm_auto_omega(x, prepared$level, keying)
    rows[[length(rows) + 1L]] <- omega$row
  }
  item_analysis <- .pm_reliability_items(x, prepared$level,
    diagnostics$covariance_matrix, diagnostics$correlation_matrix, model)
  diagnostics$omega <- if (is.null(omega)) list(status = "supplied_model") else omega$diagnostics
  diagnostics$item_deletion <- item_analysis$diagnostics
  coefficients <- do.call(rbind, rows)
  rownames(coefficients) <- NULL
  status <- if (any(coefficients$status == "ok")) "ok" else "not_evaluable"

  diagnostics$conditions <- captured_conditions
  diagnostics$semTools <- if (is.null(semtools)) list(status = "not_requested") else semtools$diagnostics
  polychoric_conditions <- captured_conditions$polychoric
  warnings_structured <- .pm_bind_condition_tables(
    omega_warning = if (is.null(omega)) NULL else omega$diagnostics$warnings,
    deletion_warning = do.call(c, lapply(item_analysis$diagnostics, function(z) z$warnings)),
    polychoric_warning = if (is.null(polychoric_conditions)) NULL else polychoric_conditions$warnings,
    semtools_warning = if (is.null(semtools)) NULL else semtools$diagnostics$warnings
  )
  errors_structured <- .pm_bind_condition_tables(
    omega_error = if (is.null(omega)) NULL else omega$diagnostics$error,
    polychoric_error = if (is.null(polychoric_conditions)) NULL else polychoric_conditions$error,
    semtools_error = if (is.null(semtools)) NULL else semtools$diagnostics$condition
  )
  requested_level <- if (inherits(data, "psychomatic_spec")) {
    data$requested_settings$level
  } else {
    level
  }
  out <- .pm_result(
    class = "psychomatic_reliability",
    call = call,
    requested_settings = list(
      items = items,
      level = requested_level,
      model = !is.null(model),
      report = report,
      language = language
    ),
    effective_settings = list(
      items = names(x),
      level = prepared$level,
      keying = keying,
      missing = list(alpha = "pairwise", omega = "listwise_per_fit",
        item_rest = "listwise_full_scale", deletion_alpha = "full_matrix_submatrix"),
      language = prepared$language
    ),
    data_summary = list(n = nrow(x), p = p, complete = sum(stats::complete.cases(x))),
    results = list(coefficients = coefficients, item_stats = item_analysis$rows),
    warnings = warnings_structured,
    errors = errors_structured,
    diagnostics = diagnostics,
    decision_log = rbind(
      .pm_decision("measurement_level", paste(requested_level, collapse = ","), prepared$level,
                   "shared measurement-level contract"),
      .pm_decision("keying", "specification", keying,
                   if (identical(keying, "none")) "no reverse keys" else "reverse keys applied")
    ),
    conclusion = if (identical(status, "ok")) NA_character_ else .pm_not_evaluable(
      "No requested reliability coefficient was evaluable."
    )
  )
  out$status <- status
  out$coefficients <- coefficients
  out$item_stats <- item_analysis$rows
  out$references <- c(
    "Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. Psychometrika, 16, 297-334. https://doi.org/10.1007/BF02310555",
    "McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence Erlbaum Associates.",
    if (prepared$level == "ordinal") c(
      "Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using structural equation modeling: An alternative to coefficient alpha. Psychometrika, 74(1), 155-167. https://doi.org/10.1007/s11336-008-9099-3",
      "Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of coefficients alpha and theta for Likert rating scales. Journal of Modern Applied Statistical Methods, 6(1), 21-29. https://doi.org/10.22237/jmasm/1177992180"))

  .pm_finish_analysis(out, output, prepared$language)
}

#' Print an automated reliability result
#'
#' @param x A `psychomatic_reliability` object.
#' @param digits Number of digits used only for display.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#' @method print psychomatic_reliability
#' @keywords internal
#' @export
print.psychomatic_reliability <- function(x, digits = NULL, ...) {
  .pm_console_report(x, digits)
}

#' Summarize an automated reliability result
#'
#' @param object A `psychomatic_reliability` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list of coefficient, setting, and diagnostic tables/lists.
#' @method summary psychomatic_reliability
#' @keywords internal
#' @export
summary.psychomatic_reliability <- function(object, ...) {
  list(
    coefficients = object$results$coefficients,
    item_stats = object$results$item_stats,
    settings = object$effective_settings,
    diagnostics = object$diagnostics,
    conclusion = object$conclusion
  )
}

#' Plot an automated reliability result
#'
#' @param x A `psychomatic_reliability` object.
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return The plotted coefficient rows, invisibly.
#' @method plot psychomatic_reliability
#' @keywords internal
#' @export
plot.psychomatic_reliability <- function(x, ...) {
  rows <- x$results$coefficients
  rows <- rows[rows$status == "ok" & is.finite(rows$estimate), , drop = FALSE]
  if (!nrow(rows)) {
    graphics::plot.new()
    graphics::title("Reliability not evaluable")
    return(invisible(rows))
  }
  labels <- ifelse(is.na(rows$factor), rows$coefficient,
                   paste(rows$factor, rows$coefficient, sep = ": "))
  has_groups <- FALSE
  if ("group" %in% names(rows)) {
    group_keys <- unique(.pm_group_key(rows$group))
    has_groups <- length(group_keys) > 1L ||
      any(!is.na(rows$group) & nzchar(rows$group))
  }
  if (has_groups) {
    labels <- paste(.pm_group_label(rows$group), labels, sep = " | ")
  }
  rows$plot_label <- labels
  graphics::barplot(rows$estimate, names.arg = rows$plot_label, las = 2L, ...)
  invisible(rows)
}
