#' Sensitivity analysis for factor-retention decisions
#'
#' Repeats parallel analysis across defensible correlation methods, reference
#' quantiles, and random seeds. The result describes how much the retained
#' factor count changes across configurations; it does not select a winner.
#' Every configuration uses the same complete-case sample. When `data` is a
#' `psychomatic_spec`, its declared item measurement levels and language are
#' honored.
#'
#' @param data Data frame, matrix, or `psychomatic_spec` containing item
#'   responses.
#' @param correlations Correlation methods to compare. The default, `"auto"`,
#'   uses Pearson and Spearman for continuous items, polychoric and Spearman for
#'   ordinal items, and tetrachoric and Pearson for binary items.
#' @param quantiles One or more reference quantiles for parallel analysis.
#' @param n_iter Number of simulated data sets per configuration.
#' @param seed One or more non-negative integer seeds.
#' @inheritParams efa_auto
#' @param language Output language: `"esp"` or `"eng"`. `NULL` inherits the
#'   language of a `psychomatic_spec`, then the package option.
#'
#' @return An object of class `psychomatic_dimensionality` containing the full
#'   configuration grid, skipped-cell reasons, the observed solution range, and
#'   descriptive agreement. Effective correlation matrices are retained in
#'   `effective_matrices` and `results$effective_matrices`, keyed by requested
#'   method.
#'
#' @examples
#' set.seed(10)
#' factor_score <- rnorm(80)
#' items <- data.frame(
#'   i1 = factor_score + rnorm(80, sd = .5),
#'   i2 = factor_score + rnorm(80, sd = .5),
#'   i3 = factor_score + rnorm(80, sd = .5),
#'   i4 = factor_score + rnorm(80, sd = .5)
#' )
#' dimensionality_sensitivity(items, n_iter = 10, seed = 10)
#'
#' @keywords internal
#' @noRd
dimensionality_sensitivity <- function(data,
                                       correlations = "auto",
                                       quantiles = 0.95,
                                       n_iter = 500,
                                       seed = 1234,
                                       report = FALSE,
                                       language = NULL,
                                       digits = 3, file_name = NULL, format = "excel",
                                       overwrite = FALSE, verbose = TRUE) {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  call <- match.call()
  specification <- if (inherits(data, "psychomatic_spec")) data else NULL
  inherited_language <- if (is.null(specification)) NULL else specification$language
  language <- .pm_normalize_language(language, inherited = inherited_language)
  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("`report` must be TRUE or FALSE.", call. = FALSE)
  }
  n_iter <- .pm_validate_scalar_number(
    n_iter,
    name = "n_iter",
    lower = 1,
    integer = TRUE
  )
  if (!is.numeric(quantiles) || length(quantiles) == 0L || anyNA(quantiles) ||
      any(!is.finite(quantiles)) || any(quantiles < 0.5) || any(quantiles > 0.999)) {
    stop("`quantiles` must contain finite values between 0.5 and 0.999.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) == 0L || anyNA(seed) ||
      any(!is.finite(seed)) || any(seed < 0) ||
      any(seed > .Machine$integer.max) || any(seed != trunc(seed))) {
    stop("`seed` must contain non-negative integers.", call. = FALSE)
  }
  quantiles <- unique(as.numeric(quantiles))
  seeds <- unique(as.integer(seed))

  if (!is.null(specification)) {
    item_data <- specification$data[specification$item_index]
    names(item_data) <- specification$items
    measurement_levels <- specification$measurement_level
    if (!is.null(names(measurement_levels))) {
      measurement_levels <- measurement_levels[names(item_data)]
    }
    requested_missing <- specification$missing
  } else {
    item_data <- .pm_as_data_frame(data)
    measurement_levels <- NULL
    requested_missing <- "unspecified"
  }
  if (nrow(item_data) < 3L || ncol(item_data) < 2L) {
    stop("`data` must contain at least three rows and two items.", call. = FALSE)
  }
  n_original <- nrow(item_data)
  complete_rows <- stats::complete.cases(item_data)
  item_data <- item_data[complete_rows, , drop = FALSE]
  if (nrow(item_data) < 3L) {
    stop("`data` must contain at least three complete cases.", call. = FALSE)
  }

  method_resolution <- .dimensionality_methods(
    item_data,
    correlations,
    measurement_levels = measurement_levels
  )
  methods <- method_resolution$methods
  roles <- method_resolution$roles
  grid <- expand.grid(
    requested_method = methods,
    quantile = quantiles,
    seed = seeds,
    stringsAsFactors = FALSE
  )
  grid$method_role <- unname(roles[grid$requested_method])
  grid$effective_method <- NA_character_
  grid$n_iter <- as.integer(n_iter)
  grid$n_factors <- NA_integer_
  grid$status <- "skipped"
  grid$reason <- NA_character_
  grid$repaired <- FALSE
  grid$pre_repair_min_eigenvalue <- NA_real_

  warning_conditions <- list()
  decision_rows <- list(.pm_decision(
    "missing_data",
    requested_missing,
    "complete_cases",
    sprintf(
      "A fixed complete-case sample removed %d of %d rows before every configuration.",
      n_original - nrow(item_data),
      n_original
    )
  ))
  correlation_cache <- list()
  for (method in methods) {
    correlation_cache[[method]] <- tryCatch(
      .efa_estimate_correlations(item_data, method = method, repair = "none"),
      error = identity
    )
  }
  effective_methods <- vapply(
    methods,
    function(method) {
      correlation <- correlation_cache[[method]]
      if (inherits(correlation, "condition")) NA_character_ else correlation$effective_method
    },
    character(1)
  )
  names(effective_methods) <- methods
  effective_matrices <- lapply(methods, function(method) {
    correlation <- correlation_cache[[method]]
    if (inherits(correlation, "condition")) NULL else correlation$matrix
  })
  names(effective_matrices) <- methods

  for (i in seq_len(nrow(grid))) {
    method <- grid$requested_method[i]
    correlation <- correlation_cache[[method]]
    if (inherits(correlation, "condition")) {
      grid$reason[i] <- conditionMessage(correlation)
      decision_rows[[length(decision_rows) + 1L]] <- .pm_decision(
        "dimensionality_configuration",
        method,
        "skipped",
        conditionMessage(correlation),
        "warning"
      )
      next
    }

    grid$effective_method[i] <- correlation$effective_method
    grid$repaired[i] <- correlation$repaired
    grid$pre_repair_min_eigenvalue[i] <- correlation$pre_repair_min_eigenvalue
    warning_conditions <- c(warning_conditions, correlation$warnings)
    min_eigen <- min(eigen(correlation$matrix, symmetric = TRUE, only.values = TRUE)$values)
    if (!is.finite(min_eigen) || min_eigen <= sqrt(.Machine$double.eps)) {
      grid$reason[i] <- paste0(
        "The correlation matrix is not positive definite (minimum eigenvalue ",
        format(min_eigen, digits = 5), "); repair was not requested."
      )
      decision_rows[[length(decision_rows) + 1L]] <- .pm_decision(
        "dimensionality_configuration",
        method,
        "skipped",
        grid$reason[i],
        "warning"
      )
      next
    }

    parallel <- tryCatch(
      .efa_parallel_analysis(
        correlation$matrix,
        n_obs = nrow(item_data),
        n_iter = n_iter,
        quantile = grid$quantile[i],
        seed = grid$seed[i]
      ),
      error = identity
    )
    if (inherits(parallel, "condition")) {
      grid$reason[i] <- conditionMessage(parallel)
      decision_rows[[length(decision_rows) + 1L]] <- .pm_decision(
        "dimensionality_configuration",
        method,
        "skipped",
        conditionMessage(parallel),
        "warning"
      )
      next
    }

    warning_conditions <- c(warning_conditions, parallel$warnings)
    grid$n_factors[i] <- parallel$n_factors
    grid$status[i] <- "ok"
    grid$reason[i] <- ""
    decision_rows[[length(decision_rows) + 1L]] <- .pm_decision(
      "dimensionality_configuration",
      paste(method, grid$quantile[i], grid$seed[i], sep = "/"),
      as.character(parallel$n_factors),
      "Parallel-analysis factor count retained for descriptive sensitivity."
    )
  }

  evaluated <- grid$n_factors[grid$status == "ok"]
  solution_range <- c(
    minimum = if (length(evaluated)) min(evaluated) else NA_integer_,
    maximum = if (length(evaluated)) max(evaluated) else NA_integer_
  )
  agreement <- data.frame(
    n_configurations = nrow(grid),
    n_evaluable = length(evaluated),
    n_skipped = sum(grid$status == "skipped"),
    n_distinct_solutions = length(unique(evaluated)),
    exact_agreement = length(evaluated) > 0L && length(unique(evaluated)) == 1L,
    observed_solutions = paste(sort(unique(evaluated)), collapse = ", "),
    stringsAsFactors = FALSE
  )
  conclusion <- if (language == "esp") {
    if (length(evaluated)) {
      sprintf(
        "Las configuraciones evaluables produjeron un rango descriptivo de %d a %d factores; interprete la sensibilidad junto con la teoria y la calidad de los items.",
        solution_range[["minimum"]], solution_range[["maximum"]]
      )
    } else {
      "Ninguna configuracion fue evaluable; revise las razones registradas por celda."
    }
  } else if (length(evaluated)) {
    sprintf(
      "Evaluable configurations produced a descriptive range of %d to %d factors; interpret sensitivity alongside theory and item quality.",
      solution_range[["minimum"]], solution_range[["maximum"]]
    )
  } else {
    "No configuration was evaluable; inspect the reason recorded for each cell."
  }
  decision_log <- if (length(decision_rows)) {
    do.call(rbind, decision_rows)
  } else {
    .pm_decision(
      "dimensionality_configuration",
      "none",
      "none",
      "No configuration was generated.",
      "warning"
    )[0, , drop = FALSE]
  }

  out <- .pm_result(
    class = "psychomatic_dimensionality",
    call = call,
    requested_settings = list(
      correlations = correlations,
      quantiles = quantiles,
      n_iter = as.integer(n_iter),
      seed = seeds,
      missing = requested_missing,
      report = report,
      language = language
    ),
    effective_settings = list(
      correlations = effective_methods,
      quantiles = quantiles,
      n_iter = as.integer(n_iter),
      seed = seeds,
      missing = "complete_cases",
      language = language
    ),
    data_summary = list(
      n = nrow(item_data),
      n_original = n_original,
      n_complete = nrow(item_data),
      p = ncol(item_data),
      measurement = method_resolution$measurement
    ),
    results = list(
      configurations = grid,
      solution_range = solution_range,
      agreement = agreement,
      effective_matrices = effective_matrices
    ),
    warnings = warning_conditions,
    errors = list(),
    diagnostics = list(
      skipped = grid[grid$status == "skipped", , drop = FALSE],
      correlations = correlation_cache
    ),
    decision_log = decision_log,
    conclusion = conclusion
  )
  out$configurations <- grid
  out$solution_range <- solution_range
  out$agreement <- agreement
  out$effective_matrices <- effective_matrices
  out$language <- language

  .pm_finish_analysis(out, output, language)
}

.dimensionality_methods <- function(data, correlations, measurement_levels = NULL) {
  supported <- c("pearson", "spearman", "polychoric", "tetrachoric")
  if (!is.character(correlations) || length(correlations) == 0L || anyNA(correlations)) {
    stop("`correlations` must be 'auto' or supported method names.", call. = FALSE)
  }
  correlations <- tolower(trimws(correlations))
  if (identical(correlations, "auto")) {
    levels <- if (is.null(measurement_levels)) {
      vapply(data, .pm_measurement_level, character(1))
    } else {
      measurement_levels <- as.character(measurement_levels)
      if (length(measurement_levels) != ncol(data) || anyNA(measurement_levels) ||
          any(!measurement_levels %in% c("ordinal", "continuous"))) {
        stop(
          "The specification must declare each item as ordinal or continuous.",
          call. = FALSE
        )
      }
      measurement_levels
    }
    category_counts <- vapply(data, function(x) length(unique(stats::na.omit(x))), integer(1))
    if (all(levels == "ordinal") && all(category_counts == 2L)) {
      methods <- c("tetrachoric", "pearson")
      roles <- c(
        tetrachoric = "primary for binary items",
        pearson = "sensitivity contrast for binary items"
      )
      measurement <- "binary"
    } else if (all(levels == "ordinal")) {
      methods <- c("polychoric", "spearman")
      roles <- c(
        polychoric = "primary for ordinal items",
        spearman = "rank-based sensitivity contrast"
      )
      measurement <- "ordinal"
    } else if (all(levels %in% c("continuous", "ordinal")) &&
               all(vapply(data, is.numeric, logical(1)))) {
      methods <- c("pearson", "spearman")
      roles <- c(
        pearson = "linear-correlation configuration",
        spearman = "rank-based sensitivity contrast"
      )
      measurement <- if (all(levels == "continuous")) "continuous" else "mixed_numeric"
    } else {
      stop(
        "Automatic correlation selection requires numeric items or explicitly ordered factors.",
        call. = FALSE
      )
    }
  } else {
    unknown <- setdiff(correlations, supported)
    if (length(unknown)) {
      stop(
        sprintf("Unsupported `correlations`: %s.", paste(unknown, collapse = ", ")),
        call. = FALSE
      )
    }
    methods <- unique(correlations)
    roles <- stats::setNames(rep("user-requested configuration", length(methods)), methods)
    measurement <- "user_specified"
  }
  list(methods = methods, roles = roles, measurement = measurement)
}

#' Print a dimensionality-sensitivity analysis
#'
#' @param x Object returned by `dimensionality_sensitivity()`.
#' @param ... Additional unused arguments.
#' @return `x`, invisibly.
#' @method print psychomatic_dimensionality
#' @keywords internal
#' @export
print.psychomatic_dimensionality <- function(x, ...) {
  language <- if (is.null(x$language)) "eng" else x$language
  requested <- paste(unique(x$configurations$requested_method), collapse = ", ")
  effective <- paste(
    paste0(
      names(x$effective_settings$correlations),
      " -> ",
      unname(x$effective_settings$correlations)
    ),
    collapse = ", "
  )
  quantiles <- paste(format(x$effective_settings$quantiles, trim = TRUE), collapse = ", ")
  seeds <- paste(x$effective_settings$seed, collapse = ", ")
  if (language == "esp") {
    cat("Sensibilidad de dimensionalidad\n")
    cat("  Configuraciones solicitadas:", requested, "\n")
    cat("  Configuraciones efectivas:", effective, "\n")
    cat("  Cuantiles:", quantiles, "\n")
    cat("  Semillas:", seeds, "\n")
    cat("  n_iter:", x$effective_settings$n_iter, "\n")
    cat(
      "  Rango descriptivo:", x$solution_range[["minimum"]], "a",
      x$solution_range[["maximum"]], "factores\n"
    )
  } else {
    cat("Dimensionality sensitivity\n")
    cat("  Requested configurations:", requested, "\n")
    cat("  Effective configurations:", effective, "\n")
    cat("  Quantiles:", quantiles, "\n")
    cat("  Seeds:", seeds, "\n")
    cat("  n_iter:", x$effective_settings$n_iter, "\n")
    cat(
      "  Descriptive range:", x$solution_range[["minimum"]], "to",
      x$solution_range[["maximum"]], "factors\n"
    )
  }
  cat(" ", x$conclusion, "\n")
  invisible(x)
}

#' Summarize a dimensionality-sensitivity analysis
#'
#' @param object Object returned by `dimensionality_sensitivity()`.
#' @param ... Additional unused arguments.
#' @return A list containing the configuration table, range, agreement, and
#'   conclusion.
#' @method summary psychomatic_dimensionality
#' @keywords internal
#' @export
summary.psychomatic_dimensionality <- function(object, ...) {
  list(
    configurations = object$configurations,
    solution_range = data.frame(
      minimum = object$solution_range[["minimum"]],
      maximum = object$solution_range[["maximum"]]
    ),
    agreement = object$agreement,
    conclusion = object$conclusion
  )
}

#' Plot a dimensionality-sensitivity analysis
#'
#' @param x Object returned by `dimensionality_sensitivity()`.
#' @param ... Additional graphical parameters passed to `plot()`.
#' @return The plotted configuration rows, invisibly.
#' @method plot psychomatic_dimensionality
#' @keywords internal
#' @export
plot.psychomatic_dimensionality <- function(x, ...) {
  plotted <- x$configurations[x$configurations$status == "ok", , drop = FALSE]
  if (nrow(plotted) == 0L) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No evaluable configurations")
    return(invisible(plotted))
  }
  methods <- unique(plotted$requested_method)
  symbols <- match(plotted$requested_method, methods)
  graphics::plot(
    seq_len(nrow(plotted)),
    plotted$n_factors,
    pch = symbols,
    xaxt = "n",
    xlab = "Configuration",
    ylab = "Retained factors",
    ...
  )
  graphics::axis(1, at = seq_len(nrow(plotted)), labels = seq_len(nrow(plotted)))
  graphics::legend("topright", legend = methods, pch = seq_along(methods), bty = "n")
  invisible(plotted)
}
