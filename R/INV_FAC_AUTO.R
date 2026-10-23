.ifa_validate_packages <- function(measurement = "continuous",
                                   checker = requireNamespace) {
  if (!is.function(checker)) {
    stop("'checker' must be a function.", call. = FALSE)
  }
  has_package <- function(package) {
    isTRUE(checker(package, quietly = TRUE))
  }
  if (!has_package("lavaan")) {
    stop(
      "Package 'lavaan' is required. Install it with install.packages('lavaan').",
      call. = FALSE
    )
  }
  if (!identical(measurement, "continuous") && !has_package("semTools")) {
    stop(
      paste0(
        "Package 'semTools' is required for identified ordinal or mixed ",
        "measurement invariance. Install it with install.packages('semTools')."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.ifa_extract_input <- function(data) {
  specification <- if (inherits(data, "psychomatic_spec")) data else NULL
  if (!is.null(specification)) {
    raw <- specification$data
    spec_levels <- specification$measurement_level
    spec_language <- specification$language
  } else {
    raw <- data
    spec_levels <- NULL
    spec_language <- NULL
  }
  raw <- .pm_as_data_frame(raw)
  list(
    data = raw,
    specification = specification,
    measurement_levels = spec_levels,
    language = spec_language
  )
}

.ifa_unique_internal_name <- function(data_names,
                                      stem = ".psychomatic_group") {
  candidate <- stem
  suffix <- 1L
  while (candidate %in% data_names) {
    suffix <- suffix + 1L
    candidate <- paste0(stem, "_", suffix)
  }
  candidate
}

.ifa_resolve_group <- function(data, group_expr, environment) {
  data_names <- names(data)
  original_indices <- seq_len(nrow(data))
  group_value <- NULL
  group_name <- NULL
  source <- NULL

  if (is.symbol(group_expr)) {
    candidate <- as.character(group_expr)
    occurrences <- sum(!is.na(data_names) & data_names == candidate)
    if (occurrences > 1L) {
      stop("The grouping column name is duplicated; use a positive column index.", call. = FALSE)
    }
    if (occurrences == 1L) {
      group_name <- candidate
      group_value <- data[[match(candidate, data_names)]]
      source <- "column_symbol"
    }
  }

  if (is.null(group_value)) {
    evaluated <- tryCatch(
      eval(group_expr, envir = environment),
      error = function(e) e
    )
    if (inherits(evaluated, "error")) {
      stop(
        paste0("The grouping variable was not found in data or the calling environment: ",
               conditionMessage(evaluated)),
        call. = FALSE
      )
    }

    if (is.character(evaluated) && length(evaluated) == 1L) {
      occurrences <- sum(!is.na(data_names) & data_names == evaluated)
      if (occurrences > 1L) {
        stop("The grouping column name is duplicated; use a positive column index.", call. = FALSE)
      }
      if (occurrences == 1L) {
        group_name <- evaluated
        group_value <- data[[match(evaluated, data_names)]]
        source <- "column_name"
      }
    }

    if (is.null(group_value) && is.numeric(evaluated) &&
        length(evaluated) == 1L && !is.na(evaluated) &&
        is.finite(evaluated) && evaluated == trunc(evaluated)) {
      if (evaluated < 1L || evaluated > ncol(data)) {
        stop("A numeric grouping-column index must be within data.", call. = FALSE)
      }
      group_position <- as.integer(evaluated)
      selected_name <- data_names[group_position]
      group_value <- data[[group_position]]
      name_is_unique <- !is.na(selected_name) && nzchar(selected_name) &&
        sum(!is.na(data_names) & data_names == selected_name) == 1L
      group_name <- if (name_is_unique) {
        selected_name
      } else {
        .ifa_unique_internal_name(data_names)
      }
      source <- "column_index"
    }

    if (is.null(group_value)) {
      if (length(evaluated) != nrow(data)) {
        stop(
          "A direct grouping vector must have exactly one value per data row.",
          call. = FALSE
        )
      }
      group_name <- .ifa_unique_internal_name(data_names)
      group_value <- evaluated
      data[[group_name]] <- group_value
      source <- "direct_vector"
    }
  }

  if (is.data.frame(group_value) || is.matrix(group_value) ||
      is.list(group_value) || is.complex(group_value)) {
    stop("The grouping variable must be an atomic vector.", call. = FALSE)
  }
  if (all(is.na(group_value))) {
    stop("The grouping variable contains only missing values.", call. = FALSE)
  }

  keep <- !is.na(group_value)
  dropped <- sum(!keep)
  data <- data[keep, , drop = FALSE]
  original_indices <- original_indices[keep]
  group_value <- group_value[keep]

  if (is.factor(group_value)) {
    group_factor <- droplevels(group_value)
  } else {
    group_labels <- unique(as.character(group_value))
    group_factor <- factor(as.character(group_value), levels = group_labels)
  }
  if (nlevels(group_factor) < 2L) {
    stop("The grouping variable must contain at least two groups with nonempty membership.",
         call. = FALSE)
  }
  data[[group_name]] <- group_factor

  list(
    data = data,
    group_name = group_name,
    group_source = source,
    original_indices = as.integer(original_indices),
    dropped_missing_group = as.integer(dropped),
    group_coding = levels(group_factor)
  )
}

.ifa_extract_indicators <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model) ||
      !nzchar(trimws(model))) {
    stop("'model' must be one non-empty lavaan model string.", call. = FALSE)
  }
  parsed <- tryCatch(
    lavaan::lavaanify(
      model = model,
      auto = TRUE,
      fixed.x = FALSE,
      model.type = "cfa"
    ),
    error = function(e) {
      stop(
        paste0("The lavaan model could not be parsed: ", conditionMessage(e)),
        call. = FALSE
      )
    }
  )
  loadings <- parsed[parsed$op == "=~", , drop = FALSE]
  if (nrow(loadings) == 0L) {
    stop("The model must contain at least one loading statement using '=~'.", call. = FALSE)
  }
  latent <- unique(loadings$lhs)
  indicators <- unique(loadings$rhs[!loadings$rhs %in% latent])
  if (length(indicators) == 0L) {
    stop("No observed indicators were detected in the model.", call. = FALSE)
  }
  indicators
}

.ifa_observed_values <- function(x) {
  x[!is.na(x)]
}

.ifa_numeric_integerish <- function(x, tolerance = 1e-8) {
  observed <- .ifa_observed_values(x)
  is.numeric(observed) && length(observed) > 0L &&
    all(is.finite(observed)) &&
    all(abs(observed - round(observed)) < tolerance)
}

.ifa_base_item_level <- function(x) {
  observed <- .ifa_observed_values(x)
  if (length(observed) == 0L) return("invalid")
  if (is.numeric(x) && any(is.nan(x))) return("invalid")
  if (is.numeric(x) && any(!is.finite(observed))) return("invalid")
  if (length(unique(observed)) < 2L) return("invalid")
  if (is.ordered(x) || is.logical(x)) return("ordinal")
  if (is.factor(x) || is.character(x)) return("ambiguous")
  if (!is.numeric(x) || inherits(x, c("Date", "POSIXt"))) return("ineligible")
  if (length(unique(observed)) <= 7L && .ifa_numeric_integerish(x)) {
    return("ordinal")
  }
  "continuous"
}

.ifa_resolve_ordered_selection <- function(ordered, data, indicators) {
  if (is.null(ordered)) return(NULL)
  if (is.logical(ordered)) {
    if (length(ordered) != 1L || is.na(ordered)) {
      stop("'ordered' must be NULL, one TRUE/FALSE value, or item names/indices.", call. = FALSE)
    }
    return(if (isTRUE(ordered)) indicators else character())
  }
  if (is.character(ordered)) {
    if (!length(ordered) || anyNA(ordered) || any(!nzchar(ordered)) ||
        anyDuplicated(ordered)) {
      stop("'ordered' item names must be unique and non-empty.", call. = FALSE)
    }
    unknown <- setdiff(ordered, indicators)
    if (length(unknown)) {
      stop(
        paste0("Unknown ordered indicators: ", paste(unknown, collapse = ", "), "."),
        call. = FALSE
      )
    }
    return(ordered)
  }
  if (is.numeric(ordered)) {
    positions <- .pm_resolve_columns(ordered, names(data), "ordered")
    selected <- names(data)[positions]
    unknown <- setdiff(selected, indicators)
    if (length(unknown)) {
      stop(
        paste0("Ordered indices must identify model indicators, not: ",
               paste(unknown, collapse = ", "), "."),
        call. = FALSE
      )
    }
    return(selected)
  }
  stop("'ordered' must be NULL, one TRUE/FALSE value, or item names/indices.", call. = FALSE)
}

.ifa_category_profile <- function(measurement, ordered, category_counts) {
  if (measurement == "continuous") return("continuous")
  if (measurement == "mixed") return("mixed")
  if (length(ordered) && all(category_counts[ordered] == 2L)) {
    return("binary")
  }
  "polytomous"
}

.ifa_prepare_items <- function(data,
                               indicators,
                               ordered = NULL,
                               specification_levels = NULL) {
  data_names <- names(data)
  missing_items <- setdiff(indicators, data_names)
  if (length(missing_items)) {
    stop(
      paste0("Indicators not found in data: ", paste(missing_items, collapse = ", "), "."),
      call. = FALSE
    )
  }
  duplicated_items <- indicators[
    vapply(indicators, function(x) sum(!is.na(data_names) & data_names == x), integer(1)) > 1L
  ]
  if (length(duplicated_items)) {
    stop("Model indicator names are duplicated in data.", call. = FALSE)
  }

  raw_levels <- vapply(data[indicators], .ifa_base_item_level, character(1))
  invalid <- indicators[raw_levels %in% c("invalid", "ineligible")]
  if (length(invalid)) {
    stop(
      paste0(
        "Indicators must contain at least two observed finite values and use ",
        "eligible encodings: ", paste(invalid, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  ambiguous <- indicators[raw_levels == "ambiguous"]
  if (length(ambiguous)) {
    stop(
      paste0(
        "Character or unordered factor indicators require explicit levels in ",
        "psychomatic_spec; values are never alphabetically ordered: ",
        paste(ambiguous, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  selected <- .ifa_resolve_ordered_selection(ordered, data, indicators)
  if (is.null(selected)) {
    if (!is.null(specification_levels)) {
      declared <- specification_levels[indicators]
      if (length(declared) == length(indicators) && !anyNA(declared)) {
        selected <- indicators[declared == "ordinal"]
      } else {
        selected <- indicators[raw_levels == "ordinal"]
      }
    } else {
      selected <- indicators[raw_levels == "ordinal"]
    }
  }

  semantic_categories <- indicators[vapply(
    data[indicators],
    function(x) is.ordered(x) || is.factor(x) || is.logical(x),
    logical(1)
  )]
  if (!identical(ordered, FALSE)) {
    inherently_ordered <- indicators[vapply(
      data[indicators],
      function(x) is.ordered(x) || is.logical(x),
      logical(1)
    )]
    selected <- unique(c(selected, inherently_ordered))
  }
  silently_continuous <- setdiff(semantic_categories, selected)
  if (length(silently_continuous)) {
    stop(
      paste0(
        "Semantically categorical indicators cannot be treated as continuous numeric codes",
        if (identical(ordered, FALSE)) " with ordered = FALSE" else "",
        ": ", paste(silently_continuous, collapse = ", "),
        ". Declare them as ordered indicators or provide explicit ordinal levels."
      ),
      call. = FALSE
    )
  }

  prepared <- data
  for (item in indicators) {
    x <- prepared[[item]]
    if (item %in% selected) {
      if (is.ordered(x)) {
        prepared[[item]] <- droplevels(x)
      } else if (is.logical(x)) {
        prepared[[item]] <- ordered(x, levels = c(FALSE, TRUE))
      } else {
        observed_levels <- sort(unique(x[!is.na(x)]))
        prepared[[item]] <- ordered(x, levels = observed_levels)
      }
    } else if (is.ordered(x) || is.factor(x)) {
      prepared[[item]] <- as.numeric(x)
    } else if (is.logical(x)) {
      prepared[[item]] <- as.numeric(x)
    } else {
      prepared[[item]] <- as.numeric(x)
    }
  }

  category_counts <- vapply(
    prepared[indicators],
    function(x) length(unique(x[!is.na(x)])),
    integer(1)
  )
  measurement <- if (!length(selected)) {
    "continuous"
  } else if (length(selected) == length(indicators)) {
    "ordinal"
  } else {
    "mixed"
  }
  category_profile <- .ifa_category_profile(
    measurement, selected, category_counts
  )

  list(
    data = prepared,
    measurement = measurement,
    ordered = selected,
    item_levels = stats::setNames(
      ifelse(indicators %in% selected, "ordinal", "continuous"),
      indicators
    ),
    category_counts = category_counts,
    category_profile = category_profile
  )
}

.ifa_refresh_analytic_items <- function(items, data, indicators) {
  if (!is.list(items) || !is.data.frame(data) ||
      !all(indicators %in% names(data))) {
    stop("Prepared item metadata and analytic data are required.", call. = FALSE)
  }
  refreshed <- data
  ordered_items <- intersect(items$ordered, indicators)
  for (item in ordered_items) {
    if (is.factor(refreshed[[item]])) {
      refreshed[[item]] <- droplevels(refreshed[[item]])
    }
  }
  category_counts <- vapply(
    refreshed[indicators],
    function(x) length(unique(x[!is.na(x)])),
    integer(1)
  )
  items$data <- refreshed
  items$category_counts <- category_counts
  items$category_profile <- .ifa_category_profile(
    items$measurement, ordered_items, category_counts
  )
  items
}

.ifa_resolve_missing <- function(missing, measurement) {
  if (!is.character(missing) || length(missing) != 1L || is.na(missing) ||
      !nzchar(trimws(missing))) {
    stop("'missing' must be 'auto', 'fiml', or 'listwise'.", call. = FALSE)
  }
  value <- tolower(trimws(missing))
  if (value == "auto") {
    value <- if (measurement == "continuous") "fiml" else "listwise"
  }
  if (!value %in% c("fiml", "listwise")) {
    stop("'missing' must be 'auto', 'fiml', or 'listwise'.", call. = FALSE)
  }
  if (measurement != "continuous" && value == "fiml") {
    stop("FIML is incompatible with categorical WLSMV/ULSMV; use listwise missing handling.", call. = FALSE)
  }
  value
}

.ifa_filter_analysis_rows <- function(data,
                                      indicators,
                                      group_name,
                                      original_indices,
                                      missing) {
  observed <- !is.na(data[indicators])
  keep <- if (missing == "listwise") {
    rowSums(observed) == length(indicators)
  } else {
    rowSums(observed) > 0L
  }
  filtered_data <- data[keep, , drop = FALSE]
  if (!group_name %in% names(filtered_data)) {
    stop("The frozen grouping column is unavailable during row filtering.", call. = FALSE)
  }
  if (is.factor(filtered_data[[group_name]])) {
    filtered_data[[group_name]] <- droplevels(filtered_data[[group_name]])
  } else {
    labels <- unique(as.character(filtered_data[[group_name]]))
    filtered_data[[group_name]] <- factor(
      as.character(filtered_data[[group_name]]),
      levels = labels
    )
  }
  list(
    data = filtered_data,
    original_indices = as.integer(original_indices[keep]),
    dropped_missing_items = as.integer(sum(!keep)),
    group_coding = levels(filtered_data[[group_name]])
  )
}

.ifa_validate_analysis_data <- function(data,
                                        indicators,
                                        group_name,
                                        min_group_n) {
  min_group_n <- .pm_validate_scalar_number(
    min_group_n,
    name = "min_group_n",
    lower = 2,
    integer = TRUE
  )
  group_sizes <- table(data[[group_name]], useNA = "no")
  group_sizes <- group_sizes[group_sizes > 0L]
  if (length(group_sizes) < 2L) {
    stop("The analytic sample must contain at least two groups with nonempty membership.",
         call. = FALSE)
  }
  if (any(group_sizes < min_group_n)) {
    counts <- paste(names(group_sizes), as.integer(group_sizes), sep = "=", collapse = ", ")
    stop(
      paste0(
        "Every analytic group must satisfy min_group_n=", as.integer(min_group_n),
        "; observed counts: ", counts, "."
      ),
      call. = FALSE
    )
  }
  if (!nrow(data)) {
    stop("No analytic cases remain after missing-data handling.", call. = FALSE)
  }

  for (item in indicators) {
    x <- data[[item]]
    observed <- x[!is.na(x)]
    if (is.numeric(observed) && any(!is.finite(observed))) {
      stop(paste0("Indicator '", item, "' contains Inf or another nonfinite value."), call. = FALSE)
    }
    if (length(unique(observed)) < 2L) {
      stop(paste0("Indicator '", item, "' is constant or has fewer than two observed values."), call. = FALSE)
    }
    by_group <- split(x, data[[group_name]], drop = TRUE)
    bad <- names(by_group)[vapply(
      by_group,
      function(values) length(unique(values[!is.na(values)])) < 2L,
      logical(1)
    )]
    if (length(bad)) {
      stop(
        paste0(
          "Indicator '", item, "' is group-degenerate in: ",
          paste(bad, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  }
  list(
    group_sizes = group_sizes,
    min_group_n = as.integer(min(group_sizes)),
    total_n = as.integer(nrow(data))
  )
}

.ifa_normality_from_groups <- function(details) {
  required <- c("skew_p", "kurtosis_p", "status")
  if (!is.data.frame(details) || !all(required %in% names(details)) ||
      nrow(details) == 0L) {
    return(list(normal = NA, reason = "missing_group_evidence"))
  }
  evaluable <- details$status == "ok" &
    is.finite(details$skew_p) & is.finite(details$kurtosis_p)
  if (!all(evaluable)) {
    return(list(normal = NA, reason = "non_evaluable_group"))
  }
  normal <- all(details$skew_p >= .05 & details$kurtosis_p >= .05)
  list(
    normal = normal,
    reason = if (normal) "all_group_p_values_at_least_.05" else "at_least_one_group_p_value_below_.05"
  )
}

.ifa_group_normality <- function(data, indicators, group_name) {
  groups <- levels(data[[group_name]])
  rows <- vector("list", length(groups))
  conditions <- list()
  for (i in seq_along(groups)) {
    group_rows <- data[[group_name]] == groups[i]
    complete <- stats::complete.cases(data[group_rows, indicators, drop = FALSE])
    x <- data[group_rows, indicators, drop = FALSE][complete, , drop = FALSE]
    captured <- .pm_capture_conditions(
      .pm_mardia_blockwise(as.matrix(x))
    )
    if (!is.null(captured$error)) {
      conditions[[length(conditions) + 1L]] <- captured$error
      rows[[i]] <- data.frame(
        group = groups[i],
        n_complete = nrow(x),
        skew_p = NA_real_,
        kurtosis_p = NA_real_,
        status = "error",
        message = conditionMessage(captured$error),
        stringsAsFactors = FALSE
      )
    } else {
      statistic <- captured$value
      rows[[i]] <- data.frame(
        group = groups[i],
        n_complete = nrow(x),
        skew_p = unname(statistic[["skew_p"]]),
        kurtosis_p = unname(statistic[["kurtosis_p"]]),
        status = "ok",
        message = "",
        stringsAsFactors = FALSE
      )
    }
    conditions <- c(conditions, captured$warnings, captured$messages)
  }
  details <- do.call(rbind, rows)
  rownames(details) <- NULL
  overall <- .ifa_normality_from_groups(details)
  list(
    overall_normal = overall$normal,
    overall_test = list(
      normal = overall$normal,
      skew_p = NA_real_,
      kurt_p = NA_real_,
      status = overall$reason,
      n_complete = sum(details$n_complete)
    ),
    by_group = details,
    conditions = conditions
  )
}

.ifa_resolve_criteria <- function(total_n,
                                  group_ratio,
                                  normal,
                                  measurement,
                                  criteria = "auto") {
  total_n <- .pm_validate_scalar_number(total_n, "total_n", lower = 2, integer = TRUE)
  group_ratio <- .pm_validate_scalar_number(group_ratio, "group_ratio", lower = 1)
  if (!is.logical(normal) || length(normal) != 1L) {
    stop("'normal' must be TRUE, FALSE, or NA.", call. = FALSE)
  }
  if (!measurement %in% c("continuous", "ordinal", "mixed")) {
    stop("'measurement' must be continuous, ordinal, or mixed.", call. = FALSE)
  }

  defaults <- c(
    configural_cfi = .90,
    configural_rmsea = .08,
    configural_srmr = .08
  )
  if (is.character(criteria) && length(criteria) == 1L &&
      !is.na(criteria) && tolower(trimws(criteria)) == "auto") {
    # Noninvariance pattern is unknown before testing. Do not pretend that
    # N, group balance or Mardia identify a Chen simulation scenario.
    cutoffs <- c(cfi = .010, rmsea = .015, srmr_loading = .030, srmr_other = .010)
    preset <- if (measurement != "continuous") {
      "categorical_descriptive_extrapolation"
    } else {
      "chen_informed_guidelines"
    }
    resolved <- as.list(c(cutoffs, defaults))
    resolved$preset <- preset
    return(resolved)
  }

  if (is.list(criteria)) {
    if (is.null(names(criteria))) {
      stop("Custom criteria must be fully named.", call. = FALSE)
    }
    criteria <- unlist(criteria, recursive = FALSE, use.names = TRUE)
  }
  if (!is.numeric(criteria) || is.null(names(criteria)) || anyNA(names(criteria)) ||
      any(!nzchar(names(criteria))) || anyDuplicated(names(criteria))) {
    stop("Custom criteria must be a fully named numeric vector with unique names.", call. = FALSE)
  }
  allowed <- c(
    "cfi", "rmsea", "srmr_loading", "srmr_other",
    "configural_cfi", "configural_rmsea", "configural_srmr"
  )
  unknown <- setdiff(names(criteria), allowed)
  if (length(unknown)) {
    stop(paste0("Custom criteria contain unknown names: ", paste(unknown, collapse = ", "), "."), call. = FALSE)
  }
  required <- c("cfi", "rmsea", "srmr_loading", "srmr_other")
  missing_required <- setdiff(required, names(criteria))
  if (length(missing_required)) {
    stop(
      paste0("Custom criteria must include: ", paste(required, collapse = ", "),
             "; missing ", paste(missing_required, collapse = ", "), "."),
      call. = FALSE
    )
  }
  if (anyNA(criteria) || any(!is.finite(criteria)) || any(criteria < 0)) {
    stop("Custom criteria must be finite, scalar per name, and nonnegative.", call. = FALSE)
  }
  completed <- c(criteria[required], defaults)
  optional <- intersect(names(criteria), names(defaults))
  completed[optional] <- criteria[optional]
  resolved <- as.list(completed)
  resolved$preset <- "custom"
  resolved
}

.ifa_resolve_estimator <- function(estimator,
                                   measurement,
                                   normal,
                                   ordered_argument = NULL,
                                   allow_numeric_wls = FALSE) {
  if (!is.character(estimator) || length(estimator) != 1L || is.na(estimator) ||
      !nzchar(trimws(estimator))) {
    stop("'estimator' must be 'auto', 'ML', 'MLR', 'ULS', 'WLSMV', or 'ULSMV'.",
         call. = FALSE)
  }
  requested <- trimws(estimator)
  upper <- toupper(requested)
  if (tolower(requested) == "auto") {
    effective <- if (measurement == "continuous") {
      if (isTRUE(normal)) "ML" else "MLR"
    } else {
      "WLSMV"
    }
    reason <- if (measurement != "continuous") {
      "categorical_measurement"
    } else if (isTRUE(normal)) {
      "all_groups_passed_mardia"
    } else {
      "nonnormal_or_unevaluable_mardia"
    }
  } else {
    allowed <- c("ML", "MLR", "ULS", "WLSMV", "ULSMV")
    if (!upper %in% allowed) {
      stop("'estimator' must be 'auto', 'ML', 'MLR', 'ULS', 'WLSMV', or 'ULSMV'.",
           call. = FALSE)
    }
    effective <- upper
    reason <- "manual_override"
  }

  categorical_estimators <- c("WLSMV", "ULSMV")
  if (measurement == "continuous" && effective %in% categorical_estimators &&
      !isTRUE(allow_numeric_wls)) {
    detail <- if (identical(ordered_argument, FALSE)) " with ordered = FALSE" else ""
    stop(
      paste0("Estimator ", effective, " is incompatible with continuous indicators", detail,
             "; declare ordered indicators or use ML/MLR/ULS."),
      call. = FALSE
    )
  }
  if (measurement != "continuous" && effective %in% c("ML", "MLR", "ULS")) {
    stop(
      paste0("Estimator ", effective,
             " is incompatible with ordered indicators; use WLSMV or ULSMV."),
      call. = FALSE
    )
  }
  list(
    requested = requested,
    estimator = effective,
    reason = reason,
    uses_ordered = measurement != "continuous",
    selection_mode = if (tolower(requested) == "auto") "auto" else "manual",
    handling = if (measurement == "continuous") "continuous" else "ordered"
  )
}

.ifa_stage_sequence <- function(measurement, category_profile) {
  if (measurement == "continuous") {
    return(list(
      configural = list(detail = "configural", constraints = character()),
      metric = list(detail = "loadings", constraints = "loadings"),
      scalar = list(
        detail = "loadings_intercepts",
        constraints = c("loadings", "intercepts")
      ),
      strict = list(
        detail = "loadings_intercepts_residuals",
        constraints = c("loadings", "intercepts", "residuals")
      )
    ))
  }
  if (identical(category_profile, "binary")) {
    return(list(
      configural = list(detail = "configural", constraints = character()),
      scalar = list(
        detail = "binary_thresholds_loadings_intercepts_joint",
        constraints = c("thresholds", "loadings", "intercepts")
      ),
      strict = list(
        detail = "binary_joint_plus_residuals",
        constraints = c("thresholds", "loadings", "intercepts", "residuals")
      )
    ))
  }
  list(
    configural = list(detail = "configural", constraints = character()),
    thresholds = list(detail = "thresholds", constraints = "thresholds"),
    metric = list(
      detail = "thresholds_loadings",
      constraints = c("thresholds", "loadings")
    ),
    scalar = list(
      detail = "thresholds_loadings_intercepts",
      constraints = c("thresholds", "loadings", "intercepts")
    ),
    strict = list(
      detail = "thresholds_loadings_intercepts_residuals",
      constraints = c("thresholds", "loadings", "intercepts", "residuals")
    )
  )
}

.ifa_normalize_constraints <- function(constraints) {
  if (is.null(constraints) || length(constraints) == 0L ||
      (length(constraints) == 1L && identical(constraints, ""))) {
    return(character())
  }
  if (!is.character(constraints) || anyNA(constraints) ||
      any(!nzchar(constraints)) || anyDuplicated(constraints)) {
    stop("Constraints must be unique non-empty character values.", call. = FALSE)
  }
  allowed <- c("thresholds", "loadings", "intercepts", "residuals")
  unknown <- setdiff(constraints, allowed)
  if (length(unknown)) {
    stop(paste0("Unknown invariance constraints: ", paste(unknown, collapse = ", "), "."),
         call. = FALSE)
  }
  constraints
}

.ifa_deteriorations <- function(previous, current) {
  previous <- as.numeric(previous[c("cfi", "rmsea", "srmr")])
  current <- as.numeric(current[c("cfi", "rmsea", "srmr")])
  names(previous) <- names(current) <- c("cfi", "rmsea", "srmr")
  c(
    delta_cfi = current[["cfi"]] - previous[["cfi"]],
    delta_rmsea = current[["rmsea"]] - previous[["rmsea"]],
    delta_srmr = current[["srmr"]] - previous[["srmr"]]
  )
}

.ifa_change_decision <- function(delta_cfi,
                                 delta_rmsea,
                                 delta_srmr,
                                 cutoffs,
                                 level) {
  values <- c(delta_cfi, delta_rmsea, delta_srmr)
  if (length(values) != 3L || any(!is.finite(values))) {
    return(list(
      decision = "not_evaluable",
      status = "not_evaluable",
      reason = "finite_fit_change_indices_are_required"
    ))
  }
  srmr_cutoff <- if (level %in% c("metric", "loadings", "thresholds_loadings")) {
    cutoffs$srmr_loading
  } else {
    cutoffs$srmr_other
  }
  failed <- delta_cfi < -cutoffs$cfi &&
    (delta_rmsea > cutoffs$rmsea || delta_srmr > srmr_cutoff)
  list(
    decision = if (failed) "not_supported" else "supported",
    status = "ok",
    reason = if (failed) {
      "fit_deterioration_exceeded_the_joint_heuristic"
    } else {
      "fit_deterioration_remained_within_the_joint_heuristic"
    }
  )
}

.ifa_configural_gate <- function(indices,
                                 criteria = list(
                                   configural_cfi = .90,
                                   configural_rmsea = .08,
                                   configural_srmr = .08
                                 )) {
  values <- as.numeric(indices[c("cfi", "rmsea", "srmr")])
  if (length(values) != 3L || any(!is.finite(values))) {
    return(list(
      passed = FALSE,
      decision = "not_evaluable",
      status = "not_evaluable",
      reason = "configural_indices_missing_or_nonfinite"
    ))
  }
  passed <- values[1L] >= criteria$configural_cfi &&
    values[2L] <= criteria$configural_rmsea &&
    values[3L] <= criteria$configural_srmr
  list(
    passed = passed,
    decision = if (passed) "supported" else "not_supported",
    status = if (passed) "ok" else "failed_configural_fit",
    reason = if (passed) {
      "configural_fit_met_all_three_descriptive_gates"
    } else {
      "configural_fit_failed_at_least_one_descriptive_gate"
    }
  )
}

.ifa_record_priority <- function(record, index) {
  if (index == "srmr") return("base")
  if (index %in% c("cfi", "rmsea") && identical(record$fit_variant, "robust")) {
    return("robust")
  }
  corrected <- any(toupper(c(record$estimator, record$estimator_orig)) %in%
    c("MLR", "MLM", "MLMV", "MLMVS", "WLSMV", "ULSMV"))
  if (corrected) "scaled" else "base"
}

.ifa_select_variant <- function(x, priority = c("scaled", "base")) {
  if (is.null(x) || !length(x)) return(list(variant = NA_character_, value = NA_real_))
  original_names <- names(x)
  x <- as.numeric(x)
  names(x) <- if (is.null(original_names)) rep("", length(x)) else original_names
  for (variant in priority) {
    position <- match(variant, names(x))
    if (!is.na(position) && is.finite(x[position])) {
      return(list(variant = variant, value = unname(x[position])))
    }
  }
  list(variant = NA_character_, value = NA_real_)
}

.ifa_common_variant <- function(previous,
                                current,
                                priority = c("scaled", "base")) {
  previous_names <- names(previous)
  current_names <- names(current)
  for (variant in priority) {
    p <- match(variant, previous_names)
    q <- match(variant, current_names)
    if (!is.na(p) && !is.na(q) && is.finite(previous[p]) && is.finite(current[q])) {
      return(list(
        variant = variant,
        previous = unname(previous[p]),
        current = unname(current[q])
      ))
    }
  }
  NULL
}

.ifa_admissibility <- function(record) {
  reasons <- character()
  if (!isTRUE(record$converged)) reasons <- c(reasons, "converged")
  if (!isTRUE(record$post_check)) reasons <- c(reasons, "post_check")
  if (!isTRUE(record$finite_parameters)) reasons <- c(reasons, "finite_parameters")
  if (isTRUE(record$negative_observed_variance)) {
    reasons <- c(reasons, "negative_observed_variance")
  }
  if (isTRUE(record$negative_latent_variance)) {
    reasons <- c(reasons, "negative_latent_variance")
  }
  for (index in c("cfi", "rmsea", "srmr", "df")) {
    selected <- .ifa_select_variant(record$measures[[index]], .ifa_record_priority(record, index))
    if (!is.finite(selected$value)) reasons <- c(reasons, paste0("nonfinite_", index))
  }
  list(
    admissible = length(reasons) == 0L,
    reason = if (length(reasons)) paste(unique(reasons), collapse = ";") else "admissible"
  )
}

.ifa_comparison_metadata <- function(record) {
  if (is.list(record$fit_metadata)) return(record$fit_metadata)
  list(
    case_indices = record$case_indices,
    variables = record$variables,
    estimator_effective = record$estimator,
    estimator_orig = record$estimator,
    estimand = record$estimand,
    ordered = record$ordered,
    missing = record$missing,
    identification = record$identification,
    parameterization = record$parameterization,
    std_lv = identical(record$identification, "std.lv"),
    group_coding = record$group_coding,
    measurement = record$measurement
  )
}

.ifa_compare_records <- function(previous, current) {
  if (!isTRUE(previous$admissible)) {
    return(list(comparable = FALSE, reason = "previous_admissible"))
  }
  if (!isTRUE(current$admissible)) {
    return(list(comparable = FALSE, reason = "current_admissible"))
  }
  previous_metadata <- .ifa_comparison_metadata(previous)
  current_metadata <- .ifa_comparison_metadata(current)
  frozen <- c(
    "case_indices", "variables", "estimator_effective", "estimator_orig",
    "estimand", "ordered", "missing", "identification", "parameterization",
    "std_lv", "group_coding", "measurement"
  )
  for (field in frozen) {
    if (!identical(previous_metadata[[field]], current_metadata[[field]])) {
      return(list(comparable = FALSE, reason = paste0("mismatch_", field)))
    }
  }
  old_constraints <- .ifa_normalize_constraints(previous$effective_constraints)
  new_constraints <- .ifa_normalize_constraints(current$effective_constraints)
  if (!all(old_constraints %in% new_constraints) ||
      length(new_constraints) <= length(old_constraints)) {
    return(list(comparable = FALSE, reason = "constraints_not_nested_cumulative"))
  }

  common <- lapply(
    c("cfi", "rmsea", "srmr", "df"),
    function(index) .ifa_common_variant(
      previous$measures[[index]], current$measures[[index]],
      intersect(.ifa_record_priority(previous, index), .ifa_record_priority(current, index))
    )
  )
  names(common) <- c("cfi", "rmsea", "srmr", "df")
  missing_variant <- names(common)[vapply(common, is.null, logical(1))]
  if (length(missing_variant)) {
    return(list(
      comparable = FALSE,
      reason = paste0("no_common_variant_", paste(missing_variant, collapse = "_"))
    ))
  }
  delta_df <- common$df$current - common$df$previous
  if (!is.finite(delta_df) || delta_df <= 0) {
    return(list(comparable = FALSE, reason = "delta_df_not_positive", delta_df = delta_df))
  }
  deltas <- .ifa_deteriorations(
    previous = c(
      cfi = common$cfi$previous,
      rmsea = common$rmsea$previous,
      srmr = common$srmr$previous
    ),
    current = c(
      cfi = common$cfi$current,
      rmsea = common$rmsea$current,
      srmr = common$srmr$current
    )
  )
  list(
    comparable = TRUE,
    reason = "comparable",
    delta_df = unname(delta_df),
    delta_cfi = unname(deltas[["delta_cfi"]]),
    delta_rmsea = unname(deltas[["delta_rmsea"]]),
    delta_srmr = unname(deltas[["delta_srmr"]]),
    variants = vapply(common, `[[`, character(1), "variant"),
    current_values = vapply(common, `[[`, numeric(1), "current")
  )
}

.ifa_normalize_group_partial <- function(group_partial) {
  if (is.null(group_partial) || length(group_partial) == 0L ||
      (length(group_partial) == 1L && identical(group_partial, ""))) {
    return("")
  }
  if (!is.character(group_partial) || anyNA(group_partial) ||
      any(!nzchar(trimws(group_partial))) || anyDuplicated(group_partial)) {
    stop("'group_partial' must contain unique non-empty lavaan parameter expressions.",
         call. = FALSE)
  }
  trimws(group_partial)
}

.ifa_build_syntax <- function(model,
                              data,
                              group_name,
                              measurement,
                              ordered,
                              estimator,
                              constraints = character(),
                              group_partial = "",
                              missing = "listwise",
                              parameterization = "theta",
                              ID.fac = "std.lv",
                              ID.cat = "Wu.Estabrook.2016",
                              native = FALSE) {
  constraints <- .ifa_normalize_constraints(constraints)
  group_partial <- .ifa_normalize_group_partial(group_partial)
  if (!is.data.frame(data) || !group_name %in% names(data)) {
    stop("Frozen analysis data and its grouping column are required.", call. = FALSE)
  }
  engine_requests <- list(
    estimator = as.character(estimator[[1L]]),
    missing = as.character(missing[[1L]]),
    parameterization = if (identical(measurement, "continuous")) {
      NA_character_
    } else {
      as.character(parameterization[[1L]])
    },
    ID.fac = as.character(ID.fac[[1L]]),
    ID.cat = if (identical(measurement, "continuous")) {
      NA_character_
    } else {
      as.character(ID.cat[[1L]])
    }
  )

  if (measurement == "continuous" || isTRUE(native)) {
    partial_for_lavaan <- if (identical(group_partial, "")) character() else group_partial
    engine_args <- list(
      model = model,
      data = data,
      group = group_name,
      estimator = estimator,
      std.lv = identical(ID.fac, "std.lv"),
      meanstructure = TRUE,
      missing = missing,
      ordered = if (measurement == "continuous") FALSE else if (isTRUE(native) &&
        setequal(ordered, .ifa_extract_indicators(model))) TRUE else ordered,
      group.equal = constraints,
      group.partial = partial_for_lavaan
    )
    captured <- .pm_capture_conditions(do.call(lavaan::cfa, engine_args))
    syntax <- model
    engine <- "lavaan::cfa"
  } else {
    engine_args <- list(
      configural.model = model,
      data = data,
      group = group_name,
      ordered = ordered,
      estimator = estimator,
      parameterization = parameterization,
      ID.fac = ID.fac,
      ID.cat = ID.cat,
      group.equal = if (length(constraints)) constraints else "",
      group.partial = group_partial,
      return.fit = TRUE
    )
    captured <- .pm_capture_conditions(do.call(semTools::measEq.syntax, engine_args))
    syntax <- if (!is.null(captured$value) && inherits(captured$value, "lavaan")) {
      as.character(stats::getCall(captured$value)$model)
    } else {
      model
    }
    engine <- "semTools::measEq.syntax"
  }

  fit <- captured$value
  parameter_table <- if (!is.null(fit) && inherits(fit, "lavaan")) {
    lavaan::parameterTable(fit)
  } else {
    data.frame()
  }
  list(
    syntax = syntax,
    fit = fit,
    engine = engine,
    engine_args = engine_args,
    engine_requests = engine_requests,
    requested_constraints = constraints,
    effective_constraints = constraints,
    group_partial = group_partial,
    parameter_table = parameter_table,
    warnings = captured$warnings,
    messages = captured$messages,
    error = captured$error
  )
}

.ifa_measure_vector <- function(all_measures, robust, scaled, base) {
  get_one <- function(name) {
    if (!name %in% names(all_measures)) return(NA_real_)
    value <- suppressWarnings(as.numeric(all_measures[[name]]))
    if (length(value) != 1L) NA_real_ else value
  }
  c(
    robust = get_one(robust),
    scaled = get_one(scaled),
    base = get_one(base)
  )
}

.ifa_extract_measures <- function(fit) {
  if (is.null(fit) || !inherits(fit, "lavaan")) {
    missing <- c(robust = NA_real_, scaled = NA_real_, base = NA_real_)
    out <- list(
      chi2 = missing, df = missing, p = missing,
      cfi = missing, rmsea = missing, srmr = missing
    )
    attr(out, "conditions") <- list()
    attr(out, "error") <- NULL
    return(out)
  }
  captured <- .pm_capture_conditions(lavaan::fitMeasures(fit))
  values <- if (is.null(captured$error)) captured$value else numeric()
  out <- list(
    chi2 = .ifa_measure_vector(values, "chisq.robust", "chisq.scaled", "chisq"),
    df = .ifa_measure_vector(values, "df.robust", "df.scaled", "df"),
    p = .ifa_measure_vector(values, "pvalue.robust", "pvalue.scaled", "pvalue"),
    cfi = .ifa_measure_vector(values, "cfi.robust", "cfi.scaled", "cfi"),
    rmsea = .ifa_measure_vector(values, "rmsea.robust", "rmsea.scaled", "rmsea"),
    srmr = c(
      robust = NA_real_,
      scaled = NA_real_,
      base = if ("srmr" %in% names(values)) as.numeric(values[["srmr"]]) else NA_real_
    )
  )
  attr(out, "conditions") <- c(captured$warnings, captured$messages)
  attr(out, "error") <- captured$error
  out
}

.ifa_effective_lavaan_settings <- function(fit) {
  if (is.null(fit) || !inherits(fit, "lavaan")) {
    return(list(
      estimator = NA_character_, estimator_orig = NA_character_,
      missing = NA_character_, parameterization = NA_character_, std.lv = NA
    ))
  }
  options <- lavaan::lavInspect(fit, "options")
  scalar_character <- function(value) {
    if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1L]])
  }
  list(
    estimator = scalar_character(options$estimator),
    estimator_orig = scalar_character(options$estimator.orig),
    missing = .ifa_canonical_missing(options$missing),
    parameterization = scalar_character(options$parameterization),
    std.lv = isTRUE(options$std.lv)
  )
}

.ifa_canonical_missing <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[[1L]])) {
    return(NA_character_)
  }
  value <- tolower(as.character(value[[1L]]))
  if (value %in% c("ml", "ml.x", "fiml")) return("fiml")
  if (value %in% c("listwise", "default")) return("listwise")
  value
}

.ifa_fit_metadata <- function(fit, blueprint) {
  if (is.null(fit) || !inherits(fit, "lavaan")) return(NULL)

  options <- lavaan::lavInspect(fit, "options")
  variables <- as.character(lavaan::lavNames(fit, type = "ov"))
  ordered <- as.character(lavaan::lavNames(fit, type = "ov.ord"))
  group_coding <- as.character(lavaan::lavInspect(fit, "group.label"))

  case_groups <- lavaan::lavInspect(fit, "case.idx")
  fit_data <- lavaan::lavInspect(fit, "data")
  if (!is.list(case_groups)) case_groups <- list(case_groups)
  if (!is.list(fit_data)) fit_data <- list(fit_data)
  if (length(case_groups) != length(fit_data)) {
    used_positions <- integer()
  } else {
    used_positions <- unlist(Map(
      function(indices, values) {
        values <- as.matrix(values)
        nonempty <- if (!nrow(values)) logical() else rowSums(!is.na(values)) > 0L
        as.integer(indices[nonempty])
      },
      case_groups,
      fit_data
    ), use.names = FALSE)
  }
  used_positions <- sort(unique(used_positions))
  source_indices <- blueprint$case_indices
  case_indices <- if (
    length(used_positions) && length(source_indices) >= max(used_positions)
  ) {
    as.integer(source_indices[used_positions])
  } else {
    integer()
  }

  scalar_character <- function(value) {
    if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1L]])
  }
  parameterization <- scalar_character(options$parameterization)
  std_lv <- isTRUE(options$std.lv)
  measurement <- if (!length(ordered)) {
    "continuous"
  } else if (setequal(ordered, variables)) {
    "ordinal"
  } else {
    "mixed"
  }
  identification <- if (std_lv) {
    "std.lv"
  } else {
    parameter_table <- lavaan::parameterTable(fit)
    latent <- as.character(lavaan::lavNames(fit, type = "lv"))
    unit_variances <- parameter_table$op == "~~" &
      parameter_table$lhs == parameter_table$rhs &
      parameter_table$lhs %in% latent &
      parameter_table$free == 0L &
      is.finite(parameter_table$ustart) &
      parameter_table$ustart == 1
    if (length(latent) && any(unit_variances)) "unit_variance" else "marker"
  }

  list(
    case_indices = case_indices,
    variables = variables,
    ordered = ordered,
    group_coding = group_coding,
    estimator_effective = scalar_character(options$estimator),
    estimator_orig = scalar_character(options$estimator.orig),
    missing = .ifa_canonical_missing(options$missing),
    parameterization = parameterization,
    std_lv = std_lv,
    identification = identification,
    measurement = measurement,
    estimand = if (measurement == "continuous") {
      "means_covariances"
    } else {
      "thresholds_polychoric_correlations"
    }
  )
}

.ifa_effective_id_fac <- function(metadata) {
  if (!is.list(metadata)) return(NA_character_)
  if (isTRUE(metadata$std_lv)) return("std.lv")
  switch(
    as.character(metadata$identification[[1L]]),
    unit_variance = "uv",
    marker = "marker",
    NA_character_
  )
}

.ifa_effective_id_cat <- function(metadata, engine_requests) {
  if (!is.list(metadata) || identical(metadata$measurement, "continuous")) {
    return(NA_character_)
  }
  requested <- engine_requests$ID.cat
  if (identical(requested, "Wu.Estabrook.2016")) "wu" else NA_character_
}

.ifa_extract_fit_record <- function(built,
                                    stage,
                                    stage_detail,
                                    blueprint,
                                    indicators,
                                    fit_variant = "scaled") {
  fit <- built$fit
  measures <- .ifa_extract_measures(fit)
  converged <- FALSE
  post_check <- FALSE
  finite_parameters <- FALSE
  negative_observed <- FALSE
  negative_latent <- FALSE
  parameter_table <- built$parameter_table

  if (!is.null(fit) && inherits(fit, "lavaan")) {
    converged_capture <- .pm_capture_conditions(lavaan::lavInspect(fit, "converged"))
    post_capture <- .pm_capture_conditions(lavaan::lavInspect(fit, "post.check"))
    estimates_capture <- .pm_capture_conditions(lavaan::parameterEstimates(fit))
    converged <- isTRUE(converged_capture$value)
    post_check <- isTRUE(post_capture$value)
    estimates <- estimates_capture$value
    if (!is.null(estimates) && nrow(estimates)) {
      finite_parameters <- all(is.finite(estimates$est)) && all(is.finite(estimates$se))
      variances <- estimates[
        estimates$op == "~~" & estimates$lhs == estimates$rhs,
        , drop = FALSE
      ]
      negative_observed <- any(
        variances$lhs %in% indicators & is.finite(variances$est) & variances$est < 0
      )
      negative_latent <- any(
        !variances$lhs %in% indicators & is.finite(variances$est) & variances$est < 0
      )
    }
  }

  effective_lavaan <- .ifa_effective_lavaan_settings(fit)
  fit_metadata <- .ifa_fit_metadata(fit, blueprint)
  engine_requests <- built$engine_requests
  if (!is.list(engine_requests)) {
    engine_requests <- list(
      estimator = if (!is.null(built$engine_args$estimator)) {
        as.character(built$engine_args$estimator[[1L]])
      } else {
        NA_character_
      },
      missing = if (!is.null(built$engine_args$missing)) {
        as.character(built$engine_args$missing[[1L]])
      } else {
        NA_character_
      },
      parameterization = if (!is.null(built$engine_args$parameterization)) {
        as.character(built$engine_args$parameterization[[1L]])
      } else {
        NA_character_
      },
      ID.fac = if (!is.null(built$engine_args$ID.fac)) {
        as.character(built$engine_args$ID.fac[[1L]])
      } else if (isTRUE(built$engine_args$std.lv)) {
        "std.lv"
      } else {
        NA_character_
      },
      ID.cat = if (!is.null(built$engine_args$ID.cat)) {
        as.character(built$engine_args$ID.cat[[1L]])
      } else {
        NA_character_
      }
    )
  }
  effective_value <- function(field, fallback) {
    if (is.list(fit_metadata) && !is.null(fit_metadata[[field]])) {
      fit_metadata[[field]]
    } else {
      fallback
    }
  }
  additional_conditions <- c(
    attr(measures, "conditions"),
    if (exists("converged_capture", inherits = FALSE)) {
      c(converged_capture$warnings, converged_capture$messages,
        post_capture$warnings, post_capture$messages,
        estimates_capture$warnings, estimates_capture$messages)
    } else {
      list()
    }
  )
  additional_errors <- Filter(Negate(is.null), list(
    attr(measures, "error"),
    if (exists("converged_capture", inherits = FALSE)) converged_capture$error else NULL,
    if (exists("post_capture", inherits = FALSE)) post_capture$error else NULL,
    if (exists("estimates_capture", inherits = FALSE)) estimates_capture$error else NULL
  ))
  record <- list(
    fit_variant = fit_variant,
    stage = stage,
    stage_detail = stage_detail,
    fit = fit,
    syntax = built$syntax,
    engine = built$engine,
    engine_args = built$engine_args,
    engine_requests = engine_requests,
    requested_constraints = built$requested_constraints,
    effective_constraints = built$effective_constraints,
    group_partial = built$group_partial,
    parameter_table = parameter_table,
    fit_metadata = fit_metadata,
    measures = measures,
    converged = converged,
    post_check = post_check,
    finite_parameters = finite_parameters,
    negative_observed_variance = negative_observed,
    negative_latent_variance = negative_latent,
    estimator_engine_request = engine_requests$estimator,
    missing_engine_request = engine_requests$missing,
    parameterization_engine_request = engine_requests$parameterization,
    ID.fac_engine_request = engine_requests$ID.fac,
    ID.cat_engine_request = engine_requests$ID.cat,
    ID.fac = engine_requests$ID.fac,
    ID.cat = engine_requests$ID.cat,
    case_indices = effective_value("case_indices", integer()),
    variables = effective_value("variables", character()),
    estimator = effective_value("estimator_effective", effective_lavaan$estimator),
    estimator_effective = effective_value(
      "estimator_effective", effective_lavaan$estimator
    ),
    estimator_orig = effective_value("estimator_orig", effective_lavaan$estimator_orig),
    estimand = effective_value("estimand", NA_character_),
    ordered = effective_value("ordered", character()),
    missing = effective_value("missing", effective_lavaan$missing),
    identification = effective_value("identification", NA_character_),
    parameterization = effective_value(
      "parameterization", effective_lavaan$parameterization
    ),
    std_lv = effective_value("std_lv", effective_lavaan$std.lv),
    ID.fac_effective = .ifa_effective_id_fac(fit_metadata),
    ID.cat_effective = .ifa_effective_id_cat(fit_metadata, engine_requests),
    group_coding = effective_value("group_coding", character()),
    measurement = effective_value("measurement", NA_character_),
    warnings = c(built$warnings, additional_conditions),
    messages = built$messages,
    error = if (!is.null(built$error)) {
      built$error
    } else if (length(additional_errors)) {
      additional_errors[[1L]]
    } else {
      NULL
    }
  )
  admissibility <- .ifa_admissibility(record)
  record$admissible <- admissibility$admissible
  record$admissibility_reason <- admissibility$reason
  record
}

.ifa_family_for_row <- function(lhs, op, rhs, indicators) {
  if (op == "=~" && rhs %in% indicators) return("loadings")
  if (op == "|" && lhs %in% indicators) return("thresholds")
  if (op == "~1" && lhs %in% indicators) return("intercepts")
  if (op == "~~" && identical(lhs, rhs) && lhs %in% indicators) return("residuals")
  NA_character_
}

.ifa_canonical_parameter <- function(lhs, op, rhs) {
  if (op == "~1") return(paste(lhs, "~ 1"))
  paste(lhs, op, rhs)
}

.ifa_constraint_families <- function(record, indicators) {
  pt <- record$parameter_table
  empty <- data.frame(
    stage = character(), family = character(), canonical = character(),
    group = integer(), group_label = character(), free = integer(),
    label = character(), plabel = character(), equality_constrained = logical(),
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(pt) || !nrow(pt) || !length(record$effective_constraints)) {
    return(empty)
  }
  family <- vapply(
    seq_len(nrow(pt)),
    function(i) .ifa_family_for_row(pt$lhs[i], pt$op[i], pt$rhs[i], indicators),
    character(1)
  )
  keep <- !is.na(family) & family %in% record$effective_constraints
  if (!any(keep)) return(empty)
  selected <- pt[keep, , drop = FALSE]
  family <- family[keep]
  canonical <- vapply(
    seq_len(nrow(selected)),
    function(i) .ifa_canonical_parameter(selected$lhs[i], selected$op[i], selected$rhs[i]),
    character(1)
  )
  labels <- as.character(selected$label)
  plabels <- if ("plabel" %in% names(selected)) as.character(selected$plabel) else rep("", nrow(selected))
  group_number <- as.integer(selected$group)
  group_label <- record$group_coding[group_number]
  equality <- logical(nrow(selected))
  keys <- paste(family, canonical, sep = "\r")
  for (key in unique(keys)) {
    idx <- which(keys == key)
    observed_groups <- sort(unique(group_label[idx]))
    actual_labels <- labels[idx]
    equality[idx] <- length(idx) == length(record$group_coding) &&
      identical(observed_groups, sort(record$group_coding)) &&
      all(nzchar(actual_labels)) && length(unique(actual_labels)) == 1L
  }
  data.frame(
    stage = rep(record$stage, nrow(selected)),
    family = family,
    canonical = canonical,
    group = group_number,
    group_label = group_label,
    free = as.integer(selected$free),
    label = labels,
    plabel = plabels,
    equality_constrained = equality,
    stringsAsFactors = FALSE
  )
}

.ifa_diagnostics_long <- function(records) {
  rows <- list()
  for (stage in names(records)) {
    record <- records[[stage]]
    for (index in names(record$measures)) {
      vector <- record$measures[[index]]
      for (variant in names(vector)) {
        rows[[length(rows) + 1L]] <- data.frame(
          level = stage,
          index = index,
          variant = variant,
          value = unname(vector[[variant]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(level = character(), index = character(),
                      variant = character(), value = numeric()))
  }
  do.call(rbind, rows)
}

.ifa_fit_table <- function(records, criteria, measurement) {
  rows <- vector("list", length(records))
  names(rows) <- names(records)
  predecessor <- NULL
  predecessor_level <- NA_character_
  chain_open <- TRUE

  for (i in seq_along(records)) {
    level <- names(records)[i]
    record <- records[[i]]
    selected <- lapply(names(record$measures), function(index)
      .ifa_select_variant(record$measures[[index]], .ifa_record_priority(record, index)))
    names(selected) <- names(record$measures)
    row <- data.frame(
      level = level,
      stage_detail = record$stage_detail,
      comparison_from = NA_character_,
      requested_constraints = paste(record$requested_constraints, collapse = "+"),
      effective_constraints = paste(record$effective_constraints, collapse = "+"),
      chi2 = selected$chi2$value,
      df = selected$df$value,
      p = selected$p$value,
      cfi = selected$cfi$value,
      rmsea = selected$rmsea$value,
      srmr = selected$srmr$value,
      chi2_variant = selected$chi2$variant,
      df_variant = selected$df$variant,
      p_variant = selected$p$variant,
      cfi_variant = selected$cfi$variant,
      rmsea_variant = selected$rmsea$variant,
      srmr_variant = selected$srmr$variant,
      delta_df = NA_real_,
      delta_cfi = NA_real_,
      delta_rmsea = NA_real_,
      delta_srmr = NA_real_,
      converged = record$converged,
      post_check = record$post_check,
      admissible = record$admissible,
      comparable = if (i == 1L) NA else FALSE,
      status = "not_evaluable",
      decision = "not_evaluable",
      reason = record$admissibility_reason,
      stringsAsFactors = FALSE
    )

    if (i == 1L) {
      if (record$admissible) {
        gate <- .ifa_configural_gate(
          c(cfi = row$cfi, rmsea = row$rmsea, srmr = row$srmr),
          criteria
        )
        row$status <- gate$status
        row$decision <- gate$decision
        row$reason <- gate$reason
        chain_open <- isTRUE(gate$passed)
        predecessor <- record
        predecessor_level <- level
      } else {
        chain_open <- FALSE
        predecessor <- record
        predecessor_level <- level
      }
    } else if (!record$admissible) {
      row$comparison_from <- predecessor_level
      row$comparable <- FALSE
      chain_open <- FALSE
      predecessor <- record
      predecessor_level <- level
    } else if (!chain_open) {
      row$comparison_from <- predecessor_level
      row$comparable <- FALSE
      row$status <- "not_evaluable"
      row$decision <- "not_evaluable"
      row$reason <- "a_less_restrictive_stage_was_not_supported"
      predecessor <- record
      predecessor_level <- level
    } else if (!is.null(predecessor)) {
      comparison <- .ifa_compare_records(predecessor, record)
      row$comparison_from <- predecessor_level
      row$comparable <- comparison$comparable
      if (!is.null(comparison$delta_df)) row$delta_df <- comparison$delta_df
      if (isTRUE(comparison$comparable)) {
        row$delta_cfi <- comparison$delta_cfi
        row$delta_rmsea <- comparison$delta_rmsea
        row$delta_srmr <- comparison$delta_srmr
        row$cfi_variant <- comparison$variants[["cfi"]]
        row$rmsea_variant <- comparison$variants[["rmsea"]]
        row$srmr_variant <- comparison$variants[["srmr"]]
        row$df_variant <- comparison$variants[["df"]]
        row$cfi <- comparison$current_values[["cfi"]]
        row$rmsea <- comparison$current_values[["rmsea"]]
        row$srmr <- comparison$current_values[["srmr"]]
        row$df <- comparison$current_values[["df"]]
        change <- .ifa_change_decision(
          comparison$delta_cfi,
          comparison$delta_rmsea,
          comparison$delta_srmr,
          criteria,
          level
        )
        if (chain_open) {
          row$status <- change$status
          row$decision <- change$decision
          row$reason <- change$reason
          if (!identical(change$decision, "supported")) chain_open <- FALSE
        } else {
          row$status <- "not_evaluable"
          row$decision <- "not_evaluable"
          row$reason <- "a_less_restrictive_stage_was_not_supported"
        }
        predecessor <- record
        predecessor_level <- level
      } else {
        optional_identification <- comparison$reason == "delta_df_not_positive" &&
          (level == "thresholds" ||
             (level == "strict" && measurement != "continuous"))
        if (optional_identification) {
          row$status <- "skipped_not_identified"
          row$decision <- "not_evaluable"
          row$reason <- comparison$reason
        } else {
          row$status <- "not_evaluable"
          row$decision <- "not_evaluable"
          row$reason <- comparison$reason
          chain_open <- FALSE
          predecessor <- record
          predecessor_level <- level
        }
      }
    }
    rows[[i]] <- row
  }
  table <- do.call(rbind, rows)
  rownames(table) <- names(records)
  table
}

.ifa_display_table <- function(fit_table, digits = 3L) {
  display <- fit_table
  numeric_columns <- vapply(display, is.numeric, logical(1))
  display[numeric_columns] <- lapply(
    display[numeric_columns],
    function(value) round(value, digits)
  )
  display
}

.ifa_conclusion_text <- function(fit_table, language) {
  supported <- fit_table$level[fit_table$decision == "supported"]
  highest <- if (length(supported)) tail(supported, 1L) else "none"
  if (language == "esp") {
    if (highest == "none") {
      return(paste0(
        "La evidencia configural no fue suficiente o no fue evaluable; ",
        "no corresponde interpretar niveles posteriores de invariancia."
      ))
    }
    return(paste0(
      "La heur\u00edstica de cambios de ajuste fue compatible secuencialmente hasta el nivel '",
      switch(highest, configural = "configural", thresholds = "umbrales",
        metric = "metrico", thresholds_loadings = "umbrales y cargas",
        scalar = "escalar", strict = "estricto", highest), "'."
    ))
  }
  if (highest == "none") {
    return(paste0(
      "Configural evidence was insufficient or not evaluable; later invariance ",
      "levels must not be interpreted."
    ))
  }
  paste0(
    "Fit-change heuristics were sequentially compatible through the '",
    if (highest == "thresholds_loadings") "equal thresholds and loadings" else highest,
    "' level."
  )
}

.ifa_condition <- function(message, code, type = c("warning", "error")) {
  type <- match.arg(type)
  structure(
    list(message = as.character(message), call = NULL, code = as.character(code)),
    class = c(paste0("psychomatic_invariance_", type), type, "condition")
  )
}

.ifa_preflight_report <- function(report, overwrite) {
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.logical(report)) {
    if (length(report) != 1L || is.na(report)) {
      stop("'report' must be FALSE, TRUE, or one output path.", call. = FALSE)
    }
    if (!report) return(NULL)
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required when an Excel report is requested.",
           call. = FALSE)
    }
    return(.pm_safe_output_path(
      .pm_report_path(NULL, "excel"),
      extension = ".xlsx",
      overwrite = overwrite
    ))
  }
  if (!is.character(report) || length(report) != 1L || is.na(report) ||
      !nzchar(trimws(report))) {
    stop("'report' must be FALSE, TRUE, or one output path.", call. = FALSE)
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required when an Excel report is requested.",
         call. = FALSE)
  }
  .pm_safe_output_path(report, extension = ".xlsx", overwrite = overwrite)
}

.ifa_conditions_for_report <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) {
    return(data.frame(
      code = character(), message = character(), class = character(),
      call = character(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    code = as.character(x$code),
    message = as.character(x$message),
    class = vapply(x$class, paste, character(1), collapse = "/"),
    call = vapply(
      x$call,
      function(value) if (is.null(value)) "" else paste(deparse(value), collapse = " "),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

.ifa_reference_table <- function() {
  data.frame(
    reference = c(
      "Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement invariance. Structural Equation Modeling: A Multidisciplinary Journal, 14(3), 464-504. https://doi.org/10.1080/10705510701301834",
      "Bowen, N. K., & Masa, R. D. (2015). Conducting measurement invariance tests with ordinal data: A guide for social work researchers. Journal of the Society for Social Work and Research, 6(2), 229-249. https://doi.org/10.1086/681607",
      "Hirschfeld, G., & von Brachel, R. (2014). Multiple-group confirmatory factor analysis in R: A tutorial in measurement invariance with continuous and ordinal indicators. Practical Assessment, Research & Evaluation, 19, Article 7. https://doi.org/10.7275/qazy-2946",
      "Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor analysis models of different levels of invariance for ordered categorical outcomes. Psychometrika, 81(4), 1014-1045. https://doi.org/10.1007/s11336-016-9506-0",
      "Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. Journal of Statistical Software, 48(2), 1-36. https://doi.org/10.18637/jss.v048.i02",
      .pm_semtools_reference()
    ),
    identifier = c(
      "doi:10.1080/10705510701301834",
      "doi:10.1086/681607",
      "Hirschfeld_vonBrachel_2014",
      "doi:10.1007/s11336-016-9506-0",
      "doi:10.18637/jss.v048.i02",
      "https://CRAN.R-project.org/package=semTools"
    ),
    stringsAsFactors = FALSE
  )
}

.ifa_report_text <- function(result) {
  lines <- c(
    if (inherits(result, "ordinal_invariance_auto")) {
      if (result$effective_settings$language == "esp") "PsychoMatic: invarianza ordinal" else
        "PsychoMatic: ordinal measurement invariance"
    } else if (result$effective_settings$language == "esp") {
      "PsychoMatic: invariancia factorial completa"
    } else {
      "PsychoMatic: complete factorial invariance"
    },
    paste0("Estimator: ", result$effective_settings$estimator),
    paste0("Measurement: ", result$measurement),
    paste0("Criteria: ", result$criteria$preset),
    result$conclusion
  )
  paste(lines, collapse = "\n")
}

#' Automated complete measurement invariance
#'
#' Fits configural, metric, scalar and strict multi-group CFA models using
#' native lavaan CFA calls and their default marker identification. WLSMV/ULSMV
#' default to ordered responses; only the four requested stages are fitted. The
#' function automates estimator selection and syntax construction while keeping
#' identification, fit variants, admissibility, and model-comparison evidence
#' available for inspection. Fit-change cutoffs are descriptive heuristics, not
#' proof of invariance.
#'
#' @inheritParams efa_auto
#' @param data A data frame, matrix, or `psychomatic_spec` containing indicators
#'   and (unless supplied as a vector) the grouping variable.
#' @param group A grouping-column name or symbol, one positive column index, or
#'   a vector with one value per original row.
#' @param model One lavaan CFA model string.
#' @param language Output language, Spanish (`"esp"`) or English (`"eng"`).
#' @param estimator `"auto"`, `"ML"`, `"MLR"`, `"ULS"`, `"WLSMV"`, or
#'   `"ULSMV"`. Explicit WLSMV/ULSMV default to `ordered = TRUE` and
#'   listwise handling. Automatic selection uses WLSMV for detected ordinal
#'   responses, otherwise ML when all groups have favourable Mardia evidence
#'   and MLR otherwise. ML/MLR/ULS analyze numeric scores as continuous.
#' @param ordered `NULL`, `TRUE`/`FALSE`, or indicator names/indices. With
#'   explicit WLSMV/ULSMV, `NULL` declares every model indicator ordered.
#'   `FALSE` explicitly requests numeric-score estimation instead. Character
#'   and unordered-factor responses require declared levels. For an identified
#'   threshold/loading invariance sequence use [ordinal_invariance_auto()].
#' @param fit_variant `"scaled"` (default) or `"robust"` for CFI and RMSEA,
#'   their changes, decisions and reports. Corrected estimators always use scaled
#'   chi-square, df and p; SRMR is uncorrected. Uncorrected ML/ULS use base
#'   indices by default. Unavailable requested robust indices remain missing;
#'   they are never silently replaced by scaled or base indices.
#' @param report `FALSE`, `TRUE` for a report in `getwd()`, or an explicit
#'   Excel path.
#' @param criteria `"auto"` or a fully named nonnegative numeric vector with
#'   `cfi`, `rmsea`, `srmr_loading`, and `srmr_other`; the three optional
#'   `configural_*` entries customize the configural gate.
#' @param missing `"auto"`, `"fiml"`, or `"listwise"`. Categorical analyses
#'   require listwise handling.
#' @param min_group_n Operational minimum analytic group size. The default of
#'   20 does not guarantee power, identification, or stable estimation.
#' @param overwrite Whether an existing explicit report path may be replaced.
#'
#' @return An object of classes `factorial_invariance_auto` and
#'   `psychomatic_result`. Raw values are not rounded. Direct legacy aliases are
#'   retained; `$warnings` is now a structured data frame and `$warnings_list`
#'   (also `$notes`) is the character compatibility view.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' model <- "F1 =~ R1 + R2 + R3 + R4"
#' factorial_invariance_auto(MY_DATA, group = "Country", model = model)
#' }
#'
#' @export
factorial_invariance_auto <- function(
  data,
  group,
  model,
  language = c("eng", "esp"),
  estimator = "auto",
  ordered = NULL,
  report = FALSE,
  criteria = "auto",
  missing = "auto",
  min_group_n = 20L,
  overwrite = FALSE,
  digits = 3, file_name = NULL, format = "excel", verbose = TRUE,
  fit_variant = c("scaled", "robust")
) {
  .ifa_auto(data, group, model, language, estimator, ordered, report, criteria,
    missing, min_group_n, overwrite, digits, file_name, format, verbose,
    group_expression = substitute(group), group_environment = parent.frame(),
    original_call = match.call(), workflow = "factorial", fit_variant = fit_variant)
}

# Shared estimation engine; public wrappers preserve the caller's group expression.
# Shared item preparation keeps complete and partial native workflows identical.
.ifa_prepare_workflow_items <- function(data, indicators, ordered, specification_levels,
                                        estimator, workflow) {
  requested_estimator <- if (is.character(estimator) && length(estimator) == 1L &&
    !is.na(estimator)) toupper(trimws(estimator)) else ""
  ordered_preparation <- if (workflow == "factorial" && is.null(ordered) &&
    requested_estimator %in% c("WLSMV", "ULSMV")) TRUE else ordered
  items <- .ifa_prepare_items(
    data,
    indicators,
    ordered = if (workflow == "factorial" && identical(ordered_preparation, FALSE)) NULL else ordered_preparation,
    specification_levels = specification_levels
  )
  source_measurement <- items$measurement
  score_coding <- list()
  if (workflow == "factorial" && (identical(ordered, FALSE) ||
      requested_estimator %in% c("ML", "MLR", "ULS"))) {
    for (item in indicators) {
      x <- items$data[[item]]
      if (is.ordered(x)) {
        labels <- suppressWarnings(as.numeric(levels(x)))
        numeric_labels <- length(labels) && all(is.finite(labels)) && !anyDuplicated(labels)
        score_coding[[item]] <- if (numeric_labels) "numeric_category_labels" else "ordered_level_ranks"
        items$data[[item]] <- if (numeric_labels) labels[as.integer(x)] else as.numeric(x)
      } else {
        score_coding[[item]] <- "original_numeric_scores"
      }
    }
    items$measurement <- "continuous"
    items$ordered <- character()
    items$item_levels[] <- "continuous"
    items$category_profile <- "continuous"
  }
  list(items = items, source_measurement = source_measurement, score_coding = score_coding)
}

.ifa_auto <- function(
  data, group, model, language, estimator, ordered, report, criteria,
  missing, min_group_n, overwrite, digits, file_name, format, verbose,
  group_expression, group_environment, original_call, workflow,
  fit_variant = c("scaled", "robust")
) {
  call <- original_call
  language <- match.arg(language, c("eng", "esp"))
  fit_variant <- match.arg(fit_variant, c("scaled", "robust"))
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  input <- .ifa_extract_input(data)
  resolved_group <- .ifa_resolve_group(
    input$data,
    group_expression,
    group_environment
  )
  group_data <- resolved_group$data
  indicators <- .ifa_extract_indicators(model)
  item_preparation <- .ifa_prepare_workflow_items(group_data, indicators, ordered,
    input$measurement_levels, estimator, workflow)
  items <- item_preparation$items
  source_measurement <- item_preparation$source_measurement
  score_coding <- item_preparation$score_coding
  if (workflow == "ordinal" && items$measurement != "ordinal") {
    stop("ordinal_invariance_auto requires all model indicators to be ordinal; declare ordered indicators explicitly when needed.", call. = FALSE)
  }
  least_squares_scores <- workflow == "factorial" && is.character(estimator) &&
    length(estimator) == 1L && !is.na(estimator) &&
    toupper(trimws(estimator)) %in% c("WLSMV", "ULSMV")
  missing_effective <- .ifa_resolve_missing(missing,
    if (least_squares_scores) "ordinal" else items$measurement)
  filtered <- .ifa_filter_analysis_rows(
    items$data,
    indicators,
    resolved_group$group_name,
    resolved_group$original_indices,
    missing_effective
  )
  items <- .ifa_refresh_analytic_items(items, filtered$data, indicators)
  analysis_data <- items$data
  sample_validation <- .ifa_validate_analysis_data(
    analysis_data,
    indicators,
    resolved_group$group_name,
    min_group_n
  )
  .ifa_validate_packages(items$measurement)

  group_sizes <- sample_validation$group_sizes
  ratio <- as.numeric(max(group_sizes) / min(group_sizes))
  if (items$measurement == "continuous") {
    normality <- .ifa_group_normality(
      analysis_data,
      indicators,
      resolved_group$group_name
    )
  } else {
    normality <- list(
      overall_normal = NA,
      overall_test = list(
        normal = NA, skew_p = NA_real_, kurt_p = NA_real_,
        status = "not_applicable_to_categorical_measurement",
        n_complete = NA_integer_
      ),
      by_group = data.frame(
        group = levels(analysis_data[[resolved_group$group_name]]),
        n_complete = as.integer(group_sizes),
        skew_p = NA_real_, kurtosis_p = NA_real_,
        status = "not_applicable", message = "",
        stringsAsFactors = FALSE
      ),
      conditions = list()
    )
  }
  criteria_effective <- .ifa_resolve_criteria(
    sample_validation$total_n,
    ratio,
    normality$overall_normal,
    items$measurement,
    criteria
  )
  estimator_info <- .ifa_resolve_estimator(
    estimator,
    items$measurement,
    normality$overall_normal,
    ordered_argument = ordered,
    allow_numeric_wls = workflow == "factorial"
  )
  parameterization <- if (workflow == "factorial") "delta" else "theta"
  parameterization_engine_request <- if (items$measurement == "continuous") {
    NA_character_
  } else {
    parameterization
  }
  identification <- if (workflow == "factorial") "marker" else "std.lv"
  id_cat <- if (items$measurement == "continuous" || workflow == "factorial") {
    NA_character_
  } else "Wu.Estabrook.2016"
  estimand <- if (items$measurement == "continuous") {
    "means_covariances"
  } else {
    "thresholds_polychoric_correlations"
  }
  blueprint <- list(
    case_indices = filtered$original_indices,
    variables = indicators,
    group_name = resolved_group$group_name,
    group_source = resolved_group$group_source,
    group_coding = filtered$group_coding,
    estimator = estimator_info$estimator,
    estimator_engine_request = estimator_info$estimator,
    estimand = estimand,
    ordered = items$ordered,
    missing_engine_request = missing_effective,
    missing = missing_effective,
    identification = identification,
    parameterization_engine_request = parameterization_engine_request,
    parameterization = parameterization,
    identification_engine_request = identification,
    ID.fac = identification,
    ID.cat = id_cat,
    ID.fac_engine_request = identification,
    ID.cat_engine_request = id_cat,
    measurement = items$measurement,
    category_profile = items$category_profile,
    model = model,
    effective_syntax = list()
  )

  stages <- .ifa_stage_sequence(if (workflow == "factorial") "continuous" else
    items$measurement, items$category_profile)
  if (workflow == "ordinal") {
    stages <- list(
      configural = list(detail = "configural", constraints = character()),
      thresholds = list(detail = "thresholds", constraints = "thresholds"),
      thresholds_loadings = list(detail = "thresholds_loadings",
        constraints = c("thresholds", "loadings"))
    )
  }
  records <- vector("list", length(stages))
  names(records) <- names(stages)
  for (stage in names(stages)) {
    built <- .ifa_build_syntax(
      model = model,
      data = analysis_data,
      group_name = resolved_group$group_name,
      measurement = items$measurement,
      ordered = items$ordered,
      estimator = estimator_info$estimator,
      constraints = stages[[stage]]$constraints,
      group_partial = "",
      missing = missing_effective,
      parameterization = parameterization,
      ID.fac = identification,
      ID.cat = if (items$measurement == "continuous") "Wu.Estabrook.2016" else id_cat,
      native = workflow == "factorial"
    )
    records[[stage]] <- .ifa_extract_fit_record(
      built,
      stage,
      stages[[stage]]$detail,
      blueprint,
      indicators,
      fit_variant = fit_variant
    )
    blueprint$effective_syntax[[stage]] <- built$syntax
  }

  actual <- .ifa_effective_lavaan_settings(records[[1L]]$fit)
  actual_metadata <- records[[1L]]$fit_metadata
  actual_identification <- if (is.list(actual_metadata)) {
    actual_metadata$identification
  } else {
    NA_character_
  }
  actual_id_fac <- .ifa_effective_id_fac(actual_metadata)
  actual_id_cat <- .ifa_effective_id_cat(
    actual_metadata, records[[1L]]$engine_requests
  )
  blueprint$estimator_effective <- if (is.list(actual_metadata)) {
    actual_metadata$estimator_effective
  } else {
    actual$estimator
  }
  blueprint$estimator_orig <- if (is.list(actual_metadata)) {
    actual_metadata$estimator_orig
  } else {
    actual$estimator_orig
  }
  blueprint$missing <- actual$missing
  blueprint$parameterization <- actual$parameterization
  blueprint$identification <- actual_identification
  blueprint$std_lv <- if (is.list(actual_metadata)) actual_metadata$std_lv else NA
  blueprint$ID.fac_effective <- actual_id_fac
  blueprint$ID.cat_effective <- actual_id_cat

  fit_table <- .ifa_fit_table(records, criteria_effective, items$measurement)
  display_table <- .ifa_display_table(fit_table)
  diagnostics <- .ifa_diagnostics_long(records)
  family_rows <- lapply(records, .ifa_constraint_families, indicators = indicators)
  constraint_families <- do.call(rbind, family_rows)
  rownames(constraint_families) <- NULL

  fit_conditions <- unlist(lapply(records, function(x) c(x$warnings, x$messages)),
                           recursive = FALSE)
  fit_errors <- Filter(Negate(is.null), lapply(records, `[[`, "error"))
  warnings <- c(normality$conditions, fit_conditions)
  if (resolved_group$dropped_missing_group > 0L) {
    warnings[[length(warnings) + 1L]] <- .ifa_condition(
      paste0(resolved_group$dropped_missing_group,
             " cases with missing group membership were removed."),
      "missing_group_rows_removed"
    )
  }
  if (filtered$dropped_missing_items > 0L) {
    warnings[[length(warnings) + 1L]] <- .ifa_condition(
      paste0(filtered$dropped_missing_items,
             " cases were removed by the item missing-data policy."),
      "missing_item_rows_removed"
    )
  }
  conclusion <- .ifa_conclusion_text(fit_table, language)

  requested_settings <- list(
    fit_variant = fit_variant,
    language = language,
    estimator = estimator,
    ordered = ordered,
    report = report,
    criteria = criteria,
    missing = missing,
    min_group_n = min_group_n,
    overwrite = overwrite
  )
  effective_settings <- list(
    fit_variant = fit_variant,
    language = language,
    estimator = estimator_info$estimator,
    estimator_engine_request = estimator_info$estimator,
    estimator_engine = actual$estimator,
    estimator_orig = actual$estimator_orig,
    ordered = items$ordered,
    missing_engine_request = missing_effective,
    missing = actual$missing,
    parameterization_engine_request = parameterization_engine_request,
    parameterization = actual$parameterization,
    identification = actual_identification,
    ID.fac = identification,
    ID.cat = id_cat,
    ID.fac_engine_request = identification,
    ID.cat_engine_request = id_cat,
    ID.fac_effective = actual_id_fac,
    ID.cat_effective = actual_id_cat,
    measurement = items$measurement
  )
  sample_info <- list(
    total_n = sample_validation$total_n,
    group_sizes = group_sizes,
    min_group_n = sample_validation$min_group_n,
    operational_minimum = as.integer(min_group_n),
    dropped_missing_group = resolved_group$dropped_missing_group,
    dropped_missing_items = filtered$dropped_missing_items
  )
  item_info <- list(
    all_indicators = indicators,
    category_counts = items$category_counts,
    ordinal_flags = indicators %in% items$ordered,
    is_ordinal = items$measurement == "ordinal",
    min_categories = min(items$category_counts),
    max_categories = max(items$category_counts),
    ordered_indicators = items$ordered,
    category_profile = items$category_profile,
    item_levels = items$item_levels
  )
  decisions <- data.frame(
    level = fit_table$level,
    decision = fit_table$decision,
    status = fit_table$status,
    reason = fit_table$reason,
    stringsAsFactors = FALSE
  )
  decision_log <- rbind(
    .pm_decision(
      "estimator",
      estimator,
      estimator_info$estimator,
      estimator_info$reason
    ),
    .pm_decision(
      "missing",
      missing,
      actual$missing,
      "resolved_for_measurement_and_estimator"
    ),
    .pm_decision(
      "criteria",
      if (is.character(criteria)) criteria else "custom",
      criteria_effective$preset,
      "chen_heuristic_scenario"
    )
  )
  if (workflow == "factorial") {
    decision_log <- rbind(decision_log, .pm_decision("measurement",
      source_measurement, items$measurement,
      "native_lavaan_four_stage_workflow"))
  }
  data_summary <- list(
    n_original = nrow(input$data),
    n_analytic = nrow(analysis_data),
    variables = indicators,
    groups = as.list(group_sizes),
    measurement = items$measurement
  )
  result <- .pm_result(
    class = if (workflow == "ordinal") {
      c("ordinal_invariance_auto", "factorial_invariance_auto")
    } else "factorial_invariance_auto",
    call = call,
    requested_settings = requested_settings,
    effective_settings = effective_settings,
    data_summary = data_summary,
    results = list(
      fit_table = fit_table,
      fitted_models = lapply(records, `[[`, "fit"),
      constraint_families = constraint_families
    ),
    warnings = warnings,
    errors = fit_errors,
    diagnostics = diagnostics,
    decision_log = decision_log,
    conclusion = conclusion
  )
  warning_messages <- result$warnings$message
  result$warnings_list <- warning_messages
  result$notes <- warning_messages
  result$estimator <- estimator_info$estimator
  result$estimator_info <- estimator_info
  result$estimator_reason <- estimator_info$reason
  result$indicators <- indicators
  result$item_info <- item_info
  result$sample_info <- sample_info
  result$normality <- normality
  result$requested_estimator <- estimator
  result$requested_ordered <- ordered
  result$chen_scenario <- criteria_effective$preset
  result$measurement <- items$measurement
  result$analysis_data <- analysis_data
  result$analysis_blueprint <- blueprint
  result$fit_records <- records
  result$fitted_models <- lapply(records, `[[`, "fit")
  result$fit_table <- fit_table
  result$display_table <- display_table
  result$decisions <- decisions
  result$criteria <- criteria_effective
  result$constraint_families <- constraint_families
  result$references <- .ifa_reference_table()
  result$workflow <- workflow
  result$effective_settings$workflow <- workflow
  result$requested_settings$workflow <- workflow
  result$effective_settings$source_measurement <- source_measurement
  result$effective_settings$score_coding <- score_coding
  if (workflow == "ordinal") {
    result$references <- rbind(result$references, data.frame(
      reference = paste0("Svetina, D., Rutkowski, L., & Rutkowski, D. (2020). ",
        "Multiple-group invariance with categorical outcomes using updated guidelines: ",
        "An illustration using Mplus and the lavaan/semTools packages. ",
        "Structural Equation Modeling: A Multidisciplinary Journal, 27(1), 111-130. ",
        "https://doi.org/10.1080/10705511.2019.1602776"),
      identifier = "doi:10.1080/10705511.2019.1602776"))
  }
  result$tables <- list(
    fit = display_table,
    decisions = decisions,
    diagnostics = diagnostics,
    criteria = data.frame(
      criterion = names(criteria_effective),
      value = vapply(criteria_effective, as.character, character(1)),
      stringsAsFactors = FALSE
    ),
    constraints = constraint_families,
    warnings = .ifa_conditions_for_report(result$warnings),
    errors = .ifa_conditions_for_report(result$errors),
    conclusion = data.frame(
      status = if (any(fit_table$decision == "supported")) "interpretable" else "not_evaluable",
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ),
    settings = data.frame(
      setting = c(
        "measurement", "estimator_api", "estimator_engine_request",
        "estimator_engine", "estimator_orig",
        "missing_engine_request", "missing",
        "parameterization_engine_request", "parameterization",
        "identification_engine_request", "identification",
        "ID.fac_engine_request", "ID.fac_effective",
        "ID.cat_engine_request", "ID.cat_effective"
      ),
      value = c(
        items$measurement, estimator_info$estimator, estimator_info$estimator,
        actual$estimator,
        actual$estimator_orig, missing_effective, actual$missing,
        parameterization_engine_request, actual$parameterization,
        identification, actual_identification,
        identification, actual_id_fac, id_cat, actual_id_cat
      ),
      stringsAsFactors = FALSE
    ),
    references = result$references
  )
  result$report_text <- .ifa_report_text(result)
  result["report_path"] <- list(NULL)

  .pm_finish_analysis(result, output, language)
}

#' Print a complete measurement-invariance result
#'
#' @param x Object returned by [factorial_invariance_auto()].
#' @param digits Display precision.
#' @param ... Additional unused arguments.
#'
#' @return `x`, invisibly.
#' @method print factorial_invariance_auto
#' @keywords internal
#' @export
print.factorial_invariance_auto <- function(x, digits = NULL, ...) {
  .pm_console_report(x, digits)
}

#' Plot complete measurement-invariance fit changes
#'
#' @param x Object returned by [factorial_invariance_auto()].
#' @param ... Additional graphical arguments (currently unused).
#'
#' @return Plotting data, invisibly.
#' @method plot factorial_invariance_auto
#' @keywords internal
#' @export
plot.factorial_invariance_auto <- function(x, ...) {
  plotting <- x$fit_table[c("level", "delta_cfi", "delta_rmsea", "delta_srmr")]
  values <- as.matrix(plotting[c("delta_cfi", "delta_rmsea", "delta_srmr")])
  if (!nrow(values) || !any(is.finite(values))) {
    graphics::plot.new()
    graphics::text(.5, .5, "No evaluable fit changes")
  } else {
    mat <- t(values)
    graphics::matplot(
      seq_len(nrow(values)), values,
      type = "b", pch = c(1, 2, 3), lty = 1,
      xaxt = "n", xlab = "Invariance stage", ylab = "Raw fit change"
    )
    labels <- plotting$level
    if (inherits(x, "ordinal_invariance_auto")) {
      es <- identical(x$effective_settings$language, "esp")
      labels <- if (es) c("Configural", "Umbrales", "Umbrales y cargas") else
        c("Configural", "Thresholds", "Thresholds and loadings")
    }
    graphics::axis(1, at = seq_len(nrow(values)), labels = labels)
    graphics::abline(h = 0, col = "grey70", lty = 2)
    graphics::legend(
      "topleft", legend = c("delta CFI", "delta RMSEA", "delta SRMR"),
      col = 1:3, pch = c(1, 2, 3), lty = 1, bty = "n"
    )
  }
  invisible(plotting)
}
