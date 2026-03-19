
.ifa_validate_packages <- function() {
  needed <- "lavaan"
  missing_pkgs <- needed[!vapply(
    needed,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]

  if (length(missing_pkgs) > 0L) {
    stop(
      sprintf(
        "Install the following package(s) before running the function: %s",
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

.ifa_resolve_group <- function(data, group) {
  if (is.character(group) && length(group) == 1L) {
    group_name <- group
  } else {
    group_expr <- substitute(group)
    group_name <- if (is.symbol(group_expr)) {
      as.character(group_expr)
    } else {
      deparse(group_expr)
    }
  }

  if (!group_name %in% names(data) && grepl("\\$", group_name)) {
    group_name <- sub("^.*\\$", "", group_name)
  }

  if (!group_name %in% names(data)) {
    stop("The grouping variable was not found in `data`.", call. = FALSE)
  }

  group_values <- data[[group_name]]
  keep_rows <- !is.na(group_values)

  if (!any(keep_rows)) {
    stop("The grouping variable contains only missing values.", call. = FALSE)
  }

  if (!all(keep_rows)) {
    data <- data[keep_rows, , drop = FALSE]
  }

  data[[group_name]] <- as.factor(data[[group_name]])

  if (nlevels(data[[group_name]]) < 2L) {
    stop("The grouping variable must contain at least two groups.", call. = FALSE)
  }

  list(
    data = data,
    group_name = group_name,
    dropped_missing_group = sum(!keep_rows)
  )
}

.ifa_extract_indicators <- function(model) {
  par_table <- lavaan::lavaanify(
    model = model,
    auto = TRUE,
    fixed.x = FALSE,
    model.type = "cfa"
  )

  loading_rows <- par_table[par_table$op == "=~", , drop = FALSE]
  if (nrow(loading_rows) == 0L) {
    stop("The model must include at least one CFA loading statement using `=~`.", call. = FALSE)
  }

  latent_names <- unique(loading_rows$lhs)
  rhs_names <- unique(loading_rows$rhs)
  indicators <- rhs_names[!rhs_names %in% latent_names]

  if (length(indicators) == 0L) {
    stop("No observed indicators were detected in the supplied model.", call. = FALSE)
  }

  indicators
}

.ifa_is_integerish <- function(x, tolerance = 1e-08) {
  x <- x[!is.na(x)]

  if (!length(x)) {
    return(FALSE)
  }

  if (is.factor(x) || is.ordered(x) || is.character(x)) {
    return(TRUE)
  }

  if (!is.numeric(x)) {
    return(FALSE)
  }

  all(abs(x - round(x)) < tolerance)
}

.ifa_prepare_numeric <- function(x) {
  if (is.ordered(x)) {
    return(as.numeric(x))
  }

  if (is.factor(x)) {
    return(as.numeric(x))
  }

  if (is.character(x)) {
    suppressWarnings(numeric_x <- as.numeric(x))
    if (all(!is.na(numeric_x[!is.na(x)]))) {
      return(numeric_x)
    }
    return(as.numeric(factor(x)))
  }

  as.numeric(x)
}

.ifa_detect_item_type <- function(data, indicators) {
  category_counts <- vapply(
    indicators,
    function(item) length(unique(data[[item]][!is.na(data[[item]])])),
    numeric(1)
  )

  ordinal_flags <- vapply(
    indicators,
    function(item) {
      item_values <- data[[item]]
      n_cat <- length(unique(item_values[!is.na(item_values)]))
      n_cat >= 2L && n_cat <= 7L && .ifa_is_integerish(item_values)
    },
    logical(1)
  )

  is_ordinal <- all(ordinal_flags)

  list(
    all_indicators = indicators,
    category_counts = category_counts,
    ordinal_flags = ordinal_flags,
    is_ordinal = is_ordinal,
    min_categories = min(category_counts),
    max_categories = max(category_counts),
    ordered_indicators = indicators[ordinal_flags]
  )
}

.ifa_mardia_test <- function(x) {
  x <- as.data.frame(x)
  x <- x[stats::complete.cases(x), , drop = FALSE]

  if (!nrow(x)) {
    return(list(
      normal = NA,
      skew_p = NA_real_,
      kurt_p = NA_real_,
      status = "no_complete_cases",
      n_complete = 0L
    ))
  }

  x[] <- lapply(x, .ifa_prepare_numeric)
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  n_obs <- nrow(x)
  n_var <- ncol(x)
  sampled <- FALSE

  if (n_obs > 2000L) {
    sample_idx <- unique(round(seq(1, n_obs, length.out = 2000L)))
    x <- x[sample_idx, , drop = FALSE]
    n_obs <- nrow(x)
    sampled <- TRUE
  }

  if (n_obs < (n_var + 5L)) {
    return(list(
      normal = NA,
      skew_p = NA_real_,
      kurt_p = NA_real_,
      status = "insufficient_n",
      n_complete = n_obs,
      sampled = sampled
    ))
  }

  covariance <- stats::cov(x)
  if (any(!is.finite(covariance)) || qr(covariance)$rank < ncol(covariance)) {
    return(list(
      normal = NA,
      skew_p = NA_real_,
      kurt_p = NA_real_,
      status = "singular_covariance",
      n_complete = n_obs,
      sampled = sampled
    ))
  }

  centered <- scale(x, center = TRUE, scale = FALSE)
  inv_cov <- solve(covariance)
  distance_matrix <- centered %*% inv_cov %*% t(centered)

  b1p <- sum(distance_matrix^3) / (n_obs^2)
  skew_chisq <- n_obs * b1p / 6
  skew_df <- n_var * (n_var + 1) * (n_var + 2) / 6
  skew_p <- stats::pchisq(skew_chisq, df = skew_df, lower.tail = FALSE)

  mahal_sq <- diag(distance_matrix)
  b2p <- mean(mahal_sq^2)
  z_kurt <- (b2p - (n_var * (n_var + 2))) / sqrt(8 * n_var * (n_var + 2) / n_obs)
  kurt_p <- 2 * stats::pnorm(abs(z_kurt), lower.tail = FALSE)

  list(
    normal = (skew_p > 0.05) && (kurt_p > 0.05),
    skew_p = skew_p,
    kurt_p = kurt_p,
    status = "ok",
    n_complete = n_obs,
    sampled = sampled
  )
}

.ifa_group_normality <- function(data, indicators, group_name) {
  split_data <- split(data[, indicators, drop = FALSE], data[[group_name]])

  by_group <- lapply(split_data, .ifa_mardia_test)
  group_names <- names(by_group)
  normal_flags <- vapply(by_group, function(x) x$normal, logical(1))
  overall_test <- .ifa_mardia_test(data[, indicators, drop = FALSE])

  if (any(!is.na(normal_flags) & !normal_flags)) {
    overall_normal <- FALSE
  } else if (all(!is.na(normal_flags)) && all(normal_flags)) {
    overall_normal <- TRUE
  } else {
    overall_normal <- overall_test$normal
  }

  details <- data.frame(
    group = group_names,
    normal = normal_flags,
    skew_p = vapply(by_group, function(x) x$skew_p, numeric(1)),
    kurt_p = vapply(by_group, function(x) x$kurt_p, numeric(1)),
    status = vapply(by_group, function(x) x$status, character(1)),
    n_complete = vapply(by_group, function(x) x$n_complete, numeric(1)),
    stringsAsFactors = FALSE
  )

  list(
    overall_normal = overall_normal,
    overall_test = overall_test,
    by_group = details
  )
}

.ifa_chen_scenario <- function(group_sizes) {
  ratio <- max(group_sizes) / min(group_sizes)
  min_group_n <- min(group_sizes)

  if (min_group_n >= 300L && ratio <= 1.50) {
    return("large_equal")
  }

  "small_unequal"
}

.ifa_select_estimator <- function(item_info, sample_info, normality_info) {
  max_cat <- item_info$max_categories

  if (item_info$is_ordinal && max_cat <= 5L) {
    return(list(
      estimator = "WLSMV",
      uses_ordered = TRUE,
      ordered = item_info$ordered_indicators,
      handling = "ordered",
      reason_key = "few_categories",
      selection_mode = "auto"
    ))
  }

  if (item_info$is_ordinal && max_cat >= 6L) {
    return(list(
      estimator = "MLR",
      uses_ordered = FALSE,
      ordered = NULL,
      handling = "continuous",
      reason_key = "six_plus_categories",
      selection_mode = "auto"
    ))
  }

  if (isTRUE(normality_info$overall_normal)) {
    return(list(
      estimator = "ML",
      uses_ordered = FALSE,
      ordered = NULL,
      handling = "continuous",
      reason_key = "continuous_normal",
      selection_mode = "auto"
    ))
  }

  list(
    estimator = "MLR",
    uses_ordered = FALSE,
    ordered = NULL,
    handling = "continuous",
    reason_key = "continuous_nonnormal_or_unknown",
    selection_mode = "auto"
  )
}

.ifa_resolve_estimator <- function(estimator, ordered, item_info, sample_info, normality_info) {
  estimator <- trimws(estimator)
  if (!nzchar(estimator)) {
    stop("`estimator` must be 'auto', 'ML', 'MLR', 'ULS', or 'WLSMV'.", call. = FALSE)
  }

  estimator_upper <- toupper(estimator)
  if (tolower(estimator) == "auto") {
    resolved <- .ifa_select_estimator(item_info, sample_info, normality_info)
  } else {
    allowed <- c("ML", "MLR", "ULS", "WLSMV")
    if (!estimator_upper %in% allowed) {
      stop("`estimator` must be 'auto', 'ML', 'MLR', 'ULS', or 'WLSMV'.", call. = FALSE)
    }

    uses_ordered <- estimator_upper %in% c("ULS", "WLSMV") && item_info$is_ordinal

    resolved <- list(
      estimator = estimator_upper,
      uses_ordered = uses_ordered,
      ordered = if (uses_ordered) item_info$ordered_indicators else NULL,
      handling = if (uses_ordered) "ordered" else "continuous",
      reason_key = "manual_override",
      selection_mode = "manual"
    )
  }

  if (!is.null(ordered)) {
    if (!is.logical(ordered) || length(ordered) != 1L || is.na(ordered)) {
      stop("`ordered` must be TRUE, FALSE, or NULL.", call. = FALSE)
    }

    if (isTRUE(ordered)) {
      if (!resolved$estimator %in% c("ULS", "WLSMV")) {
        stop("`ordered = TRUE` is only compatible with `ULS` or `WLSMV`.", call. = FALSE)
      }

      resolved$uses_ordered <- TRUE
      resolved$ordered <- item_info$all_indicators
      resolved$handling <- "ordered"
    } else {
      resolved$uses_ordered <- FALSE
      resolved$ordered <- NULL
      resolved$handling <- "continuous"
    }
  }

  resolved
}

.ifa_prepare_analysis_data <- function(data, indicators, estimator_info) {
  prepared <- data

  if (estimator_info$uses_ordered) {
    for (item in indicators) {
      if (is.character(prepared[[item]])) {
        prepared[[item]] <- ordered(prepared[[item]])
      } else if (is.factor(prepared[[item]]) && !is.ordered(prepared[[item]])) {
        prepared[[item]] <- ordered(prepared[[item]], levels = levels(prepared[[item]]))
      }
    }
    return(prepared)
  }

  prepared[indicators] <- lapply(prepared[indicators], .ifa_prepare_numeric)
  prepared
}

.ifa_constraints <- function(uses_ordered) {
  list(
    configural = character(0),
    metric = "loadings",
    scalar = c("loadings", "intercepts"),
    strict = c("loadings", "intercepts", "residuals")
  )
}

.ifa_fit_model <- function(model, data, group_name, estimator_info, group_equal = character(0)) {
  warnings_collected <- character(0)

  fit_args <- list(
    model = model,
    data = data,
    group = group_name,
    estimator = estimator_info$estimator,
    std.lv = TRUE,
    meanstructure = TRUE
  )

  if (estimator_info$uses_ordered) {
    fit_args$ordered <- TRUE
    fit_args$parameterization <- "delta"
    fit_args$missing <- "listwise"
  } else {
    fit_args$missing <- "fiml"
  }

  if (length(group_equal) > 0L) {
    fit_args$group.equal <- group_equal
  }

  fit <- withCallingHandlers(
    tryCatch(
      do.call(lavaan::cfa, fit_args),
      error = function(e) {
        structure(
          list(message = conditionMessage(e)),
          class = "ifa_fit_error"
        )
      }
    ),
    warning = function(w) {
      warnings_collected <<- c(warnings_collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    fit = fit,
    warnings = unique(warnings_collected),
    estimator = estimator_info$estimator
  )
}

.ifa_pick_measure <- function(fit, candidates) {
  values <- suppressWarnings(lavaan::fitMeasures(fit, candidates))

  for (candidate in candidates) {
    if (candidate %in% names(values) && !is.na(values[[candidate]])) {
      return(unname(values[[candidate]]))
    }
  }

  NA_real_
}

.ifa_extract_fit_info <- function(fit_result) {
  if (inherits(fit_result$fit, "ifa_fit_error")) {
    return(data.frame(
      chi2 = NA_real_,
      df = NA_real_,
      p = NA_real_,
      cfi = NA_real_,
      rmsea = NA_real_,
      srmr = NA_real_,
      converged = FALSE,
      status = "error",
      message = fit_result$fit$message,
      warnings = paste(fit_result$warnings, collapse = " | "),
      stringsAsFactors = FALSE
    ))
  }

  fit <- fit_result$fit
  estimator_name <- toupper(fit_result$estimator)
  converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  status <- if (converged) "ok" else "no_convergence"

  if (identical(estimator_name, "ML")) {
    chi2_candidates <- c("chisq")
    df_candidates <- c("df")
    p_candidates <- c("pvalue")
    cfi_candidates <- c("cfi")
    rmsea_candidates <- c("rmsea")
  } else {
    chi2_candidates <- c("chisq.scaled", "chisq")
    df_candidates <- c("df.scaled", "df")
    p_candidates <- c("pvalue.scaled", "pvalue")
    cfi_candidates <- c("cfi.scaled", "cfi")
    rmsea_candidates <- c("rmsea.scaled", "rmsea")
  }

  data.frame(
    chi2 = .ifa_pick_measure(fit, chi2_candidates),
    df = .ifa_pick_measure(fit, df_candidates),
    p = .ifa_pick_measure(fit, p_candidates),
    cfi = .ifa_pick_measure(fit, cfi_candidates),
    rmsea = .ifa_pick_measure(fit, rmsea_candidates),
    srmr = .ifa_pick_measure(fit, c("srmr")),
    converged = converged,
    status = status,
    message = "",
    warnings = paste(fit_result$warnings, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

.ifa_chen_thresholds <- function(level, scenario) {
  if (!level %in% c("metric", "scalar", "strict")) {
    stop("Chen thresholds are only defined for metric, scalar, and strict levels.", call. = FALSE)
  }

  list(
    cfi = if (scenario == "large_equal") -0.010 else -0.005,
    rmsea = 0.015
  )
}

.ifa_decide_level <- function(level, delta_cfi, delta_rmsea, scenario, is_estimable) {
  if (level == "configural") {
    return(list(decision = "baseline", thresholds = NULL))
  }

  if (!is_estimable || any(is.na(c(delta_cfi, delta_rmsea)))) {
    return(list(decision = "not_estimated", thresholds = .ifa_chen_thresholds(level, scenario)))
  }

  thresholds <- .ifa_chen_thresholds(level, scenario)

  non_invariant <- (delta_cfi <= thresholds$cfi) &&
    (delta_rmsea > thresholds$rmsea)

  list(
    decision = if (non_invariant) "not_supported" else "supported",
    thresholds = thresholds
  )
}

.ifa_estimator_reason <- function(estimator_info, item_info, sample_info, normality_info, language) {
  n_text <- paste0("min group n = ", sample_info$min_group_n)
  cat_text <- paste0(
    item_info$min_categories,
    "-",
    item_info$max_categories,
    " categories"
  )

  normality_text <- if (isTRUE(normality_info$overall_normal)) {
    "multivariate normality was not rejected"
  } else if (isFALSE(normality_info$overall_normal)) {
    "multivariate normality was rejected"
  } else {
    "multivariate normality could not be determined"
  }

  if (language == "esp") {
    normality_text <- if (isTRUE(normality_info$overall_normal)) {
      "la normalidad multivariada no fue rechazada"
    } else if (isFALSE(normality_info$overall_normal)) {
      "la normalidad multivariada fue rechazada"
    } else {
      "la normalidad multivariada no pudo determinarse"
    }

    return(switch(
      estimator_info$reason_key,
      few_categories = paste(
        "Se selecciono WLSMV porque todos los indicadores tienen 5 o menos categorias de respuesta",
        paste0("(", cat_text, ")"),
        "y se tratan como ordinales."
      ),
      six_plus_categories = paste(
        "Se selecciono MLR porque se detectaron 6 o mas opciones de respuesta",
        paste0("(", cat_text, "),"),
        "por lo que los indicadores se trataran como aproximacion continua robusta."
      ),
      ordinal_large_n_nonnormal = paste(
        "Se selecciono WLSMV porque los indicadores son ordinales",
        paste0("(", cat_text, "),"),
        normality_text,
        "y el tamano muestral es amplio",
        paste0("(", n_text, ").")
      ),
      ordinal_nonnormal = paste(
        "Se selecciono MLR porque los indicadores tienen 5-7 categorias,",
        normality_text,
        "y se analizaran como aproximacion continua robusta."
      ),
      ordinal_normal = paste(
        "Se selecciono ML porque los indicadores tienen 5-7 categorias,",
        normality_text,
        "y el tamano muestral es adecuado",
        paste0("(", n_text, ").")
      ),
      ordinal_small_n_or_unknown = paste(
        "Se selecciono MLR porque los indicadores tienen 5-7 categorias y el tamano muestral es limitado o la normalidad no fue concluyente.",
        paste0("(", n_text, ").")
      ),
      continuous_normal = paste(
        "Se selecciono ML porque los indicadores se trataron como continuos y",
        normality_text,
        "."
      ),
      continuous_nonnormal_or_unknown = paste(
        "Se selecciono MLR porque los indicadores se trataron como continuos y",
        normality_text,
        "."
      ),
      manual_override = paste(
        "El estimador fue fijado manualmente en",
        estimator_info$estimator,
        "por el usuario."
      ),
      "Heuristic estimator selection."
    ))
  }

  switch(
    estimator_info$reason_key,
    few_categories = paste(
      "WLSMV was selected because all indicators have 5 or fewer response categories",
      paste0("(", cat_text, ")"),
      "and are treated as ordinal."
    ),
    six_plus_categories = paste(
      "MLR was selected because 6 or more response options were detected",
      paste0("(", cat_text, "),"),
      "so the indicators are analyzed as approximately continuous with robust corrections."
    ),
    ordinal_large_n_nonnormal = paste(
      "WLSMV was selected because the indicators are ordinal",
      paste0("(", cat_text, "),"),
      normality_text,
      "and the sample size is large",
      paste0("(", n_text, ").")
    ),
    ordinal_nonnormal = paste(
      "MLR was selected because the indicators have 5-7 categories,",
      normality_text,
      "and they are analyzed as approximately continuous with robust corrections."
    ),
    ordinal_normal = paste(
      "ML was selected because the indicators have 5-7 categories,",
      normality_text,
      "and the sample size is adequate",
      paste0("(", n_text, ").")
    ),
    ordinal_small_n_or_unknown = paste(
      "MLR was selected because the indicators have 5-7 categories and the sample size is limited or normality was inconclusive.",
      paste0("(", n_text, ").")
    ),
    continuous_normal = paste(
      "ML was selected because the indicators were treated as continuous and",
      normality_text,
      "."
    ),
    continuous_nonnormal_or_unknown = paste(
      "MLR was selected because the indicators were treated as continuous and",
      normality_text,
      "."
    ),
    manual_override = paste(
      "The estimator was manually set to",
      estimator_info$estimator,
      "by the user."
    ),
    "Heuristic estimator selection."
  )
}

.ifa_format_p <- function(x) {
  if (is.na(x)) {
    return("NA")
  }

  formatC(x, format = "f", digits = 3)
}

.ifa_normality_label <- function(value, language) {
  if (isTRUE(value)) {
    return(if (language == "esp") "Si" else "Yes")
  }

  if (isFALSE(value)) {
    return(if (language == "esp") "No" else "No")
  }

  if (language == "esp") "No evaluable" else "Not assessable"
}

.ifa_normality_text <- function(normality_info, language) {
  overall_label <- .ifa_normality_label(normality_info$overall_normal, language)
  skew_p <- .ifa_format_p(normality_info$overall_test$skew_p)
  kurt_p <- .ifa_format_p(normality_info$overall_test$kurt_p)

  if (language == "esp") {
    return(paste0(
      overall_label,
      " (p asimetria = ", skew_p,
      "; p curtosis = ", kurt_p, ")"
    ))
  }

  paste0(
    overall_label,
    " (skewness p = ", skew_p,
    "; kurtosis p = ", kurt_p, ")"
  )
}

.ifa_scenario_label <- function(scenario, language) {
  if (language == "esp") {
    if (scenario == "large_equal") {
      return("N grande y grupos aproximadamente equilibrados")
    }
    return("N pequeno o grupos desbalanceados")
  }

  if (scenario == "large_equal") {
    return("Large N and approximately balanced groups")
  }

  "Small N or unbalanced groups"
}

.ifa_status_label <- function(status, language) {
  if (language == "esp") {
    return(switch(
      status,
      ok = "OK",
      no_convergence = "No converge",
      error = "Error",
      status
    ))
  }

  switch(
    status,
    ok = "OK",
    no_convergence = "No convergence",
    error = "Error",
    status
  )
}

.ifa_decision_label <- function(decision, language) {
  if (language == "esp") {
    return(switch(
      decision,
      baseline = "Base",
      supported = "Cumple",
      not_supported = "No cumple",
      not_estimated = "No estimado",
      decision
    ))
  }

  switch(
    decision,
    baseline = "Baseline",
    supported = "Supported",
    not_supported = "Not supported",
    not_estimated = "Not estimated",
    decision
  )
}

.ifa_handling_label <- function(handling, language) {
  if (language == "esp") {
    return(if (handling == "ordered") "Ordinal" else "Continuo")
  }

  if (handling == "ordered") "Ordinal" else "Continuous"
}

.ifa_conclusion <- function(highest_level, language) {
  if (language == "esp") {
    if (identical(highest_level, "strict")) {
      return("Se logro la invarianza a nivel estricto.")
    }

    if (identical(highest_level, "scalar")) {
      return("Se logro la invarianza hasta el nivel escalar.")
    }

    if (identical(highest_level, "metric")) {
      return("Se logro la invarianza hasta el nivel metrico.")
    }

    return("No se pudo demostrar la invarianza; se sugiere verificar invarianza parcial o invarianza por alineamiento.")
  }

  if (identical(highest_level, "strict")) {
    return("Strict invariance was achieved.")
  }

  if (identical(highest_level, "scalar")) {
    return("Invariance was achieved up to the scalar level.")
  }

  if (identical(highest_level, "metric")) {
    return("Invariance was achieved up to the metric level.")
  }

  "Invariance could not be demonstrated; partial invariance or alignment optimization should be considered."
}

.ifa_prepare_report_path <- function(report, language) {
  if (identical(report, FALSE)) {
    return(NULL)
  }

  if (identical(report, TRUE)) {
    file_name <- if (language == "esp") {
      "reporte_invarianza_factorial.xlsx"
    } else {
      "factorial_invariance_report.xlsx"
    }
    return(file.path(getwd(), file_name))
  }

  if (is.character(report) && length(report) == 1L && nzchar(report)) {
    return(report)
  }

  stop("`report` must be FALSE, TRUE, or a valid file path.", call. = FALSE)
}

.ifa_summary_sheet <- function(
  estimator_info,
  estimator_reason,
  item_info,
  sample_info,
  normality_info,
  scenario,
  conclusion,
  report_path,
  language,
  notes = character(0)
) {
  label_col <- if (language == "esp") "Elemento" else "Item"
  value_col <- if (language == "esp") "Valor" else "Value"

  base_rows <- data.frame(
    label = c(
      if (language == "esp") "Estimador seleccionado" else "Selected estimator",
      if (language == "esp") "Tratamiento de los indicadores" else "Indicator handling",
      if (language == "esp") "Justificacion" else "Rationale",
      if (language == "esp") "Tamano muestral total" else "Total sample size",
      if (language == "esp") "Tamano por grupos" else "Group sizes",
      if (language == "esp") "Opciones de respuesta" else "Response options",
      if (language == "esp") "Normalidad multivariada (Mardia)" else "Multivariate normality (Mardia)",
      if (language == "esp") "Escenario Chen (2007)" else "Chen (2007) scenario",
      if (language == "esp") "Conclusion final" else "Final conclusion",
      if (language == "esp") "Archivo Excel" else "Excel file"
    ),
    value = c(
      estimator_info$estimator,
      .ifa_handling_label(estimator_info$handling, language),
      estimator_reason,
      sample_info$total_n,
      paste(names(sample_info$group_sizes), sample_info$group_sizes, sep = "=", collapse = "; "),
      paste0(item_info$min_categories, "-", item_info$max_categories),
      .ifa_normality_text(normality_info, language),
      .ifa_scenario_label(scenario, language),
      conclusion,
      if (is.null(report_path)) {
        if (language == "esp") "No solicitado" else "Not requested"
      } else {
        normalizePath(report_path, winslash = "/", mustWork = FALSE)
      }
    ),
    stringsAsFactors = FALSE
  )

  if (length(notes) > 0L) {
    note_label <- if (language == "esp") "Nota" else "Note"
    note_rows <- data.frame(
      label = rep(note_label, length(notes)),
      value = notes,
      stringsAsFactors = FALSE
    )
    base_rows <- rbind(base_rows, note_rows)
  }

  names(base_rows) <- c(label_col, value_col)
  base_rows
}

.ifa_export_table <- function(table_df, language) {
  export_df <- table_df
  export_df <- data.frame(
    Level = rownames(export_df),
    export_df,
    row.names = NULL,
    check.names = FALSE
  )
  names(export_df)[1] <- if (language == "esp") "Nivel" else "Level"
  export_df
}

.ifa_fit_issues <- function(fit_table, language) {
  issues <- character(0)

  for (i in seq_len(nrow(fit_table))) {
    level_name <- rownames(fit_table)[i]
    level_label <- if (language == "esp") {
      switch(
        level_name,
        configural = "Configural",
        metric = "Metrica",
        scalar = "Escalar",
        strict = "Estricta",
        level_name
      )
    } else {
      switch(
        level_name,
        configural = "Configural",
        metric = "Metric",
        scalar = "Scalar",
        strict = "Strict",
        level_name
      )
    }

    if (!is.na(fit_table$message[i]) && nzchar(fit_table$message[i])) {
      issues <- c(
        issues,
        if (language == "esp") {
          paste0("Error en ", level_label, ": ", fit_table$message[i])
        } else {
          paste0("Error in ", level_label, ": ", fit_table$message[i])
        }
      )
    }

    if (!is.na(fit_table$warnings[i]) && nzchar(fit_table$warnings[i])) {
      issues <- c(
        issues,
        if (language == "esp") {
          paste0("Advertencia en ", level_label, ": ", fit_table$warnings[i])
        } else {
          paste0("Warning in ", level_label, ": ", fit_table$warnings[i])
        }
      )
    }
  }

  unique(issues)
}

.ifa_ordered_data_issues <- function(data, indicators, group_name, language) {
  issues <- character(0)

  for (item in indicators) {
    group_categories <- tapply(
      data[[item]],
      data[[group_name]],
      function(x) length(unique(x[!is.na(x)]))
    )

    bad_groups <- names(group_categories)[is.na(group_categories) | group_categories < 2L]

    if (length(bad_groups) > 0L) {
      issues <- c(
        issues,
        if (language == "esp") {
          paste0(
            "El item `", item, "` tiene menos de dos categorias observadas en: ",
            paste(bad_groups, collapse = ", "),
            ". Esto puede impedir el ajuste con WLSMV/ULS."
          )
        } else {
          paste0(
            "Item `", item, "` has fewer than two observed categories in: ",
            paste(bad_groups, collapse = ", "),
            ". This can prevent estimation with WLSMV/ULS."
          )
        }
      )
    }
  }

  unique(issues)
}

.ifa_write_excel <- function(summary_df, table_df, path) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Summary")
    openxlsx::writeData(wb, "Summary", summary_df)
    openxlsx::addWorksheet(wb, "Invariance")
    openxlsx::writeData(wb, "Invariance", table_df)
    openxlsx::saveWorkbook(wb, normalized_path, overwrite = TRUE)
    return(normalized_path)
  }

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      x = list(
        Summary = summary_df,
        Invariance = table_df
      ),
      path = normalized_path
    )
    return(normalized_path)
  }

  warning(
    "Excel export was requested, but neither `openxlsx` nor `writexl` is installed. The analysis was completed without the workbook.",
    call. = FALSE
  )

  NULL
}

.ifa_report_text <- function(
  table_df,
  estimator_info,
  estimator_reason,
  sample_info,
  item_info,
  normality_info,
  scenario,
  conclusion,
  language,
  report_path = NULL,
  ordered_note = FALSE,
  notes = character(0)
) {
  threshold_metric <- .ifa_chen_thresholds("metric", scenario)
  threshold_scalar <- .ifa_chen_thresholds("scalar", scenario)
  threshold_strict <- .ifa_chen_thresholds("strict", scenario)

  format_rule <- function(thresholds, language) {
    if (language == "esp") {
      paste0(
        "\u0394CFI <= ", formatC(thresholds$cfi, format = "f", digits = 3),
        " y \u0394RMSEA > ", formatC(thresholds$rmsea, format = "f", digits = 3)
      )
    } else {
      paste0(
        "\u0394CFI <= ", formatC(thresholds$cfi, format = "f", digits = 3),
        " and \u0394RMSEA > ", formatC(thresholds$rmsea, format = "f", digits = 3)
      )
    }
  }

  table_output <- capture.output(print(table_df))

  if (language == "esp") {
    lines <- c(
      "========================================",
      "REPORTE DE INVARIANZA FACTORIAL",
      "========================================",
      paste0("Estimador seleccionado: ", estimator_info$estimator),
      paste0("Tratamiento de los indicadores: ", .ifa_handling_label(estimator_info$handling, language)),
      paste0("Justificacion: ", estimator_reason),
      paste0("Tamano muestral total: ", sample_info$total_n),
      paste0("Tamano por grupos: ", paste(names(sample_info$group_sizes), sample_info$group_sizes, sep = "=", collapse = "; ")),
      paste0("Opciones de respuesta detectadas: ", item_info$min_categories, "-", item_info$max_categories),
      paste0("Normalidad multivariada (Mardia): ", .ifa_normality_text(normality_info, language)),
      paste0("Escenario Chen (2007): ", .ifa_scenario_label(scenario, language)),
      "Puntos de corte de Chen (2007):",
      paste0("  Metrica: no invariancia si ", format_rule(threshold_metric, language)),
      paste0("  Escalar: no invariancia si ", format_rule(threshold_scalar, language)),
      paste0("  Estricta: no invariancia si ", format_rule(threshold_strict, language)),
      if (ordered_note) {
        "Nota: con indicadores ordinales, la funcion usa parameterization = 'delta' y aplica restricciones secuenciales sobre loadings, intercepts y residuals."
      } else {
        NULL
      },
      if (length(notes) > 0L) c("Notas:", paste0("- ", notes)) else NULL,
      "",
      table_output,
      "",
      paste0("Conclusion final: ", conclusion),
      if (!is.null(report_path)) {
        paste0("Archivo Excel: ", normalizePath(report_path, winslash = "/", mustWork = FALSE))
      } else {
        NULL
      }
    )

    return(paste(lines, collapse = "\n"))
  }

  lines <- c(
    "========================================",
    "FACTORIAL INVARIANCE REPORT",
    "========================================",
    paste0("Selected estimator: ", estimator_info$estimator),
    paste0("Indicator handling: ", .ifa_handling_label(estimator_info$handling, language)),
    paste0("Rationale: ", estimator_reason),
    paste0("Total sample size: ", sample_info$total_n),
    paste0("Group sizes: ", paste(names(sample_info$group_sizes), sample_info$group_sizes, sep = "=", collapse = "; ")),
    paste0("Detected response options: ", item_info$min_categories, "-", item_info$max_categories),
    paste0("Multivariate normality (Mardia): ", .ifa_normality_text(normality_info, language)),
    paste0("Chen (2007) scenario: ", .ifa_scenario_label(scenario, language)),
    "Chen (2007) decision rules:",
    paste0("  Metric: noninvariance if ", format_rule(threshold_metric, language)),
    paste0("  Scalar: noninvariance if ", format_rule(threshold_scalar, language)),
    paste0("  Strict: noninvariance if ", format_rule(threshold_strict, language)),
    if (ordered_note) {
      "Note: with ordinal indicators, the function uses parameterization = 'delta' and applies sequential constraints on loadings, intercepts, and residuals."
    } else {
      NULL
    },
    if (length(notes) > 0L) c("Notes:", paste0("- ", notes)) else NULL,
    "",
    table_output,
    "",
    paste0("Final conclusion: ", conclusion),
    if (!is.null(report_path)) {
      paste0("Excel file: ", normalizePath(report_path, winslash = "/", mustWork = FALSE))
    } else {
      NULL
    }
  )

  paste(lines, collapse = "\n")
}

#' Invarianza factorial multigrupo automatizada
#'
#' Ajusta secuencialmente modelos configural, metrico, escalar y estricto para
#' evaluar el nivel maximo de invarianza soportado por los datos.
#'
#' @param data Data frame con los items y la variable de grupo.
#' @param group Nombre o expresion de la variable de agrupacion.
#' @param model Modelo en sintaxis `lavaan`.
#' @param language Idioma del reporte: `"esp"` o `"eng"`.
#' @param estimator Estimador solicitado o `"auto"` para seleccion automatica.
#' @param ordered Configuracion opcional para el tratamiento ordinal.
#' @param report `FALSE`, `TRUE` o ruta para exportar un reporte.
#'
#' @return Objeto de clase `factorial_invariance_auto`, devuelto invisiblemente.
#'
#' @examples
#' # modelo <- "F1 =~ i1 + i2 + i3"
#' # factorial_invariance_auto(mi_datos, group = "pais", model = modelo)
#'
#' @export
factorial_invariance_auto <- function(
  data,
  group,
  model,
  language = c("esp", "eng"),
  estimator = "auto",
  ordered = NULL,
  report = FALSE
) {
  .ifa_validate_packages()

  language <- match.arg(language)

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(model) || length(model) != 1L) {
    stop("`model` must be a single lavaan model string.", call. = FALSE)
  }

  resolved_group <- .ifa_resolve_group(data, group)
  data <- resolved_group$data
  group_name <- resolved_group$group_name

  indicators <- .ifa_extract_indicators(model)
  missing_indicators <- setdiff(indicators, names(data))
  if (length(missing_indicators) > 0L) {
    stop(
      sprintf(
        "The following indicators from the model were not found in `data`: %s",
        paste(missing_indicators, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  item_info <- .ifa_detect_item_type(data, indicators)

  group_sizes <- table(data[[group_name]])
  sample_info <- list(
    total_n = nrow(data),
    group_sizes = group_sizes,
    min_group_n = min(group_sizes),
    dropped_missing_group = resolved_group$dropped_missing_group
  )

  normality_info <- .ifa_group_normality(data, indicators, group_name)
  scenario <- .ifa_chen_scenario(group_sizes)
  estimator_info <- .ifa_resolve_estimator(estimator, ordered, item_info, sample_info, normality_info)
  estimator_reason <- .ifa_estimator_reason(
    estimator_info = estimator_info,
    item_info = item_info,
    sample_info = sample_info,
    normality_info = normality_info,
    language = language
  )

  analysis_data <- .ifa_prepare_analysis_data(data, indicators, estimator_info)
  level_constraints <- .ifa_constraints(estimator_info$uses_ordered)
  level_order <- c("configural", "metric", "scalar", "strict")

  fitted_models <- lapply(level_order, function(level) {
    .ifa_fit_model(
      model = model,
      data = analysis_data,
      group_name = group_name,
      estimator_info = estimator_info,
      group_equal = level_constraints[[level]]
    )
  })
  names(fitted_models) <- level_order

  fit_table <- do.call(
    rbind,
    lapply(fitted_models, .ifa_extract_fit_info)
  )
  rownames(fit_table) <- level_order

  fit_table$delta_cfi <- c(NA_real_, diff(fit_table$cfi))
  fit_table$delta_rmsea <- c(NA_real_, diff(fit_table$rmsea))

  decisions <- vector("list", length(level_order))
  names(decisions) <- level_order

  for (idx in seq_along(level_order)) {
    level <- level_order[idx]
    previous_ok <- if (idx == 1L) {
      TRUE
    } else {
      fit_table$converged[idx] && fit_table$converged[idx - 1L]
    }

    decisions[[level]] <- .ifa_decide_level(
      level = level,
      delta_cfi = fit_table$delta_cfi[idx],
      delta_rmsea = fit_table$delta_rmsea[idx],
      scenario = scenario,
      is_estimable = previous_ok
    )
  }

  fit_table$decision <- vapply(decisions, function(x) x$decision, character(1))

  highest_level <- "none"
  if (fit_table$converged[1]) {
    highest_level <- "configural"
  }
  if (fit_table["metric", "decision"] == "supported") {
    highest_level <- "metric"
  }
  if (
    fit_table["metric", "decision"] == "supported" &&
      fit_table["scalar", "decision"] == "supported"
  ) {
    highest_level <- "scalar"
  }
  if (
    fit_table["metric", "decision"] == "supported" &&
      fit_table["scalar", "decision"] == "supported" &&
      fit_table["strict", "decision"] == "supported"
  ) {
    highest_level <- "strict"
  }

  conclusion <- .ifa_conclusion(highest_level, language)

  display_levels <- if (language == "esp") {
    c(configural = "Configural", metric = "Metrica", scalar = "Escalar", strict = "Estricta")
  } else {
    c(configural = "Configural", metric = "Metric", scalar = "Scalar", strict = "Strict")
  }

  display_table <- fit_table[, c(
    "chi2", "df", "p", "cfi", "delta_cfi",
    "rmsea", "delta_rmsea", "srmr",
    "decision", "status"
  )]

  rownames(display_table) <- unname(display_levels[rownames(display_table)])
  display_table[] <- lapply(display_table, function(column) {
    if (is.numeric(column)) {
      round(column, 3)
    } else {
      column
    }
  })

  names(display_table) <- if (language == "esp") {
    c("chi2", "gl", "p", "CFI", "\u0394CFI", "RMSEA", "\u0394RMSEA", "SRMR", "Decision", "Estado")
  } else {
    c("chi2", "df", "p", "CFI", "\u0394CFI", "RMSEA", "\u0394RMSEA", "SRMR", "Decision", "Status")
  }

  decision_col <- if (language == "esp") "Decision" else "Decision"
  status_col <- if (language == "esp") "Estado" else "Status"
  display_table[[decision_col]] <- vapply(display_table[[decision_col]], .ifa_decision_label, character(1), language = language)
  display_table[[status_col]] <- vapply(display_table[[status_col]], .ifa_status_label, character(1), language = language)

  notes <- character(0)
  if (estimator_info$uses_ordered) {
    notes <- c(
      notes,
      if (language == "esp") {
        "Con indicadores ordinales, la funcion usa parameterization = 'delta' y aplica: metrica = loadings, escalar = loadings + intercepts, estricta = loadings + intercepts + residuals."
      } else {
        "With ordinal indicators, the function uses parameterization = 'delta' and applies: metric = loadings, scalar = loadings + intercepts, strict = loadings + intercepts + residuals."
      },
      .ifa_ordered_data_issues(analysis_data, indicators, group_name, language)
    )
  }
  if (sample_info$dropped_missing_group > 0L) {
    notes <- c(
      notes,
      if (language == "esp") {
        paste0("Se eliminaron ", sample_info$dropped_missing_group, " casos con grupo faltante antes del analisis.")
      } else {
        paste0(sample_info$dropped_missing_group, " cases with missing group membership were removed before the analysis.")
      }
    )
  }
  if (!is.null(ordered)) {
    notes <- c(
      notes,
      if (language == "esp") {
        paste0("Argumento `ordered` fijado por el usuario en ", ordered, ".")
      } else {
        paste0("`ordered` argument set by the user to ", ordered, ".")
      }
    )
  }
  notes <- c(notes, .ifa_fit_issues(fit_table, language))

  report_path <- .ifa_prepare_report_path(report, language)
  summary_df <- .ifa_summary_sheet(
    estimator_info = estimator_info,
    estimator_reason = estimator_reason,
    item_info = item_info,
    sample_info = sample_info,
    normality_info = normality_info,
    scenario = scenario,
    conclusion = conclusion,
    report_path = report_path,
    language = language,
    notes = notes
  )

  export_df <- .ifa_export_table(display_table, language)
  written_report <- if (!is.null(report_path)) {
    .ifa_write_excel(summary_df, export_df, report_path)
  } else {
    NULL
  }

  report_text <- .ifa_report_text(
    table_df = display_table,
    estimator_info = estimator_info,
    estimator_reason = estimator_reason,
    sample_info = sample_info,
    item_info = item_info,
    normality_info = normality_info,
    scenario = scenario,
    conclusion = conclusion,
    language = language,
    report_path = written_report,
    ordered_note = estimator_info$uses_ordered,
    notes = notes
  )

  cat(report_text, sep = "\n")

  invisible(structure(
    list(
      estimator = estimator_info$estimator,
      estimator_info = estimator_info,
      estimator_reason = estimator_reason,
      indicators = indicators,
      item_info = item_info,
      sample_info = sample_info,
      normality = normality_info,
      requested_estimator = estimator,
      requested_ordered = ordered,
      chen_scenario = scenario,
      fitted_models = lapply(fitted_models, function(x) x$fit),
      fit_table = fit_table,
      display_table = display_table,
      conclusion = conclusion,
      report_text = report_text,
      report_path = written_report
    ),
    class = "factorial_invariance_auto"
  ))
}
