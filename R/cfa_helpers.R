# =============================================================================
# FUNCIONES AUXILIARES (internas)
# =============================================================================

# Normalizar idioma del reporte
.normalize_language <- function(language = "eng") {
  if (is.null(language) || length(language) == 0) return("eng")

  language <- tolower(trimws(as.character(language[1])))

  if (language %in% c("esp", "es", "spa", "spanish", "espanol")) {
    return("esp")
  }

  if (language %in% c("eng", "en", "english", "ingles")) {
    return("eng")
  }

  stop("Argument 'language' must be one of: 'esp' or 'eng'.")
}

# Traducir tipo de modelo
.model_type_label <- function(model_type_code, language = "eng") {
  language <- .normalize_language(language)

  if (language == "eng") {
    return(switch(
      model_type_code,
      bifactor = "Bifactor",
      second_order = "Second-order",
      unidimensional = "Unidimensional",
      correlated = "Correlated factors",
      model_type_code
    ))
  }

  switch(
    model_type_code,
    bifactor = "Bifactor",
    second_order = "Segundo orden",
    unidimensional = "Unidimensional",
    correlated = "Factores correlacionados",
    model_type_code
  )
}

# Etiquetas de indices de ajuste
.fit_index_base_label <- function(metric, language = "eng") {
  language <- .normalize_language(language)

  if (language == "eng") {
    return(switch(
      metric,
      chisq = "Chi-square",
      df = "df",
      pvalue = "p (chi-square)",
      chisq_ratio = "Chi-square/df",
      cfi = "CFI",
      tli = "TLI",
      rmsea = "RMSEA",
      "rmsea.ci.lower" = "RMSEA lower CI",
      "rmsea.ci.upper" = "RMSEA upper CI",
      "rmsea.pvalue" = "RMSEA p-close",
      srmr = "SRMR",
      error = "Error",
      metric
    ))
  }

  switch(
    metric,
    chisq = "Chi-cuadrado",
    df = "gl",
    pvalue = "p (chi-sq)",
    chisq_ratio = "Chi-sq/gl",
    cfi = "CFI",
    tli = "TLI",
    rmsea = "RMSEA",
    "rmsea.ci.lower" = "RMSEA IC inferior",
    "rmsea.ci.upper" = "RMSEA IC superior",
    "rmsea.pvalue" = "RMSEA p-close",
    srmr = "SRMR",
    error = "Error",
    metric
  )
}

.fit_index_label <- function(metric, variant = "base", language = "eng") {
  base_label <- .fit_index_base_label(metric, language)

  if (metric %in% c("df", "chisq_ratio", "srmr", "error")) {
    return(base_label)
  }

  switch(
    variant,
    robust = paste0(base_label, " (robust)"),
    scaled = paste0(base_label, " (scaled)"),
    base_label
  )
}

# Etiquetas de impresion del reporte
.cfa_auto_labels <- function(language = "eng") {
  language <- .normalize_language(language)

  if (language == "eng") {
    return(list(
      report_title = "       CONFIRMATORY FACTOR ANALYSIS - PSYCHOMETRIC REPORT",
      warnings_title = "*** WARNINGS ***",
      model_info_title = "MODEL INFORMATION",
      model_type_label = "  Model type           :",
      factors_label = "  Factors              :",
      observed_items_label = "  Observed items       :",
      n_total_label = "  N (total)            :",
      n_complete_label = "  N (complete)         :",
      missing_suffix = function(pct) sprintf("(%.1f%% with missing data)", pct),
      categories_label = "  Response categories  :",
      estimator_label = "  Estimator            :",
      justification_label = "  Rationale            :",
      converged_label = "  Converged            :",
      yes = "Yes",
      no = "No",
      mardia_title = "MULTIVARIATE NORMALITY (Mardia)",
      skew_label = "  Skewness  : Statistic",
      kurt_label = "  Kurtosis  : Statistic",
      rejected = "[Rejected]",
      not_rejected = "[Not rejected]",
      fit_title = "FIT INDICES",
      note_label = "  Note:",
      loadings_title = "STANDARDIZED FACTOR LOADINGS",
      factor_label = "  Factor:",
      ci_label = "95% CI",
      second_order_loadings_title = "SECOND-ORDER LOADINGS",
      r2_title = "R-SQUARED BY ITEM",
      factor_corr_title = "FACTOR CORRELATIONS",
      reliability_title = "RELIABILITY",
      respec_title = "RESPECIFICATION SUGGESTIONS",
      respec_note = "  (Consider only with theoretical/substantive justification)",
      respec_none_title = "RESPECIFICATIONS",
      respec_none = "  No modification indices above the threshold were found.",
      references_title = "REFERENCES (APA 7)",
      references_note = "  (Only references relevant to the analysis performed are included)",
      footer = "lavaan object available in $fit | All indices in $all_fit_measures"
    ))
  }

  list(
    report_title = "       ANALISIS FACTORIAL CONFIRMATORIO - REPORTE PSICOMETRICO",
    warnings_title = "*** ADVERTENCIAS ***",
    model_info_title = "INFORMACION DEL MODELO",
    model_type_label = "  Tipo de modelo       :",
    factors_label = "  Factores             :",
    observed_items_label = "  N items observados   :",
    n_total_label = "  N (total)            :",
    n_complete_label = "  N (completos)        :",
    missing_suffix = function(pct) sprintf("(%.1f%% con datos faltantes)", pct),
    categories_label = "  Categorias de resp.  :",
    estimator_label = "  Estimador            :",
    justification_label = "  Justificacion        :",
    converged_label = "  Convergio            :",
    yes = "Si",
    no = "No",
    mardia_title = "NORMALIDAD MULTIVARIADA (Mardia)",
    skew_label = "  Asimetria  : Estadistico",
    kurt_label = "  Curtosis   : Estadistico",
    rejected = "[Rechazada]",
    not_rejected = "[No rechazada]",
    fit_title = "INDICES DE AJUSTE",
    note_label = "  Nota:",
    loadings_title = "CARGAS FACTORIALES ESTANDARIZADAS",
    factor_label = "  Factor:",
    ci_label = "IC 95%",
    second_order_loadings_title = "CARGAS DE SEGUNDO ORDEN",
    r2_title = "R-CUADRADO POR ITEM",
    factor_corr_title = "CORRELACIONES ENTRE FACTORES",
    reliability_title = "CONFIABILIDAD",
    respec_title = "SUGERENCIAS DE REESPECIFICACION",
    respec_note = "  (Solo considerar con justificacion teorica/sustantiva)",
    respec_none_title = "REESPECIFICACIONES",
    respec_none = "  No se encontraron indices de modificacion superiores al umbral.",
    references_title = "REFERENCIAS BIBLIOGRAFICAS (APA 7)",
    references_note = "  (Se incluyen solo las referencias pertinentes al analisis realizado)",
    footer = "Objeto lavaan accesible en $fit | Todos los indices en $all_fit_measures"
  )
}

# Evaluar indices de ajuste
.evaluar_fi <- function(valor, bueno, aceptable, direccion = ">=", language = "eng") {
  language <- .normalize_language(language)
  if (is.na(valor)) return("")

  adequate <- if (language == "eng") "Adequate" else "Adecuado"
  acceptable_label <- if (language == "eng") "Acceptable" else "Aceptable"
  poor <- if (language == "eng") "Poor" else "Inadecuado"

  if (direccion == ">=") {
    if (valor >= bueno) adequate
    else if (valor >= aceptable) acceptable_label
    else poor
  } else {
    if (valor <= bueno) adequate
    else if (valor <= aceptable) acceptable_label
    else poor
  }
}

.pm_mvn_available <- function() {
  suppressPackageStartupMessages(requireNamespace("MVN", quietly = TRUE))
}

.pm_mvn_call <- function(data, mvn_test, use_population, alpha, tidy, descriptives) {
  MVN::mvn(
    data,
    mvn_test = mvn_test,
    use_population = use_population,
    alpha = alpha,
    tidy = tidy,
    descriptives = descriptives
  )
}

.pm_bind_condition_tables <- function(...) {
  entries <- list(...)
  if (!length(entries)) return(.pm_empty_condition_table())
  codes <- names(entries)
  if (is.null(codes)) codes <- rep("condition", length(entries))
  codes[is.na(codes) | !nzchar(codes)] <- "condition"
  tables <- Map(
    function(value, code) .pm_conditions_to_table(value, code),
    entries,
    codes
  )
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}

.pm_mvn_column <- function(x, candidates) {
  normalized <- tolower(gsub("[^[:alnum:]]", "", names(x)))
  target <- tolower(gsub("[^[:alnum:]]", "", candidates))
  position <- match(target, normalized, nomatch = 0L)
  position <- position[position > 0L]
  if (!length(position)) return(NULL)
  x[[position[1L]]]
}

.pm_mvn_mardia <- function(data, alpha) {
  alpha <- .pm_validate_scalar_number(alpha, "alpha", lower = 0, upper = 1)
  if (!.pm_mvn_available()) {
    return(list(
      status = "not_available",
      reason = "mvn_not_installed",
      result = NULL,
      condition = NULL,
      normal = NA,
      skew_stat = NA_real_, kurt_stat = NA_real_,
      p_skew = NA_real_, p_kurt = NA_real_
    ))
  }

  captured <- .pm_capture_conditions(.pm_mvn_call(
    data = data,
    mvn_test = "mardia",
    use_population = TRUE,
    alpha = alpha,
    tidy = FALSE,
    descriptives = FALSE
  ))
  if (!is.null(captured$error)) {
    return(list(
      status = "not_evaluable",
      reason = "mvn_error",
      result = NULL,
      condition = captured$error,
      normal = NA,
      warnings = captured$warnings,
      messages = captured$messages,
      skew_stat = NA_real_, kurt_stat = NA_real_,
      p_skew = NA_real_, p_kurt = NA_real_
    ))
  }

  table <- captured$value[["multivariate_normality"]]
  if (!is.data.frame(table)) {
    return(list(
      status = "not_evaluable",
      reason = "mvn_output_unrecognized",
      result = captured$value,
      condition = NULL,
      normal = NA,
      warnings = captured$warnings,
      messages = captured$messages,
      skew_stat = NA_real_, kurt_stat = NA_real_,
      p_skew = NA_real_, p_kurt = NA_real_
    ))
  }
  labels <- .pm_mvn_column(table, c("Test", "test"))
  statistic <- .pm_mvn_column(table, c("Statistic", "statistic"))
  p_value <- .pm_mvn_column(table, c("p.value", "p value", "pvalue"))
  if (is.null(labels) || is.null(statistic) || is.null(p_value)) {
    return(list(
      status = "not_evaluable",
      reason = "mvn_output_unrecognized",
      result = captured$value,
      condition = NULL,
      normal = NA,
      warnings = captured$warnings,
      messages = captured$messages,
      skew_stat = NA_real_, kurt_stat = NA_real_,
      p_skew = NA_real_, p_kurt = NA_real_
    ))
  }
  normalized_labels <- tolower(trimws(as.character(labels)))
  skew <- grep("mardia.*skew|skew.*mardia", normalized_labels)[1L]
  kurt <- grep("mardia.*kurt|kurt.*mardia", normalized_labels)[1L]
  if (is.na(skew) || is.na(kurt)) {
    status <- "not_evaluable"
    reason <- "mvn_tests_unrecognized"
  } else {
    status <- "ok"
    reason <- NA_character_
  }

  skew_stat <- if (is.na(skew)) NA_real_ else as.numeric(statistic[skew])
  kurt_stat <- if (is.na(kurt)) NA_real_ else as.numeric(statistic[kurt])
  p_skew <- if (is.na(skew)) NA_real_ else as.numeric(p_value[skew])
  p_kurt <- if (is.na(kurt)) NA_real_ else as.numeric(p_value[kurt])
  normal <- if (identical(status, "ok") && is.finite(p_skew) && is.finite(p_kurt)) {
    p_skew >= alpha && p_kurt >= alpha
  } else {
    NA
  }

  list(
    status = status,
    reason = reason,
    result = captured$value,
    condition = NULL,
    normal = normal,
    warnings = captured$warnings,
    messages = captured$messages,
    skew_stat = skew_stat,
    kurt_stat = kurt_stat,
    p_skew = p_skew,
    p_kurt = p_kurt
  )
}

.pm_fit_measure_value <- function(fit, name) {
  if (!name %in% names(fit)) return(NA_real_)
  value <- suppressWarnings(as.numeric(fit[[name]])[1L])
  if (!is.finite(value)) NA_real_ else value
}

.pm_pick_fit_variant <- function(fit, base, priority = c("robust", "scaled", "base")) {
  candidates <- c(robust = paste0(base, ".robust"),
                  scaled = paste0(base, ".scaled"),
                  base = base)
  for (variant in priority) {
    value <- .pm_fit_measure_value(fit, candidates[[variant]])
    if (is.finite(value)) return(list(value = value, variant = variant))
  }
  list(value = NA_real_, variant = "base")
}

.pm_omega_from_matrices <- function(lambda, psi, theta) {
  lambda <- as.matrix(lambda)
  psi <- as.matrix(psi)
  theta <- as.matrix(theta)
  if (!is.numeric(lambda) || !is.numeric(psi) || !is.numeric(theta) ||
      ncol(lambda) != nrow(psi) || nrow(psi) != ncol(psi) ||
      nrow(lambda) != nrow(theta) || nrow(theta) != ncol(theta) ||
      any(!is.finite(lambda)) || any(!is.finite(psi)) || any(!is.finite(theta))) {
    return(NA_real_)
  }
  common <- lambda %*% psi %*% t(lambda)
  denominator <- sum(common + theta)
  if (!is.finite(denominator) || denominator <= sqrt(.Machine$double.eps)) {
    return(NA_real_)
  }
  value <- sum(common) / denominator
  if (!is.finite(value) || value < 0 || value > 1) NA_real_ else value
}

.pm_lav_matrix_sets <- function(fit, what = "std") {
  matrices <- lavaan::lavInspect(fit, what)
  if (all(c("lambda", "theta", "psi") %in% names(matrices))) {
    return(stats::setNames(list(matrices), "group1"))
  }
  keep <- vapply(
    matrices,
    function(x) is.list(x) && all(c("lambda", "theta", "psi") %in% names(x)),
    logical(1)
  )
  matrices[keep]
}

.pm_cfa_admissibility <- function(fit,
                                  standardized_solution = NULL,
                                  converged = NULL,
                                  tolerance = sqrt(.Machine$double.eps)) {
  if (is.null(converged)) {
    converged <- tryCatch(
      lavaan::lavInspect(fit, "converged"),
      error = function(e) FALSE
    )
  }
  post_check <- tryCatch(
    suppressWarnings(lavaan::lavInspect(fit, "post.check")),
    error = function(e) NA
  )
  parameter_table <- tryCatch(
    lavaan::parameterEstimates(fit),
    error = function(e) NULL
  )
  if (is.null(standardized_solution)) {
    standardized_solution <- tryCatch(
      lavaan::standardizedSolution(fit),
      error = function(e) NULL
    )
  }

  negative_variances <- if (is.data.frame(parameter_table)) {
    parameter_table[
      parameter_table$op == "~~" &
        parameter_table$lhs == parameter_table$rhs &
        is.finite(parameter_table$est) & parameter_table$est < -tolerance,
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }
  heywood_loadings <- if (is.data.frame(standardized_solution)) {
    standardized_solution[
      standardized_solution$op == "=~" &
        is.finite(standardized_solution$est.std) &
        abs(standardized_solution$est.std) > 1 + tolerance,
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }

  reason_codes <- character()
  if (!isTRUE(converged)) reason_codes <- c(reason_codes, "not_converged")
  if (!isTRUE(post_check)) reason_codes <- c(reason_codes, "post_check_failed")
  if (nrow(negative_variances)) {
    reason_codes <- c(reason_codes, "negative_variance")
  }
  if (nrow(heywood_loadings)) {
    reason_codes <- c(reason_codes, "heywood_loading")
  }
  reason_codes <- unique(reason_codes)
  admissible <- length(reason_codes) == 0L

  list(
    admissible = admissible,
    status = if (admissible) "ok" else "not_evaluable",
    reason_codes = reason_codes,
    post_check = post_check,
    negative_variances = negative_variances,
    heywood_loadings = heywood_loadings
  )
}

.pm_detect_bifactor_structure <- function(items_per_factor,
                                           observed_items,
                                           is_second_order = FALSE) {
  observed_items <- unique(as.character(observed_items))
  factor_names <- names(items_per_factor)
  if (!is.list(items_per_factor) || is.null(factor_names)) {
    stop("`items_per_factor` must be a named list.", call. = FALSE)
  }
  observed_membership <- lapply(
    items_per_factor,
    function(items) intersect(unique(as.character(items)), observed_items)
  )
  general_candidates <- factor_names[vapply(
    observed_membership,
    function(items) length(items) == length(observed_items) &&
      setequal(items, observed_items),
    logical(1)
  )]
  empty <- list(
    is_bifactor = FALSE,
    code = "ambiguous_general_structure",
    general_candidates = general_candidates,
    general_factor = NULL,
    specific_factors = NULL,
    specialized_indices_supported = FALSE
  )
  if (isTRUE(is_second_order)) {
    empty$code <- "not_applicable_second_order"
    return(empty)
  }
  if (length(general_candidates) != 1L) return(empty)

  general <- general_candidates[[1L]]
  specifics <- setdiff(factor_names, general)
  n_observed <- vapply(observed_membership[specifics], length, integer(1))
  if (length(specifics) < 2L || any(n_observed == 0L) ||
      any(n_observed >= length(observed_items))) {
    empty$code <- "unsupported_specific_structure"
    return(empty)
  }

  list(
    is_bifactor = TRUE,
    code = "ok",
    general_candidates = general_candidates,
    general_factor = general,
    specific_factors = specifics,
    specialized_indices_supported = TRUE
  )
}

.pm_validate_puc_structure <- function(membership, general, orthogonal = TRUE) {
  fail <- function(code) list(supported = FALSE, code = code, puc = NA_real_)
  if (!is.list(membership) || is.null(names(membership)) ||
      !general %in% names(membership)) {
    return(fail("missing_general_factor"))
  }
  if (!isTRUE(orthogonal)) return(fail("nonorthogonal_factors"))
  general_items <- unique(as.character(membership[[general]]))
  specifics <- membership[setdiff(names(membership), general)]
  specific_items <- lapply(specifics, function(x) unique(as.character(x)))
  if (!length(general_items) ||
      !all(unique(unlist(specific_items, use.names = FALSE)) %in% general_items)) {
    return(fail("missing_general_loadings"))
  }
  if (length(specific_items) > 1L) {
    all_specific <- unlist(specific_items, use.names = FALSE)
    if (anyDuplicated(all_specific)) return(fail("unsupported_overlap"))
  }
  total_pairs <- length(general_items) * (length(general_items) - 1) / 2
  if (total_pairs <= 0) return(fail("insufficient_items"))
  contaminated <- sum(vapply(
    specific_items,
    function(items) length(items) * (length(items) - 1) / 2,
    numeric(1)
  ))
  list(
    supported = TRUE,
    code = "ok",
    puc = (total_pairs - contaminated) / total_pairs
  )
}

.pm_structural_membership <- function(fit, factors) {
  pt <- lavaan::parTable(fit)
  rows <- pt$op == "=~" & pt$lhs %in% factors
  split(as.character(pt$rhs[rows]), factor(pt$lhs[rows], levels = factors))
}

.pm_bifactor_is_orthogonal <- function(fit, factors) {
  pt <- lavaan::parTable(fit)
  pairs <- utils::combn(factors, 2L, simplify = FALSE)
  groups <- unique(pt$group)
  for (group in groups) {
    group_pt <- pt[pt$group == group & pt$op == "~~", , drop = FALSE]
    for (pair in pairs) {
      row <- group_pt[
        (group_pt$lhs == pair[1L] & group_pt$rhs == pair[2L]) |
          (group_pt$lhs == pair[2L] & group_pt$rhs == pair[1L]),
        , drop = FALSE
      ]
      if (nrow(row) != 1L || row$free[1L] != 0L ||
          !is.finite(row$ustart[1L]) || abs(row$ustart[1L]) > sqrt(.Machine$double.eps)) {
        return(FALSE)
      }
    }
  }
  TRUE
}

.pm_semtools_long <- function(x,
                              coefficient,
                              group_labels = character(),
                              level = "L1") {
  empty <- function() {
    data.frame(
      group = character(),
      level = character(),
      factor = character(),
      coefficient = character(),
      value = numeric(),
      stringsAsFactors = FALSE
    )
  }
  coefficient <- as.character(coefficient)[1L]
  level <- as.character(level)[1L]
  group_labels <- as.character(group_labels)
  group_labels <- group_labels[!is.na(group_labels) & nzchar(group_labels)]
  make_rows <- function(group, factor, value) {
    n <- max(length(group), length(factor), length(value))
    if (n == 0L) return(empty())
    data.frame(
      group = rep_len(as.character(group), n),
      level = rep(level, n),
      factor = rep_len(as.character(factor), n),
      coefficient = rep(coefficient, n),
      value = suppressWarnings(as.numeric(rep_len(value, n))),
      stringsAsFactors = FALSE
    )
  }
  if (is.null(x) || length(x) == 0L) return(empty())

  if (is.matrix(x)) {
    frame <- as.data.frame(x, stringsAsFactors = FALSE)
    matrix_groups <- rownames(x)
    if (is.null(matrix_groups) ||
        identical(matrix_groups, as.character(seq_len(nrow(x))))) {
      matrix_groups <- if (length(group_labels) == nrow(x)) {
        group_labels
      } else {
        paste0("group", seq_len(nrow(x)))
      }
    }
    frame$group <- matrix_groups
    x <- frame
  }

  if (is.data.frame(x)) {
    normalized_names <- tolower(gsub("[^[:alnum:]]", "", names(x)))
    group_position <- match("group", normalized_names, nomatch = 0L)
    factor_position <- match("factor", normalized_names, nomatch = 0L)
    value_position <- match("value", normalized_names, nomatch = 0L)
    if (factor_position > 0L && value_position > 0L) {
      groups <- if (group_position > 0L) {
        x[[group_position]]
      } else if (length(group_labels) == nrow(x)) {
        group_labels
      } else {
        rep("group1", nrow(x))
      }
      return(make_rows(groups, x[[factor_position]], x[[value_position]]))
    }

    groups <- if (group_position > 0L) {
      as.character(x[[group_position]])
    } else if (length(group_labels) == nrow(x)) {
      group_labels
    } else if (!is.null(rownames(x)) &&
               !identical(rownames(x), as.character(seq_len(nrow(x))))) {
      rownames(x)
    } else {
      rep("group1", nrow(x))
    }
    factor_positions <- setdiff(seq_along(x), group_position)
    rows <- lapply(factor_positions, function(position) {
      make_rows(groups, names(x)[[position]], x[[position]])
    })
    if (!length(rows)) return(empty())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  if (is.list(x)) {
    list_names <- names(x)
    if (is.null(list_names)) list_names <- rep("", length(x))
    group_oriented <- length(group_labels) > 0L &&
      length(list_names) > 0L && all(list_names %in% group_labels)
    rows <- list()
    if (group_oriented) {
      for (i in seq_along(x)) {
        values <- unlist(x[[i]], use.names = TRUE)
        factors <- names(values)
        if (is.null(factors) || any(!nzchar(factors))) {
          factors <- paste0("factor", seq_along(values))
        }
        rows[[length(rows) + 1L]] <- make_rows(
          list_names[[i]], factors, values
        )
      }
    } else {
      for (i in seq_along(x)) {
        values <- unlist(x[[i]], use.names = TRUE)
        groups <- names(values)
        if (is.null(groups) || any(!nzchar(groups))) {
          groups <- if (length(group_labels) == length(values)) {
            group_labels
          } else if (length(values) == 1L) {
            if (length(group_labels) == 1L) group_labels else "group1"
          } else {
            paste0("group", seq_along(values))
          }
        }
        factor <- if (nzchar(list_names[[i]])) {
          list_names[[i]]
        } else {
          paste0("factor", i)
        }
        rows[[length(rows) + 1L]] <- make_rows(groups, factor, values)
      }
    }
    if (!length(rows)) return(empty())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  values <- suppressWarnings(as.numeric(x))
  factors <- names(x)
  if (is.null(factors) || any(!nzchar(factors))) {
    factors <- if (length(values) == 1L) "total" else paste0("factor", seq_along(values))
  }
  make_rows("group1", factors, values)
}

.pm_standard_reliability_values <- function(fit) {
  if (!.pm_semtools_available()) {
    return(list(status = "not_available", reason = "semtools_not_installed"))
  }
  cr <- .pm_capture_conditions(semTools::compRelSEM(fit))
  ave <- .pm_capture_conditions(semTools::AVE(fit))
  if (!is.null(cr$error) || !is.null(ave$error)) {
    return(list(
      status = "not_evaluable",
      reason = "semtools_error",
      condition = if (!is.null(cr$error)) cr$error else ave$error
    ))
  }
  group_labels <- tryCatch(
    lavaan::lavInspect(fit, "group.label"),
    error = function(e) character()
  )
  rows <- rbind(
    .pm_semtools_long(
      cr$value, "composite_reliability", group_labels = group_labels,
      level = "L1"
    ),
    .pm_semtools_long(
      ave$value, "ave", group_labels = group_labels, level = "L1"
    )
  )
  alpha <- .pm_cfa_alpha_rows(fit)
  rows <- rbind(rows, alpha$rows)
  valid <- is.finite(rows$value) & rows$value >= 0 & rows$value <= 1
  valid[rows$coefficient %in% c("alpha_observed", "alpha_ordinal")] <-
    is.finite(rows$value[rows$coefficient %in% c("alpha_observed", "alpha_ordinal")]) &
    rows$value[rows$coefficient %in% c("alpha_observed", "alpha_ordinal")] <= 1
  rows$status <- ifelse(valid, "ok", "not_evaluable")
  rows$reason <- ifelse(
    valid,
    NA_character_,
    ifelse(is.finite(rows$value),
           "outside_mathematical_domain", "non_finite_coefficient")
  )
  rows$value[!valid] <- NA_real_

  values <- NULL
  if (nrow(rows) && length(unique(rows$group)) == 1L) {
    factors <- unique(rows$factor)
    values <- matrix(
      NA_real_, nrow = 3L, ncol = length(factors),
      dimnames = list(c("omega", "alpha", "avevar"), factors)
    )
    omega_rows <- rows[rows$coefficient == "composite_reliability", , drop = FALSE]
    ave_rows <- rows[rows$coefficient == "ave", , drop = FALSE]
    values["omega", omega_rows$factor] <- omega_rows$value
    values["avevar", ave_rows$factor] <- ave_rows$value
    alpha_rows <- rows[rows$coefficient == "alpha_observed", , drop = FALSE]
    values["alpha", alpha_rows$factor] <- alpha_rows$value
  }
  list(
    status = if (nrow(rows)) "ok" else "not_evaluable",
    reason = if (nrow(rows)) NA_character_ else "semtools_output_unrecognized",
    values = values,
    rows = rows,
    compRelSEM = cr$value, AVE = ave$value,
    warnings = c(cr$warnings, ave$warnings, alpha$warnings),
    messages = c(cr$messages, ave$messages, alpha$messages),
    alpha_errors = alpha$errors
  )
}

# Confiabilidad estandar (unidimensional o factores correlacionados)
.compute_standard_reliability <- function(fit) {
  value <- .pm_standard_reliability_values(fit)
  c(list(type = "standard"), value)
}

.pm_second_order_composite <- function(lambda,
                                        beta,
                                        phi,
                                        psi_u,
                                        theta,
                                        weights = rep(1, nrow(lambda))) {
  lambda <- as.matrix(lambda)
  beta <- as.matrix(beta)
  phi <- as.matrix(phi)
  psi_u <- as.matrix(psi_u)
  theta <- as.matrix(theta)
  weights <- as.numeric(weights)
  valid_dimensions <- nrow(lambda) > 0L && ncol(lambda) > 0L &&
    nrow(beta) == ncol(lambda) && ncol(beta) == nrow(phi) &&
    nrow(phi) == ncol(phi) && nrow(psi_u) == ncol(lambda) &&
    ncol(psi_u) == ncol(lambda) && nrow(theta) == nrow(lambda) &&
    ncol(theta) == nrow(lambda) && length(weights) == nrow(lambda)
  valid_values <- all(vapply(
    list(lambda, beta, phi, psi_u, theta, weights),
    function(value) is.numeric(value) && all(is.finite(value)),
    logical(1)
  ))
  fail <- function(reason, raw_omega = NA_real_, raw_ave = NA_real_,
                   components = NULL) {
    list(
      status = "not_evaluable",
      reason = reason,
      estimand = "unit_weighted_observed_composite",
      weights = weights,
      omega = NA_real_,
      ave = NA_real_,
      raw_omega = raw_omega,
      raw_ave = raw_ave,
      components = components
    )
  }
  if (!valid_dimensions) return(fail("invalid_matrix_dimensions"))
  if (!valid_values) return(fail("non_finite_matrix_input"))

  higher_order <- lambda %*% beta %*% phi %*% t(beta) %*% t(lambda)
  first_order_disturbance <- lambda %*% psi_u %*% t(lambda)
  total <- higher_order + first_order_disturbance + theta
  components <- list(
    higher_order = higher_order,
    first_order_disturbance = first_order_disturbance,
    residual = theta,
    total = total
  )
  total_composite_variance <- as.numeric(crossprod(weights, total %*% weights))
  raw_omega <- if (is.finite(total_composite_variance) &&
                   total_composite_variance > sqrt(.Machine$double.eps)) {
    as.numeric(crossprod(weights, higher_order %*% weights)) /
      total_composite_variance
  } else {
    NA_real_
  }
  ave_denominator <- sum(diag(total))
  raw_ave <- if (is.finite(ave_denominator) &&
                 ave_denominator > sqrt(.Machine$double.eps)) {
    sum(diag(higher_order)) / ave_denominator
  } else {
    NA_real_
  }
  if (!is.finite(raw_omega) || !is.finite(raw_ave)) {
    return(fail("non_positive_total_variance", raw_omega, raw_ave, components))
  }
  if (raw_omega < 0 || raw_omega > 1 || raw_ave < 0 || raw_ave > 1) {
    return(fail("outside_mathematical_domain", raw_omega, raw_ave, components))
  }

  list(
    status = "ok",
    reason = NA_character_,
    estimand = "unit_weighted_observed_composite",
    weights = weights,
    omega = raw_omega,
    ave = raw_ave,
    raw_omega = raw_omega,
    raw_ave = raw_ave,
    components = components
  )
}

.pm_second_order_matrix_rows <- function(fit, so_factors) {
  sets <- .pm_lav_matrix_sets(fit, "std")
  rows <- list()
  for (group in names(sets)) {
    matrices <- sets[[group]]
    beta <- matrices$beta
    if (is.null(beta)) next
    for (factor in so_factors) {
      if (!factor %in% colnames(beta)) next
      loadings <- beta[, factor]
      first_order <- names(loadings)[
        is.finite(loadings) & abs(loadings) > sqrt(.Machine$double.eps)
      ]
      if (!length(first_order) ||
          !all(first_order %in% colnames(matrices$lambda)) ||
          !all(first_order %in% rownames(matrices$psi)) ||
          !factor %in% rownames(matrices$psi)) {
        next
      }
      lambda_all <- as.matrix(matrices$lambda[, first_order, drop = FALSE])
      downstream <- rowSums(abs(lambda_all) > sqrt(.Machine$double.eps)) > 0L
      if (!any(downstream)) next
      item_names <- rownames(lambda_all)[downstream]
      lambda <- lambda_all[downstream, , drop = FALSE]
      b <- beta[first_order, factor, drop = FALSE]
      phi <- matrices$psi[factor, factor, drop = FALSE]
      disturbances <- matrices$psi[first_order, first_order, drop = FALSE]
      theta <- matrices$theta[item_names, item_names, drop = FALSE]
      composite <- .pm_second_order_composite(
        lambda = lambda,
        beta = b,
        phi = phi,
        psi_u = disturbances,
        theta = theta,
        weights = rep(1, length(item_names))
      )
      values <- c(composite$omega, composite$ave)
      rows[[length(rows) + 1L]] <- data.frame(
        group = group,
        level = "L2",
        factor = factor,
        metric = c("composite_reliability", "ave"),
        value = values,
        status = rep(composite$status, 2L),
        reason = rep(composite$reason, 2L),
        estimand = c(
          composite$estimand,
          "mean_observed_variance_due_to_higher_order_factor"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(
      group = character(), level = character(), factor = character(),
      metric = character(), value = numeric(), status = character(),
      reason = character(), estimand = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

# Canonical group identity helpers. Missing and empty labels are distinct
# identities; neither is allowed to behave as a wildcard during row selection.
.pm_group_key <- function(group) {
  group <- as.character(group)
  if (!length(group)) return(character())
  out <- paste0("value:", group)
  out[is.na(group)] <- "missing"
  out[!is.na(group) & group == ""] <- "empty"
  out
}

.pm_group_label <- function(group) {
  group <- as.character(group)
  out <- group
  labelled <- !is.na(group) & nzchar(group)
  out[labelled] <- gsub("\\", "\\\\", group[labelled], fixed = TRUE)
  reserved <- !is.na(group) & group %in% c("<missing>", "<empty>")
  out[reserved] <- paste0("\\", out[reserved])
  out[is.na(group)] <- "<missing>"
  out[!is.na(group) & group == ""] <- "<empty>"
  out
}

.pm_group_matches <- function(values, group) {
  values <- as.character(values)
  group <- as.character(group)[1L]
  if (is.na(group)) {
    is.na(values)
  } else if (identical(group, "")) {
    !is.na(values) & values == ""
  } else {
    !is.na(values) & values == group
  }
}

# Confiabilidad para segundo orden
.compute_second_order_reliability <- function(fit, so_factors) {
  l1 <- .pm_standard_reliability_values(fit)
  l1_rows <- data.frame(
    group = character(), level = character(), factor = character(),
    metric = character(), value = numeric(), status = character(),
    reason = character(), estimand = character(), stringsAsFactors = FALSE
  )
  if (identical(l1$status, "ok") && is.data.frame(l1$rows)) {
    l1_rows <- data.frame(
      group = l1$rows$group,
      level = l1$rows$level,
      factor = l1$rows$factor,
      metric = l1$rows$coefficient,
      value = l1$rows$value,
      status = l1$rows$status,
      reason = l1$rows$reason,
      estimand = ifelse(
        l1$rows$coefficient == "ave",
        "average_variance_extracted",
        "semtools_composite"
      ),
      stringsAsFactors = FALSE
    )
  }
  l2_rows <- .pm_second_order_matrix_rows(fit, so_factors)
  fit_groups <- tryCatch(
    as.character(lavaan::lavInspect(fit, "group.label")),
    error = function(e) character()
  )
  candidate_groups <- c(
    fit_groups,
    if (nrow(l1_rows)) as.character(l1_rows$group) else character(),
    if (nrow(l2_rows)) as.character(l2_rows$group) else character()
  )
  group_names <- candidate_groups[
    !duplicated(.pm_group_key(candidate_groups))
  ]
  factor_names <- unique(as.character(so_factors))
  metric_names <- c("composite_reliability", "ave")
  metric_estimands <- c(
    composite_reliability = "unit_weighted_observed_composite",
    ave = "mean_observed_variance_due_to_higher_order_factor"
  )

  expected_rows <- list()
  if (length(group_names) && length(factor_names)) {
    for (group in group_names) {
      group_rows <- l2_rows[
        .pm_group_matches(l2_rows$group, group), , drop = FALSE
      ]
      for (factor in factor_names) {
        factor_rows <- group_rows[
          !is.na(group_rows$factor) & group_rows$factor == factor, , drop = FALSE
        ]
        for (metric in metric_names) {
          observed <- factor_rows[
            !is.na(factor_rows$metric) & factor_rows$metric == metric,
            , drop = FALSE
          ]
          if (nrow(observed)) {
            row <- observed[1L, c(
              "group", "level", "factor", "metric", "value", "status",
              "reason", "estimand"
            ), drop = FALSE]
            row$group <- group
            row$level <- "L2"
            row$factor <- factor
            row$metric <- metric
          } else {
            missing_reason <- if (!nrow(group_rows)) {
              "second_order_group_not_evaluable"
            } else if (!nrow(factor_rows)) {
              "second_order_factor_not_evaluable"
            } else {
              "second_order_metric_not_evaluable"
            }
            row <- data.frame(
              group = group,
              level = "L2",
              factor = factor,
              metric = metric,
              value = NA_real_,
              status = "not_evaluable",
              reason = missing_reason,
              estimand = unname(metric_estimands[[metric]]),
              stringsAsFactors = FALSE
            )
          }
          expected_rows[[length(expected_rows) + 1L]] <- row
        }
      }
    }
  }
  if (length(expected_rows)) {
    l2_rows <- do.call(rbind, expected_rows)
    rownames(l2_rows) <- NULL
  } else {
    l2_rows <- l2_rows[FALSE, , drop = FALSE]
  }

  l2_cell_ok <- nrow(l2_rows) > 0L &
    !is.na(l2_rows$status) & l2_rows$status == "ok" &
    is.finite(l2_rows$value)
  l2_evaluable <- any(l2_cell_ok)
  l2_complete <- nrow(l2_rows) > 0L && all(l2_cell_ok)
  l2_status <- if (l2_complete) {
    "ok"
  } else if (l2_evaluable) {
    "partially_evaluable"
  } else {
    "not_evaluable"
  }
  l2_reason <- switch(
    l2_status,
    ok = NA_character_,
    partially_evaluable = "second_order_partially_evaluable",
    not_evaluable = "second_order_not_evaluable"
  )

  l2_legacy <- list()
  if (nrow(l2_rows)) {
    by_group <- lapply(group_names, function(group) {
      group_rows <- l2_rows[
        .pm_group_matches(l2_rows$group, group), , drop = FALSE
      ]
      by_factor <- lapply(factor_names, function(factor) {
        group_rows[
          !is.na(group_rows$factor) & group_rows$factor == factor,
          c("metric", "value", "status", "reason"), drop = FALSE
        ]
      })
      names(by_factor) <- factor_names
      by_factor
    })
    names(by_group) <- make.unique(.pm_group_label(group_names), sep = " #")

    # Preserve the historical factor-indexed alias for one-group fits. A
    # multigroup fit needs the additional group level so coefficients from
    # identically named factors cannot be silently collapsed.
    l2_legacy <- if (length(by_group) == 1L) by_group[[1L]] else by_group
  }
  list(
    type = "second_order",
    status = l2_status,
    reason = l2_reason,
    first_order = if (!is.null(l1$values)) l1$values else l1$rows,
    second_order = l2_legacy,
    so_factors = so_factors,
    rows = rbind(l1_rows, l2_rows),
    compRelSEM = l1$compRelSEM,
    AVE = l1$AVE,
    warnings = l1$warnings,
    messages = l1$messages,
    condition = l1$condition
  )
}

.pm_bifactor_group_reliability <- function(matrices,
                                            general_factor,
                                            specific_factors,
                                            puc_diagnostic) {
  lambda <- as.matrix(matrices$lambda)
  psi <- as.matrix(matrices$psi)
  theta <- as.matrix(matrices$theta)
  common <- lambda %*% psi %*% t(lambda)
  denominator <- sum(common + theta)
  omega_total <- .pm_omega_from_matrices(lambda, psi, theta)

  general_lambda <- lambda[, general_factor, drop = FALSE]
  general_psi <- psi[general_factor, general_factor, drop = FALSE]
  general_common <- general_lambda %*% general_psi %*% t(general_lambda)
  omega_h <- if (is.finite(denominator) && denominator > sqrt(.Machine$double.eps)) {
    sum(general_common) / denominator
  } else {
    NA_real_
  }
  if (!is.finite(omega_h) || omega_h < 0 || omega_h > 1) omega_h <- NA_real_

  general_loadings <- lambda[, general_factor]
  specific_loadings <- lambda[, specific_factors, drop = FALSE]
  ss_general <- sum(general_loadings^2)
  ss_specific <- sum(specific_loadings^2)
  ecv <- if (ss_general + ss_specific > 0) ss_general / (ss_general + ss_specific) else NA_real_
  ecv_denominator <- general_loadings^2 + rowSums(specific_loadings^2)
  ecv_item <- ifelse(ecv_denominator > 0, general_loadings^2 / ecv_denominator, NA_real_)

  omega_hs <- vapply(specific_factors, function(factor) {
    positions <- which(abs(lambda[, factor]) > sqrt(.Machine$double.eps))
    if (!length(positions)) return(NA_real_)
    sub_lambda <- lambda[positions, , drop = FALSE]
    sub_theta <- theta[positions, positions, drop = FALSE]
    sub_common <- sub_lambda %*% psi %*% t(sub_lambda)
    sub_denominator <- sum(sub_common + sub_theta)
    factor_lambda <- lambda[positions, factor, drop = FALSE]
    factor_common <- factor_lambda %*% psi[factor, factor, drop = FALSE] %*%
      t(factor_lambda)
    if (!is.finite(sub_denominator) || sub_denominator <= 0) return(NA_real_)
    value <- sum(factor_common) / sub_denominator
    if (!is.finite(value) || value < 0 || value > 1) NA_real_ else value
  }, numeric(1))

  list(
    omega_total = omega_total,
    omega_h = omega_h,
    omega_hs = omega_hs,
    ecv = ecv,
    ecv_item = ecv_item,
    puc = puc_diagnostic$puc,
    puc_diagnostic = puc_diagnostic,
    h_general = .compute_h_index(general_loadings[general_loadings != 0]),
    h_specific = vapply(specific_factors, function(factor) {
      values <- lambda[lambda[, factor] != 0, factor]
      .compute_h_index(values)
    }, numeric(1)),
    lambda = lambda,
    psi = psi,
    theta = theta
  )
}

# Confiabilidad para bifactor
.compute_bifactor_reliability <- function(fit, general_factor, specific_factors,
                                          observed_items) {
  factors <- c(general_factor, specific_factors)
  membership <- .pm_structural_membership(fit, factors)
  structural_orthogonality <- .pm_bifactor_is_orthogonal(fit, factors)
  puc_diagnostic <- .pm_validate_puc_structure(
    membership,
    general = general_factor,
    orthogonal = structural_orthogonality
  )
  sets <- .pm_lav_matrix_sets(fit, "std")
  groups <- lapply(
    sets,
    .pm_bifactor_group_reliability,
    general_factor = general_factor,
    specific_factors = specific_factors,
    puc_diagnostic = puc_diagnostic
  )

  semtools <- if (.pm_semtools_available()) {
    .pm_capture_conditions(semTools::compRelSEM(fit))
  } else {
    list(value = NULL, error = NULL, warnings = list(), messages = list())
  }
  out <- list(
    type = "bifactor",
    status = if (length(groups)) "ok" else "not_evaluable",
    reason = if (length(groups)) NA_character_ else "matrix_extraction_failed",
    semtools_rel = semtools$value,
    semtools_condition = semtools$error,
    warnings = semtools$warnings,
    messages = semtools$messages,
    groups = groups,
    puc_diagnostic = puc_diagnostic,
    general_factor = general_factor,
    specific_factors = specific_factors,
    observed_items = observed_items
  )
  if (length(groups) == 1L) {
    for (field in names(groups[[1L]])) out[[field]] <- groups[[1L]][[field]]
  }
  out
}

# H index (construct replicability)
.compute_h_index <- function(loadings) {
  loadings <- loadings[loadings != 0]
  if (length(loadings) == 0) return(NA_real_)
  1 / (1 + 1 / sum(loadings^2 / (1 - loadings^2)))
}

# Construir referencias bibliograficas contextuales (APA 7)
.build_references <- function(x) {
  refs <- character(0)

  # --- Siempre incluir: AFC general e indices de ajuste ---
  refs <- c(refs,
    paste0("Brown, T. A. (2015). Confirmatory factor analysis for applied research ",
           "(2nd ed.). Guilford Press."),
    paste0("Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in ",
           "covariance structure analysis: Conventional criteria versus new alternatives. ",
           "Structural Equation Modeling, 6(1), 1-55. ",
           "https://doi.org/10.1080/10705519909540118"),
    paste0("Kline, R. B. (2016). Principles and practice of structural equation ",
           "modeling (4th ed.). Guilford Press.")
  )

  # --- Estimador: WLSMV / ULSMV (datos ordinales) ---
  if (x$estimator %in% c("WLSMV", "ULSMV")) {
    refs <- c(refs,
      paste0("Flora, D. B., & Curran, P. J. (2004). An empirical evaluation of ",
             "alternative methods of estimation for confirmatory factor analysis with ",
             "ordinal data. Psychological Methods, 9(4), 466-491. ",
             "https://doi.org/10.1037/1082-989X.9.4.466"),
      paste0("Li, C.-H. (2016). Confirmatory factor analysis with ordinal data: ",
             "Comparing robust maximum likelihood and diagonally weighted least squares. ",
             "Behavior Research Methods, 48(3), 936-949. ",
             "https://doi.org/10.3758/s13428-015-0619-7")
    )
  }

  # ULSMV especificamente para muestras pequenas
  if (x$estimator == "ULSMV") {
    refs <- c(refs,
      paste0("Forero, C. G., Maydeu-Olivares, A., & Gallardo-Pujol, D. (2009). ",
             "Factor analysis with ordinal indicators: A Monte Carlo study comparing ",
             "DWLS and ULS estimation. Structural Equation Modeling, 16(4), 625-641. ",
             "https://doi.org/10.1080/10705510903203573")
    )
  }

  # --- Estimador: MLR (robusto) ---
  if (x$estimator == "MLR") {
    refs <- c(refs,
      paste0("Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics ",
             "and standard errors in covariance structure analysis. In A. von Eye & ",
             "C. C. Clogg (Eds.), Latent variables analysis: Applications for ",
             "developmental research (pp. 399-419). Sage.")
    )
  }

  # --- Indices scaled (comun a WLSMV, ULSMV, MLR) ---
  if (x$estimator %in% c("WLSMV", "ULSMV", "MLR")) {
    refs <- c(refs,
      paste0("Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square ",
             "test statistic for moment structure analysis. Psychometrika, 66(4), ",
             "507-514. https://doi.org/10.1007/BF02296192")
    )
  }

  # --- Normalidad multivariada (Mardia) ---
  if (!is.na(x$mardia$p_skew)) {
    refs <- c(refs,
      paste0("Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis ",
             "with applications. Biometrika, 57(3), 519-530. ",
             "https://doi.org/10.2307/2334770")
    )
  }

  # --- Confiabilidad ---
  rel <- x$reliability

  if (rel$type != "error") {
    if (x$estimator %in% c("WLSMV", "ULSMV")) {
      refs <- c(refs,
        "Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using structural equation modeling: An alternative to coefficient alpha. Psychometrika, 74(1), 155-167. https://doi.org/10.1007/s11336-008-9099-3",
        "Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of coefficients alpha and theta for Likert rating scales. Journal of Modern Applied Statistical Methods, 6(1), 21-29. https://doi.org/10.22237/jmasm/1177992180")
    }
    # Omega siempre
    refs <- c(refs,
      paste0("McDonald, R. P. (1999). Test theory: A unified treatment. ",
             "Lawrence Erlbaum Associates."),
      paste0("Viladrich, C., Angulo-Brunet, A., & Doval, E. (2017). A journey around ",
             "alpha and omega to estimate internal consistency reliability. Anales de ",
             "Psicologia, 33(3), 755-782. https://doi.org/10.6018/analesps.33.3.268401")
    )

    # Alfa si se reporta
    if (rel$type == "standard" || rel$type == "second_order") {
      refs <- c(refs,
        paste0("Cronbach, L. J. (1951). Coefficient alpha and the internal structure ",
               "of tests. Psychometrika, 16(3), 297-334. ",
               "https://doi.org/10.1007/BF02310555")
      )
    }

    # Bifactor: indices especializados
    if (rel$type == "bifactor") {
      refs <- c(refs,
        paste0("Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013). ",
               "Multidimensionality and structural coefficient bias in structural equation ",
               "modeling: A bifactor perspective. Educational and Psychological Measurement, ",
               "73(1), 5-26. https://doi.org/10.1177/0013164412449831"),
        paste0("Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating ",
               "bifactor models: Calculating and interpreting statistical indices. ",
               "Psychological Methods, 21(2), 137-150. ",
               "https://doi.org/10.1037/met0000045")
      )
    }

    # Segundo orden
    if (rel$type == "second_order") {
      refs <- c(refs,
        paste0("Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct ",
               "reliability within latent variable systems. In R. Cudeck, S. du Toit, ",
               "& D. Sorbom (Eds.), Structural equation modeling: Present and future ",
               "(pp. 195-216). Scientific Software International.")
      )
    }
  }

  # --- Indices de modificacion ---
  if (length(x$respecifications) > 0) {
    refs <- c(refs,
      paste0("Sorbom, D. (1989). Model modification. Psychometrika, 54(3), ",
             "371-384. https://doi.org/10.1007/BF02294623")
    )
  }

  # Eliminar duplicados y ordenar alfabeticamente
  refs <- sort(unique(refs))
  return(refs)
}

# Construir sugerencias de reespecificacion
.build_respec_suggestions <- function(mi_df, items_per_factor, language = "eng") {
  language <- .normalize_language(language)
  is_english <- identical(language, "eng")
  respec <- list()
  if (is.null(mi_df) || nrow(mi_df) == 0) return(respec)

  memberships <- lapply(items_per_factor, unique)
  all_items <- unique(unlist(memberships, use.names = FALSE))
  specific_memberships <- memberships
  if (length(memberships) > 1L) {
    general_like <- vapply(
      memberships,
      function(items) length(setdiff(all_items, items)) == 0L,
      logical(1)
    )
    if (any(general_like) && any(!general_like)) {
      specific_memberships <- memberships[!general_like]
    }
  }

  for (i in seq_len(nrow(mi_df))) {
    row <- mi_df[i, ]
    lhs <- as.character(row$lhs)
    op <- as.character(row$op)
    rhs <- as.character(row$rhs)
    syntax <- paste(lhs, op, rhs)
    same_factor <- if (identical(op, "~~") && !identical(lhs, rhs)) {
      any(vapply(
        specific_memberships,
        function(items) lhs %in% items && rhs %in% items,
        logical(1)
      ))
    } else {
      FALSE
    }
    type <- switch(
      op,
      "=~" = "cross_loading",
      "~~" = if (identical(lhs, rhs)) "variance" else "residual_covariance",
      "~" = if (identical(rhs, "1")) "intercept" else "regression",
      "|" = "threshold",
      "parameter"
    )
    theory_warning <- if (is_english) {
      "Consider only with substantive and methodologically defensible theory."
    } else {
      "Considerar solo con teoria sustantiva y metodologicamente defendible."
    }
    note <- if (identical(type, "residual_covariance") && same_factor) {
      if (is_english) {
        "Same specific factor; inspect wording or method effects before freeing the covariance."
      } else {
        "Mismo factor especifico; revisar redaccion o efectos de metodo antes de liberar la covarianza."
      }
    } else if (identical(type, "residual_covariance")) {
      if (is_english) {
        "Different specific factors; the residual covariance may indicate model misspecification."
      } else {
        "Factores especificos distintos; la covarianza residual puede indicar mala especificacion."
      }
    } else {
      theory_warning
    }

    item <- list(
      lhs = lhs,
      op = op,
      rhs = rhs,
      type = type,
      syntax = syntax,
      mi = as.numeric(row$mi),
      epc = as.numeric(row$epc),
      same_factor = same_factor,
      theory_warning = theory_warning,
      note = note
    )
    item$tipo <- item$type
    item$sintaxis <- item$syntax
    item$mismo_factor <- item$same_factor
    item$nota <- item$note
    respec[[length(respec) + 1L]] <- item
  }

  return(respec)
}


