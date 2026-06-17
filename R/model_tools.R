#' Fit basic IRT models automatically
#'
#' Fits a one-dimensional Rasch, 2PL, or graded response model using `mirt`.
#' The dependency is optional and is loaded only when this function is called.
#'
#' @param data Data frame or matrix containing item responses.
#' @param model IRT model: `"auto"`, `"rasch"`, `"2pl"`, or `"grm"`.
#' @param factors Number of latent dimensions passed to `mirt::mirt()`.
#' @param itemtype Optional `mirt` item type. If `NULL`, it is selected from
#'   `model`.
#' @param technical Optional technical arguments passed to `mirt::mirt()`.
#' @param language Interpretation language: `"eng"` or `"esp"`.
#' @param digits Number of decimal places used in interpretation tables.
#'
#' @return A list of class `psychomatic_irt` with fitted model, item
#'   parameters, item-level interpretation, a model-level summary, and
#'   supporting APA 7 references.
#'
#' @examples
#' if (identical(Sys.getenv("PSYCHOMATIC_RUN_SLOW_EXAMPLES"), "true") &&
#'     requireNamespace("mirt", quietly = TRUE)) {
#'   data(psychomatic_ordinal)
#'   irt_auto(psychomatic_ordinal[1:60, 1:4], model = "grm")
#' }
#'
#' @export
irt_auto <- function(data,
                     model = c("auto", "rasch", "2pl", "grm"),
                     factors = 1,
                     itemtype = NULL,
                     technical = list(),
                     language = c("eng", "esp"),
                     digits = 3) {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package `mirt` is required for `irt_auto()`.", call. = FALSE)
  }
  items <- .pm_prepare_numeric_items(data)
  model <- match.arg(model)
  language <- match.arg(language)
  if (!is.numeric(factors) || length(factors) != 1L || factors < 1) {
    stop("`factors` must be a positive integer.", call. = FALSE)
  }

  n_categories <- vapply(items, function(x) length(unique(stats::na.omit(x))), integer(1))
  if (model == "auto") {
    model <- if (max(n_categories) <= 2L) "rasch" else "grm"
  }
  if (is.null(itemtype)) {
    itemtype <- switch(
      model,
      rasch = "Rasch",
      `2pl` = "2PL",
      grm = "graded"
    )
  }

  fit <- mirt::mirt(
    data = items,
    model = as.integer(factors),
    itemtype = itemtype,
    verbose = FALSE,
    technical = technical
  )
  coefficients <- mirt::coef(fit, IRTpars = TRUE, simplify = TRUE)
  item_parameters <- .pm_mirt_item_parameters(coefficients)
  interpretation <- .pm_irt_interpretation(
    item_parameters,
    model = model,
    itemtype = itemtype,
    language = language,
    digits = digits
  )
  model_summary <- .pm_irt_model_summary(
    interpretation,
    model = model,
    itemtype = itemtype,
    n = nrow(items),
    language = language,
    digits = digits
  )

  result <- list(
    fit = fit,
    coefficients = coefficients,
    item_parameters = item_parameters,
    interpretation = interpretation,
    model_summary = model_summary,
    model = model,
    itemtype = itemtype,
    factors = as.integer(factors),
    language = language,
    n = nrow(items),
    n_items = ncol(items),
    references = .pm_irt_references(),
    call = match.call()
  )
  class(result) <- "psychomatic_irt"
  result
}

#' Screen differential item functioning by group
#'
#' Performs a lightweight DIF screening using logistic regression. Polytomous
#' items are dichotomized at their median for this screening, so flagged results
#' should be followed up with a dedicated IRT or ordinal DIF analysis.
#'
#' @param data Data frame or matrix containing item responses. If `group` is a
#'   column name in `data`, that column is used as the grouping variable.
#' @param group Grouping vector or the name of the grouping column in `data`.
#' @param items Optional item names. Defaults to all numeric item columns except
#'   the grouping column.
#' @param alpha Significance level used for DIF flags.
#' @param language Interpretation language: `"eng"` or `"esp"`.
#' @param digits Number of decimal places used in interpretation tables.
#'
#' @return A list of class `psychomatic_dif` with screening results,
#'   item-level interpretation (`dif_type`, decision, severity, and
#'   recommendation), an overall DIF summary, and supporting APA 7 references.
#'
#' @examples
#' set.seed(1)
#' group <- rep(c("A", "B"), each = 20)
#' theta <- rnorm(40)
#' items <- data.frame(
#'   item1 = as.integer(theta + rnorm(40) > 0),
#'   item2 = as.integer(theta + rnorm(40) > 0)
#' )
#' dif_auto(items, group = group)
#'
#' @export
dif_auto <- function(data,
                     group,
                     items = NULL,
                     alpha = 0.01,
                     language = c("eng", "esp"),
                     digits = 3) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }
  language <- match.arg(language)
  data <- as.data.frame(data)
  group_name <- NULL
  if (is.character(group) && length(group) == 1L && group %in% names(data)) {
    group_name <- group
    group_vector <- data[[group]]
  } else {
    group_vector <- group
  }
  group_factor <- factor(group_vector)
  if (length(group_factor) != nrow(data)) {
    stop("`group` must have the same length as the number of rows in `data`.",
         call. = FALSE)
  }
  if (nlevels(group_factor) < 2L) {
    stop("`group` must contain at least two groups.", call. = FALSE)
  }

  if (is.null(items)) {
    numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
    items <- setdiff(numeric_cols, group_name)
  }
  if (length(items) == 0L) {
    stop("No numeric items were available for DIF screening.", call. = FALSE)
  }
  missing_items <- setdiff(items, names(data))
  if (length(missing_items) > 0L) {
    stop("Items not found in `data`: ", paste(missing_items, collapse = ", "), call. = FALSE)
  }

  score_data <- data[items]
  total_score <- rowMeans(score_data, na.rm = TRUE)
  total_score[!is.finite(total_score)] <- NA_real_

  rows <- lapply(items, function(item) {
    x <- data[[item]]
    outcome_info <- .pm_dif_outcome(x)
    complete <- stats::complete.cases(outcome_info$y, total_score, group_factor)
    analysis_data <- data.frame(
      y = outcome_info$y[complete],
      score = total_score[complete],
      group = droplevels(group_factor[complete])
    )

    if (nrow(analysis_data) < 10L || length(unique(analysis_data$y)) < 2L ||
        nlevels(analysis_data$group) < 2L) {
      return(data.frame(
        item = item,
        p_uniform = NA_real_,
        p_nonuniform = NA_real_,
        flagged = FALSE,
        method = outcome_info$method,
        note = "Insufficient variation for logistic DIF screening",
        stringsAsFactors = FALSE
      ))
    }

    p_values <- tryCatch({
      base <- stats::glm(y ~ score, data = analysis_data, family = stats::binomial())
      uniform <- stats::glm(y ~ score + group, data = analysis_data, family = stats::binomial())
      nonuniform <- stats::glm(y ~ score * group, data = analysis_data, family = stats::binomial())
      list(
        uniform = stats::anova(base, uniform, test = "Chisq")$`Pr(>Chi)`[2],
        nonuniform = stats::anova(uniform, nonuniform, test = "Chisq")$`Pr(>Chi)`[2]
      )
    }, error = function(e) list(uniform = NA_real_, nonuniform = NA_real_))

    data.frame(
      item = item,
      p_uniform = p_values$uniform,
      p_nonuniform = p_values$nonuniform,
      flagged = isTRUE(!is.na(p_values$uniform) && p_values$uniform < alpha) ||
        isTRUE(!is.na(p_values$nonuniform) && p_values$nonuniform < alpha),
      method = outcome_info$method,
      note = outcome_info$note,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  interpretation <- .pm_dif_interpretation(
    results,
    alpha = alpha,
    language = language,
    digits = digits
  )
  summary <- .pm_dif_summary(
    interpretation,
    alpha = alpha,
    group_levels = levels(group_factor),
    language = language,
    digits = digits
  )
  out <- list(
    results = results,
    interpretation = interpretation,
    summary = summary,
    alpha = alpha,
    group_levels = levels(group_factor),
    language = language,
    references = .pm_dif_references(),
    call = match.call()
  )
  class(out) <- "psychomatic_dif"
  out
}

#' Compare alternative CFA models
#'
#' Fits a named set of lavaan CFA models and returns a compact comparison table.
#' The recommendation is descriptive only and favors converged models with
#' stronger comparative fit indices. Information criteria are reported when
#' available, but they are not required and are not used for WLSMV-style
#' categorical estimators.
#'
#' @param data Data frame containing observed variables.
#' @param models Named list or named character vector of lavaan model strings.
#' @param ordered Optional ordered variables passed to `lavaan::cfa()`. If
#'   `NULL` and a categorical estimator such as WLSMV is used, observed items
#'   referenced in `models` are treated as ordered.
#' @param estimator Estimator passed to `lavaan::cfa()`.
#' @param language Recommendation language: `"eng"` or `"esp"`.
#'
#' @return A list of class `psychomatic_model_comparison`.
#'
#' @examples
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   models <- list(one = "visual =~ x1 + x2 + x3")
#'   compare_models_auto(lavaan::HolzingerSwineford1939, models, estimator = "ML")
#' }
#'
#' @export
compare_models_auto <- function(data,
                                models,
                                ordered = NULL,
                                estimator = "ML",
                                language = c("eng", "esp")) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package `lavaan` is required for `compare_models_auto()`.", call. = FALSE)
  }
  language <- match.arg(language)
  if (is.character(models)) {
    models <- as.list(models)
  }
  if (!is.list(models) || length(models) == 0L || is.null(names(models)) ||
      any(!nzchar(names(models)))) {
    stop("`models` must be a named list or named character vector.", call. = FALSE)
  }
  estimator <- toupper(estimator)
  data <- as.data.frame(data)
  if (is.null(ordered) && estimator %in% c("WLSMV", "DWLS", "ULSMV", "ULS")) {
    ordered <- .pm_model_observed_variables(models, data)
  }

  fits <- list()
  rows <- lapply(names(models), function(model_name) {
    model_syntax <- models[[model_name]]
    fit_args <- list(model = model_syntax, data = data, estimator = estimator)
    if (!is.null(ordered)) fit_args$ordered <- ordered
    fit <- tryCatch(do.call(lavaan::cfa, fit_args), error = function(e) e)
    fits[[model_name]] <<- fit

    if (inherits(fit, "error")) {
      return(data.frame(
        model = model_name,
        converged = FALSE,
        estimator = estimator,
        cfi = NA_real_,
        tli = NA_real_,
        rmsea = NA_real_,
        srmr = NA_real_,
        cfi_scaled = NA_real_,
        tli_scaled = NA_real_,
        rmsea_scaled = NA_real_,
        cfi_robust = NA_real_,
        tli_robust = NA_real_,
        rmsea_robust = NA_real_,
        comparison_cfi = NA_real_,
        comparison_tli = NA_real_,
        comparison_rmsea = NA_real_,
        comparison_srmr = NA_real_,
        fit_score = NA_real_,
        aic = NA_real_,
        bic = NA_real_,
        error = conditionMessage(fit),
        stringsAsFactors = FALSE
      ))
    }

    measure_names <- c(
      "cfi", "tli", "rmsea", "srmr", "aic", "bic",
      "cfi.scaled", "tli.scaled", "rmsea.scaled",
      "cfi.robust", "tli.robust", "rmsea.robust"
    )
    measures <- tryCatch(
      lavaan::fitMeasures(fit, measure_names),
      error = function(e) stats::setNames(rep(NA_real_, length(measure_names)), measure_names)
    )
    comparison_cfi <- .pm_first_finite(measures[c("cfi.scaled", "cfi.robust", "cfi")])
    comparison_tli <- .pm_first_finite(measures[c("tli.scaled", "tli.robust", "tli")])
    comparison_rmsea <- .pm_first_finite(measures[c("rmsea.scaled", "rmsea.robust", "rmsea")])
    comparison_srmr <- .pm_first_finite(measures["srmr"])

    data.frame(
      model = model_name,
      converged = isTRUE(lavaan::lavInspect(fit, "converged")),
      estimator = estimator,
      cfi = unname(measures["cfi"]),
      tli = unname(measures["tli"]),
      rmsea = unname(measures["rmsea"]),
      srmr = unname(measures["srmr"]),
      cfi_scaled = unname(measures["cfi.scaled"]),
      tli_scaled = unname(measures["tli.scaled"]),
      rmsea_scaled = unname(measures["rmsea.scaled"]),
      cfi_robust = unname(measures["cfi.robust"]),
      tli_robust = unname(measures["tli.robust"]),
      rmsea_robust = unname(measures["rmsea.robust"]),
      comparison_cfi = comparison_cfi,
      comparison_tli = comparison_tli,
      comparison_rmsea = comparison_rmsea,
      comparison_srmr = comparison_srmr,
      fit_score = .pm_fit_score(comparison_cfi, comparison_tli, comparison_rmsea, comparison_srmr),
      aic = .pm_information_criterion(measures["aic"], estimator),
      bic = .pm_information_criterion(measures["bic"], estimator),
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  fit_table <- do.call(rbind, rows)
  rownames(fit_table) <- NULL
  recommendation <- .pm_model_recommendation(fit_table, language)

  out <- list(
    fit = fit_table,
    recommendation = recommendation,
    fitted_models = fits,
    estimator = estimator,
    ordered = ordered,
    call = match.call()
  )
  class(out) <- "psychomatic_model_comparison"
  out
}

.pm_mirt_item_parameters <- function(coefficients) {
  items <- if (is.list(coefficients) && "items" %in% names(coefficients)) {
    coefficients$items
  } else {
    coefficients
  }
  out <- as.data.frame(items)
  out$item <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("item", setdiff(names(out), "item")), drop = FALSE]
  out
}

.pm_irt_interpretation <- function(item_parameters,
                                   model,
                                   itemtype,
                                   language = "eng",
                                   digits = 3) {
  if (!is.data.frame(item_parameters) || nrow(item_parameters) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(item_parameters)), function(i) {
    row <- item_parameters[i, , drop = FALSE]
    discrimination <- .pm_item_discrimination(row)
    thresholds <- .pm_item_thresholds(row)
    difficulty_center <- if (length(thresholds) > 0L) {
      mean(thresholds, na.rm = TRUE)
    } else {
      .pm_item_difficulty(row)
    }
    threshold_min <- if (length(thresholds) > 0L) min(thresholds, na.rm = TRUE) else difficulty_center
    threshold_max <- if (length(thresholds) > 0L) max(thresholds, na.rm = TRUE) else difficulty_center
    threshold_ordered <- if (length(thresholds) <= 1L) TRUE else all(diff(thresholds) > 0)

    data.frame(
      item = as.character(row$item),
      discrimination = round(discrimination, digits),
      discrimination_level = .pm_discrimination_level(discrimination, language),
      difficulty_center = round(difficulty_center, digits),
      threshold_min = round(threshold_min, digits),
      threshold_max = round(threshold_max, digits),
      threshold_ordered = threshold_ordered,
      recommendation = .pm_irt_recommendation(
        discrimination = discrimination,
        difficulty_center = difficulty_center,
        threshold_min = threshold_min,
        threshold_max = threshold_max,
        threshold_ordered = threshold_ordered,
        language = language
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.pm_irt_model_summary <- function(interpretation,
                                  model,
                                  itemtype,
                                  n,
                                  language = "eng",
                                  digits = 3) {
  if (!is.data.frame(interpretation) || nrow(interpretation) == 0L) {
    return(data.frame(
      model = model,
      itemtype = itemtype,
      n = n,
      n_items = 0L,
      items_to_review = 0L,
      mean_discrimination = NA_real_,
      recommendation = if (language == "esp") {
        "No hay parametros interpretables disponibles."
      } else {
        "No interpretable item parameters were available."
      },
      stringsAsFactors = FALSE
    ))
  }

  review_flags <- !grepl(
    if (language == "esp") "^Sin alerta principal" else "^No major flag",
    interpretation$recommendation
  )
  items_to_review <- sum(review_flags, na.rm = TRUE)
  mean_discrimination <- mean(interpretation$discrimination, na.rm = TRUE)
  recommendation <- .pm_irt_overall_recommendation(
    items_to_review = items_to_review,
    n_items = nrow(interpretation),
    mean_discrimination = mean_discrimination,
    language = language
  )

  data.frame(
    model = model,
    itemtype = itemtype,
    n = n,
    n_items = nrow(interpretation),
    items_to_review = items_to_review,
    mean_discrimination = round(mean_discrimination, digits),
    recommendation = recommendation,
    stringsAsFactors = FALSE
  )
}

.pm_irt_overall_recommendation <- function(items_to_review,
                                           n_items,
                                           mean_discrimination,
                                           language) {
  review_pct <- if (n_items > 0L) items_to_review / n_items else 1
  if (language == "esp") {
    if (items_to_review == 0L && is.finite(mean_discrimination) && mean_discrimination >= 0.65) {
      return("Sin alertas principales; revisar contenido e informacion del test antes de concluir.")
    }
    if (review_pct <= 0.25) {
      return("Pocos items requieren revision; priorice contenido, discriminacion y umbrales.")
    }
    return("Varios items requieren revision; use estos resultados como cribado antes de decisiones definitivas.")
  }
  if (items_to_review == 0L && is.finite(mean_discrimination) && mean_discrimination >= 0.65) {
    return("No major item flags; review item content and test information before concluding.")
  }
  if (review_pct <= 0.25) {
    return("Few items require review; prioritize content, discrimination, and thresholds.")
  }
  "Several items require review; treat these results as screening evidence before final decisions."
}

.pm_item_discrimination <- function(row) {
  candidates <- intersect(c("a", "a1", "a2", "a3"), names(row))
  values <- as.numeric(row[1, candidates, drop = TRUE])
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  sqrt(sum(values^2))
}

.pm_item_difficulty <- function(row) {
  candidates <- intersect(c("b", "b1"), names(row))
  values <- as.numeric(row[1, candidates, drop = TRUE])
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  values[1L]
}

.pm_item_thresholds <- function(row) {
  threshold_names <- grep("^b[0-9]+$", names(row), value = TRUE)
  values <- as.numeric(row[1, threshold_names, drop = TRUE])
  values[is.finite(values)]
}

.pm_discrimination_level <- function(discrimination, language) {
  if (!is.finite(discrimination)) {
    return(if (language == "esp") "No disponible" else "Not available")
  }
  if (discrimination < 0.65) {
    if (language == "esp") "Baja" else "Low"
  } else if (discrimination < 1.35) {
    if (language == "esp") "Moderada" else "Moderate"
  } else if (discrimination < 1.70) {
    if (language == "esp") "Alta" else "High"
  } else {
    if (language == "esp") "Muy alta" else "Very high"
  }
}

.pm_irt_recommendation <- function(discrimination,
                                   difficulty_center,
                                   threshold_min,
                                   threshold_max,
                                   threshold_ordered,
                                   language) {
  issues <- character(0)
  if (!is.finite(discrimination) || discrimination < 0.65) {
    issues <- c(issues, if (language == "esp") "discriminacion baja" else "low discrimination")
  }
  if (is.finite(difficulty_center) && abs(difficulty_center) > 3) {
    issues <- c(issues, if (language == "esp") "dificultad extrema" else "extreme difficulty")
  }
  if (is.finite(threshold_min) && is.finite(threshold_max) &&
      (threshold_min < -4 || threshold_max > 4)) {
    issues <- c(issues, if (language == "esp") "umbrales extremos" else "extreme thresholds")
  }
  if (!isTRUE(threshold_ordered)) {
    issues <- c(issues, if (language == "esp") "umbrales desordenados" else "disordered thresholds")
  }
  if (length(issues) == 0L) {
    if (language == "esp") {
      "Sin alerta principal; revisar junto con contenido y precision del test."
    } else {
      "No major flag; review together with item content and test information."
    }
  } else if (language == "esp") {
    paste("Revisar item por", paste(unique(issues), collapse = ", "), ".")
  } else {
    paste("Review item for", paste(unique(issues), collapse = ", "), ".")
  }
}

.pm_dif_outcome <- function(x) {
  observed <- sort(unique(stats::na.omit(x)))
  if (length(observed) <= 2L) {
    high <- max(observed)
    return(list(
      y = as.integer(x == high),
      method = "binary logistic regression",
      note = ""
    ))
  }
  cutoff <- stats::median(x, na.rm = TRUE)
  list(
    y = as.integer(x > cutoff),
    method = "median-dichotomized logistic regression",
    note = "Polytomous item dichotomized at the median"
  )
}

.pm_dif_interpretation <- function(results,
                                   alpha,
                                   language = "eng",
                                   digits = 3) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(results)), function(i) {
    row <- results[i, , drop = FALSE]
    p_uniform <- as.numeric(row$p_uniform)
    p_nonuniform <- as.numeric(row$p_nonuniform)
    uniform_flag <- isTRUE(is.finite(p_uniform) && p_uniform < alpha)
    nonuniform_flag <- isTRUE(is.finite(p_nonuniform) && p_nonuniform < alpha)
    assessable <- is.finite(p_uniform) || is.finite(p_nonuniform)
    dif_type <- .pm_dif_type(uniform_flag, nonuniform_flag, assessable)
    severity <- .pm_dif_severity(dif_type)
    data.frame(
      item = row$item,
      dif_type = dif_type,
      p_uniform = round(p_uniform, digits),
      p_nonuniform = round(p_nonuniform, digits),
      decision = .pm_dif_decision(dif_type, language),
      severity = severity,
      recommendation = .pm_dif_recommendation(dif_type, row$method, row$note, language),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.pm_dif_summary <- function(interpretation,
                            alpha,
                            group_levels,
                            language = "eng",
                            digits = 3) {
  n_items <- nrow(interpretation)
  uniform_flags <- sum(interpretation$dif_type %in% c("uniform", "uniform_and_nonuniform"), na.rm = TRUE)
  nonuniform_flags <- sum(interpretation$dif_type %in% c("nonuniform", "uniform_and_nonuniform"), na.rm = TRUE)
  flagged_items <- sum(interpretation$severity == "review", na.rm = TRUE)
  flagged_pct <- if (n_items > 0L) flagged_items / n_items else NA_real_
  data.frame(
    n_items = n_items,
    flagged_items = flagged_items,
    uniform_flags = uniform_flags,
    nonuniform_flags = nonuniform_flags,
    flagged_pct = round(flagged_pct, digits),
    alpha = alpha,
    groups = paste(group_levels, collapse = ", "),
    recommendation = .pm_dif_overall_recommendation(flagged_items, n_items, language),
    stringsAsFactors = FALSE
  )
}

.pm_dif_type <- function(uniform_flag, nonuniform_flag, assessable) {
  if (!isTRUE(assessable)) return("not_evaluable")
  if (uniform_flag && nonuniform_flag) return("uniform_and_nonuniform")
  if (uniform_flag) return("uniform")
  if (nonuniform_flag) return("nonuniform")
  "none"
}

.pm_dif_severity <- function(dif_type) {
  if (dif_type == "not_evaluable") return("not_evaluable")
  if (dif_type == "none") return("no_flag")
  "review"
}

.pm_dif_decision <- function(dif_type, language) {
  if (language == "esp") {
    return(switch(
      dif_type,
      none = "Sin evidencia estadistica de DIF con este cribado.",
      uniform = "Posible DIF uniforme.",
      nonuniform = "Posible DIF no uniforme.",
      uniform_and_nonuniform = "Posible DIF uniforme y no uniforme.",
      not_evaluable = "No evaluable con este cribado."
    ))
  }
  switch(
    dif_type,
    none = "No statistical DIF evidence in this screening.",
    uniform = "Possible uniform DIF.",
    nonuniform = "Possible nonuniform DIF.",
    uniform_and_nonuniform = "Possible uniform and nonuniform DIF.",
    not_evaluable = "Not assessable with this screening."
  )
}

.pm_dif_recommendation <- function(dif_type, method, note, language) {
  extra_note <- if (!is.na(note) && nzchar(note)) paste0(" ", note) else ""
  if (language == "esp") {
    return(switch(
      dif_type,
      none = paste0("Conservar salvo evidencia teorica o clinica adicional.", extra_note),
      uniform = paste0("Revisar contenido del item y diferencias sistematicas entre grupos.", extra_note),
      nonuniform = paste0("Revisar interaccion item-grupo y considerar analisis DIF ordinal/IRT.", extra_note),
      uniform_and_nonuniform = paste0("Priorizar revision del item y confirmar con un metodo DIF dedicado.", extra_note),
      not_evaluable = paste0("Aumentar variabilidad o tamano muestral antes de interpretar.", extra_note)
    ))
  }
  switch(
    dif_type,
    none = paste0("Retain unless theory or substantive evidence suggests review.", extra_note),
    uniform = paste0("Review item content and systematic group differences.", extra_note),
    nonuniform = paste0("Review item-by-group interaction and consider ordinal or IRT DIF follow-up.", extra_note),
    uniform_and_nonuniform = paste0("Prioritize item review and confirm with a dedicated DIF method.", extra_note),
    not_evaluable = paste0("Increase variation or sample size before interpreting.", extra_note)
  )
}

.pm_dif_overall_recommendation <- function(flagged_items, n_items, language) {
  flagged_pct <- if (n_items > 0L) flagged_items / n_items else 1
  if (language == "esp") {
    if (flagged_items == 0L) {
      return("No se detectaron items con DIF en este cribado; interpretar junto con teoria y evidencia externa.")
    }
    if (flagged_pct <= 0.25) {
      return("Se detectaron pocos items con posible DIF; revisar esos items antes de decidir eliminacion o ajuste.")
    }
    return("Se detectaron varios items con posible DIF; confirmar con analisis DIF dedicado antes de decisiones finales.")
  }
  if (flagged_items == 0L) {
    return("No DIF-flagged items in this screening; interpret together with theory and external evidence.")
  }
  if (flagged_pct <= 0.25) {
    return("Few items show possible DIF; review those items before deletion or adjustment decisions.")
  }
  "Several items show possible DIF; confirm with a dedicated DIF analysis before final decisions."
}

.pm_irt_references <- function() {
  data.frame(
    method = c(
      "Ordered categories and Rasch-family response models",
      "Nominal and polytomous item response models",
      "Marginal maximum likelihood and EM estimation",
      "Generalized partial credit and polytomous IRT",
      "Rasch, 2PL, graded response, and multidimensional IRT implementation"
    ),
    reference = c(
      "Andrich, D. (1978). A rating formulation for ordered response categories. Psychometrika, 43, 561-573. https://doi.org/10.1007/BF02293814",
      "Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or more nominal categories. Psychometrika, 37(1), 29-51. https://doi.org/10.1007/BF02291411",
      "Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm. Psychometrika, 46, 443-459. https://doi.org/10.1007/BF02293801",
      "Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. Applied Psychological Measurement, 16(2), 159-176. https://doi.org/10.1177/014662169201600206",
      "Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment. Journal of Statistical Software, 48(6), 1-29. https://doi.org/10.18637/jss.v048.i06"
    ),
    stringsAsFactors = FALSE
  )
}

.pm_dif_references <- function() {
  data.frame(
    method = c(
      "Contingency-table and logit item-bias models",
      "Logistic regression DIF screening",
      "Logistic regression DIF effect-size evaluation",
      "IRT-based DIF effect-size reporting",
      "Ordinal logistic regression and IRT DIF implementation"
    ),
    reference = c(
      "Mellenbergh, G. J. (1982). Contingency table models for assessing item bias. Journal of Educational Statistics, 7(2), 105-118. https://doi.org/10.3102/10769986007002105",
      "Swaminathan, H., & Rogers, H. J. (1990). Detecting differential item functioning using logistic regression procedures. Journal of Educational Measurement, 27(4), 361-370. https://doi.org/10.1111/j.1745-3984.1990.tb00754.x",
      "Jodoin, M. G., & Gierl, M. J. (2001). Evaluating Type I error and power rates using an effect size measure with the logistic regression procedure for DIF detection. Applied Measurement in Education, 14(4), 329-349. https://doi.org/10.1207/S15324818AME1404_2",
      "Steinberg, L., & Thissen, D. (2006). Using effect sizes for research reporting: Examples using item response theory to analyze differential item functioning. Psychological Methods, 11(4), 402-415. https://doi.org/10.1037/1082-989X.11.4.402",
      "Choi, S. W., Gibbons, L. E., & Crane, P. K. (2011). lordif: An R package for detecting differential item functioning using iterative hybrid ordinal logistic regression/item response theory and Monte Carlo simulations. Journal of Statistical Software, 39(8), 1-30. https://doi.org/10.18637/jss.v039.i08"
    ),
    stringsAsFactors = FALSE
  )
}

.pm_model_recommendation <- function(fit_table, language) {
  candidates <- fit_table[fit_table$converged & is.finite(fit_table$fit_score), , drop = FALSE]
  if (nrow(candidates) == 0L) {
    message <- if (language == "esp") {
      "Ningun modelo convergente con indices de ajuste comparables; revise especificaciones y datos antes de decidir."
    } else {
      "No converged model with comparable fit indices was found; inspect specifications and data before deciding."
    }
    return(list(model = NA_character_, basis = "none", message = message))
  }
  best <- candidates[.pm_rank_fit_models(candidates)[1L], , drop = FALSE]
  message <- if (language == "esp") {
    sprintf(
      "El modelo `%s` presenta la combinacion descriptivamente mas favorable de CFI/TLI altos y RMSEA/SRMR bajos; interprete esta recomendacion junto con teoria, cargas y residuos.",
      best$model
    )
  } else {
    sprintf(
      "Model `%s` has the most favorable descriptive combination of higher CFI/TLI and lower RMSEA/SRMR; interpret this recommendation together with theory, loadings, and residuals.",
      best$model
    )
  }
  list(model = best$model, basis = "fit_indices", message = message)
}

.pm_model_observed_variables <- function(models, data) {
  vars <- unique(unlist(lapply(models, function(model) {
    pt <- tryCatch(lavaan::lavaanify(model), error = function(e) NULL)
    if (is.null(pt)) return(character(0))
    unique(pt$rhs[pt$op == "=~"])
  }), use.names = FALSE))
  intersect(vars, names(data))
}

.pm_first_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else x[1L]
}

.pm_information_criterion <- function(x, estimator) {
  if (estimator %in% c("WLSMV", "DWLS", "ULSMV", "ULS")) return(NA_real_)
  value <- as.numeric(x)
  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

.pm_fit_score <- function(cfi, tli, rmsea, srmr) {
  score <- 0
  if (is.finite(cfi)) score <- score + ifelse(cfi >= .95, 2, ifelse(cfi >= .90, 1, 0))
  if (is.finite(tli)) score <- score + ifelse(tli >= .95, 2, ifelse(tli >= .90, 1, 0))
  if (is.finite(rmsea)) score <- score + ifelse(rmsea <= .06, 2, ifelse(rmsea <= .08, 1, 0))
  if (is.finite(srmr)) score <- score + ifelse(srmr <= .08, 2, ifelse(srmr <= .10, 1, 0))
  score
}

.pm_rank_fit_models <- function(candidates) {
  order(
    -candidates$fit_score,
    -ifelse(is.finite(candidates$comparison_cfi), candidates$comparison_cfi, -Inf),
    -ifelse(is.finite(candidates$comparison_tli), candidates$comparison_tli, -Inf),
    ifelse(is.finite(candidates$comparison_rmsea), candidates$comparison_rmsea, Inf),
    ifelse(is.finite(candidates$comparison_srmr), candidates$comparison_srmr, Inf)
  )
}
