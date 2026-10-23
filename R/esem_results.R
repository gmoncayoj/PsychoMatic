.pm_esem_indices <- function(fit, variant) {
  metrics <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea",
               "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
  cap <- .pm_capture_conditions(lavaan::fitMeasures(fit))
  fm <- cap$value
  get <- function(key) {
    x <- fm[key]
    if (length(x) != 1L || !is.finite(x)) NA_real_ else unname(x)
  }
  rows <- lapply(metrics, function(metric) {
    base <- get(metric)
    scaled <- get(paste0(metric, ".scaled"))
    robust <- get(paste0(metric, ".robust"))
    effective <- if (metric %in% c("chisq", "df", "pvalue")) {
      if (variant == "base") "base" else "scaled"
    } else if (metric == "srmr") "base" else variant
    value <- switch(effective, base = base, scaled = scaled, robust = robust)
    data.frame(metric = metric, value = value, variant = effective,
      base = base, scaled = scaled, robust = robust,
      reason = if (is.na(value)) paste("Requested", effective, "index unavailable",
        if (!is.null(cap$error)) paste0(": ", conditionMessage(cap$error)) else "") else "",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  attr(out, "conditions") <- cap[c("warnings", "error")]
  out
}

.pm_esem_result <- function(prepared, fitted, call, fit_variant, language) {
  fit <- fitted$fit
  conditions <- fitted$warnings
  errors <- if (is.null(fitted$error)) list() else list(fitted$error)
  settings <- prepared$settings
  results <- list(fit = fit, loadings = data.frame(), item_diagnostics = data.frame(),
    factor_correlations = data.frame(), fit_indices = data.frame(), residuals = data.frame(),
    target = prepared$target, item_map = prepared$item_map,
    factor_map = prepared$factor_map, syntax = prepared$syntax,
    preparation = prepared$preparation_table)
  diagnostics <- list(converged = FALSE, admissible = FALSE, finite_se = FALSE,
    case_indices = integer(), saturated = FALSE, underidentified = FALSE)
  reason <- if (length(errors)) conditionMessage(errors[[1]]) else "No fitted model was returned."
  used <- 0L
  if (!is.null(fit)) {
    captures <- list(
      convergence = .pm_capture_conditions(lavaan::lavInspect(fit, "converged")),
      post_check = .pm_capture_conditions(lavaan::lavInspect(fit, "post.check")),
      standardized = .pm_capture_conditions(lavaan::standardizedSolution(fit)),
      parameters = .pm_capture_conditions(lavaan::parameterEstimates(fit)),
      cases = .pm_capture_conditions(lavaan::lavInspect(fit, "case.idx")),
      nobs = .pm_capture_conditions(lavaan::lavInspect(fit, "nobs")),
      vcov = .pm_capture_conditions(lavaan::vcov(fit)),
      options = .pm_capture_conditions(lavaan::lavInspect(fit, "options")))
    for (cap in captures) {
      conditions <- c(conditions, cap$warnings)
      if (!is.null(cap$error)) errors <- c(errors, list(cap$error))
    }
    diagnostics$converged <- isTRUE(captures$convergence$value)
    diagnostics$admissible <- isTRUE(captures$post_check$value)
    diagnostics$case_indices <- as.integer(unlist(captures$cases$value, use.names = FALSE))
    used <- as.integer(sum(unlist(captures$nobs$value)))
    ss <- captures$standardized$value
    pe <- captures$parameters$value
    if (!is.null(ss) && !is.null(pe)) {
      ld <- ss[ss$op == "=~", , drop = FALSE]
      raw <- pe[match(paste(ld$lhs, ld$rhs), paste(pe$lhs, pe$rhs)), , drop = FALSE]
      item_map <- prepared$item_map
      factor_map <- prepared$factor_map
      results$loadings <- data.frame(
        factor = factor_map$label[match(ld$lhs, factor_map$factor)],
        factor_id = ld$lhs, item = item_map$item[match(ld$rhs, item_map$alias)],
        estimate = raw$est, raw_se = raw$se,
        raw_ci_lower = raw$ci.lower, raw_ci_upper = raw$ci.upper,
        standardized = ld$est.std, se = ld$se, z = ld$z, p_value = ld$pvalue,
        ci_lower = ld$ci.lower, ci_upper = ld$ci.upper, stringsAsFactors = FALSE)
      diagnostics$finite_se <- nrow(ld) == nrow(item_map) * nrow(factor_map) &&
        all(is.finite(ld$se)) && all(is.finite(ld$est.std)) &&
        all(is.finite(raw$est)) && all(is.finite(raw$se))
      diagnostics$parameters_finite <- all(is.finite(pe$est))
      results$item_diagnostics <- do.call(rbind, lapply(item_map$item, function(item) {
        rows <- results$loadings[results$loadings$item == item, , drop = FALSE]
        rank <- order(abs(rows$standardized), decreasing = TRUE, na.last = TRUE)
        a <- rows[rank[1L], ]; b <- rows[rank[2L], ]
        expected <- if (settings$rotation == "target") {
          which(vapply(prepared$factors, function(v) item %in% v, logical(1)))
        } else integer()
        expected_label <- if (length(expected)) factor_map$label[expected] else NA_character_
        expected_loading <- if (length(expected)) rows$standardized[rows$factor_id == factor_map$factor[expected]] else NA_real_
        data.frame(item = item, largest_factor = a$factor,
          largest_loading = a$standardized, second_factor = b$factor,
          second_loading = b$standardized,
          difference = abs(a$standardized) - abs(b$standardized),
          expected_factor = expected_label, expected_loading = expected_loading,
          expected_is_largest = if (length(expected)) expected_label == a$factor else NA,
          stringsAsFactors = FALSE)
      }))
      fc <- ss[ss$op == "~~" & ss$lhs %in% factor_map$factor &
                 ss$rhs %in% factor_map$factor & ss$lhs != ss$rhs, , drop = FALSE]
      results$factor_correlations <- data.frame(
        factor1 = factor_map$label[match(fc$lhs, factor_map$factor)],
        factor2 = factor_map$label[match(fc$rhs, factor_map$factor)],
        correlation = fc$est.std, se = fc$se,
        ci_lower = fc$ci.lower, ci_upper = fc$ci.upper, stringsAsFactors = FALSE)
    }
    vc <- captures$vcov$value
    diagnostics$vcov_finite <- !is.null(vc) && all(is.finite(vc))
    opt <- captures$options$value
    if (!is.null(opt)) {
      settings$native_estimator <- opt$estimator
      settings$native_missing <- opt$missing
    }
    warn_text <- vapply(conditions, conditionMessage, character(1))
    diagnostics$rotation_warning <- any(grepl("rotat|GPA", warn_text, ignore.case = TRUE) &
      grepl("converg|fail|singular|maximum|iteration", warn_text, ignore.case = TRUE))
    # Robust contrast/baseline failures concern fit indices, not rotated SEs.
    parameter_warnings <- warn_text[!grepl("lav_test_|robust test statistic|baseline model",
                                          warn_text, ignore.case = TRUE)]
    diagnostics$inference_warning <- any(
      grepl("vcov|information matrix|standard errors|identif", parameter_warnings, ignore.case = TRUE) &
      grepl("not.*positive definite|could not|not.*identif|singular|invert", parameter_warnings, ignore.case = TRUE))
    reason <- if (!diagnostics$converged) "The ESEM optimizer did not converge." else
      if (!diagnostics$admissible) "The native post-check rejected the solution." else
      if (!diagnostics$finite_se || !isTRUE(diagnostics$parameters_finite) ||
          !diagnostics$vcov_finite) "Rotated parameter inference is not finite or unavailable." else
      if (diagnostics$rotation_warning) "A native rotation warning prevents reliable interpretation." else
      if (diagnostics$inference_warning) "A native identification or parameter-inference warning prevents reliable interpretation." else
      if (length(errors)) conditionMessage(errors[[1L]]) else NA_character_
    if (diagnostics$converged) {
      results$fit_indices <- .pm_esem_indices(fit, fit_variant)
      index_conditions <- attr(results$fit_indices, "conditions")
      conditions <- c(conditions, index_conditions$warnings)
      diagnostics$fit_indices_error <- if (is.null(index_conditions$error)) NULL else
        conditionMessage(index_conditions$error)
      if (!is.null(index_conditions$error)) errors <- c(errors, list(index_conditions$error))
      df <- results$fit_indices$base[results$fit_indices$metric == "df"]
      diagnostics$saturated <- length(df) == 1L && is.finite(df) && df == 0
      diagnostics$underidentified <- length(df) == 1L && is.finite(df) && df < 0
      if (diagnostics$underidentified && is.na(reason)) {
        reason <- "The ESEM model has negative degrees of freedom and is underidentified."
      }
      residual <- .pm_capture_conditions(lavaan::residuals(fit, type = "raw"))
      conditions <- c(conditions, residual$warnings)
      if (!is.null(residual$value$cov)) {
        r <- residual$value$cov
        rownames(r) <- prepared$item_map$item[match(rownames(r), prepared$item_map$alias)]
        colnames(r) <- prepared$item_map$item[match(colnames(r), prepared$item_map$alias)]
        results$residuals <- .pm_matrix_table(r)
      }
    }
  }
  out <- .pm_result(class = "esem_auto", call = call,
    requested_settings = list(fit_variant = fit_variant, language = language),
    effective_settings = c(settings, list(fit_variant = fit_variant)),
    data_summary = list(n_total = nrow(prepared$analysis_data), n = used,
      n_complete = sum(stats::complete.cases(prepared$analysis_data)),
      n_items = ncol(prepared$analysis_data), n_factors = settings$n_factors,
      missing_cells_pct = 100 * mean(is.na(prepared$analysis_data))),
    results = results, warnings = conditions, errors = errors,
    diagnostics = diagnostics, decision_log = prepared$decisions,
    conclusion = list(status = if (is.na(reason)) "evaluable" else "not_evaluable", reason = reason))
  out$fit <- fit
  out$references <- data.frame(reference = c(
    "Asparouhov, T., & Muth\u00e9n, B. (2009). Exploratory structural equation modeling. Structural Equation Modeling, 16(3), 397\u2013438. https://doi.org/10.1080/10705510903008204",
    "Marsh, H., & Alamer, A. (2024). When and how to use set-exploratory structural equation modelling to test structural models: A tutorial using the R package lavaan. British Journal of Mathematical and Statistical Psychology, 77(3), 459\u2013476. https://doi.org/10.1111/bmsp.12336"),
    stringsAsFactors = FALSE)
  out$versions <- list(PsychoMatic = as.character(utils::packageVersion("PsychoMatic")),
                        lavaan = as.character(utils::packageVersion("lavaan")),
                        R = as.character(getRversion()))
  out$manuscript_text <- .pm_esem_manuscript(out, language)
  out$results$manuscript <- data.frame(Text = out$manuscript_text, stringsAsFactors = FALSE)
  out
}

.pm_esem_compare_cfa <- function(prepared, esem_fit, seed, fit_variant) {
  args <- prepared$native_args
  args$model <- prepared$cfa_syntax
  args$rotation <- NULL
  cap <- .efa_with_seed(seed, .pm_capture_conditions(do.call(lavaan::cfa, args)))
  fit <- cap$value
  reason <- if (!is.null(cap$error)) conditionMessage(cap$error) else NA_character_
  if (!is.null(fit)) {
    check <- .pm_capture_conditions({
      isTRUE(lavaan::lavInspect(fit, "converged")) &&
        isTRUE(lavaan::lavInspect(fit, "post.check")) &&
        isTRUE(lavaan::lavInspect(esem_fit, "converged")) &&
        isTRUE(lavaan::lavInspect(esem_fit, "post.check")) &&
        identical(lavaan::lavInspect(fit, "case.idx"), lavaan::lavInspect(esem_fit, "case.idx")) &&
        identical(lavaan::lavNames(fit, "ov.ord"), lavaan::lavNames(esem_fit, "ov.ord"))
    })
    cap$warnings <- c(cap$warnings, check$warnings)
    if (!isTRUE(check$value)) reason <- "The reference fit or common-sample comparability check failed."
  }
  table <- data.frame()
  if (is.na(reason) && !is.null(esem_fit)) {
    a <- .pm_esem_indices(esem_fit, fit_variant)
    b <- .pm_esem_indices(fit, fit_variant)
    for (indices in list(a, b)) {
      index_conditions <- attr(indices, "conditions")
      cap$warnings <- c(cap$warnings, index_conditions$warnings)
      if (!is.null(index_conditions$error)) cap$warnings <- c(cap$warnings,
        list(simpleWarning(paste("CFA comparison index calculation:", conditionMessage(index_conditions$error)))))
    }
    table <- data.frame(metric = a$metric, variant = a$variant,
                         esem = a$value, cfa = b$value, delta = a$value - b$value,
                         reason = ifelse(is.na(a$value) | is.na(b$value), "Index unavailable in at least one model", ""),
                         stringsAsFactors = FALSE)
  }
  list(status = if (is.na(reason)) "evaluable" else "not_evaluable", reason = reason,
       fit = fit, table = table, warnings = cap$warnings,
       errors = if (is.null(cap$error)) list() else list(cap$error))
}
