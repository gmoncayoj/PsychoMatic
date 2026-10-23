# =============================================================================
# FUNCION PRINCIPAL
# =============================================================================

#' Automated CFA with psychometric reporting
#'
#' Runs a confirmatory factor analysis with heuristic estimator selection,
#' multivariate normality diagnostics, fit indices, reliability estimates, and
#' model respecification suggestions.
#' `lavaan` performs the model fit. When installed, `MVN` supplies the raw
#' Mardia diagnostic used by automatic estimator selection and `semTools`
#' supplies composite reliability and AVE. Their absence is recorded explicitly
#' and does not prevent the CFA from being fitted; when normality cannot be
#' evaluated, automatic selection conservatively uses MLR for continuous data.
#' Fit indices retain their base, scaled, and robust values, and the effective
#' variant is identified separately for every index.
#'
#' @inheritParams efa_auto
#' @param fit_variant Preferred fit-index variant: `"scaled"` (default), `"robust"`,
#'   or `"base"`. Unavailable requested indices fall back to scaled then base.
#'   Base, scaled and robust estimates remain available in the result.
#' @param data Data frame containing the observed variables.
#' @param model Model syntax written for `lavaan`.
#' @param ordered Optional vector indicating which items should be treated as ordinal.
#' @param estimator Optional estimator. If `NULL`, the function selects one automatically.
#' @param std.lv Whether latent variables should be standardized.
#' @param mi_threshold Minimum threshold for reporting modification indices.
#' @param n_mi Maximum number of modification indices to display.
#' @param alpha_norm Alpha level used when evaluating multivariate normality.
#' @param language Report language: `"esp"` or `"eng"`.
#' @param bifactor_orthogonal Logical scalar used only for detected bifactor
#'   models. The default preserves the historical orthogonal bifactor fit; set
#'   `FALSE` only when latent-factor correlations are substantively intended.
#'
#' @return An object of classes `cfa_auto` and `psychomatic_result`. The fitted
#'   `lavaan` object is available in `fit`; reusable tidy tables are available in
#'   `tables`, `fit_indices`, `loadings_table`, and `reliability_table`.
#'   `warnings` is a structured data frame with condition metadata;
#'   `warnings_list` is the transitional character-vector view retained for
#'   code written for PsychoMatic 0.3.x.
#'
#' @section Compatibility:
#' Since PsychoMatic 1.0.0, use `result$warnings` for programmatic inspection
#' of warning sources and messages. Existing code that expects a character
#' vector can migrate to `result$warnings_list`; this compatibility alias will
#' be retained through a deprecation cycle. For example,
#' `warning_messages <- result$warnings_list` replaces code that previously
#' treated `result$warnings` as character data.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' model <- "F1 =~ R1 + R2 + R3 + R4"
#' cfa_auto(MY_DATA, model = model)
#' }
#'
#' @export
cfa_auto <- function(data,
                       model,
                       ordered = NULL,
                       estimator = NULL,
                       std.lv = FALSE,
                       mi_threshold = 10,
                       n_mi = 10,
                       alpha_norm = 0.05,
                       language = "eng",
                       bifactor_orthogonal = TRUE,
                       digits = 3, report = FALSE, file_name = NULL,
                       format = "excel", overwrite = FALSE, verbose = TRUE,
                       fit_variant = c("scaled", "robust", "base")) {

  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  fit_variant <- match.arg(fit_variant)
  pick_variant <- function(values, metric) {
    .pm_pick_fit_variant(values, metric, priority = unique(c(fit_variant, "scaled", "base")))
  }
  language <- .normalize_language(language)
  is_english <- identical(language, "eng")

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    if (is_english) {
      stop("Package 'lavaan' is required.", call. = FALSE)
    }
    stop("Se requiere el paquete 'lavaan'.", call. = FALSE)
  }
  if (!is.logical(bifactor_orthogonal) || length(bifactor_orthogonal) != 1L ||
      is.na(bifactor_orthogonal)) {
    stop("`bifactor_orthogonal` must be TRUE or FALSE.", call. = FALSE)
  }
  alpha_norm <- .pm_validate_scalar_number(
    alpha_norm, "alpha_norm", lower = 0, upper = 1
  )

  # =========================================================================
  # 1. PARSEAR MODELO
  # =========================================================================
  pt <- lavaanify(model)
  loading_rows <- pt[pt$op == "=~", ]
  factors <- unique(loading_rows$lhs)

  items_per_factor <- lapply(setNames(factors, factors), function(f) {
    loading_rows$rhs[loading_rows$lhs == f]
  })

  all_rhs <- unique(loading_rows$rhs)
  latent_vars <- factors
  observed_items <- setdiff(all_rhs, latent_vars)

  # Verificar que los items existen en los datos
  missing_items <- setdiff(observed_items, names(data))
  if (length(missing_items) > 0) {
    if (is_english) {
      stop("Items not found in the data: ", paste(missing_items, collapse = ", "))
    }
    stop("Items no encontrados en los datos: ", paste(missing_items, collapse = ", "))
  }

  # =========================================================================
  # 2. DETECTAR TIPO DE MODELO
  # =========================================================================
  items_per_factor_obs <- lapply(items_per_factor, function(x) intersect(x, observed_items))
  n_obs_per_factor <- sapply(items_per_factor_obs, length)

  # Segundo orden: al menos un factor definido enteramente por otros factores
  second_order_factors <- character(0)
  first_order_factors <- character(0)

  for (f in factors) {
    indicadores <- items_per_factor[[f]]
    if (all(indicadores %in% latent_vars)) {
      second_order_factors <- c(second_order_factors, f)
    } else if (any(indicadores %in% observed_items)) {
      first_order_factors <- c(first_order_factors, f)
    }
  }

  is_second_order <- length(second_order_factors) > 0

  # Bifactor: exactamente un factor general carga en todos los items y los
  # factores especificos cargan en subconjuntos. Las estructuras ambiguas no
  # activan indices especializados.
  bifactor_detection <- .pm_detect_bifactor_structure(
    items_per_factor_obs,
    observed_items = observed_items,
    is_second_order = is_second_order
  )
  is_bifactor <- bifactor_detection$is_bifactor
  general_factor <- bifactor_detection$general_factor
  specific_factors <- bifactor_detection$specific_factors

  is_unidimensional <- length(factors) == 1 && !is_second_order

  model_type_code <- if (is_bifactor) "bifactor"
    else if (is_second_order) "second_order"
    else if (is_unidimensional) "unidimensional"
    else "correlated"
  model_type <- .model_type_label(model_type_code, language)

  # Para segundo orden, identificar los factores de primer orden reales
  if (is_second_order) {
    fo_factors <- setdiff(factors, second_order_factors)
    # Re-armar items_per_factor solo con primer orden para confiabilidad
    items_fo <- lapply(setNames(fo_factors, fo_factors), function(f) {
      items_per_factor[[f]]
    })
  }

  # =========================================================================
  # 3. DIAGNOSTICOS DE DATOS
  # =========================================================================
  data_items <- data[, observed_items, drop = FALSE]
  n_total <- nrow(data_items)
  n_complete <- sum(complete.cases(data_items))
  n_missing_rows <- n_total - n_complete
  pct_missing <- round(100 * n_missing_rows / n_total, 1)

  # Conteo de categorias por item
  n_categories <- sapply(data_items, function(x) length(unique(na.omit(x))))
  max_cat <- max(n_categories)
  min_cat <- min(n_categories)

  # Proporcion de missing por item
  pct_na_item <- sapply(data_items, function(x) round(100 * mean(is.na(x)), 1))

  # =========================================================================
  # 4. SELECCION AUTOMATICA DE ESTIMADOR
  # =========================================================================
  mvn_result <- NULL
  mardia_skew_stat <- NA
  mardia_kurt_stat <- NA
  mardia_p_skew <- NA
  mardia_p_kurt <- NA
  mardia_data <- as.data.frame(lapply(data_items, function(x) {
    if (is.ordered(x)) as.integer(x) else x
  }))
  mvn_diagnostic <- .pm_mvn_mardia(
    mardia_data[complete.cases(mardia_data), , drop = FALSE], alpha = alpha_norm)
  mvn_result <- mvn_diagnostic$result
  mardia_p_skew <- mvn_diagnostic$p_skew
  mardia_p_kurt <- mvn_diagnostic$p_kurt
  mardia_skew_stat <- mvn_diagnostic$skew_stat
  mardia_kurt_stat <- mvn_diagnostic$kurt_stat

  selection <- .pm_select_cfa_estimator(max_cat, min_cat, n_complete,
    mvn_diagnostic, estimator, ordered, language)
  estimator <- selection$estimator
  ordered <- selection$ordered
  estimator_reason <- selection$reason

  # =========================================================================
  # 5. AJUSTAR MODELO
  # =========================================================================
  fit_args <- list(model = model, data = data, estimator = estimator, std.lv = std.lv)

  if (isTRUE(ordered)) {
    fit_args$ordered <- TRUE
  } else if (is.character(ordered)) {
    fit_args$ordered <- ordered
  }

  bifactor_orthogonal_effective <- if (is_bifactor) {
    isTRUE(bifactor_orthogonal)
  } else {
    NA
  }

  # Bifactor: la ortogonalidad se controla explicitamente.
  if (is_bifactor) {
    fit_args$orthogonal <- bifactor_orthogonal_effective
  }

  fitted_conditions <- .pm_capture_conditions(do.call(lavaan::cfa, fit_args))
  if (!is.null(fitted_conditions$error)) stop(fitted_conditions$error)
  fit <- fitted_conditions$value
  vcov_diagnostic <- tryCatch({
    values <- eigen(lavaan::lavInspect(fit, "vcov"), symmetric = TRUE, only.values = TRUE)$values
    list(minimum_eigenvalue = min(values),
         tolerance = sqrt(.Machine$double.eps) * max(1, max(abs(values))))
  }, error = function(e) list(minimum_eigenvalue = NA_real_, tolerance = NA_real_, reason = conditionMessage(e)))
  converged <- lavInspect(fit, "converged")
  bifactor_orthogonal_fitted <- if (is_bifactor) {
    .pm_bifactor_is_orthogonal(fit, c(general_factor, specific_factors))
  } else {
    NA
  }

  # =========================================================================
  # 6. INDICES DE AJUSTE
  # =========================================================================
  all_fi <- tryCatch(fitMeasures(fit), error = function(e) NULL)

  if (!is.null(all_fi)) {
    # Extraer cada indice
    p_chisq  <- pick_variant(all_fi, "chisq")
    p_df     <- pick_variant(all_fi, "df")
    p_pval   <- pick_variant(all_fi, "pvalue")
    p_cfi    <- pick_variant(all_fi, "cfi")
    p_tli    <- pick_variant(all_fi, "tli")
    p_rmsea  <- pick_variant(all_fi, "rmsea")
    p_rmsea_lo <- pick_variant(all_fi, "rmsea.ci.lower")
    p_rmsea_hi <- pick_variant(all_fi, "rmsea.ci.upper")
    p_rmsea_p  <- pick_variant(all_fi, "rmsea.pvalue")

    # SRMR: no tiene version scaled/robust tipicamente; probar srmr_bentler tambien
    p_srmr <- pick_variant(all_fi, "srmr")
    if (is.na(p_srmr$value)) {
      # Intentar srmr_bentler (comun en WLSMV)
      if ("srmr_bentler" %in% names(all_fi) && !is.na(all_fi["srmr_bentler"])) {
        p_srmr <- list(value = unname(all_fi["srmr_bentler"]), variant = "base")
      }
    }

    # Chi-sq/gl
    chisq_ratio <- if (!is.na(p_chisq$value) && !is.na(p_df$value) && p_df$value > 0) {
      p_chisq$value / p_df$value
    } else NA_real_
    ratio_variants <- c(
      chisq_ratio = .pm_fit_measure_value(all_fi, "chisq") /
        .pm_fit_measure_value(all_fi, "df"),
      chisq_ratio.scaled = .pm_fit_measure_value(all_fi, "chisq.scaled") /
        .pm_fit_measure_value(all_fi, "df.scaled"),
      chisq_ratio.robust = .pm_fit_measure_value(all_fi, "chisq.robust") /
        .pm_fit_measure_value(all_fi, "df.robust")
    )
    ratio_variants[!is.finite(ratio_variants)] <- NA_real_
    p_ratio <- pick_variant(ratio_variants, "chisq_ratio")
    chisq_ratio <- p_ratio$value

    # Determinar variante dominante para la nota informativa
    variants_used <- unique(c(p_cfi$variant, p_tli$variant, p_rmsea$variant))
    variants_used <- setdiff(variants_used, "base")
    fi_variant_note <- if (length(variants_used) > 0) {
      if (is_english) {
        paste0("Reported corrected indices: ", paste(variants_used, collapse = "/"))
      } else {
        paste0("Indices corregidos reportados: ", paste(variants_used, collapse = "/"))
      }
    } else {
      if (is_english) {
        "Uncorrected indices reported (classical ML estimator)"
      } else {
        "Indices sin correccion (estimador ML clasico)"
      }
    }

    fit_summary <- data.frame(
      Metric = c(
        "chisq", "df", "pvalue", "chisq_ratio", "cfi", "tli",
        "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "srmr"
      ),
      Index = c(
        .fit_index_label("chisq", p_chisq$variant, language),
        .fit_index_base_label("df", language),
        .fit_index_label("pvalue", p_pval$variant, language),
        .fit_index_base_label("chisq_ratio", language),
        .fit_index_label("cfi", p_cfi$variant, language),
        .fit_index_label("tli", p_tli$variant, language),
        .fit_index_label("rmsea", p_rmsea$variant, language),
        .fit_index_label("rmsea.ci.lower", p_rmsea_lo$variant, language),
        .fit_index_label("rmsea.ci.upper", p_rmsea_hi$variant, language),
        .fit_index_label("rmsea.pvalue", p_rmsea_p$variant, language),
        .fit_index_base_label("srmr", language)
      ),
      Value = c(
        p_chisq$value, p_df$value, p_pval$value, chisq_ratio,
        p_cfi$value, p_tli$value,
        p_rmsea$value, p_rmsea_lo$value, p_rmsea_hi$value, p_rmsea_p$value,
        p_srmr$value
      ),
      Base = c(
        vapply(c("chisq", "df", "pvalue"), function(x) .pm_fit_measure_value(all_fi, x), numeric(1)),
        ratio_variants[["chisq_ratio"]],
        vapply(c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "srmr"),
               function(x) .pm_fit_measure_value(all_fi, x), numeric(1))
      ),
      Scaled = c(
        vapply(c("chisq", "df", "pvalue"), function(x) .pm_fit_measure_value(all_fi, paste0(x, ".scaled")), numeric(1)),
        ratio_variants[["chisq_ratio.scaled"]],
        vapply(c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "srmr"),
               function(x) .pm_fit_measure_value(all_fi, paste0(x, ".scaled")), numeric(1))
      ),
      Robust = c(
        vapply(c("chisq", "df", "pvalue"), function(x) .pm_fit_measure_value(all_fi, paste0(x, ".robust")), numeric(1)),
        ratio_variants[["chisq_ratio.robust"]],
        vapply(c("cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue", "srmr"),
               function(x) .pm_fit_measure_value(all_fi, paste0(x, ".robust")), numeric(1))
      ),
      Variant = c(
        p_chisq$variant, p_df$variant, p_pval$variant, p_ratio$variant,
        p_cfi$variant, p_tli$variant, p_rmsea$variant, p_rmsea_lo$variant,
        p_rmsea_hi$variant, p_rmsea_p$variant, p_srmr$variant
      ),
      stringsAsFactors = FALSE
    )

    # Evaluar ajuste
    fit_summary$Evaluation <- ""
    fit_summary$Evaluation[5] <-
      .evaluar_fi(fit_summary$Value[5], 0.95, 0.90, ">=", language)
    fit_summary$Evaluation[6] <-
      .evaluar_fi(fit_summary$Value[6], 0.95, 0.90, ">=", language)
    fit_summary$Evaluation[7] <-
      .evaluar_fi(fit_summary$Value[7], 0.06, 0.08, "<=", language)
    fit_summary$Evaluation[11] <-
      .evaluar_fi(fit_summary$Value[11], 0.06, 0.08, "<=", language)
    fit_summary$Evaluation[4] <-
      .evaluar_fi(fit_summary$Value[4], 3, 5, "<=", language)
    fit_summary$Indice <- fit_summary$Index
    fit_summary$Valor <- fit_summary$Value
    fit_summary$Evaluacion <- fit_summary$Evaluation

  } else {
    fit_summary <- data.frame(
      Metric = "error",
      Index = .fit_index_base_label("error", language),
      Value = NA,
      Base = NA_real_,
      Scaled = NA_real_,
      Robust = NA_real_,
      Variant = "base",
      Evaluation = if (is_english) "Not available" else "No disponible",
      stringsAsFactors = FALSE
    )
    fit_summary$Indice <- fit_summary$Index
    fit_summary$Valor <- fit_summary$Value
    fit_summary$Evaluacion <- fit_summary$Evaluation
    fi_variant_note <- if (is_english) "Not available" else "No disponible"
  }

  # =========================================================================
  # 7. CARGAS FACTORIALES ESTANDARIZADAS
  # =========================================================================
  std_sol <- standardizedSolution(fit)
  loadings_df <- std_sol[std_sol$op == "=~",
                         c("lhs", "rhs", "est.std", "se", "z", "pvalue",
                           "ci.lower", "ci.upper")]
  names(loadings_df) <- c("Factor", "Item", "Loading", "SE", "z", "p_value",
                          "CI_lower", "CI_upper")

  # Flags de advertencia
  loadings_df$Alert <- ""
  loadings_df$Alert[abs(loadings_df$Loading) > 1.0]  <- "*** HEYWOOD"
  loadings_df$Alert[abs(loadings_df$Loading) < 0.20 &
                    loadings_df$Alert == ""]         <- if (is_english) "** Very low (<.20)" else "** Muy baja (<.20)"
  loadings_df$Alert[abs(loadings_df$Loading) < 0.30 &
                    loadings_df$Alert == ""]         <- if (is_english) "* Low (<.30)" else "* Baja (<.30)"
  loadings_df$Lambda <- loadings_df$Loading
  loadings_df$p <- loadings_df$p_value
  loadings_df$IC_inf <- loadings_df$CI_lower
  loadings_df$IC_sup <- loadings_df$CI_upper
  loadings_df$Alerta <- loadings_df$Alert

  admissibility <- .pm_cfa_admissibility(
    fit,
    standardized_solution = std_sol,
    converged = converged
  )
  admissible <- isTRUE(admissibility$admissible)
  if (!admissible) {
    evaluable_metrics <- fit_summary$Metric %in%
      c("chisq_ratio", "cfi", "tli", "rmsea", "srmr")
    fit_summary$Evaluation[evaluable_metrics] <- if (is_english) {
      "Not evaluable"
    } else {
      "No evaluable"
    }
    fit_summary$Evaluacion <- fit_summary$Evaluation
  }

  # Solo filas con items observados (excluir segundo orden -> primer orden)
  loadings_obs <- loadings_df[loadings_df$Item %in% observed_items, ]
  loadings_ho  <- loadings_df[!loadings_df$Item %in% observed_items, ]

  # =========================================================================
  # 8. R-CUADRADO POR ITEM
  # =========================================================================
  r2 <- tryCatch(lavInspect(fit, "rsquare"), error = function(e) NULL)

  # =========================================================================
  # 9. CORRELACIONES ENTRE FACTORES
  # =========================================================================
  factor_corr_df <- NULL

  if (model_type_code == "correlated") {
    corr_rows <- std_sol[std_sol$op == "~~" &
                         std_sol$lhs != std_sol$rhs &
                         std_sol$lhs %in% factors &
                         std_sol$rhs %in% factors,
                         c("lhs", "rhs", "est.std", "se", "z", "pvalue")]
    if (nrow(corr_rows) > 0) {
      names(corr_rows) <- c("Factor1", "Factor2", "Correlation", "SE", "z", "p_value")
      corr_rows$Alert <- ""
      corr_rows$Alert[abs(corr_rows$Correlation) > 0.85] <- if (is_english) {
        "Discriminant validity?"
      } else {
        "Validez discriminante?"
      }
      corr_rows$r <- corr_rows$Correlation
      corr_rows$p <- corr_rows$p_value
      corr_rows$Alerta <- corr_rows$Alert
      factor_corr_df <- corr_rows
    }
  }

  if (is_second_order) {
    so_rows <- std_sol[std_sol$op == "=~" &
                       std_sol$lhs %in% second_order_factors,
                       c("lhs", "rhs", "est.std", "se", "z", "pvalue")]
    if (nrow(so_rows) > 0) {
      names(so_rows) <- c("HigherOrderFactor", "LowerOrderFactor", "Loading", "SE", "z", "p_value")
      so_rows$Factor_SO <- so_rows$HigherOrderFactor
      so_rows$Factor_PO <- so_rows$LowerOrderFactor
      so_rows$Lambda <- so_rows$Loading
      so_rows$p <- so_rows$p_value
    }
  } else {
    so_rows <- NULL
  }

  # =========================================================================
  # 10. CONFIABILIDAD
  # =========================================================================
  reliability_results <- if (!admissible) {
    list(
      type = if (is_bifactor) {
        "bifactor"
      } else if (is_second_order) {
        "second_order"
      } else {
        "standard"
      },
      status = "not_evaluable",
      reason = "inadmissible_cfa_solution",
      rows = data.frame(
        group = character(), level = character(), factor = character(),
        coefficient = character(), metric = character(), value = numeric(),
        status = character(), reason = character(), stringsAsFactors = FALSE
      )
    )
  } else {
    tryCatch({
      if (is_bifactor) {
        .compute_bifactor_reliability(fit, general_factor, specific_factors,
                                      observed_items)
      } else if (is_second_order) {
        .compute_second_order_reliability(fit, second_order_factors)
      } else {
        .compute_standard_reliability(fit)
      }
    }, error = function(e) {
      list(
        type = "error",
        status = "not_evaluable",
        reason = "reliability_error",
        message = conditionMessage(e),
        condition = e
      )
    })
  }

  # =========================================================================
  # 11. INDICES DE MODIFICACION Y REESPECIFICACIONES
  # =========================================================================
  mi_df <- tryCatch({
    mi <- modificationIndices(fit, sort. = TRUE, minimum.value = mi_threshold)
    if (nrow(mi) > 0) utils::head(mi, n_mi) else mi
  }, error = function(e) data.frame())

  respec <- .build_respec_suggestions(mi_df, items_per_factor, language)

  # =========================================================================
  # 12. DETECCION DE PROBLEMAS
  # =========================================================================
  warnings_list <- character(0)

  if (!converged) {
    warnings_list <- c(warnings_list,
      if (is_english) {
        "The model did NOT converge. Results are not trustworthy."
      } else {
        "El modelo NO convergio. Los resultados no son confiables."
      })
  }

  if (isTRUE(converged) && !admissible) {
    warnings_list <- c(
      warnings_list,
      if (is_english) {
        paste0(
          "The model converged, but the solution is inadmissible (",
          paste(admissibility$reason_codes, collapse = ", "),
          "). Fit interpretation and reliability are not evaluable."
        )
      } else {
        paste0(
          "El modelo convergio, pero la solucion es inadmisible (",
          paste(admissibility$reason_codes, collapse = ", "),
          "). El ajuste y la confiabilidad no son evaluables."
        )
      }
    )
  }

  # Varianzas negativas (Heywood)
  var_rows <- std_sol[std_sol$op == "~~" & std_sol$lhs == std_sol$rhs, ]
  neg_var <- var_rows[var_rows$est.std < 0, ]
  if (nrow(neg_var) > 0) {
    warnings_list <- c(warnings_list,
      if (is_english) {
        paste0("Heywood case(s) detected - negative variance in: ",
               paste(neg_var$lhs, collapse = ", "))
      } else {
        paste0("Caso(s) Heywood detectado(s) - varianza negativa en: ",
               paste(neg_var$lhs, collapse = ", "))
      })
  }

  # Cargas > 1
  heywood_loads <- loadings_df[abs(loadings_df$Loading) > 1.0, ]
  if (nrow(heywood_loads) > 0) {
    warnings_list <- c(warnings_list,
      if (is_english) {
        paste0("Factor loading(s) > 1.0 detected: ",
               paste(heywood_loads$Item, collapse = ", "))
      } else {
        paste0("Carga(s) factorial(es) > 1.0 detectada(s): ",
               paste(heywood_loads$Item, collapse = ", "))
      })
  }

  # Cargas muy bajas
  low_loads <- loadings_obs[abs(loadings_obs$Loading) < 0.30, ]
  if (nrow(low_loads) > 0) {
    warnings_list <- c(warnings_list,
      if (is_english) {
        paste0("Item(s) with loadings < .30: ",
               paste(low_loads$Item, collapse = ", "),
               " - consider removing them")
      } else {
        paste0("Item(s) con cargas < .30: ",
               paste(low_loads$Item, collapse = ", "),
               " - considerar eliminacion")
      })
  }

  if (is_bifactor && !isTRUE(bifactor_orthogonal_fitted)) {
    warnings_list <- c(
      warnings_list,
      if (is_english) {
        "The fitted bifactor is not structurally orthogonal; PUC is not evaluated."
      } else {
        "El bifactor ajustado no es estructuralmente ortogonal; no se evalua PUC."
      }
    )
  }

  warnings_structured <- .pm_bind_condition_tables(
    lavaan_warning = fitted_conditions$warnings,
    cfa_warning = warnings_list,
    mvn_warning = mvn_diagnostic$warnings,
    reliability_warning = reliability_results$warnings
  )
  errors_structured <- .pm_bind_condition_tables(
    mvn_error = mvn_diagnostic$condition,
    reliability_error = if (!is.null(reliability_results$condition)) {
      reliability_results$condition
    } else {
      reliability_results$semtools_condition
    }
  )
  warnings_list <- unique(warnings_structured$message)
  estimator_call <- match.call()$estimator
  estimator_requested <- if (is.null(estimator_call)) {
    "auto"
  } else {
    paste(deparse(estimator_call), collapse = " ")
  }
  decision_log <- rbind(
    .pm_decision(
      "estimator",
      estimator_requested,
      estimator,
      estimator_reason
    ),
    .pm_decision(
      "bifactor_orthogonality",
      as.character(bifactor_orthogonal),
      if (is_bifactor) as.character(bifactor_orthogonal_effective) else "not_applicable",
      if (is_bifactor) "explicit bifactor setting" else "model is not bifactor"
    ),
    .pm_decision(
      "solution_admissibility",
      "admissible",
      if (admissible) "admissible" else "not_evaluable",
      if (admissible) {
        "lavaan post-check and parameter-domain checks passed"
      } else {
        paste(admissibility$reason_codes, collapse = ", ")
      },
      if (admissible) "info" else "warning"
    )
  )

  # =========================================================================
  # 13. COMPILAR RESULTADOS
  # =========================================================================
  results <- list(
    # Objeto lavaan
    fit              = fit,
    converged        = converged,
    admissible       = admissible,
    admissibility    = admissibility,
    language         = language,
    # Modelo
    model_type       = model_type,
    model_type_code  = model_type_code,
    factors          = factors,
    items_per_factor = items_per_factor,
    observed_items   = observed_items,
    general_factor   = general_factor,
    specific_factors = specific_factors,
    second_order_factors = if (is_second_order) second_order_factors else NULL,
    so_loadings      = so_rows,
    # Datos
    n                = n_total,
    n_complete       = n_complete,
    pct_missing      = pct_missing,
    n_items          = length(observed_items),
    n_categories     = n_categories,
    pct_na_item      = pct_na_item,
    # Estimador
    estimator        = estimator,
    estimator_reason = estimator_reason,
    ordered          = ordered,
    bifactor_orthogonal_requested = bifactor_orthogonal,
    bifactor_orthogonal_effective = bifactor_orthogonal_effective,
    bifactor_orthogonal_fitted = bifactor_orthogonal_fitted,
    mvn_result       = mvn_result,
    mardia           = list(skew_stat = mardia_skew_stat, kurt_stat = mardia_kurt_stat,
                            p_skew = mardia_p_skew, p_kurt = mardia_p_kurt,
                            n = n_complete,
                            ordinal_codes = isTRUE(ordered) || (is.character(ordered) && length(ordered) > 0L),
                            status = mvn_diagnostic$status,
                            reason = mvn_diagnostic$reason),
    # Resultados
    fit_indices      = fit_summary,
    fit_summary      = fit_summary,
    fi_variant_note  = fi_variant_note,
    all_fit_measures = all_fi,
    loadings         = loadings_obs,
    factor_loadings  = loadings_obs,
    loadings_ho      = loadings_ho,
    higher_order_loadings = loadings_ho,
    r_squared        = r2,
    factor_correlations = factor_corr_df,
    factor_correlation_table = factor_corr_df,
    reliability      = reliability_results,
    modification_indices = mi_df,
    respecifications = respec,
    respecification_suggestions = respec,
    warnings         = warnings_structured,
    warnings_list    = warnings_list,
    errors           = errors_structured,
    requested_settings = list(
      ordered = match.call()$ordered,
      estimator = match.call()$estimator,
      fit_variant = fit_variant,
      std.lv = std.lv,
      mi_threshold = mi_threshold,
      n_mi = n_mi,
      alpha_norm = alpha_norm,
      language = match.call()$language,
      bifactor_orthogonal = bifactor_orthogonal
    ),
    effective_settings = list(
      ordered = ordered,
      estimator = estimator,
      std.lv = std.lv,
      fit_variant = fit_variant,
      language = language,
      alpha_norm = alpha_norm,
      bifactor_orthogonal = bifactor_orthogonal_effective
    ),
    data_summary = list(
      n = n_total,
      p = length(observed_items),
      complete = n_complete,
      missing_rows = n_missing_rows
    ),
    decision_log = decision_log,
    diagnostics = list(
      mvn = mvn_diagnostic,
      fit_variants = fit_summary[, c("Metric", "Base", "Scaled", "Robust", "Variant"), drop = FALSE],
      admissibility = admissibility,
      vcov = vcov_diagnostic,
      bifactor = list(
        detection_code = bifactor_detection$code,
        general_candidates = bifactor_detection$general_candidates,
        specialized_indices_supported = bifactor_detection$specialized_indices_supported,
        requested_orthogonal = bifactor_orthogonal,
        effective_orthogonal = bifactor_orthogonal_effective,
        fitted_orthogonal = bifactor_orthogonal_fitted
      )
    ),
    conclusion = if (!isTRUE(converged)) {
      .pm_not_evaluable("The CFA model did not converge.")
    } else if (!admissible) {
      .pm_not_evaluable(paste0(
        "The CFA solution is inadmissible: ",
        paste(admissibility$reason_codes, collapse = ", "),
        "."
      ))
    } else {
      NA_character_
    }
  )

  results <- .cfa_add_tidy_fields(results, call = match.call())
  results$results <- list(
    fit = results$fit,
    fit_indices = results$fit_indices,
    loadings = results$loadings,
    tidy_loadings = results$loadings_table,
    higher_order_loadings = results$higher_order_loadings,
    reliability = results$reliability,
    reliability_table = results$reliability_table,
    modification_indices = results$modification_indices,
    respecification_suggestions = results$respecification_suggestions
  )
  class(results) <- c("cfa_auto", "psychomatic_result")
  .pm_finish_analysis(results, output, language)
}
