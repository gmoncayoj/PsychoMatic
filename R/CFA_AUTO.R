# =============================================================================
# FUNCION PRINCIPAL
# =============================================================================

#' Automated CFA with psychometric reporting
#'
#' Runs a confirmatory factor analysis with heuristic estimator selection,
#' multivariate normality diagnostics, fit indices, reliability estimates, and
#' model respecification suggestions.
#'
#' @param data Data frame containing the observed variables.
#' @param model Model syntax written for `lavaan`.
#' @param ordered Optional vector indicating which items should be treated as ordinal.
#' @param estimator Optional estimator. If `NULL`, the function selects one automatically.
#' @param std.lv Whether latent variables should be standardized.
#' @param mi_threshold Minimum threshold for reporting modification indices.
#' @param n_mi Maximum number of modification indices to display.
#' @param alpha_norm Alpha level used when evaluating multivariate normality.
#' @param language Report language: `"esp"` or `"eng"`.
#'
#' @return Object of class `cfa_auto`.
#'
#' @examples
#' # model <- "F1 =~ i1 + i2 + i3"
#' # result <- cfa_auto(my_data, model)
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
                       language = "esp") {

  language <- .normalize_language(language)
  is_english <- identical(language, "eng")

  # --- Verificar paquetes ---
  pkgs <- c("lavaan", "semTools", "MVN")
  faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(faltantes) > 0) {
    if (is_english) {
      stop("Required packages are not installed: ",
           paste(faltantes, collapse = ", "),
           "\nInstall with: install.packages(c(",
           paste0('"', faltantes, '"', collapse = ", "), "))")
    }
    stop("Paquetes requeridos no instalados: ", paste(faltantes, collapse = ", "),
         "\nInstalar con: install.packages(c(",
         paste0('"', faltantes, '"', collapse = ", "), "))")
  }

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

  # Bifactor: un factor carga en TODOS los items observados + factores grupales en subconjuntos
  is_bifactor <- FALSE
  general_factor <- NULL
  specific_factors <- NULL

  if (!is_second_order && length(factors) > 2) {
    potential_general <- names(which(n_obs_per_factor == length(observed_items)))
    if (length(potential_general) >= 1) {
      remaining <- setdiff(names(n_obs_per_factor), potential_general)
      if (length(remaining) >= 2 && all(n_obs_per_factor[remaining] < length(observed_items))) {
        is_bifactor <- TRUE
        general_factor <- potential_general[1]
        specific_factors <- remaining
      }
    }
  }

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

  if (is.null(estimator)) {
    if (max_cat <= 5) {
      # ------ Items ordinales ------
      if (is.null(ordered)) ordered <- TRUE

      if (n_complete < 200) {
        estimator <- "ULSMV"
        estimator_reason <- if (is_english) {
          sprintf(
            "Ordinal items (%d-%d cat.) + n=%d (<200) -> ULSMV (recommended for small ordinal samples)",
            min_cat, max_cat, n_complete
          )
        } else {
          sprintf(
            "Items ordinales (%d-%d cat.) + n=%d (<200) -> ULSMV (recomendado para muestras pequenas ordinales)",
            min_cat, max_cat, n_complete
          )
        }
      } else {
        estimator <- "WLSMV"
        estimator_reason <- if (is_english) {
          sprintf(
            "Ordinal items (%d-%d cat.) + n=%d -> WLSMV (robust WLS for ordinal data)",
            min_cat, max_cat, n_complete
          )
        } else {
          sprintf(
            "Items ordinales (%d-%d cat.) + n=%d -> WLSMV (WLS robusto para datos ordinales)",
            min_cat, max_cat, n_complete
          )
        }
      }

    } else if (max_cat <= 7 && min_cat <= 5) {
      # ------ Mezcla ------
      if (is.null(ordered)) ordered <- TRUE
      estimator <- "WLSMV"
      estimator_reason <- if (is_english) {
        sprintf(
          "Mixed item categories (%d-%d) -> WLSMV (ordinal treatment recommended)",
          min_cat, max_cat
        )
      } else {
        sprintf(
          "Mezcla de items con %d-%d categorias -> WLSMV (tratamiento ordinal recomendado)",
          min_cat, max_cat
        )
      }

    } else {
      # ------ Items continuos (>= 6 categorias) ------
      if (is.null(ordered)) ordered <- FALSE

      mvn_result <- tryCatch({
        MVN::mvn(data_items[complete.cases(data_items), ], mvnTest = "mardia")
      }, error = function(e) NULL)

      if (!is.null(mvn_result)) {
        mardia_p_skew <- as.numeric(mvn_result$multivariateNormality$`p value`[1])
        mardia_p_kurt <- as.numeric(mvn_result$multivariateNormality$`p value`[2])
        mardia_skew_stat <- as.numeric(mvn_result$multivariateNormality$Statistic[1])
        mardia_kurt_stat <- as.numeric(mvn_result$multivariateNormality$Statistic[2])
        is_mvn <- !is.na(mardia_p_skew) && !is.na(mardia_p_kurt) &&
                  (mardia_p_skew > alpha_norm) && (mardia_p_kurt > alpha_norm)

        if (is_mvn) {
          estimator <- "ML"
          estimator_reason <- if (is_english) {
            sprintf(
              "Continuous items (%d-%d cat.), multivariate normality NOT rejected (Mardia: p_skew=%.4f, p_kurt=%.4f) -> ML",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          } else {
            sprintf(
              "Items continuos (%d-%d cat.), normalidad multivariada NO rechazada (Mardia: p_asim=%.4f, p_curt=%.4f) -> ML",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          }
        } else {
          estimator <- "MLR"
          estimator_reason <- if (is_english) {
            sprintf(
              "Continuous items (%d-%d cat.), multivariate normality REJECTED (Mardia: p_skew=%.4f, p_kurt=%.4f) -> MLR",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          } else {
            sprintf(
              "Items continuos (%d-%d cat.), normalidad multivariada RECHAZADA (Mardia: p_asim=%.4f, p_curt=%.4f) -> MLR",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          }
        }
      } else {
        estimator <- "MLR"
        estimator_reason <- if (is_english) {
          sprintf(
            "Continuous items (%d-%d cat.), normality test unavailable -> MLR as a precaution",
            min_cat, max_cat
          )
        } else {
          sprintf(
            "Items continuos (%d-%d cat.), prueba de normalidad no disponible -> MLR como precaucion",
            min_cat, max_cat
          )
        }
      }
    }
  } else {
    estimator_reason <- if (is_english) {
      "Estimator specified by the user"
    } else {
      "Estimador especificado por el usuario"
    }
    if (is.null(ordered)) {
      ordered <- if (max_cat <= 5) TRUE else FALSE
    }
  }

  # =========================================================================
  # 5. AJUSTAR MODELO
  # =========================================================================
  fit_args <- list(model = model, data = data, estimator = estimator, std.lv = std.lv)

  if (isTRUE(ordered)) {
    fit_args$ordered <- TRUE
  } else if (is.character(ordered)) {
    fit_args$ordered <- ordered
  }

  # Bifactor: factores ortogonales (general y especificos no correlacionados)
  if (is_bifactor) {
    fit_args$orthogonal <- TRUE
  }

  fit <- do.call(lavaan::cfa, fit_args)
  converged <- lavInspect(fit, "converged")

  # =========================================================================
  # 6. INDICES DE AJUSTE
  # =========================================================================
  all_fi <- tryCatch(fitMeasures(fit), error = function(e) NULL)

  if (!is.null(all_fi)) {
    # Seleccionar la mejor version disponible: scaled > robust > base
    # Prioriza scaled (recomendado en la literatura para CFI, TLI, RMSEA)
    # Devuelve lista con $value y $variant ("scaled", "robust", o "base")
    .pick <- function(base, fi_vec) {
      candidates <- list(
        list(name = paste0(base, ".scaled"), variant = "scaled"),
        list(name = paste0(base, ".robust"), variant = "robust"),
        list(name = base,                   variant = "base")
      )
      for (cand in candidates) {
        if (cand$name %in% names(fi_vec) && !is.na(fi_vec[cand$name])) {
          return(list(value = unname(fi_vec[cand$name]), variant = cand$variant))
        }
      }
      return(list(value = NA_real_, variant = "base"))
    }

    # Etiqueta con sufijo segun variante
    .label <- function(base_label, variant) {
      switch(variant,
             robust  = paste0(base_label, " (robust)"),
             scaled  = paste0(base_label, " (scaled)"),
             base_label)
    }

    # Extraer cada indice
    p_chisq  <- .pick("chisq", all_fi)
    p_df     <- .pick("df", all_fi)
    p_pval   <- .pick("pvalue", all_fi)
    p_cfi    <- .pick("cfi", all_fi)
    p_tli    <- .pick("tli", all_fi)
    p_rmsea  <- .pick("rmsea", all_fi)
    p_rmsea_lo <- .pick("rmsea.ci.lower", all_fi)
    p_rmsea_hi <- .pick("rmsea.ci.upper", all_fi)
    p_rmsea_p  <- .pick("rmsea.pvalue", all_fi)

    # SRMR: no tiene version scaled/robust tipicamente; probar srmr_bentler tambien
    p_srmr <- .pick("srmr", all_fi)
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
  reliability_results <- tryCatch({

    if (is_bifactor) {
      .compute_bifactor_reliability(fit, general_factor, specific_factors,
                                    observed_items)
    } else if (is_second_order) {
      .compute_second_order_reliability(fit, second_order_factors)
    } else {
      .compute_standard_reliability(fit)
    }

  }, error = function(e) {
    list(type = "error", message = conditionMessage(e))
  })

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

  # =========================================================================
  # 13. COMPILAR RESULTADOS
  # =========================================================================
  results <- list(
    # Objeto lavaan
    fit              = fit,
    converged        = converged,
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
    mvn_result       = mvn_result,
    mardia           = list(skew_stat = mardia_skew_stat, kurt_stat = mardia_kurt_stat,
                            p_skew = mardia_p_skew, p_kurt = mardia_p_kurt),
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
    warnings         = warnings_list,
    warnings_list    = warnings_list
  )

  class(results) <- "cfa_auto"
  return(results)
}


# =============================================================================
# FUNCIONES AUXILIARES (internas)
# =============================================================================

# Normalizar idioma del reporte
.normalize_language <- function(language = "esp") {
  if (is.null(language) || length(language) == 0) return("esp")

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
.model_type_label <- function(model_type_code, language = "esp") {
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
.fit_index_base_label <- function(metric, language = "esp") {
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

.fit_index_label <- function(metric, variant = "base", language = "esp") {
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
.cfa_auto_labels <- function(language = "esp") {
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
      ci_label = "CI",
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
    ci_label = "IC",
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
.evaluar_fi <- function(valor, bueno, aceptable, direccion = ">=", language = "esp") {
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

# Confiabilidad estandar (unidimensional o factores correlacionados)
.compute_standard_reliability <- function(fit) {
  rel <- semTools::reliability(fit)
  list(
    type   = "standard",
    values = rel
  )
}

# Confiabilidad para segundo orden
.compute_second_order_reliability <- function(fit, so_factors) {
  rel_l1 <- semTools::reliability(fit)
  # Calcular reliabilityL2 para cada factor de segundo orden
  rel_l2 <- lapply(setNames(so_factors, so_factors), function(sof) {
    tryCatch(
      semTools::reliabilityL2(fit, secondFactor = sof),
      error = function(e) NULL
    )
  })
  list(
    type   = "second_order",
    first_order = rel_l1,
    second_order = rel_l2,
    so_factors = so_factors
  )
}

# Confiabilidad para bifactor
.compute_bifactor_reliability <- function(fit, general_factor, specific_factors,
                                          observed_items) {
  rel <- semTools::reliability(fit)

  # Extraer parametros estandarizados
  std_est <- lavaan::inspect(fit, what = "std")
  lambda  <- std_est$lambda
  theta   <- diag(std_est$theta)

  # Cargas del factor general
  gen_loads <- lambda[, general_factor]

  # --- Omega total ---
  sum_gen  <- sum(gen_loads)
  sum_spec <- sapply(specific_factors, function(sf) sum(lambda[, sf]))
  sum_theta <- sum(theta)
  total_var <- sum_gen^2 + sum(sum_spec^2) + sum_theta

  omega_total <- (sum_gen^2 + sum(sum_spec^2)) / total_var

  # --- Omega jerarquico (general) ---
  omega_h <- sum_gen^2 / total_var

  # --- ECV (Explained Common Variance) ---
  ss_gen  <- sum(gen_loads^2)
  ss_spec <- sum(sapply(specific_factors, function(sf) sum(lambda[, sf]^2)))
  ecv <- ss_gen / (ss_gen + ss_spec)

  # --- ECV por item ---
  ecv_item <- gen_loads^2 / (gen_loads^2 + rowSums(lambda[, specific_factors, drop = FALSE]^2))

  # --- PUC (Percent of Uncontaminated Correlations) ---
  p <- length(observed_items)
  n_spec <- sapply(specific_factors, function(sf) sum(lambda[, sf] != 0))
  puc <- (p * (p - 1) / 2 - sum(n_spec * (n_spec - 1) / 2)) / (p * (p - 1) / 2)

  # --- Omega_hs por factor especifico ---
  omega_hs <- sapply(specific_factors, function(sf) {
    items_sf <- which(lambda[, sf] != 0)
    gl <- gen_loads[items_sf]
    sl <- lambda[items_sf, sf]
    th <- theta[items_sf]
    var_sub <- sum(gl)^2 + sum(sl)^2 + sum(th)
    sum(sl)^2 / var_sub
  })

  # --- H index (construct replicability) ---
  h_general <- .compute_h_index(gen_loads[gen_loads != 0])
  h_specific <- sapply(specific_factors, function(sf) {
    sl <- lambda[lambda[, sf] != 0, sf]
    .compute_h_index(sl)
  })

  # --- FD (Factor Determinacy) ---
  fd <- tryCatch({
    scores_info <- lavInspect(fit, "fscoefficient")
    # Correlacion entre factor y sus scores
    NULL  # Se puede obtener de otra forma
  }, error = function(e) NULL)

  list(
    type             = "bifactor",
    semtools_rel     = rel,
    omega_total      = omega_total,
    omega_h          = omega_h,
    omega_hs         = omega_hs,
    ecv              = ecv,
    ecv_item         = ecv_item,
    puc              = puc,
    h_general        = h_general,
    h_specific       = h_specific,
    general_factor   = general_factor,
    specific_factors = specific_factors
  )
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
.build_respec_suggestions <- function(mi_df, items_per_factor, language = "esp") {
  language <- .normalize_language(language)
  is_english <- identical(language, "eng")
  respec <- list()
  if (is.null(mi_df) || nrow(mi_df) == 0) return(respec)

  for (i in seq_len(nrow(mi_df))) {
    row <- mi_df[i, ]

    if (row$op == "~~" && row$lhs != row$rhs) {
      # Covarianza residual
      same_f <- any(sapply(items_per_factor, function(items) {
        row$lhs %in% items && row$rhs %in% items
      }))
      respec[[length(respec) + 1]] <- list(
        type    = if (is_english) "Residual covariance" else "Covarianza residual",
        syntax = paste0(row$lhs, " ~~ ", row$rhs),
        mi      = round(row$mi, 2),
        epc     = round(row$epc, 3),
        same_factor = same_f,
        note    = if (same_f) {
          if (is_english) {
            "Same factor: consider only if there is theoretical justification (similar wording, method effects, etc.)"
          } else {
            "Mismo factor: considerar si hay justificacion teorica (items similares en redaccion, metodo, etc.)"
          }
        } else {
          if (is_english) {
            "Different factor: evaluate cautiously, it may indicate an unmodeled latent factor"
          } else {
            "Diferente factor: evaluar con cautela, puede indicar factor latente no modelado"
          }
        }
      )
      respec[[length(respec)]]$tipo <- respec[[length(respec)]]$type
      respec[[length(respec)]]$sintaxis <- respec[[length(respec)]]$syntax
      respec[[length(respec)]]$mismo_factor <- respec[[length(respec)]]$same_factor
      respec[[length(respec)]]$nota <- respec[[length(respec)]]$note

    } else if (row$op == "=~") {
      respec[[length(respec) + 1]] <- list(
        type    = if (is_english) "Cross-loading" else "Carga cruzada",
        syntax = paste0(row$lhs, " =~ ", row$rhs),
        mi      = round(row$mi, 2),
        epc     = round(row$epc, 3),
        same_factor = FALSE,
        note    = if (is_english) {
          "Potential cross-loading: review whether the item content justifies belonging to this factor"
        } else {
          "Carga cruzada potencial: revisar si el contenido del item justifica pertenencia a este factor"
        }
      )
      respec[[length(respec)]]$tipo <- respec[[length(respec)]]$type
      respec[[length(respec)]]$sintaxis <- respec[[length(respec)]]$syntax
      respec[[length(respec)]]$mismo_factor <- respec[[length(respec)]]$same_factor
      respec[[length(respec)]]$nota <- respec[[length(respec)]]$note
    }
  }

  return(respec)
}


# =============================================================================
# METODO PRINT
# =============================================================================

.print_cfa_auto_legacy <- function(x, digits = 3, ...) {

  sep  <- paste(rep("=", 70), collapse = "")
  sep2 <- paste(rep("-", 70), collapse = "")

  cat("\n", sep, "\n", sep = "")
  cat("       ANALISIS FACTORIAL CONFIRMATORIO - REPORTE PSICOMETRICO\n")
  cat(sep, "\n\n", sep = "")

  # --- ADVERTENCIAS ---
  if (length(x$warnings) > 0) {
    cat("*** ADVERTENCIAS ***\n")
    cat(sep2, "\n", sep = "")
    for (w in x$warnings) {
      cat("  [!] ", w, "\n", sep = "")
    }
    cat("\n")
  }

  # --- INFORMACION DEL MODELO ---
  cat("INFORMACION DEL MODELO\n")
  cat(sep2, "\n", sep = "")
  cat("  Tipo de modelo       :", x$model_type, "\n")
  cat("  Factores             :", paste(x$factors, collapse = ", "), "\n")
  cat("  N items observados   :", x$n_items, "\n")
  cat("  N (total)            :", x$n, "\n")
  if (x$pct_missing > 0) {
    cat("  N (completos)        :", x$n_complete,
        sprintf("(%.1f%% con datos faltantes)\n", x$pct_missing))
  }
  cat("  Categorias de resp.  :", paste0(min(x$n_categories), "-",
                                          max(x$n_categories)), "\n")
  cat("  Estimador            :", x$estimator, "\n")
  cat("  Justificacion        :", x$estimator_reason, "\n")
  cat("  Convergio            :",
      ifelse(x$converged, "Si", "*** NO ***"), "\n")
  cat("\n")

  # --- NORMALIDAD MULTIVARIADA ---
  if (!is.na(x$mardia$p_skew)) {
    cat("NORMALIDAD MULTIVARIADA (Mardia)\n")
    cat(sep2, "\n", sep = "")
    cat(sprintf("  Asimetria  : Estadistico = %.2f, p = %.4f  %s\n",
                x$mardia$skew_stat, x$mardia$p_skew,
                ifelse(x$mardia$p_skew < 0.05, "[Rechazada]", "[No rechazada]")))
    cat(sprintf("  Curtosis   : Estadistico = %.2f, p = %.4f  %s\n",
                x$mardia$kurt_stat, x$mardia$p_kurt,
                ifelse(x$mardia$p_kurt < 0.05, "[Rechazada]", "[No rechazada]")))
    cat("\n")
  }

  # --- INDICES DE AJUSTE ---
  cat("INDICES DE AJUSTE\n")
  cat(sep2, "\n", sep = "")
  cat("  Nota:", x$fi_variant_note, "\n\n")
  fi <- x$fit_indices
  for (i in seq_len(nrow(fi))) {
    val <- fi$Valor[i]
    ev  <- fi$Evaluacion[i]
    if (is.na(val)) next
    eval_str <- if (nchar(ev) > 0) paste0("  [", ev, "]") else ""
    if (fi$Indice[i] %in% c("gl")) {
      cat(sprintf("  %-22s: %d%s\n", fi$Indice[i], as.integer(val), eval_str))
    } else {
      cat(sprintf("  %-22s: %.3f%s\n", fi$Indice[i], val, eval_str))
    }
  }
  cat("\n")

  # --- CARGAS FACTORIALES ---
  cat("CARGAS FACTORIALES ESTANDARIZADAS\n")
  cat(sep2, "\n", sep = "")
  ld <- x$loadings
  for (f in unique(ld$Factor)) {
    cat(sprintf("\n  Factor: %s\n", f))
    sub <- ld[ld$Factor == f, ]
    for (j in seq_len(nrow(sub))) {
      alerta <- if (nchar(sub$Alerta[j]) > 0) paste0("  ", sub$Alerta[j]) else ""
      cat(sprintf("    %-18s %6.3f  (SE=%.3f, z=%.2f, p=%.4f)  IC[%.3f, %.3f]%s\n",
                  sub$Item[j], sub$Lambda[j], sub$SE[j], sub$z[j], sub$p[j],
                  sub$IC_inf[j], sub$IC_sup[j], alerta))
    }
  }
  cat("\n")

  # --- Cargas segundo orden (si aplica) ---
  if (!is.null(x$so_loadings) && nrow(x$so_loadings) > 0) {
    cat("CARGAS DE SEGUNDO ORDEN\n")
    cat(sep2, "\n", sep = "")
    so <- x$so_loadings
    for (j in seq_len(nrow(so))) {
      cat(sprintf("  %s -> %s : %.3f  (SE=%.3f, z=%.2f, p=%.4f)\n",
                  so$Factor_SO[j], so$Factor_PO[j], so$Lambda[j],
                  so$SE[j], so$z[j], so$p[j]))
    }
    cat("\n")
  }

  # --- R-CUADRADO ---
  if (!is.null(x$r_squared)) {
    cat("R-CUADRADO POR ITEM\n")
    cat(sep2, "\n", sep = "")
    r2 <- x$r_squared
    if (is.list(r2)) {
      for (nm in names(r2)) {
        r2_vec <- r2[[nm]]
        for (j in seq_along(r2_vec)) {
          cat(sprintf("    %-18s %.3f\n", names(r2_vec)[j], r2_vec[j]))
        }
      }
    } else {
      for (j in seq_along(r2)) {
        cat(sprintf("    %-18s %.3f\n", names(r2)[j], r2[j]))
      }
    }
    cat("\n")
  }

  # --- CORRELACIONES ENTRE FACTORES ---
  if (!is.null(x$factor_correlations)) {
    cat("CORRELACIONES ENTRE FACTORES\n")
    cat(sep2, "\n", sep = "")
    fc <- x$factor_correlations
    for (j in seq_len(nrow(fc))) {
      alerta <- if (nchar(fc$Alerta[j]) > 0) paste0("  *** ", fc$Alerta[j]) else ""
      cat(sprintf("  %-12s <-> %-12s : %6.3f  (SE=%.3f, p=%.4f)%s\n",
                  fc$Factor1[j], fc$Factor2[j], fc$r[j], fc$SE[j], fc$p[j], alerta))
    }
    cat("\n")
  }

  # --- CONFIABILIDAD ---
  cat("CONFIABILIDAD\n")
  cat(sep2, "\n", sep = "")

  rel <- x$reliability

  if (rel$type == "error") {
    cat("  Error al calcular: ", rel$message, "\n")

  } else if (rel$type == "standard") {
    cat("  Metodo: Omega de McDonald y Alfa de Cronbach (semTools::reliability)\n\n")
    rv <- rel$values
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    etiquetas <- c(omega = "Omega", alpha = "Alfa", avevar = "AVE")
    for (f in colnames(rv)) {
      if (f == "total" && ncol(rv) > 2) next
      cat(sprintf("  Factor: %s\n", f))
      for (m in show_rows) {
        cat(sprintf("    %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
      cat("\n")
    }

  } else if (rel$type == "second_order") {
    cat("  Modelo de Segundo Orden\n\n")
    cat("  -- Primer Orden (por dimension) --\n")
    rv <- rel$first_order
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    etiquetas <- c(omega = "Omega", alpha = "Alfa", avevar = "AVE")
    for (f in colnames(rv)) {
      if (f == "total") next
      cat(sprintf("    Factor: %s\n", f))
      for (m in show_rows) {
        cat(sprintf("      %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
    }
    cat("\n  -- Segundo Orden --\n")
    if (!is.null(rel$second_order)) {
      for (sof in rel$so_factors) {
        l2 <- rel$second_order[[sof]]
        if (!is.null(l2)) {
          cat(sprintf("    Factor: %s\n", sof))
          for (nm in names(l2)) {
            cat(sprintf("      %-25s: %.3f\n", nm, l2[nm]))
          }
          cat("\n")
        } else {
          cat(sprintf("    Factor: %s - No se pudo calcular\n", sof))
        }
      }
    } else {
      cat("    No se pudo calcular (revisar especificacion del modelo)\n")
    }
    cat("\n")

  } else if (rel$type == "bifactor") {
    cat("  Modelo Bifactor - Indices Especializados\n\n")
    cat(sprintf("  Omega total             : %.3f\n", rel$omega_total))
    cat(sprintf("  Omega jerarquico (w_h)  : %.3f  [%s]\n",
                rel$omega_h, rel$general_factor))
    cat(sprintf("  ECV                     : %.3f\n", rel$ecv))
    cat(sprintf("  PUC                     : %.3f\n", rel$puc))
    cat(sprintf("  H (factor general)      : %.3f\n", rel$h_general))

    cat("\n  Por factor especifico:\n")
    cat(sprintf("  %-20s  %8s  %8s\n", "Factor", "Omega_hs", "H"))
    for (i in seq_along(rel$omega_hs)) {
      cat(sprintf("  %-20s  %8.3f  %8.3f\n",
                  names(rel$omega_hs)[i], rel$omega_hs[i], rel$h_specific[i]))
    }

    cat("\n  Interpretacion:\n")
    if (rel$omega_h >= 0.80) {
      cat("    Omega_h >= .80: Factor general FUERTE.\n")
      cat("    -> La puntuacion total refleja predominantemente el factor general.\n")
    } else if (rel$omega_h >= 0.50) {
      cat("    Omega_h entre .50-.80: Factor general MODERADO.\n")
      cat("    -> Evaluar si las subescalas aportan informacion incremental.\n")
    } else {
      cat("    Omega_h < .50: Factor general DEBIL.\n")
      cat("    -> Las puntuaciones por subescala pueden ser mas informativas.\n")
    }
    # Interpretación conjunta ECV/PUC basada en Rodriguez et al. (2016) y Reise et al. (2013)
    if (rel$puc > 0.80 && rel$ecv > 0.60 && rel$omega_h >= 0.70) {
      cat("    PUC > .80, ECV > .60 y Omega_h >= .70: Evidencia favorable a la\n")
      cat("    unidimensionalidad (Reise et al., 2013).\n")
    } else if (rel$ecv > 0.70 && rel$puc > 0.70) {
      cat("    ECV > .70 y PUC > .70: Se podria concluir a favor de la\n")
      cat("    unidimensionalidad (Rodriguez et al., 2016).\n")
    } else if (rel$ecv < 0.60) {
      cat("    ECV < .60: Multidimensionalidad sustancial presente.\n")
    } else {
      cat(sprintf("    ECV = %.3f, PUC = %.3f: No se cumplen los umbrales para\n", rel$ecv, rel$puc))
      cat("    concluir a favor de la unidimensionalidad. Interpretar subescalas\n")
      cat("    con cautela.\n")
    }
    for (i in seq_along(rel$h_specific)) {
      if (rel$h_specific[i] < 0.50) {
        cat(sprintf("    H(%s) < .50: Factor especifico poco replicable.\n",
                    names(rel$h_specific)[i]))
      }
    }
    cat("\n")
  }

  # --- REESPECIFICACIONES ---
  if (length(x$respecifications) > 0) {
    cat("SUGERENCIAS DE REESPECIFICACION\n")
    cat(sep2, "\n", sep = "")
    cat("  (Solo considerar con justificacion teorica/sustantiva)\n\n")
    for (i in seq_along(x$respecifications)) {
      r <- x$respecifications[[i]]
      cat(sprintf("  %d. [%s] %s\n", i, r$tipo, r$sintaxis))
      cat(sprintf("     IM = %.2f, EPC = %.3f\n", r$mi, r$epc))
      cat(sprintf("     %s\n\n", r$nota))
    }
  } else {
    cat("REESPECIFICACIONES\n")
    cat(sep2, "\n", sep = "")
    cat("  No se encontraron indices de modificacion superiores al umbral.\n\n")
  }

  # --- REFERENCIAS BIBLIOGRAFICAS ---
  cat("REFERENCIAS BIBLIOGRAFICAS (APA 7)\n")
  cat(sep2, "\n", sep = "")
  cat("  (Se incluyen solo las referencias pertinentes al analisis realizado)\n\n")
  refs <- .build_references(x)
  for (i in seq_along(refs)) {
    cat(sprintf("  [%d] %s\n\n", i, refs[i]))
  }

  cat(sep, "\n", sep = "")
  cat("Objeto lavaan accesible en $fit | Todos los indices en $all_fit_measures\n")
  cat(sep, "\n\n", sep = "")

  invisible(x)
}

#' Print a `cfa_auto` object
#'
#' @param x Object returned by `cfa_auto()`.
#' @param digits Number of decimal places shown in the printed report.
#' @param ... Additional unused arguments.
#'
#' @rdname cfa_auto
#' @method print cfa_auto
#' @export
print.cfa_auto <- function(x, digits = 3, ...) {
  language <- if (!is.null(x$language)) .normalize_language(x$language) else "esp"
  is_english <- identical(language, "eng")
  txt <- .cfa_auto_labels(language)

  sep  <- paste(rep("=", 70), collapse = "")
  sep2 <- paste(rep("-", 70), collapse = "")

  cat("\n", sep, "\n", sep = "")
  cat(txt$report_title, "\n")
  cat(sep, "\n\n", sep = "")

  if (length(x$warnings) > 0) {
    cat(txt$warnings_title, "\n")
    cat(sep2, "\n", sep = "")
    for (w in x$warnings) {
      cat("  [!] ", w, "\n", sep = "")
    }
    cat("\n")
  }

  cat(txt$model_info_title, "\n")
  cat(sep2, "\n", sep = "")
  cat(txt$model_type_label, x$model_type, "\n")
  cat(txt$factors_label, paste(x$factors, collapse = ", "), "\n")
  cat(txt$observed_items_label, x$n_items, "\n")
  cat(txt$n_total_label, x$n, "\n")
  if (x$pct_missing > 0) {
    cat(txt$n_complete_label, x$n_complete, txt$missing_suffix(x$pct_missing), "\n")
  }
  cat(txt$categories_label, paste0(min(x$n_categories), "-", max(x$n_categories)), "\n")
  cat(txt$estimator_label, x$estimator, "\n")
  cat(txt$justification_label, x$estimator_reason, "\n")
  cat(txt$converged_label, ifelse(x$converged, txt$yes, txt$no), "\n")
  cat("\n")

  if (!is.na(x$mardia$p_skew)) {
    cat(txt$mardia_title, "\n")
    cat(sep2, "\n", sep = "")
    cat(sprintf("%s = %.2f, p = %.4f  %s\n",
                txt$skew_label,
                x$mardia$skew_stat, x$mardia$p_skew,
                ifelse(x$mardia$p_skew < 0.05, txt$rejected, txt$not_rejected)))
    cat(sprintf("%s = %.2f, p = %.4f  %s\n",
                txt$kurt_label,
                x$mardia$kurt_stat, x$mardia$p_kurt,
                ifelse(x$mardia$p_kurt < 0.05, txt$rejected, txt$not_rejected)))
    cat("\n")
  }

  cat(txt$fit_title, "\n")
  cat(sep2, "\n", sep = "")
  cat(txt$note_label, x$fi_variant_note, "\n\n")
  fi <- x$fit_indices
  metric_col <- if ("Metric" %in% names(fi)) fi$Metric else rep(NA_character_, nrow(fi))
  for (i in seq_len(nrow(fi))) {
    val <- fi$Valor[i]
    ev  <- fi$Evaluacion[i]
    if (is.na(val)) next
    eval_str <- if (nchar(ev) > 0) paste0("  [", ev, "]") else ""
    if ((!is.na(metric_col[i]) && metric_col[i] == "df") || fi$Indice[i] %in% c("gl", "df")) {
      cat(sprintf("  %-22s: %d%s\n", fi$Indice[i], as.integer(val), eval_str))
    } else {
      cat(sprintf("  %-22s: %.3f%s\n", fi$Indice[i], val, eval_str))
    }
  }
  cat("\n")

  cat(txt$loadings_title, "\n")
  cat(sep2, "\n", sep = "")
  ld <- x$loadings
  for (f in unique(ld$Factor)) {
    cat(sprintf("\n%s %s\n", txt$factor_label, f))
    sub <- ld[ld$Factor == f, ]
    for (j in seq_len(nrow(sub))) {
      alerta <- if (nchar(sub$Alerta[j]) > 0) paste0("  ", sub$Alerta[j]) else ""
      cat(sprintf("    %-18s %6.3f  (SE=%.3f, z=%.2f, p=%.4f)  %s[%.3f, %.3f]%s\n",
                  sub$Item[j], sub$Lambda[j], sub$SE[j], sub$z[j], sub$p[j],
                  txt$ci_label, sub$IC_inf[j], sub$IC_sup[j], alerta))
    }
  }
  cat("\n")

  if (!is.null(x$so_loadings) && nrow(x$so_loadings) > 0) {
    cat(txt$second_order_loadings_title, "\n")
    cat(sep2, "\n", sep = "")
    so <- x$so_loadings
    for (j in seq_len(nrow(so))) {
      cat(sprintf("  %s -> %s : %.3f  (SE=%.3f, z=%.2f, p=%.4f)\n",
                  so$Factor_SO[j], so$Factor_PO[j], so$Lambda[j],
                  so$SE[j], so$z[j], so$p[j]))
    }
    cat("\n")
  }

  if (!is.null(x$r_squared)) {
    cat(txt$r2_title, "\n")
    cat(sep2, "\n", sep = "")
    r2 <- x$r_squared
    if (is.list(r2)) {
      for (nm in names(r2)) {
        r2_vec <- r2[[nm]]
        for (j in seq_along(r2_vec)) {
          cat(sprintf("    %-18s %.3f\n", names(r2_vec)[j], r2_vec[j]))
        }
      }
    } else {
      for (j in seq_along(r2)) {
        cat(sprintf("    %-18s %.3f\n", names(r2)[j], r2[j]))
      }
    }
    cat("\n")
  }

  if (!is.null(x$factor_correlations)) {
    cat(txt$factor_corr_title, "\n")
    cat(sep2, "\n", sep = "")
    fc <- x$factor_correlations
    for (j in seq_len(nrow(fc))) {
      alerta <- if (nchar(fc$Alerta[j]) > 0) paste0("  *** ", fc$Alerta[j]) else ""
      cat(sprintf("  %-12s <-> %-12s : %6.3f  (SE=%.3f, p=%.4f)%s\n",
                  fc$Factor1[j], fc$Factor2[j], fc$r[j], fc$SE[j], fc$p[j], alerta))
    }
    cat("\n")
  }

  cat(txt$reliability_title, "\n")
  cat(sep2, "\n", sep = "")

  rel <- x$reliability

  if (rel$type == "error") {
    if (is_english) {
      cat("  Error while computing: ", rel$message, "\n")
    } else {
      cat("  Error al calcular: ", rel$message, "\n")
    }

  } else if (rel$type == "standard") {
    if (is_english) {
      cat("  Method: McDonald's Omega and Cronbach's Alpha (semTools::reliability)\n\n")
      etiquetas <- c(omega = "Omega", alpha = "Alpha", avevar = "AVE")
    } else {
      cat("  Metodo: Omega de McDonald y Alfa de Cronbach (semTools::reliability)\n\n")
      etiquetas <- c(omega = "Omega", alpha = "Alfa", avevar = "AVE")
    }
    rv <- rel$values
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    for (f in colnames(rv)) {
      if (f == "total" && ncol(rv) > 2) next
      cat(sprintf("%s %s\n", txt$factor_label, f))
      for (m in show_rows) {
        cat(sprintf("    %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
      cat("\n")
    }

  } else if (rel$type == "second_order") {
    if (is_english) {
      cat("  Second-order model\n\n")
      cat("  -- First-order (by dimension) --\n")
      etiquetas <- c(omega = "Omega", alpha = "Alpha", avevar = "AVE")
    } else {
      cat("  Modelo de Segundo Orden\n\n")
      cat("  -- Primer Orden (por dimension) --\n")
      etiquetas <- c(omega = "Omega", alpha = "Alfa", avevar = "AVE")
    }
    rv <- rel$first_order
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    for (f in colnames(rv)) {
      if (f == "total") next
      cat(sprintf("    Factor: %s\n", f))
      for (m in show_rows) {
        cat(sprintf("      %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
    }
    if (is_english) {
      cat("\n  -- Second-order --\n")
    } else {
      cat("\n  -- Segundo Orden --\n")
    }
    if (!is.null(rel$second_order)) {
      for (sof in rel$so_factors) {
        l2 <- rel$second_order[[sof]]
        if (!is.null(l2)) {
          cat(sprintf("    Factor: %s\n", sof))
          for (nm in names(l2)) {
            cat(sprintf("      %-25s: %.3f\n", nm, l2[nm]))
          }
          cat("\n")
        } else if (is_english) {
          cat(sprintf("    Factor: %s - Could not be computed\n", sof))
        } else {
          cat(sprintf("    Factor: %s - No se pudo calcular\n", sof))
        }
      }
    } else if (is_english) {
      cat("    Could not be computed (check model specification)\n")
    } else {
      cat("    No se pudo calcular (revisar especificacion del modelo)\n")
    }
    cat("\n")

  } else if (rel$type == "bifactor") {
    if (is_english) {
      cat("  Bifactor model - Specialized indices\n\n")
      cat(sprintf("  Total omega             : %.3f\n", rel$omega_total))
      cat(sprintf("  Hierarchical omega (w_h): %.3f  [%s]\n",
                  rel$omega_h, rel$general_factor))
      cat(sprintf("  ECV                     : %.3f\n", rel$ecv))
      cat(sprintf("  PUC                     : %.3f\n", rel$puc))
      cat(sprintf("  H (general factor)      : %.3f\n", rel$h_general))
      cat("\n  By specific factor:\n")
    } else {
      cat("  Modelo Bifactor - Indices Especializados\n\n")
      cat(sprintf("  Omega total             : %.3f\n", rel$omega_total))
      cat(sprintf("  Omega jerarquico (w_h)  : %.3f  [%s]\n",
                  rel$omega_h, rel$general_factor))
      cat(sprintf("  ECV                     : %.3f\n", rel$ecv))
      cat(sprintf("  PUC                     : %.3f\n", rel$puc))
      cat(sprintf("  H (factor general)      : %.3f\n", rel$h_general))
      cat("\n  Por factor especifico:\n")
    }
    cat(sprintf("  %-20s  %8s  %8s\n", "Factor", "Omega_hs", "H"))
    for (i in seq_along(rel$omega_hs)) {
      cat(sprintf("  %-20s  %8.3f  %8.3f\n",
                  names(rel$omega_hs)[i], rel$omega_hs[i], rel$h_specific[i]))
    }

    if (is_english) {
      cat("\n  Interpretation:\n")
    } else {
      cat("\n  Interpretacion:\n")
    }
    if (rel$omega_h >= 0.80) {
      if (is_english) {
        cat("    Omega_h >= .80: STRONG general factor.\n")
        cat("    -> The total score mainly reflects the general factor.\n")
      } else {
        cat("    Omega_h >= .80: Factor general FUERTE.\n")
        cat("    -> La puntuacion total refleja predominantemente el factor general.\n")
      }
    } else if (rel$omega_h >= 0.50) {
      if (is_english) {
        cat("    Omega_h between .50-.80: MODERATE general factor.\n")
        cat("    -> Evaluate whether the subscales add incremental information.\n")
      } else {
        cat("    Omega_h entre .50-.80: Factor general MODERADO.\n")
        cat("    -> Evaluar si las subescalas aportan informacion incremental.\n")
      }
    } else if (is_english) {
      cat("    Omega_h < .50: WEAK general factor.\n")
      cat("    -> Subscale scores may be more informative.\n")
    } else {
      cat("    Omega_h < .50: Factor general DEBIL.\n")
      cat("    -> Las puntuaciones por subescala pueden ser mas informativas.\n")
    }

    if (rel$puc > 0.80 && rel$ecv > 0.60 && rel$omega_h >= 0.70) {
      if (is_english) {
        cat("    PUC > .80, ECV > .60, and Omega_h >= .70: Evidence favorable to\n")
        cat("    unidimensionality (Reise et al., 2013).\n")
      } else {
        cat("    PUC > .80, ECV > .60 y Omega_h >= .70: Evidencia favorable a la\n")
        cat("    unidimensionalidad (Reise et al., 2013).\n")
      }
    } else if (rel$ecv > 0.70 && rel$puc > 0.70) {
      if (is_english) {
        cat("    ECV > .70 and PUC > .70: One could conclude in favor of\n")
        cat("    unidimensionality (Rodriguez et al., 2016).\n")
      } else {
        cat("    ECV > .70 y PUC > .70: Se podria concluir a favor de la\n")
        cat("    unidimensionalidad (Rodriguez et al., 2016).\n")
      }
    } else if (rel$ecv < 0.60) {
      if (is_english) {
        cat("    ECV < .60: Substantial multidimensionality is present.\n")
      } else {
        cat("    ECV < .60: Multidimensionalidad sustancial presente.\n")
      }
    } else if (is_english) {
      cat(sprintf("    ECV = %.3f, PUC = %.3f: The thresholds for concluding\n", rel$ecv, rel$puc))
      cat("    in favor of unidimensionality are not met. Interpret subscales\n")
      cat("    with caution.\n")
    } else {
      cat(sprintf("    ECV = %.3f, PUC = %.3f: No se cumplen los umbrales para\n", rel$ecv, rel$puc))
      cat("    concluir a favor de la unidimensionalidad. Interpretar subescalas\n")
      cat("    con cautela.\n")
    }
    for (i in seq_along(rel$h_specific)) {
      if (rel$h_specific[i] < 0.50) {
        if (is_english) {
          cat(sprintf("    H(%s) < .50: Specific factor shows low replicability.\n",
                      names(rel$h_specific)[i]))
        } else {
          cat(sprintf("    H(%s) < .50: Factor especifico poco replicable.\n",
                      names(rel$h_specific)[i]))
        }
      }
    }
    cat("\n")
  }

  if (length(x$respecifications) > 0) {
    cat(txt$respec_title, "\n")
    cat(sep2, "\n", sep = "")
    cat(txt$respec_note, "\n\n")
    for (i in seq_along(x$respecifications)) {
      r <- x$respecifications[[i]]
      cat(sprintf("  %d. [%s] %s\n", i, r$tipo, r$sintaxis))
      cat(sprintf("     %s = %.2f, EPC = %.3f\n",
                  if (is_english) "MI" else "IM", r$mi, r$epc))
      cat(sprintf("     %s\n\n", r$nota))
    }
  } else {
    cat(txt$respec_none_title, "\n")
    cat(sep2, "\n", sep = "")
    cat(txt$respec_none, "\n\n")
  }

  cat(txt$references_title, "\n")
  cat(sep2, "\n", sep = "")
  cat(txt$references_note, "\n\n")
  refs <- .build_references(x)
  for (i in seq_along(refs)) {
    cat(sprintf("  [%d] %s\n\n", i, refs[i]))
  }

  cat(sep, "\n", sep = "")
  cat(txt$footer, "\n")
  cat(sep, "\n\n", sep = "")

  invisible(x)
}

.cfa_respecifications_to_df <- function(x) {
  respec <- if (!is.null(x$respecification_suggestions)) {
    x$respecification_suggestions
  } else {
    x$respecifications
  }

  if (is.null(respec) || length(respec) == 0) {
    return(data.frame())
  }

  do.call(
    rbind,
    lapply(respec, function(item) {
      data.frame(
        Type = if (!is.null(item$type)) item$type else item$tipo,
        Syntax = if (!is.null(item$syntax)) item$syntax else item$sintaxis,
        MI = item$mi,
        EPC = item$epc,
        SameFactor = if (!is.null(item$same_factor)) item$same_factor else item$mismo_factor,
        Note = if (!is.null(item$note)) item$note else item$nota,
        stringsAsFactors = FALSE
      )
    })
  )
}

.cfa_export_excel <- function(result, file_name) {
  wb <- createWorkbook()
  fit_df <- if (!is.null(result$fit_summary)) result$fit_summary else result$fit_indices
  if (!"Index" %in% names(fit_df) && "Indice" %in% names(fit_df)) fit_df$Index <- fit_df$Indice
  if (!"Value" %in% names(fit_df) && "Valor" %in% names(fit_df)) fit_df$Value <- fit_df$Valor
  if (!"Evaluation" %in% names(fit_df) && "Evaluacion" %in% names(fit_df)) fit_df$Evaluation <- fit_df$Evaluacion

  loadings_df <- if (!is.null(result$factor_loadings)) result$factor_loadings else result$loadings
  if (!"Loading" %in% names(loadings_df) && "Lambda" %in% names(loadings_df)) loadings_df$Loading <- loadings_df$Lambda
  if (!"p_value" %in% names(loadings_df) && "p" %in% names(loadings_df)) loadings_df$p_value <- loadings_df$p
  if (!"CI_lower" %in% names(loadings_df) && "IC_inf" %in% names(loadings_df)) loadings_df$CI_lower <- loadings_df$IC_inf
  if (!"CI_upper" %in% names(loadings_df) && "IC_sup" %in% names(loadings_df)) loadings_df$CI_upper <- loadings_df$IC_sup
  if (!"Alert" %in% names(loadings_df) && "Alerta" %in% names(loadings_df)) loadings_df$Alert <- loadings_df$Alerta

  higher_order_df <- if (!is.null(result$higher_order_loadings)) result$higher_order_loadings else result$loadings_ho
  if (!is.null(higher_order_df) && nrow(higher_order_df) > 0) {
    if (!"HigherOrderFactor" %in% names(higher_order_df) && "Factor_SO" %in% names(higher_order_df)) higher_order_df$HigherOrderFactor <- higher_order_df$Factor_SO
    if (!"LowerOrderFactor" %in% names(higher_order_df) && "Factor_PO" %in% names(higher_order_df)) higher_order_df$LowerOrderFactor <- higher_order_df$Factor_PO
    if (!"Loading" %in% names(higher_order_df) && "Lambda" %in% names(higher_order_df)) higher_order_df$Loading <- higher_order_df$Lambda
    if (!"p_value" %in% names(higher_order_df) && "p" %in% names(higher_order_df)) higher_order_df$p_value <- higher_order_df$p
  }

  corr_df <- if (!is.null(result$factor_correlation_table)) result$factor_correlation_table else result$factor_correlations
  if (!is.null(corr_df) && nrow(corr_df) > 0) {
    if (!"Correlation" %in% names(corr_df) && "r" %in% names(corr_df)) corr_df$Correlation <- corr_df$r
    if (!"p_value" %in% names(corr_df) && "p" %in% names(corr_df)) corr_df$p_value <- corr_df$p
    if (!"Alert" %in% names(corr_df) && "Alerta" %in% names(corr_df)) corr_df$Alert <- corr_df$Alerta
  }

  warnings_vec <- if (!is.null(result$warnings_list)) result$warnings_list else result$warnings

  summary_df <- data.frame(
    Element = c(
      "Model type", "Factors", "Observed items", "Sample size",
      "Complete cases", "Missing rows (%)", "Estimator",
      "Estimator rationale", "Converged", "Language"
    ),
    Value = c(
      result$model_type,
      paste(result$factors, collapse = ", "),
      result$n_items,
      result$n,
      result$n_complete,
      result$pct_missing,
      result$estimator,
      result$estimator_reason,
      result$converged,
      result$language
    ),
    stringsAsFactors = FALSE
  )

  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", summary_df)

  addWorksheet(wb, "Fit Indices")
  writeData(wb, "Fit Indices", fit_df[, c("Metric", "Index", "Value", "Evaluation"), drop = FALSE])

  addWorksheet(wb, "Loadings")
  writeData(
    wb,
    "Loadings",
    loadings_df[, c("Factor", "Item", "Loading", "SE", "z", "p_value", "CI_lower", "CI_upper", "Alert"), drop = FALSE]
  )

  if (!is.null(higher_order_df) && nrow(higher_order_df) > 0) {
    addWorksheet(wb, "Higher Order")
    writeData(
      wb,
      "Higher Order",
      higher_order_df[, c("HigherOrderFactor", "LowerOrderFactor", "Loading", "SE", "z", "p_value"), drop = FALSE]
    )
  }

  if (!is.null(corr_df) && nrow(corr_df) > 0) {
    addWorksheet(wb, "Factor Corr")
    writeData(
      wb,
      "Factor Corr",
      corr_df[, c("Factor1", "Factor2", "Correlation", "SE", "z", "p_value", "Alert"), drop = FALSE]
    )
  }

  if (!is.null(result$modification_indices) && nrow(result$modification_indices) > 0) {
    addWorksheet(wb, "Mod Indices")
    writeData(wb, "Mod Indices", result$modification_indices)
  }

  respec_df <- .cfa_respecifications_to_df(result)
  if (nrow(respec_df) > 0) {
    addWorksheet(wb, "Respec")
    writeData(wb, "Respec", respec_df)
  }

  if (length(warnings_vec) > 0) {
    addWorksheet(wb, "Warnings")
    writeData(wb, "Warnings", data.frame(Warning = warnings_vec, stringsAsFactors = FALSE))
  }

  if (!grepl("\\.xlsx$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".xlsx")
  }

  saveWorkbook(wb, file_name, overwrite = TRUE)
  invisible(normalizePath(file_name))
}

.cfa_export_word <- function(result, file_name, digits = 3) {
  doc <- read_docx()
  report_lines <- capture.output(print(result, digits = digits))

  for (line in report_lines) {
    doc <- body_add_par(doc, if (nzchar(line)) line else " ", style = "Normal")
  }

  if (!grepl("\\.docx$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".docx")
  }

  print(doc, target = file_name)
  invisible(normalizePath(file_name))
}

#' Export CFA results to Excel or Word
#'
#' @param result Object of class `"cfa_auto"` returned by `cfa_auto()`.
#' @param format `"excel"` or `"word"` (you can provide both).
#' @param file_name Optional base file name without extension. If `NULL`, a
#'   timestamped name is generated automatically.
#' @param digits Number of decimal places used when exporting the Word report.
#'
#' @return The input `result`, invisibly.
#'
#' @examples
#' # export_cfa(result, format = "excel")
#' # export_cfa(result, format = "word")
#'
#' @export
export_cfa <- function(result, format = c("excel", "word"), file_name = NULL, digits = 3) {
  if (!inherits(result, "cfa_auto")) {
    stop("'result' must be an object of class 'cfa_auto'.", call. = FALSE)
  }

  format <- match.arg(format, c("excel", "word"), several.ok = TRUE)

  if (is.null(file_name)) {
    file_name <- paste0("CFA_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  for (fmt in format) {
    if (fmt == "excel") {
      .cfa_export_excel(result, file_name)
    } else {
      .cfa_export_word(result, file_name, digits = digits)
    }
  }

  invisible(result)
}

# Internal legacy alias for `export_cfa()`
#'
#' @param resultado Deprecated. Use `result`.
#' @param formato Deprecated. Use `format`.
#' @param archivo Deprecated. Use `file_name`.
#' @param digits Number of decimal places used when exporting the Word report.
#'
#' @return The input `resultado`, invisibly.
#' @noRd
exportar_cfa <- function(resultado, formato = c("excel", "word"), archivo = NULL, digits = 3) {
  warning("`exportar_cfa()` is deprecated; use `export_cfa()` instead.", call. = FALSE)
  export_cfa(result = resultado, format = formato, file_name = archivo, digits = digits)
}
