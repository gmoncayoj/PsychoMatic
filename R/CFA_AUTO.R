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
#' if (identical(Sys.getenv("PSYCHOMATIC_RUN_SLOW_EXAMPLES"), "true") &&
#'     requireNamespace("semTools", quietly = TRUE) &&
#'     requireNamespace("MVN", quietly = TRUE)) {
#'   model <- "visual =~ x1 + x2 + x3"
#'   items <- lavaan::HolzingerSwineford1939[1:80, c("x1", "x2", "x3")]
#'   result <- cfa_auto(items, model = model, estimator = "ML", language = "eng")
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

  results <- .cfa_add_tidy_fields(results, call = match.call())
  class(results) <- "cfa_auto"
  return(results)
}
