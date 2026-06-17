# ============================================================
#        ANALISIS FACTORIAL EXPLORATORIO OPTIMO (AFE)
#   Seleccion automatica de estimador, rotacion y refinamiento
#   iterativo con eliminacion de items problematicos
# ============================================================
# Paquetes requeridos: psych, GPArotation
# ============================================================

# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

# Test de Mardia para normalidad multivariada
.efa_mardia <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Submuestrear si n es muy grande para evitar problemas de memoria
  if (n > 5000) {
    set.seed(7492)
    X <- X[sample(n, 5000), ]
    n <- 5000L
  }

  X_c <- scale(X, scale = FALSE)
  S <- crossprod(X_c) / n

  S_inv <- tryCatch(solve(S), error = function(e) {
    ev <- eigen(S, symmetric = TRUE)
    vals <- ev$values
    vals[vals < 1e-10] <- Inf
    ev$vectors %*% diag(1 / vals) %*% t(ev$vectors)
  })

  G <- X_c %*% S_inv
  D <- tcrossprod(G, X_c)

  # Asimetria multivariada de Mardia
  b1p <- sum(D^3) / n^2
  chi_skew <- n * b1p / 6
  df_skew <- p * (p + 1) * (p + 2) / 6
  p_skew <- pchisq(chi_skew, df_skew, lower.tail = FALSE)

  # Curtosis multivariada de Mardia
  b2p <- sum(diag(D)^2) / n
  z_kurt <- (b2p - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  p_kurt <- 2 * pnorm(abs(z_kurt), lower.tail = FALSE)

  list(
    skewness = b1p, skewness_chi = chi_skew,
    skewness_df = df_skew, skewness_p = p_skew,
    kurtosis = b2p, kurtosis_z = z_kurt, kurtosis_p = p_kurt,
    normal = (p_skew > 0.05) & (p_kurt > 0.05)
  )
}

# Normalizar idioma del reporte
.efa_normalize_language <- function(language = "esp") {
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

.efa_warn_deprecated_arg <- function(old_name, new_name) {
  warning(
    sprintf("Argument `%s` is deprecated; use `%s` instead.", old_name, new_name),
    call. = FALSE
  )
}

.efa_normalize_rotation <- function(rotation = "oblique") {
  if (is.null(rotation) || length(rotation) == 0) return("oblique")

  rotation <- tolower(trimws(as.character(rotation[1])))

  if (rotation %in% c("oblique", "oblimin", "oblicua")) {
    return("oblique")
  }

  if (rotation %in% c("orthogonal", "varimax", "ortogonal")) {
    return("orthogonal")
  }

  stop("Argument 'rotation' must be one of: 'oblique' or 'orthogonal'.")
}

# Etiquetas y textos del reporte por consola
.efa_labels <- function(language = "esp") {
  language <- .efa_normalize_language(language)

  if (language == "eng") {
    return(list(
      required_psych = "Package 'psych' is required. Install with: install.packages('psych')",
      required_gpa = "Package 'GPArotation' is required. Install with: install.packages('GPArotation')",
      invalid_data = "'data' must be a data frame or a matrix.",
      need_numeric = "At least 3 numeric variables are required.",
      excluded_non_numeric = "NOTE: Excluded",
      non_numeric_cols = "non-numeric column(s).",
      removed_rows = "NOTE: Removed",
      rows_with_na = "row(s) with NA.",
      final_sample = "Final sample: n =",
      analysis_title = "          OPTIMAL EXPLORATORY FACTOR ANALYSIS",
      data_title = "--- DATA CHARACTERISTICS ---",
      sample_size = "  Sample size:            n =",
      n_items = "  Number of items:        p =",
      median_categories = "  Categories (median):   ",
      data_type = "  Data type:             ",
      data_type_ordinal = "Ordinal (<=5 categories)",
      data_type_continuous = "Continuous/quasi-continuous (>5 categories)",
      categories_by_item = "  Categories by item:",
      normality_title = "--- MULTIVARIATE NORMALITY (Mardia test) ---",
      normality_conclusion = "  Conclusion:",
      normality_ok = "Multivariate normality IS MET (p > .05)\n",
      normality_bad = "Multivariate normality is NOT met (p <= .05)\n",
      estimator_title = "--- ESTIMATOR SELECTION ---",
      estimator_label = "  Estimator:      ",
      justification_label = "  Rationale:      ",
      correlation_label = "  Correlation:    ",
      calculating_polychoric = "\n  Computing polychoric correlations...\n",
      polychoric_warning = "Error in polychoric correlations. Using Pearson as fallback: ",
      non_pd_matrix = "  NOTE: Matrix is not positive definite. Applying smoothing.\n",
      adequacy_title = "--- SAMPLING ADEQUACY (original items) ---",
      kmo_warning_items = "\n  WARNING - Items with individual KMO < 0.50:\n",
      kmo_warning_global = "KMO < 0.50: The data may not be adequate for factor analysis.",
      parallel_title = "--- PARALLEL ANALYSIS ---",
      optimal_factors = "  Optimal number of factors:",
      no_rotation = "  With a single factor, rotation is not applied.\n",
      rotation_label = "  Rotation: %s (%s)\n",
      oblique = "oblique",
      orthogonal = "orthogonal",
      refinement_title = "          ITERATIVE REFINEMENT PROCESS",
      factors_adjusted = "\n  >> Number of factors adjusted: %d -> %d\n",
      insufficient_items = "  >> Too few items. Reducing to",
      factor_suffix = "factor(s).\n",
      estimator_error_note = "  Note: Error with",
      switching_minres = "- switching to 'minres'.\n",
      iteration_title = "\n--- Iteration %d (%d items, %d factor(s)) ---\n",
      no_problem_items = "  No problematic items.\n",
      problematic_items = "  Problematic items detected:\n",
      removed_item = "  >> REMOVED: %s | %s\n",
      fewer_than_three = "Fewer than 3 items remain. Process stopped.",
      no_model = "No factor model could be fitted.",
      final_results_title = "          FINAL RESULTS",
      fit_indices_title = "--- FIT INDICES ---",
      chi_square = "  Chi-square:         %.3f\n",
      degrees_freedom = "  Degrees of freedom: %d\n",
      p_value = "  p-value:            %s\n",
      tli = "  TLI:                %.3f\n",
      rmsea = "  RMSEA:              %.3f [90%% CI: %.3f - %.3f]\n",
      bic = "  BIC:                %.3f\n",
      final_kmo = "\n  KMO (final solution): %.3f - %s\n",
      final_bartlett = "  Bartlett (final sol.): chi2 = %.3f, df = %.0f, p = %s\n",
      variance_title = "--- EXPLAINED VARIANCE ---",
      factor_variance = "  Factor %d: %.1f%%\n",
      total_variance = "  Cumulative total: %.1f%%\n",
      loadings_title = "--- FACTOR LOADING MATRIX ---",
      loadings_note_left = "  (Loadings with |value| <",
      loadings_note_right = "omitted for readability)\n",
      distribution_title = "--- ITEM DISTRIBUTION BY DIMENSION ---",
      distribution_factor = "\n  %s (%d items):\n",
      distribution_item = "    %-25s loading = %6.3f    h2 = %.3f\n",
      factor_corr_title = "--- FACTOR CORRELATIONS ---",
      removed_items_title = "--- REMOVED ITEMS ---",
      removed_items_line = "  %2d. %-20s %s (iter. %d)\n",
      removed_items_total = "\n  Total removed: %d of %d original items\n",
      no_items_removed = "\n  No items were removed.\n",
      reliability_title = "--- RELIABILITY BY DIMENSION ---",
      alpha_improve = "    Items whose removal would improve alpha:\n",
      alpha_if_removed = "      %-20s alpha if removed = %.3f (+%.3f)\n",
      reliability_error = "  %s: Error while computing reliability.\n",
      spearman_brown = "\n  %s (2 items): Spearman-Brown = %.3f (r = %.3f)\n",
      reliability_na = "\n  %s: %d item - reliability not computable.\n",
      optimal_solution_title = "\n  >>> OPTIMAL SOLUTION <<<\n\n",
      optimal_solution_1 = "  The factor solution does not require additional\n",
      optimal_solution_2 = "  respecifications. All quality criteria are met.\n\n",
      respec_title = "  RESPECIFICATION SUGGESTIONS\n",
      refs_title = "\n--- METHODOLOGICAL REFERENCES ---\n\n",
      refs_intro = "  The decisions in this analysis are based on:\n\n",
      refs_estimator = "  ULS/MINRES estimator:\n",
      refs_polychoric = "  Polychoric correlations for ordinal data:\n",
      refs_parallel = "  Parallel analysis for factor retention:\n",
      refs_normality = "  Multivariate normality (Mardia test):\n",
      refs_kmo = "  KMO and Bartlett's test:\n",
      refs_elimination = "  Item removal criteria and EFA best practices:\n",
      refs_reliability = "  Reliability (alpha and omega):\n",
      refs_software = "  Software:\n"
    ))
  }

  list(
    required_psych = "Paquete 'psych' requerido. Instale con: install.packages('psych')",
    required_gpa = "Paquete 'GPArotation' requerido. Instale con: install.packages('GPArotation')",
    invalid_data = "'data' debe ser un data frame o una matriz.",
    need_numeric = "Se necesitan al menos 3 variables numericas.",
    excluded_non_numeric = "NOTA: Se excluyeron",
    non_numeric_cols = "columna(s) no numerica(s).",
    removed_rows = "NOTA: Se eliminaron",
    rows_with_na = "fila(s) con NA.",
    final_sample = "Muestra final: n =",
    analysis_title = "          ANALISIS FACTORIAL EXPLORATORIO OPTIMO",
    data_title = "--- CARACTERISTICAS DE LOS DATOS ---",
    sample_size = "  Tamano de muestra:        n =",
    n_items = "  Numero de items:          p =",
    median_categories = "  Categorias (mediana):    ",
    data_type = "  Tipo de datos:           ",
    data_type_ordinal = "Ordinales (<=5 categorias)",
    data_type_continuous = "Continuos/cuasi-continuos (>5 categorias)",
    categories_by_item = "  Categorias por item:",
    normality_title = "--- NORMALIDAD MULTIVARIADA (Test de Mardia) ---",
    normality_conclusion = "  Conclusion:",
    normality_ok = "Se CUMPLE normalidad multivariada (p > .05)\n",
    normality_bad = "NO se cumple normalidad multivariada (p <= .05)\n",
    estimator_title = "--- SELECCION DEL ESTIMADOR ---",
    estimator_label = "  Estimador:     ",
    justification_label = "  Justificacion: ",
    correlation_label = "  Correlacion:   ",
    calculating_polychoric = "\n  Calculando correlaciones policoricas...\n",
    polychoric_warning = "Error en policoricas. Usando Pearson como alternativa: ",
    non_pd_matrix = "  NOTA: Matriz no definida positiva. Aplicando suavizado.\n",
    adequacy_title = "--- ADECUACION MUESTRAL (items originales) ---",
    kmo_warning_items = "\n  ADVERTENCIA - Items con KMO individual < 0.50:\n",
    kmo_warning_global = "KMO < 0.50: Los datos pueden no ser adecuados para analisis factorial.",
    parallel_title = "--- ANALISIS PARALELO ---",
    optimal_factors = "  Numero optimo de factores:",
    no_rotation = "  Con 1 solo factor, no se aplica rotacion.\n",
    rotation_label = "  Rotacion: %s (%s)\n",
    oblique = "oblicua",
    orthogonal = "ortogonal",
    refinement_title = "          PROCESO ITERATIVO DE REFINAMIENTO",
    factors_adjusted = "\n  >> Numero de factores ajustado: %d -> %d\n",
    insufficient_items = "  >> Items insuficientes. Reduciendo a",
    factor_suffix = "factor(es).\n",
    estimator_error_note = "  Nota: Error con",
    switching_minres = "- cambiando a 'minres'.\n",
    iteration_title = "\n--- Iteracion %d (%d items, %d factor(es)) ---\n",
    no_problem_items = "  Sin items problematicos.\n",
    problematic_items = "  Items problematicos detectados:\n",
    removed_item = "  >> ELIMINADO: %s | %s\n",
    fewer_than_three = "Quedan menos de 3 items. Proceso detenido.",
    no_model = "No se pudo ajustar ningun modelo factorial.",
    final_results_title = "          RESULTADOS FINALES",
    fit_indices_title = "--- INDICES DE AJUSTE ---",
    chi_square = "  Chi-cuadrado:       %.3f\n",
    degrees_freedom = "  Grados de libertad: %d\n",
    p_value = "  p-valor:            %s\n",
    tli = "  TLI:                %.3f\n",
    rmsea = "  RMSEA:              %.3f [IC 90%%: %.3f - %.3f]\n",
    bic = "  BIC:                %.3f\n",
    final_kmo = "\n  KMO (solucion final):  %.3f - %s\n",
    final_bartlett = "  Bartlett (sol. final): chi2 = %.3f, gl = %.0f, p = %s\n",
    variance_title = "--- VARIANZA EXPLICADA ---",
    factor_variance = "  Factor %d: %.1f%%\n",
    total_variance = "  Total acumulada: %.1f%%\n",
    loadings_title = "--- MATRIZ DE CARGAS FACTORIALES ---",
    loadings_note_left = "  (Cargas con |valor| <",
    loadings_note_right = "omitidas para legibilidad)\n",
    distribution_title = "--- DISTRIBUCION DE ITEMS POR DIMENSION ---",
    distribution_factor = "\n  %s (%d items):\n",
    distribution_item = "    %-25s carga = %6.3f    h2 = %.3f\n",
    factor_corr_title = "--- CORRELACIONES ENTRE FACTORES ---",
    removed_items_title = "--- ITEMS ELIMINADOS ---",
    removed_items_line = "  %2d. %-20s %s (iter. %d)\n",
    removed_items_total = "\n  Total eliminados: %d de %d items originales\n",
    no_items_removed = "\n  No se eliminaron items.\n",
    reliability_title = "--- CONFIABILIDAD POR DIMENSION ---",
    alpha_improve = "    Items cuya eliminacion mejoraria alpha:\n",
    alpha_if_removed = "      %-20s alpha si elimina = %.3f (+%.3f)\n",
    reliability_error = "  %s: Error al calcular confiabilidad.\n",
    spearman_brown = "\n  %s (2 items): Spearman-Brown = %.3f (r = %.3f)\n",
    reliability_na = "\n  %s: %d item - confiabilidad no calculable.\n",
    optimal_solution_title = "\n  >>> SOLUCION OPTIMA <<<\n\n",
    optimal_solution_1 = "  La solucion factorial no requiere reespecificaciones\n",
    optimal_solution_2 = "  adicionales. Todos los criterios de calidad se cumplen.\n\n",
    respec_title = "  SUGERENCIAS DE REESPECIFICACION\n",
    refs_title = "\n--- REFERENCIAS METODOLOGICAS ---\n\n",
    refs_intro = "  Las decisiones de este analisis se fundamentan en:\n\n",
    refs_estimator = "  Estimador ULS/MINRES:\n",
    refs_polychoric = "  Correlaciones policoricas para datos ordinales:\n",
    refs_parallel = "  Analisis paralelo para retencion de factores:\n",
    refs_normality = "  Normalidad multivariada (Test de Mardia):\n",
    refs_kmo = "  KMO y Test de Bartlett:\n",
    refs_elimination = "  Criterios de eliminacion de items y buenas practicas en AFE:\n",
    refs_reliability = "  Confiabilidad (alfa y omega):\n",
    refs_software = "  Software:\n"
  )
}

# Interpretacion del KMO
.efa_interpretar_kmo <- function(kmo, language = "esp") {
  language <- .efa_normalize_language(language)

  if (language == "eng") {
    if (kmo >= 0.90) return("Excellent")
    if (kmo >= 0.80) return("Good")
    if (kmo >= 0.70) return("Acceptable")
    if (kmo >= 0.60) return("Mediocre")
    if (kmo >= 0.50) return("Poor")
    return("Unacceptable")
  }

  if (kmo >= 0.90) return("Excelente")
  if (kmo >= 0.80) return("Bueno")
  if (kmo >= 0.70) return("Aceptable")
  if (kmo >= 0.60) return("Mediocre")
  if (kmo >= 0.50) return("Malo")
  "Inaceptable"
}

# Seleccion automatica del estimador
# ULS/MINRES es el estimador recomendado por defecto en AFE
# (Lloret-Segura et al., 2014; Ferrando & Lorenzo-Seva, 2017).
# No requiere supuestos distribucionales y es robusto ante
# violaciones de normalidad. Se prioriza siempre, combinandolo
# con la matriz de correlacion apropiada segun el tipo de datos.
.efa_seleccionar_estimador <- function(es_ordinal, normalidad_ok, n, p,
                                      language = "esp") {
  language <- .efa_normalize_language(language)

  if (es_ordinal) {
    list(
      fm = "minres",
      nombre = if (language == "eng") {
        "ULS/MINRES (Minimum Residuals / Unweighted Least Squares)"
      } else {
        "ULS/MINRES (Minimos Residuales / Unweighted Least Squares)"
      },
      justificacion = if (language == "eng") {
        paste0(
          "ULS is the recommended estimator for EFA (Lloret-Segura et al., 2014). ",
          "Ordinal data (<=5 categories): it is used with polychoric correlations. ",
          "It does not require normality assumptions. n=", n, ", p=", p, "."
        )
      } else {
        paste0(
          "ULS es el estimador recomendado para AFE (Lloret-Segura et al., 2014). ",
          "Datos ordinales (<=5 categorias): se usa con correlaciones policoricas. ",
          "No requiere supuestos de normalidad. n=", n, ", p=", p, "."
        )
      },
      tipo_cor = if (language == "eng") "Polychoric" else "Policorica",
      usar_policorica = TRUE
    )
  } else {
    list(
      fm = "minres",
      nombre = if (language == "eng") {
        "ULS/MINRES (Minimum Residuals / Unweighted Least Squares)"
      } else {
        "ULS/MINRES (Minimos Residuales / Unweighted Least Squares)"
      },
      justificacion = if (language == "eng") {
        paste0(
          "ULS is the recommended estimator for EFA (Lloret-Segura et al., 2014). ",
          "Continuous data (>5 categories): it is used with Pearson correlations. ",
          "It is robust to violations of normality. n=", n, ", p=", p,
          ifelse(normalidad_ok,
                 ". Multivariate normality: YES.",
                 ". Multivariate normality: NO (ULS does not require it).")
        )
      } else {
        paste0(
          "ULS es el estimador recomendado para AFE (Lloret-Segura et al., 2014). ",
          "Datos continuos (>5 categorias): se usa con correlaciones de Pearson. ",
          "Robusto ante violaciones de normalidad. n=", n, ", p=", p,
          ifelse(normalidad_ok,
                 ". Normalidad multivariada: SI.",
                 ". Normalidad multivariada: NO (ULS no la requiere).")
        )
      },
      tipo_cor = "Pearson",
      usar_policorica = FALSE
    )
  }
}

# Evaluacion de items problematicos
.efa_evaluar_items <- function(cargas, comunalidades, carga_min,
                               comunalidad_min, dif_cargas_cruzadas,
                               n_factores, language = "esp") {
  language <- .efa_normalize_language(language)
  problemas <- data.frame(
    Item = character(), Razon = character(),
    Prioridad = numeric(), stringsAsFactors = FALSE
  )

  for (item in rownames(cargas)) {
    cargas_abs <- abs(as.numeric(cargas[item, ]))

    # Prioridad 3: Comunalidad baja
    if (comunalidades[item] < comunalidad_min) {
      problemas <- rbind(problemas, data.frame(
        Item = item,
        Razon = if (language == "eng") {
          sprintf("Low communality (h2 = %.3f < %.2f)",
                  comunalidades[item], comunalidad_min)
        } else {
          sprintf("Comunalidad baja (h2 = %.3f < %.2f)",
                  comunalidades[item], comunalidad_min)
        },
        Prioridad = 3 + (comunalidad_min - comunalidades[item]),
        stringsAsFactors = FALSE
      ))
      next
    }

    # Prioridad 2: Carga maxima insuficiente
    if (max(cargas_abs) < carga_min) {
      problemas <- rbind(problemas, data.frame(
        Item = item,
        Razon = if (language == "eng") {
          sprintf("Insufficient maximum loading (max = %.3f < %.2f)",
                  max(cargas_abs), carga_min)
        } else {
          sprintf("Carga maxima insuficiente (max = %.3f < %.2f)",
                  max(cargas_abs), carga_min)
        },
        Prioridad = 2 + (carga_min - max(cargas_abs)),
        stringsAsFactors = FALSE
      ))
      next
    }

    # Prioridad 1: Cargas cruzadas
    if (n_factores > 1) {
      sorted_idx <- order(cargas_abs, decreasing = TRUE)
      sorted_vals <- cargas_abs[sorted_idx]
      if (length(sorted_vals) >= 2 && sorted_vals[2] >= carga_min) {
        diff_top2 <- sorted_vals[1] - sorted_vals[2]
        if (diff_top2 < dif_cargas_cruzadas) {
          f1 <- colnames(cargas)[sorted_idx[1]]
          f2 <- colnames(cargas)[sorted_idx[2]]
          problemas <- rbind(problemas, data.frame(
            Item = item,
            Razon = if (language == "eng") {
              sprintf(
                "Cross-loading (%s=%.3f, %s=%.3f, diff=%.3f < %.2f)",
                f1, sorted_vals[1], f2, sorted_vals[2],
                diff_top2, dif_cargas_cruzadas
              )
            } else {
              sprintf(
                "Carga cruzada (%s=%.3f, %s=%.3f, dif=%.3f < %.2f)",
                f1, sorted_vals[1], f2, sorted_vals[2],
                diff_top2, dif_cargas_cruzadas
              )
            },
            Prioridad = 1 + (dif_cargas_cruzadas - diff_top2),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  if (nrow(problemas) > 0) {
    problemas <- problemas[order(-problemas$Prioridad), ]
  }
  problemas
}

# Generacion de sugerencias post-analisis
.efa_sugerencias <- function(efa_obj, distribucion, fiabilidad, kmo_final,
                             bartlett_final, n_factores, n, items_eliminados,
                             carga_min, var_explicada, p_original,
                             language = "esp") {
  language <- .efa_normalize_language(language)
  sug <- character()

  # KMO bajo
  if (kmo_final$MSA < 0.60) {
    sug <- c(sug, sprintf(
      if (language == "eng") {
        "Final KMO (%.3f) is low. Review whether the items are suitable for factor analysis."
      } else {
        "KMO final (%.3f) es bajo. Revise si los items son adecuados para AF."
      },
      kmo_final$MSA
    ))
  }

  # Bartlett no significativo
  if (bartlett_final$p.value > 0.05) {
    sug <- c(sug, paste0(
      if (language == "eng") "Bartlett not significant (p = " else "Bartlett no significativo (p = ",
      sprintf("%.4f", bartlett_final$p.value),
      if (language == "eng") "). The matrix may be an identity matrix." else "). La matriz podria ser una identidad."
    ))
  }

  # Factores con pocos items
  for (f in sort(unique(distribucion$Dimension))) {
    n_items_f <- sum(distribucion$Dimension == f)
    if (n_items_f < 3) {
      sug <- c(sug, sprintf(
        if (language == "eng") {
          "%s has only %d item(s). A minimum of 3 items per dimension is recommended."
        } else {
          "%s tiene solo %d item(s). Se recomienda minimo 3 items por dimension."
        },
        f, n_items_f
      ))
    }
  }

  # Fiabilidad baja
  for (f in names(fiabilidad)) {
    a <- fiabilidad[[f]]$alpha
    if (!is.na(a) && a < 0.60) {
      sug <- c(sug, sprintf(
        if (language == "eng") {
          "%s: alpha = %.3f (< .60). The dimension shows low internal consistency."
        } else {
          "%s: alpha = %.3f (< .60). La dimension tiene baja consistencia interna."
        },
        f, a
      ))
    } else if (!is.na(a) && a < 0.70) {
      sug <- c(sug, sprintf(
        if (language == "eng") {
          "%s: alpha = %.3f (< .70). Consider reviewing the items in this dimension."
        } else {
          "%s: alpha = %.3f (< .70). Considere revisar los items de esta dimension."
        },
        f, a
      ))
    }
  }

  # Varianza total explicada
  if (is.matrix(var_explicada)) {
    total_var <- NULL
    if ("Cumulative Var" %in% rownames(var_explicada)) {
      total_var <- var_explicada["Cumulative Var", ncol(var_explicada)]
    } else if ("Cumulative Proportion" %in% rownames(var_explicada)) {
      total_var <- var_explicada["Cumulative Proportion", ncol(var_explicada)]
    }
    if (!is.null(total_var) && total_var < 0.50) {
      sug <- c(sug, sprintf(
        if (language == "eng") {
          "Total explained variance (%.1f%%) < 50%%. Consider reviewing the structure."
        } else {
          "Varianza total explicada (%.1f%%) < 50%%. Considere revisar la estructura."
        },
        total_var * 100
      ))
    }
  }

  # Correlaciones muy altas entre factores (rotacion oblicua)
  if (!is.null(efa_obj$Phi)) {
    phi <- efa_obj$Phi
    diag(phi) <- 0
    high_cor <- which(abs(phi) > 0.85, arr.ind = TRUE)
    if (nrow(high_cor) > 0) {
      pares_vistos <- character()
      for (k in 1:nrow(high_cor)) {
        par <- paste(sort(c(high_cor[k, 1], high_cor[k, 2])), collapse = "-")
        if (!(par %in% pares_vistos)) {
          pares_vistos <- c(pares_vistos, par)
          sug <- c(sug, sprintf(
            if (language == "eng") {
              "Correlation F%d-F%d = %.3f (> .85). They could be merged into a single factor."
            } else {
              "Correlacion F%d-F%d = %.3f (> .85). Podrian fusionarse en un factor."
            },
            high_cor[k, 1], high_cor[k, 2],
            phi[high_cor[k, 1], high_cor[k, 2]]
          ))
        }
      }
    }
  }

  # Porcentaje alto de eliminacion
  if (nrow(items_eliminados) > 0) {
    pct <- nrow(items_eliminados) / p_original * 100
    if (pct > 30) {
      sug <- c(sug, sprintf(
        if (language == "eng") {
          "%.0f%% of the original items were removed. Review the quality of the instrument."
        } else {
          "Se eliminaron %.0f%% de los items originales. Revise la calidad del instrumento."
        },
        pct
      ))
    }
  }

  # Casos Heywood
  if (any(efa_obj$communality > 0.999)) {
    items_hey <- names(which(efa_obj$communality > 0.999))
    sug <- c(sug, sprintf(
      if (language == "eng") {
        "Heywood case(s): %s (communality ~1). The solution may be unstable."
      } else {
        "Caso(s) Heywood: %s (comunalidad ~1). La solucion puede ser inestable."
      },
      paste(items_hey, collapse = ", ")
    ))
  }

  sug
}


