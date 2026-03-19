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

# Etiquetas y textos del reporte por consola
.efa_labels <- function(language = "esp") {
  language <- .efa_normalize_language(language)

  if (language == "eng") {
    return(list(
      required_psych = "Package 'psych' is required. Install with: install.packages('psych')",
      required_gpa = "Package 'GPArotation' is required. Install with: install.packages('GPArotation')",
      invalid_data = "'datos' must be a data frame or a matrix.",
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
    invalid_data = "'datos' debe ser un data frame o una matriz.",
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


# ============================================================
#                   FUNCION PRINCIPAL
# ============================================================

#' Analisis Factorial Exploratorio Optimo
#'
#' Realiza un AFE automatizado con seleccion de estimador, analisis paralelo,
#' rotacion y eliminacion iterativa de items problematicos.
#'
#' @param datos Data frame o matriz con los items (variables numericas).
#' @param rotacion Tipo de rotacion: "oblicua" (oblimin) u "ortogonal" (varimax).
#' @param carga_min Carga factorial minima aceptable (default: 0.30).
#' @param comunalidad_min Comunalidad minima aceptable (default: 0.30).
#' @param dif_cargas_cruzadas Diferencia minima entre las dos cargas mas altas
#'        de un item para no considerarlo carga cruzada (default: 0.15).
#' @param max_iter Maximo de iteraciones de eliminacion (default: 20).
#' @param verbose Imprimir reporte en consola (default: TRUE).
#' @param exportar Formato de exportacion opcional usado por la funcion.
#' @param archivo Nombre base opcional para archivos exportados.
#' @param language Idioma del reporte: "esp" o "eng" (default: "esp").
#'
#' @return Lista (clase "efa_auto") con todos los resultados del analisis.
#'
#' @examples
#' # resultado <- efa_auto(mi_datos, rotacion = "oblicua")
#' # resultado <- efa_auto(mi_datos, rotacion = "ortogonal")
#' # resultado <- efa_auto(mi_datos, rotacion = "oblicua", language = "eng")
#'
#' @export

efa_auto <- function(datos,
                       rotacion = c("oblicua", "ortogonal"),
                       carga_min = 0.30,
                       comunalidad_min = 0.30,
                       dif_cargas_cruzadas = 0.15,
                       max_iter = 20,
                       verbose = TRUE,
                       exportar = NULL,
                       archivo = NULL,
                       language = "esp") {

  language <- .efa_normalize_language(language)
  txt <- .efa_labels(language)

  if (!requireNamespace("psych", quietly = TRUE))
    stop(txt$required_psych)
  if (!requireNamespace("GPArotation", quietly = TRUE))
    stop(txt$required_gpa)

  rotacion <- match.arg(rotacion)
  sep <- paste0(rep("=", 65), collapse = "")
  rotacion_label <- ifelse(rotacion == "oblicua", txt$oblique, txt$orthogonal)

  # ================================================================
  # 1. VALIDACION Y PREPARACION DE DATOS
  # ================================================================
  if (!is.data.frame(datos) && !is.matrix(datos))
    stop(txt$invalid_data)

  datos <- as.data.frame(datos)
  cols_num <- sapply(datos, is.numeric)

  if (sum(cols_num) < 3)
    stop(txt$need_numeric)

  if (sum(!cols_num) > 0 && verbose)
    cat(txt$excluded_non_numeric, sum(!cols_num), txt$non_numeric_cols, "\n")

  datos <- datos[, cols_num, drop = FALSE]

  n_original <- nrow(datos)
  datos <- na.omit(datos)
  n <- nrow(datos)
  p <- ncol(datos)

  if (verbose && n_original != n)
    cat(txt$removed_rows, n_original - n, txt$rows_with_na,
        txt$final_sample, n, "\n\n")

  # ================================================================
  # 2. CARACTERISTICAS DE LOS DATOS
  # ================================================================
  n_categorias <- sapply(datos, function(x) length(unique(x)))
  mediana_cat <- median(n_categorias)
  es_ordinal <- mediana_cat <= 5

  if (verbose) {
    cat(sep, "\n")
    cat(txt$analysis_title, "\n")
    cat(sep, "\n\n")
    cat(txt$data_title, "\n")
    cat(txt$sample_size, n, "\n")
    cat(txt$n_items, p, "\n")
    cat(txt$median_categories, mediana_cat, "\n")
    cat(txt$data_type,
        ifelse(es_ordinal, txt$data_type_ordinal, txt$data_type_continuous), "\n")
    cat(txt$categories_by_item, "\n")
    cat_tbl <- data.frame(
      Item = names(n_categorias),
      Categorias = as.integer(n_categorias)
    )
    for (i in seq_len(nrow(cat_tbl)))
      cat(sprintf("    %-25s %d\n", cat_tbl$Item[i], cat_tbl$Categorias[i]))
  }

  # ================================================================
  # 3. NORMALIDAD MULTIVARIADA (Test de Mardia)
  # ================================================================
  normalidad <- .efa_mardia(datos)

  if (verbose) {
    cat("\n", txt$normality_title, "\n", sep = "")
    if (language == "eng") {
      cat(sprintf("  Skewness: chi2 = %.3f, df = %.0f, p = %.4f\n",
                  normalidad$skewness_chi, normalidad$skewness_df,
                  normalidad$skewness_p))
      cat(sprintf("  Kurtosis:  z = %.3f, p = %.4f\n",
                  normalidad$kurtosis_z, normalidad$kurtosis_p))
    } else {
      cat(sprintf("  Asimetria: chi2 = %.3f, gl = %.0f, p = %.4f\n",
                  normalidad$skewness_chi, normalidad$skewness_df,
                  normalidad$skewness_p))
      cat(sprintf("  Curtosis:  z = %.3f, p = %.4f\n",
                  normalidad$kurtosis_z, normalidad$kurtosis_p))
    }
    cat(txt$normality_conclusion,
        ifelse(normalidad$normal, txt$normality_ok, txt$normality_bad))
  }

  # ================================================================
  # 4. SELECCION DEL ESTIMADOR
  # ================================================================
  est <- .efa_seleccionar_estimador(es_ordinal, normalidad$normal, n, p, language)

  if (verbose) {
    cat("\n", txt$estimator_title, "\n", sep = "")
    cat(txt$estimator_label, est$nombre, "\n")
    cat(txt$justification_label, est$justificacion, "\n")
    cat(txt$correlation_label, est$tipo_cor, "\n")
  }

  # ================================================================
  # 5. MATRIZ DE CORRELACIONES
  # ================================================================
  if (est$usar_policorica) {
    if (verbose) cat(txt$calculating_polychoric)
    cor_result <- tryCatch(
      polychoric(datos),
      error = function(e) {
        warning(txt$polychoric_warning, e$message)
        list(rho = cor(datos))
      }
    )
    R <- cor_result$rho
  } else {
    R <- cor(datos)
  }

  # Verificar que R es definida positiva
  eigvals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0)) {
    if (verbose) cat(txt$non_pd_matrix)
    R <- as.matrix(nearPD(R, corr = TRUE)$mat)
  }

  # ================================================================
  # 6. KMO Y TEST DE BARTLETT (sobre items originales)
  # ================================================================
  kmo_original <- KMO(R)
  bartlett_original <- cortest.bartlett(R, n = n)

  if (verbose) {
    cat("\n", txt$adequacy_title, "\n", sep = "")
    if (language == "eng") {
      cat(sprintf("  Overall KMO: %.3f - %s\n",
                  kmo_original$MSA, .efa_interpretar_kmo(kmo_original$MSA, language)))
      cat(sprintf("  Bartlett: chi2 = %.3f, df = %.0f, p = %s\n",
                  bartlett_original$chisq, bartlett_original$df,
                  formatC(bartlett_original$p.value, format = "e", digits = 4)))
    } else {
      cat(sprintf("  KMO global: %.3f - %s\n",
                  kmo_original$MSA, .efa_interpretar_kmo(kmo_original$MSA, language)))
      cat(sprintf("  Bartlett: chi2 = %.3f, gl = %.0f, p = %s\n",
                  bartlett_original$chisq, bartlett_original$df,
                  formatC(bartlett_original$p.value, format = "e", digits = 4)))
    }

    kmo_bajo <- kmo_original$MSAi[kmo_original$MSAi < 0.50]
    if (length(kmo_bajo) > 0) {
      cat(txt$kmo_warning_items)
      for (i in seq_along(kmo_bajo))
        cat(sprintf("    %-25s %.3f\n", names(kmo_bajo)[i], kmo_bajo[i]))
    }
  }

  if (kmo_original$MSA < 0.50)
    warning(txt$kmo_warning_global)

  # ================================================================
  # 7. ANALISIS PARALELO
  # ================================================================
  if (verbose) cat("\n", txt$parallel_title, "\n", sep = "")

  ap <- suppressWarnings(
    fa.parallel(R, n.obs = n, fa = "fa", fm = est$fm, plot = FALSE)
  )
  n_factores <- max(1, ap$nfact)

  if (verbose) cat(txt$optimal_factors, n_factores, "\n")

  # ================================================================
  # 8. SELECCION DE ROTACION
  # ================================================================
  if (n_factores == 1) {
    tipo_rot <- "none"
    if (verbose) cat(txt$no_rotation)
  } else {
    tipo_rot <- ifelse(rotacion == "oblicua", "oblimin", "varimax")
    if (verbose)
      cat(sprintf(txt$rotation_label, tipo_rot, rotacion_label))
  }

  # ================================================================
  # 9. PROCESO ITERATIVO DE REFINAMIENTO
  # ================================================================
  if (verbose) {
    cat("\n", sep, "\n", sep = "")
    cat(txt$refinement_title, "\n")
    cat(sep, "\n")
  }

  datos_actual <- datos
  R_actual <- R
  items_eliminados <- data.frame(
    Item = character(), Razon = character(),
    Iteracion = integer(), stringsAsFactors = FALSE
  )
  iteracion <- 0
  efa <- NULL

  while (iteracion < max_iter) {
    iteracion <- iteracion + 1
    p_actual <- ncol(datos_actual)

    # Recalcular correlaciones y numero de factores tras eliminacion
    if (iteracion > 1) {
      if (est$usar_policorica) {
        R_actual <- tryCatch(
          polychoric(datos_actual)$rho,
          error = function(e) cor(datos_actual)
        )
      } else {
        R_actual <- cor(datos_actual)
      }

      # Asegurar definida positiva
      eigvals_act <- eigen(R_actual, symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals_act <= 0))
        R_actual <- as.matrix(nearPD(R_actual, corr = TRUE)$mat)

      ap_new <- suppressWarnings(
        fa.parallel(R_actual, n.obs = n, fa = "fa", fm = est$fm, plot = FALSE)
      )
      nf_new <- max(1, ap_new$nfact)

      if (nf_new != n_factores) {
        if (verbose)
          cat(sprintf(txt$factors_adjusted, n_factores, nf_new))
        n_factores <- nf_new
        tipo_rot <- if (n_factores == 1) "none" else
          ifelse(rotacion == "oblicua", "oblimin", "varimax")
      }
    }

    # Verificar items suficientes
    if (p_actual < n_factores * 3) {
      if (n_factores > 1) {
        n_factores <- max(1, floor(p_actual / 3))
        tipo_rot <- if (n_factores == 1) "none" else
          ifelse(rotacion == "oblicua", "oblimin", "varimax")
        if (verbose)
          cat(txt$insufficient_items, n_factores, txt$factor_suffix)
      } else {
        break
      }
    }

    # Ejecutar AFE
    efa <- tryCatch(
      suppressWarnings(
        fa(R_actual, nfactors = n_factores, n.obs = n,
           rotate = tipo_rot, fm = est$fm)
      ),
      error = function(e) {
        if (est$fm != "minres") {
          if (verbose)
            cat(txt$estimator_error_note, est$fm, txt$switching_minres)
          suppressWarnings(
            fa(R_actual, nfactors = n_factores, n.obs = n,
               rotate = tipo_rot, fm = "minres")
          )
        } else {
          stop(if (language == "eng") "Error in EFA: " else "Error en AFE: ",
               e$message)
        }
      }
    )

    cargas <- as.data.frame(unclass(efa$loadings))
    comunalidades <- efa$communality

    if (verbose)
      cat(sprintf(txt$iteration_title, iteracion, p_actual, n_factores))

    # Evaluar items
    problemas <- .efa_evaluar_items(
      cargas, comunalidades, carga_min,
      comunalidad_min, dif_cargas_cruzadas, n_factores, language
    )

    if (nrow(problemas) == 0) {
      if (verbose) cat(txt$no_problem_items)
      break
    }

    if (verbose) {
      cat(txt$problematic_items)
      for (i in 1:nrow(problemas))
        cat(sprintf("    - %-20s %s\n",
                    problemas$Item[i], problemas$Razon[i]))
    }

    # Eliminar el item mas problematico
    item_out <- problemas$Item[1]
    razon_out <- problemas$Razon[1]

    items_eliminados <- rbind(items_eliminados, data.frame(
      Item = item_out, Razon = razon_out,
      Iteracion = iteracion, stringsAsFactors = FALSE
    ))

    datos_actual <- datos_actual[, colnames(datos_actual) != item_out,
                                 drop = FALSE]

    if (verbose)
      cat(sprintf(txt$removed_item, item_out, razon_out))

    if (ncol(datos_actual) < 3) {
      warning(txt$fewer_than_three)
      break
    }
  }

  if (is.null(efa))
    stop(txt$no_model)

  # ================================================================
  # 10. RESULTADOS FINALES
  # ================================================================
  cargas_final <- as.data.frame(unclass(efa$loadings))
  asignacion <- apply(abs(cargas_final), 1, which.max)

  distribucion <- data.frame(
    Item = rownames(cargas_final),
    Dimension = paste0("F", asignacion),
    Carga_Principal = sapply(seq_along(asignacion), function(i)
      cargas_final[i, asignacion[i]]),
    Comunalidad = efa$communality[rownames(cargas_final)],
    stringsAsFactors = FALSE
  )
  distribucion <- distribucion[order(distribucion$Dimension,
                                     -abs(distribucion$Carga_Principal)), ]
  rownames(distribucion) <- NULL

  var_explicada <- efa$Vaccounted
  chi2 <- if (!is.null(efa$STATISTIC)) efa$STATISTIC else NA
  gl   <- if (!is.null(efa$dof)) efa$dof else NA
  p_chi2 <- if (!is.null(efa$PVAL)) efa$PVAL else NA

  # Recalcular KMO y Bartlett sobre items finales
  R_final <- R_actual
  kmo_final <- KMO(R_final)
  bartlett_final <- cortest.bartlett(R_final, n = n)

  if (verbose) {
    cat("\n", sep, "\n", sep = "")
    cat(txt$final_results_title, "\n")
    cat(sep, "\n\n")

    cat(txt$fit_indices_title, "\n")
    if (!is.na(chi2))
      cat(sprintf(txt$chi_square, chi2))
    if (!is.na(gl))
      cat(sprintf(txt$degrees_freedom, as.integer(gl)))
    if (!is.na(p_chi2))
      cat(sprintf(txt$p_value, formatC(p_chi2, format = "e", digits = 4)))
    if (!is.null(efa$TLI))
      cat(sprintf(txt$tli, efa$TLI))
    if (!is.null(efa$RMSEA))
      cat(sprintf(txt$rmsea, efa$RMSEA[1], efa$RMSEA[2], efa$RMSEA[3]))
    if (!is.null(efa$BIC))
      cat(sprintf(txt$bic, efa$BIC))

    cat(sprintf(txt$final_kmo,
                kmo_final$MSA, .efa_interpretar_kmo(kmo_final$MSA, language)))
    cat(sprintf(txt$final_bartlett,
                bartlett_final$chisq, bartlett_final$df,
                formatC(bartlett_final$p.value, format = "e", digits = 4)))

    # --- Varianza explicada ---
    cat("\n", txt$variance_title, "\n", sep = "")
    if (is.matrix(var_explicada)) {
      for (j in 1:ncol(var_explicada)) {
        if ("Proportion Var" %in% rownames(var_explicada))
          cat(sprintf(txt$factor_variance, j,
                      var_explicada["Proportion Var", j] * 100))
      }
      if ("Cumulative Var" %in% rownames(var_explicada))
        cat(sprintf(txt$total_variance,
                    var_explicada["Cumulative Var",
                                  ncol(var_explicada)] * 100))
    }

    # --- Matriz de cargas factoriales ---
    cat("\n", txt$loadings_title, "\n", sep = "")
    cargas_print <- cargas_final
    for (i in 1:nrow(cargas_print))
      for (j in 1:ncol(cargas_print))
        cargas_print[i, j] <- ifelse(
          abs(cargas_final[i, j]) < carga_min,
          NA_real_,
          round(cargas_final[i, j], 3)
        )
    cargas_print$h2 <- round(efa$communality[rownames(cargas_final)], 3)
    print(cargas_print, na.print = "   ")
    cat(txt$loadings_note_left, carga_min, txt$loadings_note_right)

    # --- Distribucion de items ---
    cat("\n", txt$distribution_title, "\n", sep = "")
    for (f in sort(unique(distribucion$Dimension))) {
      items_f <- distribucion[distribucion$Dimension == f, ]
      cat(sprintf(txt$distribution_factor, f, nrow(items_f)))
      for (i in 1:nrow(items_f))
        cat(sprintf(txt$distribution_item,
                    items_f$Item[i], items_f$Carga_Principal[i],
                    items_f$Comunalidad[i]))
    }

    # --- Correlaciones entre factores ---
    if (rotacion == "oblicua" && n_factores > 1 && !is.null(efa$Phi)) {
      cat("\n", txt$factor_corr_title, "\n", sep = "")
      print(round(efa$Phi, 3))
    }

    # --- Items eliminados ---
    if (nrow(items_eliminados) > 0) {
      cat("\n", txt$removed_items_title, "\n", sep = "")
      for (i in 1:nrow(items_eliminados))
        cat(sprintf(txt$removed_items_line,
                    i, items_eliminados$Item[i],
                    items_eliminados$Razon[i],
                    items_eliminados$Iteracion[i]))
      cat(sprintf(txt$removed_items_total,
                  nrow(items_eliminados), p))
    } else {
      cat(txt$no_items_removed)
    }
  }

  # ================================================================
  # 11. CONFIABILIDAD POR DIMENSION
  # ================================================================
  if (verbose) cat("\n", txt$reliability_title, "\n", sep = "")

  fiabilidad <- list()
  for (f in sort(unique(distribucion$Dimension))) {
    items_f <- distribucion$Item[distribucion$Dimension == f]

    if (length(items_f) >= 3) {
      alpha_f <- tryCatch(
        alpha(datos_actual[, items_f, drop = FALSE], check.keys = TRUE),
        error = function(e) NULL
      )

      omega_val <- tryCatch({
        invisible(capture.output(
          om <- omega(datos_actual[, items_f, drop = FALSE],
                      nfactors = 1, plot = FALSE)
        ))
        om$omega.tot
      }, error = function(e) NA_real_)

      if (!is.null(alpha_f)) {
        fiabilidad[[f]] <- list(
          alpha = alpha_f$total$raw_alpha,
          omega = omega_val,
          items = items_f,
          n_items = length(items_f),
          alpha_drop = alpha_f$alpha.drop
        )

        if (verbose) {
          if (language == "eng") {
            cat(sprintf("\n  %s (%d items): alpha = %.3f",
                        f, length(items_f), alpha_f$total$raw_alpha))
          } else {
            cat(sprintf("\n  %s (%d items): alpha = %.3f",
                        f, length(items_f), alpha_f$total$raw_alpha))
          }
          if (!is.na(omega_val))
            cat(sprintf(", omega = %.3f", omega_val))
          cat("\n")

          # Items cuya eliminacion mejoraria alpha
          if (length(items_f) > 3 && !is.null(alpha_f$alpha.drop)) {
            mejora <- alpha_f$alpha.drop[
              alpha_f$alpha.drop$raw_alpha > alpha_f$total$raw_alpha,
              , drop = FALSE
            ]
            if (nrow(mejora) > 0) {
              cat(txt$alpha_improve)
              for (r in rownames(mejora))
                cat(sprintf(txt$alpha_if_removed,
                            r, mejora[r, "raw_alpha"],
                            mejora[r, "raw_alpha"] - alpha_f$total$raw_alpha))
            }
          }
        }
      } else {
        fiabilidad[[f]] <- list(
          alpha = NA, omega = NA,
          items = items_f, n_items = length(items_f)
        )
        if (verbose)
          cat(sprintf(txt$reliability_error, f))
      }

    } else if (length(items_f) == 2) {
      # Con 2 items, solo Spearman-Brown
      cor_2 <- cor(datos_actual[, items_f[1]], datos_actual[, items_f[2]],
                   use = "complete.obs")
      sb <- (2 * abs(cor_2)) / (1 + abs(cor_2))
      fiabilidad[[f]] <- list(
        alpha = NA, omega = NA, spearman_brown = sb,
        items = items_f, n_items = 2
      )
      if (verbose)
        cat(sprintf(txt$spearman_brown, f, sb, cor_2))
    } else {
      fiabilidad[[f]] <- list(
        alpha = NA, omega = NA,
        items = items_f, n_items = length(items_f)
      )
      if (verbose)
        cat(sprintf(txt$reliability_na, f, length(items_f)))
    }
  }

  # ================================================================
  # 12. VEREDICTO FINAL
  # ================================================================
  sugerencias <- .efa_sugerencias(
    efa, distribucion, fiabilidad, kmo_final, bartlett_final,
    n_factores, n, items_eliminados, carga_min, var_explicada, p, language
  )

  if (verbose) {
    cat("\n", sep, "\n", sep = "")
    if (length(sugerencias) == 0) {
      cat(txt$optimal_solution_title)
      cat(txt$optimal_solution_1)
      cat(txt$optimal_solution_2)
    } else {
      cat(txt$respec_title)
      cat(sep, "\n\n")
      for (i in seq_along(sugerencias))
        cat(sprintf("  %d. %s\n", i, sugerencias[[i]]))
      cat("\n")
    }
    cat(sep, "\n")

    # --- Referencias bibliograficas ---
    cat(txt$refs_title)
    cat(txt$refs_intro)

    cat(txt$refs_estimator)
    cat("    Lloret-Segura, S., Ferreres-Traver, A., Hernandez-Baeza, A., &\n")
    cat("      Tomas-Marco, I. (2014). El analisis factorial exploratorio de los\n")
    cat("      items: una guia practica, revisada y actualizada. Anales de\n")
    cat("      Psicologia, 30(3), 1151-1169.\n")
    cat("      https://doi.org/10.6018/analesps.30.3.199361\n\n")

    cat("    Ferrando, P. J., & Lorenzo-Seva, U. (2017). Program FACTOR at 10:\n")
    cat("      Origins, development and future directions. Psicothema, 29(2),\n")
    cat("      236-240. https://doi.org/10.7334/psicothema2016.304\n\n")

    cat(txt$refs_polychoric)
    cat("    Holgado-Tello, F. P., Chacon-Moscoso, S., Barbero-Garcia, I., &\n")
    cat("      Vila-Abad, E. (2010). Polychoric versus Pearson correlations in\n")
    cat("      exploratory and confirmatory factor analysis of ordinal variables.\n")
    cat("      Quality & Quantity, 44(1), 153-166.\n")
    cat("      https://doi.org/10.1007/s11135-008-9190-y\n\n")

    cat(txt$refs_parallel)
    cat("    Horn, J. L. (1965). A rationale and test for the number of factors\n")
    cat("      in factor analysis. Psychometrika, 30(2), 179-185.\n")
    cat("      https://doi.org/10.1007/BF02289447\n\n")

    cat("    Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality\n")
    cat("      assessment of ordered polytomous items with parallel analysis.\n")
    cat("      Psychological Methods, 16(2), 209-220.\n")
    cat("      https://doi.org/10.1037/a0023353\n\n")

    cat(txt$refs_normality)
    cat("    Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis\n")
    cat("      with applications. Biometrika, 57(3), 519-530.\n")
    cat("      https://doi.org/10.1093/biomet/57.3.519\n\n")

    cat(txt$refs_kmo)
    cat("    Kaiser, H. F. (1974). An index of factorial simplicity. Psychometrika,\n")
    cat("      39(1), 31-36. https://doi.org/10.1007/BF02291575\n\n")

    cat(txt$refs_elimination)
    cat("    Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019).\n")
    cat("      Multivariate data analysis (8th ed.). Cengage.\n\n")

    cat("    Watkins, M. W. (2018). Exploratory factor analysis: A guide to best\n")
    cat("      practice. Journal of Black Psychology, 44(3), 219-246.\n")
    cat("      https://doi.org/10.1177/0095798418771807\n\n")

    cat("    Howard, M. C. (2016). A review of exploratory factor analysis\n")
    cat("      decisions and overview of current practices. International Journal\n")
    cat("      of Human-Computer Interaction, 32(1), 51-62.\n")
    cat("      https://doi.org/10.1080/10447318.2015.1087664\n\n")

    cat("    Costello, A. B., & Osborne, J. W. (2005). Best practices in\n")
    cat("      exploratory factor analysis: Four recommendations for getting the\n")
    cat("      most from your analysis. Practical Assessment, Research, and\n")
    cat("      Evaluation, 10(7), 1-9. https://doi.org/10.7275/jyj1-4868\n\n")

    cat(txt$refs_reliability)
    cat("    McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence\n")
    cat("      Erlbaum Associates.\n\n")

    cat("    Dunn, T. J., Baguley, T., & Brunsden, V. (2014). From alpha to\n")
    cat("      omega: A practical solution to the pervasive problem of internal\n")
    cat("      consistency estimation. British Journal of Psychology, 105(3),\n")
    cat("      399-412. https://doi.org/10.1111/bjop.12046\n\n")

    cat(txt$refs_software)
    cat("    Revelle, W. (2024). psych: Procedures for psychological,\n")
    cat("      psychometric, and personality research. R package.\n")
    cat("      https://CRAN.R-project.org/package=psych\n\n")

    cat(sep, "\n")
  }

  # ================================================================
  # OBJETO DE RETORNO
  # ================================================================
  resultado <- list(
    efa                  = efa,
    cargas               = cargas_final,
    comunalidades        = efa$communality,
    distribucion         = distribucion,
    language             = language,
    n_factores           = n_factores,
    estimador            = est,
    rotacion             = tipo_rot,
    kmo_original         = kmo_original,
    kmo_final            = kmo_final,
    bartlett_original    = bartlett_original,
    bartlett_final       = bartlett_final,
    chi_cuadrado         = chi2,
    grados_libertad      = gl,
    p_valor              = p_chi2,
    varianza_explicada   = var_explicada,
    normalidad           = normalidad,
    fiabilidad           = fiabilidad,
    items_eliminados     = items_eliminados,
    correlacion_factores = if (!is.null(efa$Phi)) efa$Phi else NULL,
    sugerencias          = sugerencias,
    solucion_optima      = length(sugerencias) == 0,
    datos_finales        = datos_actual,
    n                    = n,
    p_original           = p,
    p_final              = ncol(datos_actual)
  )

  class(resultado) <- "efa_auto"

  # Exportar si se solicito
  if (!is.null(exportar)) {
    exportar_efa(resultado, formato = exportar, archivo = archivo)
  }

  invisible(resultado)
}

# ============================================================
#                FUNCIONES DE EXPORTACION
# ============================================================

# --- Exportacion a Excel ---
.efa_exportar_excel <- function(res, archivo) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Paquete 'openxlsx' requerido. Instale con: install.packages('openxlsx')")

  wb <- createWorkbook()

  # Estilos reutilizables
  estilo_titulo <- createStyle(
    fontSize = 14, fontColour = "#FFFFFF", fgFill = "#2C3E50",
    halign = "center", textDecoration = "bold",
    border = "Bottom", borderColour = "#2C3E50"
  )
  estilo_encabezado <- createStyle(
    fontSize = 11, fontColour = "#FFFFFF", fgFill = "#3498DB",
    halign = "center", textDecoration = "bold",
    border = "TopBottomLeftRight", borderColour = "#BDC3C7"
  )
  estilo_seccion <- createStyle(
    fontSize = 11, fgFill = "#ECF0F1", textDecoration = "bold",
    border = "Bottom", borderColour = "#BDC3C7"
  )
  estilo_celda <- createStyle(
    border = "TopBottomLeftRight", borderColour = "#ECF0F1",
    halign = "center"
  )
  estilo_celda_izq <- createStyle(
    border = "TopBottomLeftRight", borderColour = "#ECF0F1"
  )
  estilo_numero <- createStyle(
    border = "TopBottomLeftRight", borderColour = "#ECF0F1",
    halign = "center", numFmt = "0.000"
  )
  estilo_bueno <- createStyle(fgFill = "#D5F5E3")
  estilo_alerta <- createStyle(fgFill = "#FADBD8")

  # ============================================
  # Hoja 1: RESUMEN GENERAL
  # ============================================
  addWorksheet(wb, "Resumen")

  fila <- 1
  mergeCells(wb, "Resumen", cols = 1:4, rows = fila)
  writeData(wb, "Resumen", "ANALISIS FACTORIAL EXPLORATORIO - RESUMEN",
            startRow = fila, startCol = 1)
  addStyle(wb, "Resumen", estilo_titulo, rows = fila, cols = 1:4)
  fila <- fila + 2

  # Info basica
  info <- data.frame(
    Parametro = c(
      "Tamano de muestra (n)",
      "Items originales",
      "Items finales",
      "Items eliminados",
      "Factores retenidos",
      "Estimador",
      "Tipo de correlacion",
      "Rotacion",
      "",
      "ADECUACION MUESTRAL",
      "KMO global (solucion final)",
      "Interpretacion KMO",
      "Bartlett chi-cuadrado",
      "Bartlett gl",
      "Bartlett p-valor",
      "",
      "NORMALIDAD MULTIVARIADA (Mardia)",
      "Asimetria chi2",
      "Asimetria p-valor",
      "Curtosis z",
      "Curtosis p-valor",
      "Conclusion normalidad",
      "",
      "INDICES DE AJUSTE DEL MODELO",
      "Chi-cuadrado",
      "Grados de libertad",
      "p-valor",
      "TLI",
      "RMSEA",
      "RMSEA IC 90% inferior",
      "RMSEA IC 90% superior",
      "BIC"
    ),
    Valor = c(
      as.character(res$n),
      as.character(res$p_original),
      as.character(res$p_final),
      as.character(nrow(res$items_eliminados)),
      as.character(res$n_factores),
      res$estimador$nombre,
      res$estimador$tipo_cor,
      res$rotacion,
      "",
      "",
      sprintf("%.3f", res$kmo_final$MSA),
      .efa_interpretar_kmo(res$kmo_final$MSA),
      sprintf("%.3f", res$bartlett_final$chisq),
      sprintf("%.0f", res$bartlett_final$df),
      formatC(res$bartlett_final$p.value, format = "e", digits = 4),
      "",
      "",
      sprintf("%.3f", res$normalidad$skewness_chi),
      sprintf("%.4f", res$normalidad$skewness_p),
      sprintf("%.3f", res$normalidad$kurtosis_z),
      sprintf("%.4f", res$normalidad$kurtosis_p),
      ifelse(res$normalidad$normal, "Se cumple", "No se cumple"),
      "",
      "",
      ifelse(is.na(res$chi_cuadrado), "N/D", sprintf("%.3f", res$chi_cuadrado)),
      ifelse(is.na(res$grados_libertad), "N/D",
             as.character(as.integer(res$grados_libertad))),
      ifelse(is.na(res$p_valor), "N/D",
             formatC(res$p_valor, format = "e", digits = 4)),
      ifelse(is.null(res$efa$TLI), "N/D", sprintf("%.3f", res$efa$TLI)),
      ifelse(is.null(res$efa$RMSEA), "N/D", sprintf("%.3f", res$efa$RMSEA[1])),
      ifelse(is.null(res$efa$RMSEA), "N/D", sprintf("%.3f", res$efa$RMSEA[2])),
      ifelse(is.null(res$efa$RMSEA), "N/D", sprintf("%.3f", res$efa$RMSEA[3])),
      ifelse(is.null(res$efa$BIC), "N/D", sprintf("%.3f", res$efa$BIC))
    ),
    stringsAsFactors = FALSE
  )

  writeData(wb, "Resumen", info, startRow = fila, headerStyle = estilo_encabezado)
  addStyle(wb, "Resumen", estilo_celda_izq,
           rows = (fila + 1):(fila + nrow(info)), cols = 1, gridExpand = TRUE)
  addStyle(wb, "Resumen", estilo_celda,
           rows = (fila + 1):(fila + nrow(info)), cols = 2, gridExpand = TRUE)

  # Resaltar filas de seccion
  filas_seccion <- which(info$Parametro %in%
    c("ADECUACION MUESTRAL", "NORMALIDAD MULTIVARIADA (Mardia)",
      "INDICES DE AJUSTE DEL MODELO"))
  for (fs in filas_seccion)
    addStyle(wb, "Resumen", estilo_seccion,
             rows = fila + fs, cols = 1:2, gridExpand = TRUE)

  setColWidths(wb, "Resumen", cols = 1:2, widths = c(40, 35))

  # Veredicto
  fila_v <- fila + nrow(info) + 2
  mergeCells(wb, "Resumen", cols = 1:2, rows = fila_v)
  if (res$solucion_optima) {
    writeData(wb, "Resumen", "SOLUCION OPTIMA", startRow = fila_v)
    addStyle(wb, "Resumen", createStyle(
      fontSize = 13, textDecoration = "bold", fontColour = "#27AE60",
      halign = "center", fgFill = "#D5F5E3"
    ), rows = fila_v, cols = 1:2, gridExpand = TRUE)
  } else {
    writeData(wb, "Resumen", "SE REQUIEREN REESPECIFICACIONES",
              startRow = fila_v)
    addStyle(wb, "Resumen", createStyle(
      fontSize = 13, textDecoration = "bold", fontColour = "#E74C3C",
      halign = "center", fgFill = "#FADBD8"
    ), rows = fila_v, cols = 1:2, gridExpand = TRUE)
  }

  # ============================================
  # Hoja 2: CARGAS FACTORIALES
  # ============================================
  addWorksheet(wb, "Cargas Factoriales")

  cargas_exp <- round(res$cargas, 3)
  cargas_exp$h2 <- round(res$comunalidades[rownames(cargas_exp)], 3)

  fila <- 1
  mergeCells(wb, "Cargas Factoriales",
             cols = 1:(ncol(cargas_exp) + 1), rows = fila)
  writeData(wb, "Cargas Factoriales",
            "MATRIZ DE CARGAS FACTORIALES", startRow = fila)
  addStyle(wb, "Cargas Factoriales", estilo_titulo,
           rows = fila, cols = 1:(ncol(cargas_exp) + 1))
  fila <- fila + 2

  # Escribir nombres de items en col 1
  items_names <- data.frame(Item = rownames(cargas_exp), stringsAsFactors = FALSE)
  cargas_tabla <- cbind(items_names, cargas_exp)
  rownames(cargas_tabla) <- NULL

  writeData(wb, "Cargas Factoriales", cargas_tabla,
            startRow = fila, headerStyle = estilo_encabezado)
  addStyle(wb, "Cargas Factoriales", estilo_celda_izq,
           rows = (fila + 1):(fila + nrow(cargas_tabla)), cols = 1,
           gridExpand = TRUE)
  addStyle(wb, "Cargas Factoriales", estilo_numero,
           rows = (fila + 1):(fila + nrow(cargas_tabla)),
           cols = 2:ncol(cargas_tabla), gridExpand = TRUE)

  # Resaltar cargas >= 0.30 en verde, y la columna h2
  for (i in 1:nrow(cargas_exp)) {
    for (j in 1:(ncol(cargas_exp) - 1)) {
      if (abs(cargas_exp[i, j]) >= 0.30) {
        addStyle(wb, "Cargas Factoriales",
                 createStyle(
                   fgFill = "#D5F5E3", halign = "center", numFmt = "0.000",
                   border = "TopBottomLeftRight", borderColour = "#ECF0F1"
                 ),
                 rows = fila + i, cols = j + 1)
      }
    }
  }

  setColWidths(wb, "Cargas Factoriales",
               cols = 1:ncol(cargas_tabla),
               widths = c(20, rep(12, ncol(cargas_tabla) - 1)))

  # ============================================
  # Hoja 3: DISTRIBUCION POR DIMENSION
  # ============================================
  addWorksheet(wb, "Distribucion")

  fila <- 1
  mergeCells(wb, "Distribucion", cols = 1:4, rows = fila)
  writeData(wb, "Distribucion",
            "DISTRIBUCION DE ITEMS POR DIMENSION", startRow = fila)
  addStyle(wb, "Distribucion", estilo_titulo, rows = fila, cols = 1:4)
  fila <- fila + 2

  dist_exp <- res$distribucion
  dist_exp$Carga_Principal <- round(dist_exp$Carga_Principal, 3)
  dist_exp$Comunalidad <- round(dist_exp$Comunalidad, 3)

  writeData(wb, "Distribucion", dist_exp,
            startRow = fila, headerStyle = estilo_encabezado)
  addStyle(wb, "Distribucion", estilo_celda_izq,
           rows = (fila + 1):(fila + nrow(dist_exp)), cols = 1:2,
           gridExpand = TRUE)
  addStyle(wb, "Distribucion", estilo_numero,
           rows = (fila + 1):(fila + nrow(dist_exp)), cols = 3:4,
           gridExpand = TRUE)

  # Colorear por dimension
  colores_dim <- c("#D6EAF8", "#D5F5E3", "#FDEBD0", "#F5B7B1",
                   "#D2B4DE", "#A9DFBF", "#F9E79F", "#AED6F1")
  dims_unicas <- sort(unique(dist_exp$Dimension))
  for (d in seq_along(dims_unicas)) {
    filas_d <- which(dist_exp$Dimension == dims_unicas[d])
    col_d <- colores_dim[((d - 1) %% length(colores_dim)) + 1]
    for (fd in filas_d)
      addStyle(wb, "Distribucion",
               createStyle(fgFill = col_d, border = "TopBottomLeftRight",
                           borderColour = "#ECF0F1"),
               rows = fila + fd, cols = 1:4, gridExpand = TRUE)
  }

  setColWidths(wb, "Distribucion", cols = 1:4, widths = c(20, 12, 18, 14))

  # ============================================
  # Hoja 4: VARIANZA EXPLICADA
  # ============================================
  addWorksheet(wb, "Varianza Explicada")

  fila <- 1
  mergeCells(wb, "Varianza Explicada",
             cols = 1:(res$n_factores + 1), rows = fila)
  writeData(wb, "Varianza Explicada",
            "VARIANZA EXPLICADA", startRow = fila)
  addStyle(wb, "Varianza Explicada", estilo_titulo,
           rows = fila, cols = 1:(res$n_factores + 1))
  fila <- fila + 2

  if (is.matrix(res$varianza_explicada)) {
    var_df <- as.data.frame(round(res$varianza_explicada, 4))
    var_df <- cbind(Indicador = rownames(var_df), var_df)
    rownames(var_df) <- NULL
    writeData(wb, "Varianza Explicada", var_df,
              startRow = fila, headerStyle = estilo_encabezado)
    addStyle(wb, "Varianza Explicada", estilo_celda_izq,
             rows = (fila + 1):(fila + nrow(var_df)), cols = 1,
             gridExpand = TRUE)
    addStyle(wb, "Varianza Explicada", estilo_numero,
             rows = (fila + 1):(fila + nrow(var_df)),
             cols = 2:ncol(var_df), gridExpand = TRUE)
    setColWidths(wb, "Varianza Explicada",
                 cols = 1:ncol(var_df),
                 widths = c(25, rep(14, ncol(var_df) - 1)))
  }

  # ============================================
  # Hoja 5: CONFIABILIDAD
  # ============================================
  addWorksheet(wb, "Confiabilidad")

  fila <- 1
  mergeCells(wb, "Confiabilidad", cols = 1:5, rows = fila)
  writeData(wb, "Confiabilidad",
            "CONFIABILIDAD POR DIMENSION", startRow = fila)
  addStyle(wb, "Confiabilidad", estilo_titulo, rows = fila, cols = 1:5)
  fila <- fila + 2

  fiab_df <- do.call(rbind, lapply(names(res$fiabilidad), function(f) {
    fi <- res$fiabilidad[[f]]
    data.frame(
      Dimension = f,
      N_Items = fi$n_items,
      Alpha = ifelse(is.na(fi$alpha), NA, round(fi$alpha, 3)),
      Omega = ifelse(is.na(fi$omega), NA, round(fi$omega, 3)),
      Items = paste(fi$items, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))

  writeData(wb, "Confiabilidad", fiab_df,
            startRow = fila, headerStyle = estilo_encabezado)
  addStyle(wb, "Confiabilidad", estilo_celda,
           rows = (fila + 1):(fila + nrow(fiab_df)), cols = 1:4,
           gridExpand = TRUE)
  addStyle(wb, "Confiabilidad", estilo_celda_izq,
           rows = (fila + 1):(fila + nrow(fiab_df)), cols = 5,
           gridExpand = TRUE)

  # Colorear alpha segun valor
  for (i in 1:nrow(fiab_df)) {
    if (!is.na(fiab_df$Alpha[i])) {
      color_a <- if (fiab_df$Alpha[i] >= 0.70) "#D5F5E3" else
                 if (fiab_df$Alpha[i] >= 0.60) "#FDEBD0" else "#FADBD8"
      addStyle(wb, "Confiabilidad",
               createStyle(fgFill = color_a, halign = "center",
                           border = "TopBottomLeftRight",
                           borderColour = "#ECF0F1"),
               rows = fila + i, cols = 3)
    }
  }

  setColWidths(wb, "Confiabilidad",
               cols = 1:5, widths = c(12, 10, 10, 10, 50))

  # ============================================
  # Hoja 6: ITEMS ELIMINADOS (si aplica)
  # ============================================
  if (nrow(res$items_eliminados) > 0) {
    addWorksheet(wb, "Items Eliminados")

    fila <- 1
    mergeCells(wb, "Items Eliminados", cols = 1:3, rows = fila)
    writeData(wb, "Items Eliminados",
              "ITEMS ELIMINADOS", startRow = fila)
    addStyle(wb, "Items Eliminados", estilo_titulo,
             rows = fila, cols = 1:3)
    fila <- fila + 2

    writeData(wb, "Items Eliminados", res$items_eliminados,
              startRow = fila, headerStyle = estilo_encabezado)
    addStyle(wb, "Items Eliminados", estilo_celda_izq,
             rows = (fila + 1):(fila + nrow(res$items_eliminados)),
             cols = 1:3, gridExpand = TRUE)
    setColWidths(wb, "Items Eliminados",
                 cols = 1:3, widths = c(20, 55, 12))
  }

  # ============================================
  # Hoja 7: CORRELACIONES ENTRE FACTORES (si oblicua)
  # ============================================
  if (!is.null(res$correlacion_factores)) {
    addWorksheet(wb, "Corr Factores")

    fila <- 1
    nf <- ncol(res$correlacion_factores)
    mergeCells(wb, "Corr Factores", cols = 1:(nf + 1), rows = fila)
    writeData(wb, "Corr Factores",
              "CORRELACIONES ENTRE FACTORES", startRow = fila)
    addStyle(wb, "Corr Factores", estilo_titulo,
             rows = fila, cols = 1:(nf + 1))
    fila <- fila + 2

    phi_df <- as.data.frame(round(res$correlacion_factores, 3))
    phi_df <- cbind(Factor = rownames(phi_df), phi_df)
    rownames(phi_df) <- NULL

    writeData(wb, "Corr Factores", phi_df,
              startRow = fila, headerStyle = estilo_encabezado)
    addStyle(wb, "Corr Factores", estilo_celda,
             rows = (fila + 1):(fila + nrow(phi_df)),
             cols = 1:ncol(phi_df), gridExpand = TRUE)
    setColWidths(wb, "Corr Factores",
                 cols = 1:ncol(phi_df),
                 widths = c(12, rep(12, ncol(phi_df) - 1)))
  }

  # ============================================
  # Hoja 8: SUGERENCIAS (si las hay)
  # ============================================
  if (length(res$sugerencias) > 0) {
    addWorksheet(wb, "Sugerencias")

    fila <- 1
    mergeCells(wb, "Sugerencias", cols = 1:2, rows = fila)
    writeData(wb, "Sugerencias",
              "SUGERENCIAS DE REESPECIFICACION", startRow = fila)
    addStyle(wb, "Sugerencias", estilo_titulo, rows = fila, cols = 1:2)
    fila <- fila + 2

    sug_df <- data.frame(
      N = seq_along(res$sugerencias),
      Sugerencia = unlist(res$sugerencias),
      stringsAsFactors = FALSE
    )
    writeData(wb, "Sugerencias", sug_df,
              startRow = fila, headerStyle = estilo_encabezado)
    addStyle(wb, "Sugerencias", estilo_celda_izq,
             rows = (fila + 1):(fila + nrow(sug_df)),
             cols = 1:2, gridExpand = TRUE)
    setColWidths(wb, "Sugerencias", cols = 1:2, widths = c(5, 80))
  }

  # ============================================
  # Hoja 9: REFERENCIAS
  # ============================================
  addWorksheet(wb, "Referencias")

  fila <- 1
  mergeCells(wb, "Referencias", cols = 1:2, rows = fila)
  writeData(wb, "Referencias",
            "REFERENCIAS METODOLOGICAS", startRow = fila)
  addStyle(wb, "Referencias", estilo_titulo, rows = fila, cols = 1:2)
  fila <- fila + 2

  refs <- data.frame(
    Aspecto = c(
      "Estimador ULS/MINRES",
      "Estimador ULS/MINRES",
      "Correlaciones policoricas",
      "Analisis paralelo",
      "Analisis paralelo",
      "Normalidad multivariada",
      "KMO",
      "Buenas practicas AFE",
      "Buenas practicas AFE",
      "Buenas practicas AFE",
      "Criterios de eliminacion",
      "Confiabilidad (omega)",
      "Confiabilidad (alpha/omega)",
      "Software"
    ),
    Referencia = c(
      "Lloret-Segura, S., Ferreres-Traver, A., Hernandez-Baeza, A., & Tomas-Marco, I. (2014). El analisis factorial exploratorio de los items: una guia practica, revisada y actualizada. Anales de Psicologia, 30(3), 1151-1169. https://doi.org/10.6018/analesps.30.3.199361",
      "Ferrando, P. J., & Lorenzo-Seva, U. (2017). Program FACTOR at 10: Origins, development and future directions. Psicothema, 29(2), 236-240. https://doi.org/10.7334/psicothema2016.304",
      "Holgado-Tello, F. P., Chacon-Moscoso, S., Barbero-Garcia, I., & Vila-Abad, E. (2010). Polychoric versus Pearson correlations in exploratory and confirmatory factor analysis of ordinal variables. Quality & Quantity, 44(1), 153-166. https://doi.org/10.1007/s11135-008-9190-y",
      "Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30(2), 179-185. https://doi.org/10.1007/BF02289447",
      "Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. Psychological Methods, 16(2), 209-220. https://doi.org/10.1037/a0023353",
      "Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. Biometrika, 57(3), 519-530. https://doi.org/10.1093/biomet/57.3.519",
      "Kaiser, H. F. (1974). An index of factorial simplicity. Psychometrika, 39(1), 31-36. https://doi.org/10.1007/BF02291575",
      "Watkins, M. W. (2018). Exploratory factor analysis: A guide to best practice. Journal of Black Psychology, 44(3), 219-246. https://doi.org/10.1177/0095798418771807",
      "Howard, M. C. (2016). A review of exploratory factor analysis decisions and overview of current practices. International Journal of Human-Computer Interaction, 32(1), 51-62. https://doi.org/10.1080/10447318.2015.1087664",
      "Costello, A. B., & Osborne, J. W. (2005). Best practices in exploratory factor analysis: Four recommendations for getting the most from your analysis. Practical Assessment, Research, and Evaluation, 10(7), 1-9. https://doi.org/10.7275/jyj1-4868",
      "Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019). Multivariate data analysis (8th ed.). Cengage.",
      "McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence Erlbaum Associates.",
      "Dunn, T. J., Baguley, T., & Brunsden, V. (2014). From alpha to omega: A practical solution to the pervasive problem of internal consistency estimation. British Journal of Psychology, 105(3), 399-412. https://doi.org/10.1111/bjop.12046",
      "Revelle, W. (2024). psych: Procedures for psychological, psychometric, and personality research. R package. https://CRAN.R-project.org/package=psych"
    ),
    stringsAsFactors = FALSE
  )

  writeData(wb, "Referencias", refs,
            startRow = fila, headerStyle = estilo_encabezado)
  addStyle(wb, "Referencias", estilo_celda_izq,
           rows = (fila + 1):(fila + nrow(refs)), cols = 1:2,
           gridExpand = TRUE)
  setColWidths(wb, "Referencias", cols = 1:2, widths = c(30, 100))

  # Guardar
  if (!grepl("\\.xlsx$", archivo)) archivo <- paste0(archivo, ".xlsx")
  saveWorkbook(wb, archivo, overwrite = TRUE)
  cat(sprintf("  Archivo Excel exportado: %s\n", normalizePath(archivo)))
}


# --- Exportacion a Word ---
.efa_exportar_word <- function(res, archivo) {
  if (!requireNamespace("officer", quietly = TRUE))
    stop("Paquete 'officer' requerido. Instale con: install.packages('officer')")
  if (!requireNamespace("flextable", quietly = TRUE))
    stop("Paquete 'flextable' requerido. Instale con: install.packages('flextable')")

  doc <- read_docx()

  # --- Titulo ---
  doc <- body_add_par(doc,
    "ANALISIS FACTORIAL EXPLORATORIO - REPORTE DE RESULTADOS",
    style = "heading 1")
  doc <- body_add_par(doc,
    paste("Fecha:", format(Sys.time(), "%d/%m/%Y %H:%M")),
    style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Caracteristicas de los datos ---
  doc <- body_add_par(doc, "1. Caracteristicas de los datos", style = "heading 2")

  info_datos <- data.frame(
    Parametro = c("Tamano de muestra", "Items originales", "Items finales",
                  "Items eliminados", "Factores retenidos", "Estimador",
                  "Correlacion", "Rotacion"),
    Valor = c(
      as.character(res$n), as.character(res$p_original),
      as.character(res$p_final), as.character(nrow(res$items_eliminados)),
      as.character(res$n_factores), res$estimador$nombre,
      res$estimador$tipo_cor, res$rotacion
    ),
    stringsAsFactors = FALSE
  )

  ft <- flextable(info_datos) |>
    set_header_labels(Parametro = "Parametro", Valor = "Valor") |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    width(j = 1, width = 2.5) |>
    width(j = 2, width = 4) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")

  doc <- body_add_par(doc,
    paste("Justificacion del estimador:", res$estimador$justificacion),
    style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Normalidad ---
  doc <- body_add_par(doc,
    "2. Normalidad multivariada (Test de Mardia)", style = "heading 2")

  norm_df <- data.frame(
    Estadistico = c("Asimetria (chi2)", "Asimetria (p)",
                    "Curtosis (z)", "Curtosis (p)", "Conclusion"),
    Valor = c(
      sprintf("%.3f", res$normalidad$skewness_chi),
      sprintf("%.4f", res$normalidad$skewness_p),
      sprintf("%.3f", res$normalidad$kurtosis_z),
      sprintf("%.4f", res$normalidad$kurtosis_p),
      ifelse(res$normalidad$normal, "Se cumple", "No se cumple")
    ),
    stringsAsFactors = FALSE
  )
  ft <- flextable(norm_df) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    width(j = 1, width = 2.5) |>
    width(j = 2, width = 2) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Adecuacion muestral ---
  doc <- body_add_par(doc, "3. Adecuacion muestral", style = "heading 2")

  adc_df <- data.frame(
    Indicador = c("KMO global", "Interpretacion KMO",
                  "Bartlett chi2", "Bartlett gl", "Bartlett p"),
    Valor = c(
      sprintf("%.3f", res$kmo_final$MSA),
      .efa_interpretar_kmo(res$kmo_final$MSA),
      sprintf("%.3f", res$bartlett_final$chisq),
      sprintf("%.0f", res$bartlett_final$df),
      formatC(res$bartlett_final$p.value, format = "e", digits = 4)
    ),
    stringsAsFactors = FALSE
  )
  ft <- flextable(adc_df) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    width(j = 1, width = 2.5) |>
    width(j = 2, width = 2) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Indices de ajuste ---
  doc <- body_add_par(doc, "4. Indices de ajuste", style = "heading 2")

  ajuste_params <- c("Chi-cuadrado", "Grados de libertad", "p-valor",
                     "TLI", "RMSEA", "RMSEA IC inferior", "RMSEA IC superior",
                     "BIC")
  ajuste_vals <- c(
    ifelse(is.na(res$chi_cuadrado), "N/D", sprintf("%.3f", res$chi_cuadrado)),
    ifelse(is.na(res$grados_libertad), "N/D",
           as.character(as.integer(res$grados_libertad))),
    ifelse(is.na(res$p_valor), "N/D",
           formatC(res$p_valor, format = "e", digits = 4)),
    ifelse(is.null(res$efa$TLI), "N/D", sprintf("%.3f", res$efa$TLI)),
    ifelse(is.null(res$efa$RMSEA), "N/D", sprintf("%.3f", res$efa$RMSEA[1])),
    ifelse(is.null(res$efa$RMSEA), "N/D", sprintf("%.3f", res$efa$RMSEA[2])),
    ifelse(is.null(res$efa$RMSEA), "N/D", sprintf("%.3f", res$efa$RMSEA[3])),
    ifelse(is.null(res$efa$BIC), "N/D", sprintf("%.3f", res$efa$BIC))
  )
  ajuste_df <- data.frame(Indicador = ajuste_params, Valor = ajuste_vals,
                           stringsAsFactors = FALSE)
  ft <- flextable(ajuste_df) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    width(j = 1, width = 2.5) |>
    width(j = 2, width = 2) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Varianza explicada ---
  doc <- body_add_par(doc, "5. Varianza explicada", style = "heading 2")

  if (is.matrix(res$varianza_explicada)) {
    var_df <- as.data.frame(round(res$varianza_explicada, 4))
    var_df <- cbind(Indicador = rownames(var_df), var_df)
    rownames(var_df) <- NULL
    ft <- flextable(var_df) |>
      bold(part = "header") |>
      bg(part = "header", bg = "#3498DB") |>
      color(part = "header", color = "white") |>
      colformat_double(digits = 4) |>
      width(j = 1, width = 2.5) |>
      theme_box()
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "", style = "Normal")
  }

  # --- Cargas factoriales ---
  doc <- body_add_par(doc, "6. Matriz de cargas factoriales", style = "heading 2")

  cargas_w <- round(res$cargas, 3)
  cargas_w$h2 <- round(res$comunalidades[rownames(cargas_w)], 3)
  cargas_w <- cbind(Item = rownames(cargas_w), cargas_w)
  rownames(cargas_w) <- NULL

  ft <- flextable(cargas_w) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    colformat_double(digits = 3) |>
    width(j = 1, width = 1.5) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Distribucion ---
  doc <- body_add_par(doc,
    "7. Distribucion de items por dimension", style = "heading 2")

  dist_w <- res$distribucion
  dist_w$Carga_Principal <- round(dist_w$Carga_Principal, 3)
  dist_w$Comunalidad <- round(dist_w$Comunalidad, 3)

  ft <- flextable(dist_w) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    colformat_double(digits = 3) |>
    width(j = 1, width = 1.5) |>
    width(j = 2, width = 1.2) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")

  # --- Correlaciones entre factores ---
  if (!is.null(res$correlacion_factores)) {
    doc <- body_add_par(doc,
      "8. Correlaciones entre factores", style = "heading 2")
    phi_w <- as.data.frame(round(res$correlacion_factores, 3))
    phi_w <- cbind(Factor = rownames(phi_w), phi_w)
    rownames(phi_w) <- NULL
    ft <- flextable(phi_w) |>
      bold(part = "header") |>
      bg(part = "header", bg = "#3498DB") |>
      color(part = "header", color = "white") |>
      colformat_double(digits = 3) |>
      theme_box()
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "", style = "Normal")
  }

  # --- Items eliminados ---
  sec_num <- ifelse(is.null(res$correlacion_factores), 8, 9)
  if (nrow(res$items_eliminados) > 0) {
    doc <- body_add_par(doc,
      paste0(sec_num, ". Items eliminados"), style = "heading 2")
    ft <- flextable(res$items_eliminados) |>
      bold(part = "header") |>
      bg(part = "header", bg = "#E74C3C") |>
      color(part = "header", color = "white") |>
      width(j = 1, width = 1.5) |>
      width(j = 2, width = 4) |>
      width(j = 3, width = 1) |>
      theme_box()
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "", style = "Normal")
    sec_num <- sec_num + 1
  }

  # --- Confiabilidad ---
  doc <- body_add_par(doc,
    paste0(sec_num, ". Confiabilidad por dimension"), style = "heading 2")

  fiab_w <- do.call(rbind, lapply(names(res$fiabilidad), function(f) {
    fi <- res$fiabilidad[[f]]
    data.frame(
      Dimension = f,
      N_Items = fi$n_items,
      Alpha = ifelse(is.na(fi$alpha), NA, round(fi$alpha, 3)),
      Omega = ifelse(is.na(fi$omega), NA, round(fi$omega, 3)),
      Items = paste(fi$items, collapse = ", "),
      stringsAsFactors = FALSE
    )
  }))

  ft <- flextable(fiab_w) |>
    bold(part = "header") |>
    bg(part = "header", bg = "#3498DB") |>
    color(part = "header", color = "white") |>
    colformat_double(digits = 3) |>
    width(j = 1, width = 1) |>
    width(j = 5, width = 3) |>
    theme_box()
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")
  sec_num <- sec_num + 1

  # --- Veredicto ---
  doc <- body_add_par(doc,
    paste0(sec_num, ". Conclusion"), style = "heading 2")

  if (res$solucion_optima) {
    doc <- body_add_par(doc,
      "SOLUCION OPTIMA: La solucion factorial no requiere reespecificaciones adicionales. Todos los criterios de calidad se cumplen.",
      style = "Normal")
  } else {
    doc <- body_add_par(doc,
      "SUGERENCIAS DE REESPECIFICACION:", style = "Normal")
    for (s in res$sugerencias)
      doc <- body_add_par(doc, paste("  -", s), style = "Normal")
  }
  doc <- body_add_par(doc, "", style = "Normal")
  sec_num <- sec_num + 1

  # --- Referencias ---
  doc <- body_add_par(doc,
    paste0(sec_num, ". Referencias metodologicas"), style = "heading 2")

  refs_texto <- c(
    "Costello, A. B., & Osborne, J. W. (2005). Best practices in exploratory factor analysis. Practical Assessment, Research, and Evaluation, 10(7), 1-9.",
    "Dunn, T. J., Baguley, T., & Brunsden, V. (2014). From alpha to omega: A practical solution to the pervasive problem of internal consistency estimation. British Journal of Psychology, 105(3), 399-412.",
    "Ferrando, P. J., & Lorenzo-Seva, U. (2017). Program FACTOR at 10: Origins, development and future directions. Psicothema, 29(2), 236-240.",
    "Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019). Multivariate data analysis (8th ed.). Cengage.",
    "Holgado-Tello, F. P., Chacon-Moscoso, S., Barbero-Garcia, I., & Vila-Abad, E. (2010). Polychoric versus Pearson correlations in exploratory and confirmatory factor analysis of ordinal variables. Quality & Quantity, 44(1), 153-166.",
    "Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30(2), 179-185.",
    "Howard, M. C. (2016). A review of exploratory factor analysis decisions and overview of current practices. International Journal of Human-Computer Interaction, 32(1), 51-62.",
    "Kaiser, H. F. (1974). An index of factorial simplicity. Psychometrika, 39(1), 31-36.",
    "Lloret-Segura, S., Ferreres-Traver, A., Hernandez-Baeza, A., & Tomas-Marco, I. (2014). El analisis factorial exploratorio de los items: una guia practica, revisada y actualizada. Anales de Psicologia, 30(3), 1151-1169.",
    "Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. Biometrika, 57(3), 519-530.",
    "McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence Erlbaum Associates.",
    "Revelle, W. (2024). psych: Procedures for psychological, psychometric, and personality research. R package.",
    "Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. Psychological Methods, 16(2), 209-220.",
    "Watkins, M. W. (2018). Exploratory factor analysis: A guide to best practice. Journal of Black Psychology, 44(3), 219-246."
  )
  for (ref in refs_texto)
    doc <- body_add_par(doc, ref, style = "Normal")

  # Guardar
  if (!grepl("\\.docx$", archivo)) archivo <- paste0(archivo, ".docx")
  print(doc, target = archivo)
  cat(sprintf("  Archivo Word exportado: %s\n", normalizePath(archivo)))
}


# ============================================================
#        FUNCION PUBLICA DE EXPORTACION
# ============================================================

#' Exportar resultados del AFE a Excel o Word
#'
#' @param resultado Objeto de clase "efa_auto" retornado por efa_auto().
#' @param formato "excel" o "word" (puede ser un vector con ambos).
#' @param archivo Nombre del archivo (sin extension). Si es NULL, se genera
#'        automaticamente con fecha y hora.
#'
#' @return El objeto `resultado`, invisiblemente.
#'
#' @examples
#' # exportar_efa(resultado, formato = "excel")
#' # exportar_efa(resultado, formato = "word")
#' # exportar_efa(resultado, formato = c("excel", "word"))
#'
#' @export

exportar_efa <- function(resultado, formato = c("excel", "word"), archivo = NULL) {
  if (!inherits(resultado, "efa_auto"))
    stop("'resultado' debe ser un objeto de clase 'efa_auto'.")

  formato <- match.arg(formato, c("excel", "word"), several.ok = TRUE)

  if (is.null(archivo))
    archivo <- paste0("AFE_Resultados_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  for (fmt in formato) {
    if (fmt == "excel") {
      .efa_exportar_excel(resultado, archivo)
    } else {
      .efa_exportar_word(resultado, archivo)
    }
  }

  invisible(resultado)
}

