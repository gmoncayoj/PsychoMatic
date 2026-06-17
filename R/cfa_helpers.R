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


