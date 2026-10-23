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

# Mardia statistics accumulated in row blocks. The largest temporary object is
# `block_size x n`; an n x n matrix is never retained.
.pm_mardia_block_ranges <- function(n, block_size) {
  effective_size <- min(as.integer(block_size), 256L, as.integer(n) - 1L)
  starts <- seq.int(1L, as.integer(n), by = effective_size)
  lapply(
    starts,
    function(start) seq.int(start, min(as.integer(n), start + effective_size - 1L))
  )
}

.pm_mardia_blockwise <- function(X, block_size = 256L) {
  X <- as.matrix(X)
  if (!is.numeric(X) || length(dim(X)) != 2L || nrow(X) < 2L || ncol(X) < 1L) {
    stop("`X` must be a numeric matrix with at least two rows.", call. = FALSE)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain complete finite values.", call. = FALSE)
  }
  block_size <- .pm_validate_scalar_number(
    block_size,
    name = "block_size",
    lower = 1,
    integer = TRUE
  )
  block_size <- as.integer(block_size)

  n <- nrow(X)
  p <- ncol(X)
  Xc <- sweep(X, 2L, colMeans(X), FUN = "-")
  S <- crossprod(Xc) / n
  S_inv <- tryCatch(
    solve(S),
    error = function(e) {
      condition <- structure(
        list(
          message = paste0(
            "Mardia statistics are not evaluable because the population ",
            "covariance matrix is singular: ", conditionMessage(e)
          ),
          call = sys.call(-1L),
          parent = e,
          code = "mardia_singular_covariance"
        ),
        class = c("psychomatic_mardia_singular", "error", "condition")
      )
      stop(condition)
    }
  )
  G <- Xc %*% S_inv

  cubic_sum <- 0
  block_ranges <- .pm_mardia_block_ranges(n, block_size)
  for (rows in block_ranges) {
    block <- G[rows, , drop = FALSE] %*% t(Xc)
    cubic_sum <- cubic_sum + sum(block^3)
  }

  mahalanobis_sq <- rowSums(G * Xc)
  b1p <- cubic_sum / n^2
  b2p <- sum(mahalanobis_sq^2) / n
  skew_df <- p * (p + 1) * (p + 2) / 6
  correction <- 1
  if (n < 20L) {
    denominator <- n * ((n + 1) * (p + 1) - 6)
    if (denominator > 0) {
      correction <- ((p + 1) * (n + 1) * (n + 3)) / denominator
    } else {
      correction <- NA_real_
    }
  }
  skew_statistic <- n * correction * b1p / 6
  skew_p <- stats::pchisq(skew_statistic, skew_df, lower.tail = FALSE)
  kurtosis_z <- (b2p - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  kurtosis_p <- 2 * stats::pnorm(abs(kurtosis_z), lower.tail = FALSE)

  c(
    b1p = b1p,
    skew_statistic = skew_statistic,
    skew_df = skew_df,
    skew_p = skew_p,
    b2p = b2p,
    kurtosis_z = kurtosis_z,
    kurtosis_p = kurtosis_p
  )
}

# Test de Mardia para normalidad multivariada
.efa_mardia <- function(X, block_size = 256L) {
  statistics <- .pm_mardia_blockwise(X, block_size = block_size)
  list(
    skewness = unname(statistics[["b1p"]]),
    skewness_chi = unname(statistics[["skew_statistic"]]),
    skewness_df = unname(statistics[["skew_df"]]),
    skewness_p = unname(statistics[["skew_p"]]),
    kurtosis = unname(statistics[["b2p"]]),
    kurtosis_z = unname(statistics[["kurtosis_z"]]),
    kurtosis_p = unname(statistics[["kurtosis_p"]]),
    normal = isTRUE(statistics[["skew_p"]] > 0.05) &&
      isTRUE(statistics[["kurtosis_p"]] > 0.05),
    block_size = as.integer(block_size),
    covariance_divisor = "n"
  )
}

.pm_parallel_count <- function(observed, simulated) {
  observed <- as.numeric(observed)
  simulated <- as.numeric(simulated)
  length_common <- min(length(observed), length(simulated))
  if (length_common == 0L) return(0L)
  keep <- seq_len(length_common)
  usable <- is.finite(observed[keep]) & is.finite(simulated[keep])
  comparisons <- usable & observed[keep] > simulated[keep]
  first_failure <- which(!comparisons)[1L]
  if (is.na(first_failure)) as.integer(length_common) else as.integer(first_failure - 1L)
}

.efa_with_seed <- function(seed, code) {
  if (is.null(seed)) return(eval.parent(substitute(code)))
  seed <- .pm_validate_scalar_number(
    seed,
    name = "seed",
    lower = 0,
    upper = .Machine$integer.max,
    integer = TRUE
  )
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  eval.parent(substitute(code))
}

.efa_parallel_analysis <- function(R,
                                   n_obs,
                                   fm = "minres",
                                   n_iter = 500L,
                                   quantile = 0.95,
                                   seed = NULL) {
  n_iter <- .pm_validate_scalar_number(
    n_iter,
    name = "parallel_iter",
    lower = 1,
    integer = TRUE
  )
  quantile <- .pm_validate_scalar_number(
    quantile,
    name = "parallel_quantile",
    lower = 0.5,
    upper = 0.999
  )
  n_obs <- .pm_validate_scalar_number(n_obs, name = "n_obs", lower = 2, integer = TRUE)
  R <- as.matrix(R)
  if (!is.numeric(R) || nrow(R) != ncol(R) || ncol(R) < 2L ||
      anyNA(R) || any(!is.finite(R))) {
    stop("`R` must be a complete finite correlation matrix.", call. = FALSE)
  }

  observed_reduced <- R
  diag(observed_reduced) <- .efa_smc_quiet(observed_reduced)
  observed <- eigen(observed_reduced, symmetric = TRUE, only.values = TRUE)$values
  p <- ncol(R)
  factor_simulations <- .efa_with_seed(seed, {
    simulations <- matrix(NA_real_, nrow = as.integer(n_iter), ncol = p)
    for (iteration in seq_len(as.integer(n_iter))) {
      simulated_data <- matrix(stats::rnorm(as.integer(n_obs) * p), ncol = p)
      simulated_correlation <- stats::cor(simulated_data)
      diag(simulated_correlation) <- .efa_smc_quiet(simulated_correlation)
      simulations[iteration, ] <- eigen(
        simulated_correlation,
        symmetric = TRUE,
        only.values = TRUE
      )$values
    }
    simulations
  })
  simulated_quantiles <- apply(
    factor_simulations,
    2L,
    stats::quantile,
    probs = quantile,
    na.rm = TRUE,
    names = FALSE
  )
  list(
    n_factors = .pm_parallel_count(observed, simulated_quantiles),
    observed = observed,
    simulated = simulated_quantiles,
    simulated_eigenvalues = factor_simulations,
    n_iter = as.integer(n_iter),
    quantile = quantile,
    seed = if (is.null(seed)) NULL else as.integer(seed),
    warnings = list(),
    engine = "PsychoMatic factor-eigenvalue simulation with psych::smc",
    extraction = fm,
    target = "factor_eigenvalues"
  )
}

.efa_smc_quiet <- function(R) {
  invisible(utils::capture.output(
    values <- suppressMessages(psych::smc(R))
  ))
  values
}

.efa_numeric_ordered_codes <- function(data) {
  is.data.frame(data) && ncol(data) > 0L && all(vapply(
    data,
    function(x) is.numeric(x) && identical(.pm_measurement_level(x), "ordinal"),
    logical(1)
  ))
}

.efa_estimate_correlations <- function(data,
                                       method = c("pearson", "spearman", "polychoric", "tetrachoric"),
                                       repair = c("none", "nearPD"),
                                       polychoric_fun = psych::polychoric,
                                       tetrachoric_fun = psych::tetrachoric) {
  method <- match.arg(method)
  repair <- match.arg(repair)
  data <- .pm_as_data_frame(data)
  if (ncol(data) < 2L) {
    stop("At least two items are required to estimate correlations.", call. = FALSE)
  }

  type <- switch(
    method,
    pearson = "pearson",
    spearman = "poly",
    polychoric = "poly",
    tetrachoric = "tetra"
  )
  prepared <- .cormat_prepare_data(data, type)
  validation <- check_items(prepared, language = "eng")
  .cormat_validate_observed_variation(prepared, validation$issues)
  pairwise_n <- validation$diagnostics$pairwise_n
  warnings <- list()
  decisions <- .pm_decision(
    "correlation_method",
    method,
    method,
    "The requested correlation method was attempted."
  )

  captured <- .pm_capture_conditions(switch(
      method,
      pearson = stats::cor(prepared, use = "pairwise.complete.obs", method = "pearson"),
      spearman = stats::cor(prepared, use = "pairwise.complete.obs", method = "spearman"),
      polychoric = .cormat_extract_rho(
        polychoric_fun(prepared, smooth = FALSE, progress = FALSE)
      ),
      tetrachoric = .cormat_extract_rho(
        tetrachoric_fun(prepared, smooth = FALSE)
      )
    ))
  warnings <- captured$warnings
  estimate <- if (is.null(captured$error)) captured$value else captured$error

  effective_method <- method
  if (inherits(estimate, "condition")) {
    original_condition <- estimate
    if (!method %in% c("polychoric", "tetrachoric")) stop(original_condition)
    if (!.efa_numeric_ordered_codes(data)) {
      stop(
        paste0(
          "Correlation estimation failed and Pearson fallback is allowed only ",
          "for numeric ordered codes. Original error: ", conditionMessage(original_condition)
        ),
        call. = FALSE
      )
    }
    warnings <- c(list(original_condition), warnings)
    estimate <- stats::cor(prepared, use = "pairwise.complete.obs", method = "pearson")
    effective_method <- "pearson"
    decisions <- rbind(
      decisions,
      .pm_decision(
        "correlation_fallback",
        method,
        "pearson",
        paste0("Original condition retained: ", conditionMessage(original_condition)),
        "warning"
      )
    )
  }

  estimate <- as.matrix(estimate)
  dimnames(estimate) <- list(names(prepared), names(prepared))
  eigenvalues <- eigen(estimate, symmetric = TRUE, only.values = TRUE)$values
  pre_repair_min <- min(eigenvalues)
  tolerance <- sqrt(.Machine$double.eps)
  needs_repair <- !is.finite(pre_repair_min) || pre_repair_min <= tolerance
  repaired <- FALSE
  if (needs_repair && identical(repair, "nearPD")) {
    estimate <- as.matrix(Matrix::nearPD(estimate, corr = TRUE)$mat)
    dimnames(estimate) <- list(names(prepared), names(prepared))
    repaired <- TRUE
    decisions <- rbind(
      decisions,
      .pm_decision(
        "correlation_repair",
        "nearPD",
        "nearPD",
        sprintf("Pre-repair minimum eigenvalue: %.8g.", pre_repair_min),
        "warning"
      )
    )
  } else if (needs_repair) {
    decisions <- rbind(
      decisions,
      .pm_decision(
        "correlation_repair",
        "none",
        "none",
        sprintf("Matrix retained without repair; minimum eigenvalue: %.8g.", pre_repair_min),
        "warning"
      )
    )
  }

  list(
    matrix = estimate,
    requested_method = method,
    effective_method = effective_method,
    repaired = repaired,
    repair_requested = repair,
    pre_repair_min_eigenvalue = pre_repair_min,
    warnings = warnings,
    pairwise_n = pairwise_n,
    diagnostics = validation$issues,
    decision_log = decisions
  )
}

# Normalizar idioma del reporte
.efa_normalize_language <- function(language = "eng") {
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
.efa_labels <- function(language = "eng") {
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
      analysis_title = "          AUTOMATED EXPLORATORY FACTOR ANALYSIS",
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
      optimal_factors = "  Parallel-analysis factor count:",
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
      optimal_solution_title = "\n  >>> SCREENING CRITERIA PASSED <<<\n\n",
      optimal_solution_1 = "  The automated screen generated no additional\n",
      optimal_solution_2 = "  flags; substantive interpretation remains necessary.\n\n",
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
    analysis_title = "          ANALISIS FACTORIAL EXPLORATORIO AUTOMATIZADO",
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
    optimal_factors = "  Factores segun analisis paralelo:",
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
    optimal_solution_title = "\n  >>> CRITERIOS DE CRIBADO CUMPLIDOS <<<\n\n",
    optimal_solution_1 = "  El cribado automatizado no genero alertas adicionales;\n",
    optimal_solution_2 = "  sigue siendo necesaria la interpretacion sustantiva.\n\n",
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
.efa_interpretar_kmo <- function(kmo, language = "eng") {
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
                                      language = "eng") {
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
                               n_factores, language = "eng") {
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
                             language = "eng") {
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


