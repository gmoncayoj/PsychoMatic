#' Automated exploratory factor analysis
#'
#' Runs a compact exploratory factor-analysis workflow: item checks,
#' correlation selection, Mardia diagnostics, parallel analysis, factor fitting,
#' item screening, and reliability summaries. Screening can suggest changes or
#' remove one item at a time. Every mutation is followed by a complete refit
#' that reclassifies the current item set and resolves its correlation method.
#' Omega uses the same independent one-factor model per dimension as
#' [reliability_auto()], not the loadings of the exploratory solution.
#'
#' @param data Data frame, matrix, or `psychomatic_spec` containing item
#'   responses.
#' @param rotation Rotation family: `"oblique"` or `"orthogonal"`.
#' @param min_loading Minimum acceptable absolute factor loading.
#' @param min_communality Minimum acceptable communality.
#' @param min_cross_loading_diff Minimum difference between the two largest
#'   loadings used by the cross-loading screen.
#' @param max_iter Maximum number of item-removal steps.
#' @param verbose If `TRUE`, print results directly, including when assigned.
#' @param digits Decimal places in console and reports; stored values keep full precision.
#' @param report If `TRUE`, also save a reader-facing report.
#' @param format Report format: `"excel"`, `"word"`, or both.
#' @param overwrite If `TRUE`, allow replacement of an existing report.
#' @param export_format Optional export format passed to `export_efa()`.
#' @param file_name Optional base path for exported files. When `NULL`, reports
#'   use a unique timestamped name in the working directory (`getwd()`).
#' @param language Output language: `"esp"` or `"eng"`.
#' @param seed Optional non-negative integer seed for parallel analysis.
#' @param refine `"remove"` preserves the historical automated-removal
#'   workflow; `"suggest"` retains every item and reports screening flags.
#' @param parallel_iter Number of simulations used by parallel analysis.
#' @param parallel_quantile Reference quantile used by parallel analysis.
#' @param repair Correlation-matrix repair. `"none"` never changes the matrix;
#'   `"nearPD"` explicitly opts into nearest-positive-definite repair.
#' @param ... Deprecated Spanish argument aliases.
#'
#' @return An object of classes `efa_auto` and `psychomatic_result`. Legacy
#'   fields remain available during the 1.x compatibility cycle.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' efa_auto(MY_DATA[, 3:14])
#' }
#'
#' @export
efa_auto <- function(data,
                     rotation = c("oblique", "orthogonal"),
                     min_loading = 0.30,
                     min_communality = 0.30,
                     min_cross_loading_diff = 0.15,
                     max_iter = 20,
                     verbose = TRUE,
                     export_format = NULL,
                     file_name = NULL,
                     language = "eng",
                     seed = NULL,
                     refine = c("remove", "suggest"),
                     parallel_iter = 500L,
                     parallel_quantile = 0.95,
                     repair = c("none", "nearPD"),
                     digits = 3, report = FALSE, format = "excel",
                     overwrite = FALSE,
                     ...) {
  call <- match.call()
  supplied <- c(
    data = !missing(data), rotation = !missing(rotation),
    min_loading = !missing(min_loading),
    min_communality = !missing(min_communality),
    min_cross_loading_diff = !missing(min_cross_loading_diff),
    export_format = !missing(export_format), file_name = !missing(file_name)
  )
  values <- list(
    data = if (missing(data)) NULL else data,
    rotation = rotation,
    min_loading = min_loading,
    min_communality = min_communality,
    min_cross_loading_diff = min_cross_loading_diff,
    export_format = export_format,
    file_name = file_name
  )
  legacy <- list(...)
  if (length(legacy) && (is.null(names(legacy)) || any(!nzchar(names(legacy))))) {
    stop("All additional arguments passed to `efa_auto()` must be named.", call. = FALSE)
  }
  aliases <- c(
    datos = "data", rotacion = "rotation", carga_min = "min_loading",
    comunalidad_min = "min_communality",
    dif_cargas_cruzadas = "min_cross_loading_diff",
    exportar = "export_format", archivo = "file_name"
  )
  for (old in intersect(names(aliases), names(legacy))) {
    new <- aliases[[old]]
    if (isTRUE(supplied[[new]])) {
      stop(sprintf("Use only one of `%s` or deprecated `%s`.", new, old), call. = FALSE)
    }
    values[[new]] <- legacy[[old]]
    legacy[[old]] <- NULL
    .efa_warn_deprecated_arg(old, new)
  }
  if (length(legacy)) {
    stop(
      sprintf("Unused argument(s): %s", paste(names(legacy), collapse = ", ")),
      call. = FALSE
    )
  }
  data <- values$data
  rotation <- values$rotation
  min_loading <- values$min_loading
  min_communality <- values$min_communality
  min_cross_loading_diff <- values$min_cross_loading_diff
  export_format <- values$export_format
  file_name <- values$file_name
  if (is.null(data)) stop("Argument `data` is required.", call. = FALSE)
  if (!is.null(export_format)) { report <- TRUE; format <- export_format }
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)

  language <- .efa_normalize_language(language)
  rotation <- .efa_normalize_rotation(rotation)
  refine <- match.arg(refine)
  repair <- match.arg(repair)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  min_loading <- .pm_validate_scalar_number(min_loading, "min_loading", 0, 1)
  min_communality <- .pm_validate_scalar_number(min_communality, "min_communality", 0, 1)
  min_cross_loading_diff <- .pm_validate_scalar_number(
    min_cross_loading_diff, "min_cross_loading_diff", 0, 1
  )
  max_iter <- as.integer(.pm_validate_scalar_number(
    max_iter, "max_iter", lower = 0, integer = TRUE
  ))
  parallel_iter <- as.integer(.pm_validate_scalar_number(
    parallel_iter, "parallel_iter", lower = 1, integer = TRUE
  ))
  parallel_quantile <- .pm_validate_scalar_number(
    parallel_quantile, "parallel_quantile", 0.5, 0.999
  )
  if (!is.null(seed)) {
    seed <- as.integer(.pm_validate_scalar_number(
      seed, "seed", 0, .Machine$integer.max, integer = TRUE
    ))
  }
  if (!requireNamespace("psych", quietly = TRUE)) stop("Package 'psych' is required.", call. = FALSE)
  if (!requireNamespace("GPArotation", quietly = TRUE)) stop("Package 'GPArotation' is required.", call. = FALSE)

  if (inherits(data, "psychomatic_spec")) {
    specification <- data
    data <- specification$data[specification$item_index]
    names(data) <- specification$items
  } else {
    data <- .pm_as_data_frame(data)
  }
  numeric_columns <- vapply(data, is.numeric, logical(1))
  if (sum(numeric_columns) < 3L) stop("At least 3 numeric variables are required.", call. = FALSE)
  excluded_columns <- names(data)[!numeric_columns]
  data <- data[numeric_columns]
  validation <- check_items(data, language = language)
  fatal_codes <- c("all_missing", "zero_variance", "non_finite", "duplicate_name")
  fatal <- validation$issues[validation$issues$code %in% fatal_codes, , drop = FALSE]
  if (nrow(fatal)) {
    stop(
      paste("Item validation failed:", paste(fatal$message, collapse = " | ")),
      call. = FALSE
    )
  }

  n_original <- nrow(data)
  complete_rows <- stats::complete.cases(data)
  data <- data[complete_rows, , drop = FALSE]
  if (nrow(data) < 3L) stop("At least three complete cases are required.", call. = FALSE)
  n <- nrow(data)
  p_original <- ncol(data)

  mardia_capture <- .pm_capture_conditions(.efa_mardia(data))
  mardia_condition <- mardia_capture$error
  normality <- if (is.null(mardia_condition)) {
    mardia_capture$value
  } else {
    list(
      skewness = NA_real_, skewness_chi = NA_real_, skewness_df = NA_real_,
      skewness_p = NA_real_, kurtosis = NA_real_, kurtosis_z = NA_real_,
      kurtosis_p = NA_real_, normal = FALSE, status = "not_evaluable",
      reason = conditionMessage(mardia_condition), condition = mardia_condition,
      covariance_divisor = "n"
    )
  }
  initial_settings <- .efa_resolve_current_settings(
    data,
    normality_ok = normality$normal,
    language = language
  )
  categories <- initial_settings$categories
  item_levels <- initial_settings$item_levels
  measurement <- initial_settings$measurement
  estimator <- initial_settings$estimator

  fit_current <- function(current_data) {
    current_settings <- .efa_resolve_current_settings(
      current_data,
      normality_ok = normality$normal,
      language = language
    )
    state <- .efa_fit_current(
      current_data, current_settings$requested_method, repair, parallel_iter,
      parallel_quantile, seed, rotation, current_settings$estimator$fm
    )
    state$current_settings <- current_settings
    state
  }
  current_data <- data
  removed <- data.frame(
    Item = character(), Razon = character(), Iteracion = integer(),
    stringsAsFactors = FALSE
  )
  state <- fit_current(current_data)
  last_screen <- .efa_empty_screening()
  iterations <- 0L

  if (identical(state$status, "ok")) {
    repeat {
      last_screen <- .efa_evaluar_items(
        as.data.frame(unclass(state$fit$loadings)), state$fit$communality,
        min_loading, min_communality, min_cross_loading_diff,
        state$n_factors, language
      )
      if (identical(refine, "suggest") || nrow(last_screen) == 0L ||
          iterations >= max_iter) break

      iterations <- iterations + 1L
      removed <- rbind(removed, data.frame(
        Item = last_screen$Item[1L], Razon = last_screen$Razon[1L],
        Iteracion = iterations, stringsAsFactors = FALSE
      ))
      current_data <- current_data[
        , names(current_data) != last_screen$Item[1L], drop = FALSE
      ]
      # Refit immediately after every mutation, including the last allowed one.
      state <- fit_current(current_data)
      if (!identical(state$status, "ok")) {
        last_screen <- .efa_empty_screening()
        break
      }
    }
  }

  result <- .efa_build_result(
    state, current_data, data, removed, last_screen, validation, normality,
    mardia_condition, estimator, rotation, min_loading, n_original,
    excluded_columns, categories, item_levels, measurement, refine, repair,
    parallel_iter, parallel_quantile, seed, language, call
  )
  .pm_finish_analysis(result, output, language)
}

.efa_empty_screening <- function() {
  data.frame(
    Item = character(), Razon = character(), Prioridad = numeric(),
    stringsAsFactors = FALSE
  )
}

.efa_resolve_current_settings <- function(current_data,
                                          normality_ok,
                                          language = "eng") {
  categories <- vapply(
    current_data,
    function(x) length(unique(x)),
    integer(1)
  )
  item_levels <- vapply(current_data, .pm_measurement_level, character(1))
  measurement <- if (all(item_levels == "ordinal")) {
    "ordinal"
  } else if (all(item_levels == "continuous")) {
    "continuous"
  } else {
    "mixed_numeric"
  }
  ordinal <- identical(measurement, "ordinal")

  list(
    categories = categories,
    item_levels = item_levels,
    measurement = measurement,
    requested_method = if (ordinal) "polychoric" else "pearson",
    estimator = .efa_seleccionar_estimador(
      ordinal,
      normality_ok,
      nrow(current_data),
      ncol(current_data),
      language
    )
  )
}

.efa_fit_current <- function(data,
                             requested_method,
                             repair,
                             n_iter,
                             quantile,
                             seed,
                             rotation,
                             fm) {
  if (ncol(data) < 3L) {
    return(.efa_not_evaluable_state("Fewer than three items remain after screening."))
  }
  correlation <- tryCatch(
    .efa_estimate_correlations(data, requested_method, repair),
    error = identity
  )
  if (inherits(correlation, "condition")) {
    return(.efa_not_evaluable_state(conditionMessage(correlation)))
  }
  minimum_eigenvalue <- min(eigen(
    correlation$matrix, symmetric = TRUE, only.values = TRUE
  )$values)
  if (!is.finite(minimum_eigenvalue) ||
      minimum_eigenvalue <= sqrt(.Machine$double.eps)) {
    return(.efa_not_evaluable_state(
      paste0(
        "The correlation matrix is not positive definite (minimum eigenvalue ",
        format(minimum_eigenvalue, digits = 5),
        "). Use `repair = \"nearPD\"` only when substantively justified."
      ),
      correlation
    ))
  }
  parallel <- tryCatch(
    .efa_parallel_analysis(
      correlation$matrix, nrow(data), fm, n_iter, quantile, seed
    ),
    error = identity
  )
  if (inherits(parallel, "condition")) {
    return(.efa_not_evaluable_state(
      paste0("Parallel analysis failed: ", conditionMessage(parallel)),
      correlation
    ))
  }
  if (parallel$n_factors == 0L) {
    decisions <- rbind(
      correlation$decision_log,
      .pm_decision(
        "factor_retention", "parallel_analysis", "0",
        "No observed factor eigenvalue exceeded its simulated reference quantile."
      )
    )
    return(list(
      status = "not_evaluable", reason = "Parallel analysis retained zero common factors.",
      fit = NULL, n_factors = 0L, correlation = correlation,
      parallel = parallel, rotation = "none", decision_log = decisions
    ))
  }

  retained <- parallel$n_factors
  capacity <- max(1L, floor(ncol(data) / 3L))
  retained <- min(retained, capacity)
  effective_rotation <- if (retained == 1L) {
    "none"
  } else if (rotation == "oblique") {
    "oblimin"
  } else {
    "varimax"
  }
  fit <- tryCatch(
    suppressWarnings(psych::fa(
      correlation$matrix, nfactors = retained, n.obs = nrow(data),
      rotate = effective_rotation, fm = fm
    )),
    error = identity
  )
  if (inherits(fit, "condition")) {
    return(.efa_not_evaluable_state(
      paste0("Factor fitting failed: ", conditionMessage(fit)),
      correlation, parallel
    ))
  }
  list(
    status = "ok", reason = "", fit = fit, n_factors = as.integer(retained),
    correlation = correlation, parallel = parallel, rotation = effective_rotation,
    decision_log = rbind(
      correlation$decision_log,
      .pm_decision(
        "factor_retention", as.character(parallel$n_factors), as.character(retained),
        if (retained == parallel$n_factors) "Parallel-analysis count retained."
        else "Count reduced to preserve at least three items per fitted factor."
      ),
      .pm_decision(
        "rotation", rotation, effective_rotation,
        if (retained == 1L) "Rotation is unnecessary for one factor."
        else "Requested rotation family applied."
      )
    )
  )
}

.efa_not_evaluable_state <- function(reason, correlation = NULL, parallel = NULL) {
  decisions <- if (is.null(correlation)) {
    .pm_decision("factor_fit", "fit", "not_evaluable", reason, "warning")
  } else {
    rbind(
      correlation$decision_log,
      .pm_decision("factor_fit", "fit", "not_evaluable", reason, "warning")
    )
  }
  list(
    status = "not_evaluable", reason = reason, fit = NULL,
    n_factors = if (is.null(parallel)) NA_integer_ else parallel$n_factors,
    correlation = correlation, parallel = parallel, rotation = "none",
    decision_log = decisions
  )
}

.efa_safe_adequacy <- function(R, n) {
  empty <- list(
    kmo = list(MSA = NA_real_, MSAi = numeric()),
    bartlett = list(chisq = NA_real_, df = NA_real_, p.value = NA_real_)
  )
  if (is.null(R) || ncol(R) < 2L) return(empty)
  empty$kmo <- tryCatch(
    psych::KMO(R),
    error = function(e) empty$kmo
  )
  empty$bartlett <- tryCatch(
    psych::cortest.bartlett(R, n = n),
    error = function(e) empty$bartlett
  )
  empty
}

.efa_reliability <- function(data, distribution) {
  reliability <- list()
  if (!nrow(distribution)) return(reliability)
  for (factor in sort(unique(distribution$Dimension))) {
    items <- distribution$Item[distribution$Dimension == factor]
    if (length(items) >= 3L) {
      alpha_fit <- tryCatch({
        invisible(utils::capture.output(
          fitted_alpha <- suppressWarnings(suppressMessages(psych::alpha(
            data[items], check.keys = TRUE
          )))
        ))
        fitted_alpha
      }, error = function(e) NULL)
      prepared <- .pm_prepare_reliability_input(data[items], NULL, "auto", "eng")
      omega <- .pm_auto_omega(prepared$data, prepared$level, prepared$keying)
      reliability[[factor]] <- list(
        alpha = if (is.null(alpha_fit)) NA_real_ else alpha_fit$total$raw_alpha,
        omega = omega$row$estimate,
        omega_status = omega$row$status,
        omega_reason = omega$row$reason,
        omega_diagnostics = omega$diagnostics,
        items = items,
        n_items = length(items),
        alpha_drop = if (is.null(alpha_fit)) NULL else alpha_fit$alpha.drop
      )
    } else if (length(items) == 2L) {
      item_correlation <- stats::cor(
        data[[items[1L]]], data[[items[2L]]], use = "complete.obs"
      )
      reliability[[factor]] <- list(
        alpha = NA_real_, omega = NA_real_,
        spearman_brown = .pm_spearman_brown(item_correlation),
        items = items, n_items = 2L
      )
    } else {
      reliability[[factor]] <- list(
        alpha = NA_real_, omega = NA_real_, items = items,
        n_items = length(items)
      )
    }
  }
  reliability
}

.efa_build_result <- function(state,
                              current_data,
                              original_data,
                              removed,
                              screening,
                              validation,
                              normality,
                              mardia_condition,
                              estimator,
                              rotation_requested,
                              min_loading,
                              n_original,
                              excluded_columns,
                              categories,
                              item_levels,
                              measurement,
                              refine,
                              repair,
                              parallel_iter,
                              parallel_quantile,
                              seed,
                              language,
                              call) {
  n <- nrow(current_data)
  p_original <- ncol(original_data)
  final_settings <- state$current_settings
  final_measurement <- final_settings$measurement
  correlation_method <- if (isTRUE(estimator$usar_policorica)) "polychoric" else "pearson"
  original_correlation <- tryCatch(
    .efa_estimate_correlations(original_data, correlation_method, repair),
    error = function(e) NULL
  )
  original_adequacy <- .efa_safe_adequacy(
    if (is.null(original_correlation)) NULL else original_correlation$matrix,
    nrow(original_data)
  )

  if (identical(state$status, "ok")) {
    fit <- state$fit
    loadings <- as.data.frame(unclass(fit$loadings))
    assignment <- apply(abs(loadings), 1L, which.max)
    distribution <- data.frame(
      Item = rownames(loadings),
      Dimension = paste0("F", assignment),
      Carga_Principal = vapply(
        seq_along(assignment),
        function(i) loadings[i, assignment[i]],
        numeric(1)
      ),
      Comunalidad = fit$communality[rownames(loadings)],
      stringsAsFactors = FALSE
    )
    distribution <- distribution[
      order(distribution$Dimension, -abs(distribution$Carga_Principal)),
      , drop = FALSE
    ]
    rownames(distribution) <- NULL
    reliability <- .efa_reliability(current_data, distribution)
    adequacy <- .efa_safe_adequacy(state$correlation$matrix, n)
    variance <- fit$Vaccounted
    methodological <- .efa_sugerencias(
      fit, distribution, reliability, adequacy$kmo, adequacy$bartlett,
      state$n_factors, n, removed, min_loading, variance, p_original, language
    )
    chi_square <- if (is.null(fit$STATISTIC)) NA_real_ else fit$STATISTIC
    degrees_freedom <- if (is.null(fit$dof)) NA_real_ else fit$dof
    p_value <- if (is.null(fit$PVAL)) NA_real_ else fit$PVAL
  } else {
    fit <- NULL
    loadings <- data.frame(row.names = names(current_data))
    distribution <- data.frame(
      Item = character(), Dimension = character(), Carga_Principal = numeric(),
      Comunalidad = numeric(), stringsAsFactors = FALSE
    )
    reliability <- list()
    adequacy <- .efa_safe_adequacy(
      if (is.null(state$correlation)) NULL else state$correlation$matrix,
      n
    )
    variance <- NULL
    methodological <- state$reason
    chi_square <- degrees_freedom <- p_value <- NA_real_
  }
  suggestions <- unique(c(screening$Razon, methodological))
  suggestions <- suggestions[!is.na(suggestions) & nzchar(suggestions)]
  screening_passed <- identical(state$status, "ok") && !length(suggestions)
  conclusion <- .efa_conclusion(state, screening_passed, language)

  analysis_conditions <- if (is.null(state$correlation)) list() else state$correlation$warnings
  if (!is.null(mardia_condition)) analysis_conditions <- c(analysis_conditions, list(mardia_condition))
  measurement_reason <- switch(
    measurement,
    ordinal = "All items were classified as ordinal numeric codes (2-7 categories).",
    continuous = "All items were classified as continuous.",
    mixed_numeric = paste0(
      "Mixed ordinal/continuous numeric items were detected; Pearson correlations ",
      "were used instead of treating the full set as ordinal."
    )
  )
  final_measurement_reason <- switch(
    final_measurement,
    ordinal = "All final items were classified as ordinal numeric codes (2-7 categories).",
    continuous = "All final items were classified as continuous.",
    mixed_numeric = paste0(
      "The final item set contains mixed ordinal/continuous numeric items; Pearson ",
      "correlations were requested."
    )
  )
  measurement_transition <- if (identical(measurement, final_measurement)) {
    NULL
  } else {
    .pm_decision(
      "measurement_transition", measurement, final_measurement,
      "Measurement classification was recomputed after item screening."
    )
  }
  decisions <- rbind(
    .pm_decision(
      "missing_data", "unspecified", "complete_cases",
      sprintf("A fixed complete-case sample removed %d of %d rows before all iterations.",
              n_original - nrow(original_data), n_original)
    ),
    .pm_decision("measurement_classification", "auto", measurement, measurement_reason),
    .pm_decision(
      "measurement_classification_final", "auto", final_measurement,
      final_measurement_reason
    ),
    .pm_decision(
      "correlation_policy", "adaptive_by_current_items",
      "adaptive_by_current_items",
      "Measurement classification and requested correlations were resolved for each current item set."
    ),
    measurement_transition,
    state$decision_log
  )

  result <- list(
    call = call,
    status = state$status,
    reason = state$reason,
    requested_settings = list(
      rotation = rotation_requested, refine = refine, repair = repair,
      missing = "complete_cases",
      correlation_policy = "adaptive_by_current_items",
      parallel_analysis = list(
        n_iter = as.integer(parallel_iter), quantile = parallel_quantile,
        seed = seed, target = "factor_eigenvalues"
      )
    ),
    effective_settings = list(
      rotation = state$rotation, refine = refine,
      repair = if (!is.null(state$correlation) && state$correlation$repaired) "nearPD" else "none",
      correlation = if (is.null(state$correlation)) NA_character_ else state$correlation$effective_method,
      measurement = final_measurement,
      missing = "complete_cases",
      parallel_analysis = list(
        n_iter = as.integer(parallel_iter), quantile = parallel_quantile,
        seed = seed, target = "factor_eigenvalues"
      )
    ),
    data_summary = list(
      n_original = n_original, n_complete = nrow(original_data),
      p_original = p_original, p_final = ncol(current_data),
      excluded_columns = excluded_columns, categories = categories,
      item_levels = item_levels, measurement = final_measurement,
      measurement_initial = measurement,
      measurement_final = final_measurement,
      categories_final = final_settings$categories,
      item_levels_final = final_settings$item_levels
    ),
    results = list(
      fit = fit, loadings = loadings, distribution = distribution,
      reliability = reliability, screening = screening
    ),
    warnings = .pm_conditions_to_table(analysis_conditions, "warning"),
    errors = .pm_empty_condition_table(),
    diagnostics = list(
      item_validation = validation$issues, normality = normality,
      correlation = state$correlation, parallel_analysis = state$parallel
    ),
    decision_log = decisions,
    conclusion = conclusion,
    efa = fit,
    cargas = loadings,
    comunalidades = if (is.null(fit)) stats::setNames(numeric(), character()) else fit$communality,
    distribucion = distribution,
    language = language,
    seed = seed,
    n_factores = state$n_factors,
    estimador = final_settings$estimator,
    rotation = state$rotation,
    rotacion = state$rotation,
    kmo_original = original_adequacy$kmo,
    kmo_final = adequacy$kmo,
    bartlett_original = original_adequacy$bartlett,
    bartlett_final = adequacy$bartlett,
    chi_cuadrado = chi_square,
    grados_libertad = degrees_freedom,
    p_valor = p_value,
    varianza_explicada = variance,
    normalidad = normality,
    fiabilidad = reliability,
    items_eliminados = removed,
    correlacion_factores = if (is.null(fit) || is.null(fit$Phi)) NULL else fit$Phi,
    sugerencias = suggestions,
    screening_suggestions = screening,
    screening_passed = screening_passed,
    solucion_optima = screening_passed,
    final_data = current_data,
    datos_finales = current_data,
    correlation = state$correlation,
    parallel_analysis = state$parallel,
    n = n,
    p_original = p_original,
    p_final = ncol(current_data),
    compatibility = list(aliases = list(
      solucion_optima = list(
        replacement = "screening_passed", deprecated = TRUE,
        removal = "after the 1.x compatibility cycle"
      ),
      datos_finales = list(
        replacement = "final_data", deprecated = FALSE,
        removal = NA_character_
      )
    ))
  )
  result <- .efa_add_tidy_fields(result, call)
  class(result) <- c("efa_auto", "psychomatic_result")
  result
}

.efa_conclusion <- function(state, screening_passed, language) {
  if (identical(state$status, "not_evaluable")) {
    if (language == "esp") paste("El analisis factorial no fue evaluable:", state$reason)
    else paste("The factor analysis was not evaluable:", state$reason)
  } else if (screening_passed) {
    if (language == "esp") {
      "Se cumplieron los criterios automatizados de cribado; la interpretacion sustantiva sigue siendo necesaria."
    } else {
      "Automated screening criteria were met; substantive interpretation remains necessary."
    }
  } else if (language == "esp") {
    "El cribado genero observaciones que requieren revision sustantiva."
  } else {
    "Screening generated observations that require substantive review."
  }
}

.efa_legacy <- function(x, name = c("solucion_optima", "datos_finales")) {
  if (!inherits(x, "efa_auto")) stop("`x` must be an `efa_auto` result.", call. = FALSE)
  name <- match.arg(name)
  metadata <- x$compatibility$aliases[[name]]
  if (isTRUE(metadata$deprecated)) {
    warning(
      sprintf("`%s` is deprecated; use `%s` instead.", name, metadata$replacement),
      call. = FALSE
    )
  }
  x[[name]]
}

#' Print an automated EFA result
#'
#' @param x Object returned by `efa_auto()`.
#' @param digits Decimal places for display. `NULL` uses the analysis setting.
#' @param ... Additional unused arguments.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.efa_auto <- function(x, digits = NULL, ...) {
  .pm_console_report(x, digits)
}
