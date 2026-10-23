# Reader-facing presentation shared by console and reports. Never changes raw data.
.pm_output_options <- function(report = FALSE, format = "excel", file_name = NULL,
                               overwrite = FALSE, digits = 3, verbose = TRUE) {
  digits <- .pm_validate_scalar_number(digits, "digits", 0, 15, integer = TRUE)
  for (name in c("overwrite", "verbose")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    }
  }
  if (is.character(report) && length(report) == 1L && !is.na(report)) {
    if (!is.null(file_name)) stop("Use `file_name` or a path in `report`, not both.", call. = FALSE)
    file_name <- report
    report <- TRUE
  }
  if (!is.logical(report) || length(report) != 1L || is.na(report)) {
    stop("`report` must be TRUE, FALSE, or one output path.", call. = FALSE)
  }
  format <- match.arg(format, c("excel", "word"), several.ok = TRUE)
  if (anyDuplicated(format)) stop("Duplicate report formats.", call. = FALSE)
  paths <- NULL
  if (report) {
    .pm_check_report_dependencies(format)
    if (is.null(file_name)) file_name <- .pm_report_path(NULL, format[[1L]])
    paths <- vapply(format, function(fmt) {
      .pm_safe_output_path(.pm_report_path(file_name, fmt),
                           if (fmt == "excel") ".xlsx" else ".docx", overwrite)
    }, character(1))
  }
  list(report = report, format = format, paths = paths, overwrite = overwrite,
       digits = digits, verbose = verbose)
}

.pm_finish_analysis <- function(x, output, language = NULL) {
  attr(x, "display_digits") <- output$digits
  if (!is.null(language)) attr(x, "display_language") <- language
  if (output$report) {
    sections <- list(.pm_analysis_sections(x, language, destination = "report"))
    paths <- lapply(seq_along(output$format), function(i) {
      if (output$format[i] == "excel") {
        .pm_write_analysis_excel(sections, output$paths[i], output$digits, output$overwrite)
      } else {
        .pm_write_analysis_word(sections, output$paths[i], output$digits, output$overwrite)
      }
    })
    names(paths) <- output$format
    attr(x, "report_paths") <- vapply(paths, as.character, character(1))
    if (is.list(x) && !is.data.frame(x)) x$report_path <- unname(attr(x, "report_paths"))
    message(if (identical(language, "esp")) "Reporte guardado: " else "Report saved: ",
            paste(attr(x, "report_paths"), collapse = "; "))
  }
  if (output$verbose) .pm_console_report(x, output$digits, language)
  invisible(x)
}

.pm_matrix_table <- function(x, label = "Item") {
  if (is.null(x) || !length(x)) return(NULL)
  x <- as.data.frame(x, check.names = FALSE)
  labels <- rownames(x)
  if (is.null(labels)) labels <- as.character(seq_len(nrow(x)))
  out <- data.frame(labels, x, check.names = FALSE)
  names(out)[1] <- label
  rownames(out) <- NULL
  out
}

.pm_semtools_reference <- function() {
  year <- if (requireNamespace("semTools", quietly = TRUE)) {
    tryCatch(suppressWarnings(as.character(utils::citation("semTools")[[1L]]$year)),
      error = function(e) character())
  } else character()
  if (length(year) != 1L || is.na(year) || !nzchar(year)) year <- "n.d."
  version <- if (requireNamespace("semTools", quietly = TRUE)) paste0(" (Version ", utils::packageVersion("semTools"), ")") else ""
  paste0("Jorgensen, T. D., Pornprasertmanit, S., Schoemann, A. M., & Rosseel, Y. (", year,
         "). semTools: Useful tools for structural equation modeling", version,
         " [R package]. https://CRAN.R-project.org/package=semTools")
}

.pm_reference_text <- function(x) {
  if (is.null(x)) return(character())
  if (is.data.frame(x)) {
    key <- intersect(c("reference", "Reference", "Referencia"), names(x))
    if (!length(key)) return(character())
    x <- x[[key[1]]]
  }
  x <- as.character(x)
  sort(unique(x[!is.na(x) & nzchar(trimws(x))]))
}

# Suppress only redundant context columns, never the corresponding raw fields.
.pm_cfa_context <- function(table) {
  if (!is.data.frame(table)) return(table)
  for (key in intersect(c("group", "level"), names(table))) {
    values <- unique(table[[key]][!is.na(table[[key]])])
    if (length(values) <= 1L) table[[key]] <- NULL
  }
  if ("level" %in% names(table)) names(table)[names(table) == "level"] <- "FactorLevel"
  table
}

.pm_cfa_conclusion <- function(x, language) {
  text <- if (is.list(x$conclusion)) x$conclusion$reason else x$conclusion
  if (!identical(language, "esp") || !length(text) || all(is.na(text))) return(text)
  codes <- x$diagnostics$admissibility$reason_codes
  labels <- c(not_converged = "el modelo no convergio",
    post_check_failed = "no supero la comprobacion de admisibilidad",
    negative_variance = "se encontro una varianza negativa",
    heywood_loading = "se encontro una carga estandarizada fuera del rango admisible")
  if (length(codes)) {
    reasons <- unname(labels[codes])
    reasons[is.na(reasons)] <- "otra comprobacion diagnostica fallo"
    return(paste0("La solucion CFA no es interpretable: ", paste(reasons, collapse = "; "), "."))
  }
  if (identical(text, "The CFA model did not converge.")) return("El modelo CFA no convergio.")
  text
}

# Reader interpretation is separate from computational status in the raw result.
.pm_coefficient_evaluation <- function(table, value, metric, language) {
  if (!nrow(table)) return(table)
  es <- identical(language, "esp")
  cutoff <- ifelse(table[[metric]] == "ave", .50, .70)
  usable <- !is.na(table$status) & table$status == "ok" & is.finite(table[[value]])
  adequate <- table[[value]] >= cutoff
  table$status[usable] <- ifelse(adequate[usable],
    if (es) "Adecuado" else "Adequate", if (es) "No adecuado" else "Inadequate")
  table$reason[usable] <- paste(ifelse(adequate[usable], ">=", "<"),
    formatC(cutoff[usable], format="f", digits=2))
  empty <- !usable & (is.na(table$reason) | !nzchar(table$reason))
  table$reason[empty] <- if (es) "No se obtuvo una estimacion interpretable." else "No interpretable estimate was obtained."
  table
}

.pm_analysis_sections <- function(x, language = NULL, destination = "console") {
  if (inherits(x, "esem_auto")) return(.pm_esem_sections(x, language, destination))
  if (inherits(x, "psychomatic_network")) return(.pm_network_sections(x, language, destination))
  if (is.null(language) && is.list(x) && !is.data.frame(x)) {
    language <- .pm_null_default(x$language, x$effective_settings$language)
  }
  if (is.null(language)) language <- attr(x, "display_language", exact = TRUE)
  language <- .pm_normalize_language(language)
  es <- identical(language, "esp")
  tr <- function(esp, eng) if (es) esp else eng
  sections <- list()
  add <- function(title, value) {
    if (is.null(value) || !length(value)) return(invisible(NULL))
    if (is.character(value) && !is.data.frame(value)) {
      value <- value[!is.na(value) & nzchar(value)]
      if (!length(value)) return(invisible(NULL))
      value <- data.frame(Text = value, stringsAsFactors = FALSE)
    }
    if (is.matrix(value)) value <- .pm_matrix_table(value)
    if (is.data.frame(value) && nrow(value)) sections[[title]] <<- value
  }
  pick <- function(value, columns) {
    if (!is.data.frame(value)) return(NULL)
    value[, intersect(columns, names(value)), drop = FALSE]
  }
  meta <- function(...) {
    values <- list(...)
    data.frame(Information = names(values), Value = vapply(values, function(v) {
      if (is.null(v) || !length(v)) return("-")
      paste(v, collapse = ", ")
    }, character(1)), check.names = FALSE)
  }
  title <- tr("Resultados", "Results")
  refs <- if (is.list(x) && !is.data.frame(x)) x$references else NULL
  if (inherits(x, "efa_auto")) {
    title <- tr("Analisis factorial exploratorio", "Exploratory factor analysis")
    add(tr("Datos y decisiones", "Data and decisions"), meta(
      N = x$n, Items = paste(x$p_original, "->", x$p_final),
      Factors = x$n_factores, Correlation = x$effective_settings$correlation,
      Estimator = x$estimador$nombre, Rotation = x$effective_settings$rotation,
      Missing = x$effective_settings$missing,
      Excluded = x$data_summary$n_original - x$n))
    add(tr("Justificacion", "Rationale"), x$estimador$justificacion)
    norm <- x$normalidad
    add(tr("Normalidad multivariada Mardia", "Mardia multivariate normality"), data.frame(
      Test = c("Skewness", "Kurtosis"), Statistic = c(norm$skewness_chi, norm$kurtosis_z),
      df = c(norm$skewness_df, NA_real_), p = c(norm$skewness_p, norm$kurtosis_p)))
    add(tr("Adecuacion muestral", "Sampling adequacy"), data.frame(
      Sample = c(tr("Inicial", "Initial"), tr("Final", "Final")),
      KMO = c(x$kmo_original$MSA, x$kmo_final$MSA),
      Bartlett = c(x$bartlett_original$chisq, x$bartlett_final$chisq),
      df = c(x$bartlett_original$df, x$bartlett_final$df),
      p = c(x$bartlett_original$p.value, x$bartlett_final$p.value)))
    pa <- x$effective_settings$parallel_analysis
    add(tr("Analisis paralelo", "Parallel analysis"), meta(
      Factors = x$n_factores, Iterations = pa$n_iter, Quantile = pa$quantile, Seed = pa$seed))
    add(tr("Criterio del analisis paralelo", "Parallel-analysis criterion"), tr(
      "Se comparan los autovalores observados con el percentil de referencia de los datos aleatorios. La semilla permite reproducir la simulacion; '-' indica que no se fijo desde la funcion.",
      "Observed eigenvalues are compared with the reference percentile from random data. The seed makes the simulation reproducible; '-' means no seed was set by the function."))
    add(tr("Indices de ajuste", "Fit indices"), x$fit_table)
    variance <- x$varianza_explicada
    if (is.matrix(variance) && "Proportion Var" %in% rownames(variance)) {
      variance_table <- data.frame(Factor = colnames(variance),
        TotalPercent = unname(100 * variance["Proportion Var", ]),
        CumulativePercent = unname(100 * cumsum(variance["Proportion Var", ])))
      add(tr("Varianza explicada", "Explained variance"), variance_table)
      add(tr("Interpretacion de la varianza", "Variance interpretation"), tr(
        "Los porcentajes usan como denominador la varianza total de los items (respuestas latentes si se usan policoricas). En rotacion oblicua, la contribucion por factor no es varianza exclusiva: los factores estan correlacionados.",
        "Percentages use total item variance as the denominator (latent responses for polychoric correlations). With oblique rotation, per-factor contributions are not unique variance: factors are correlated."))
    }
    loads <- .pm_matrix_table(x$loadings)
    if (!is.null(loads)) loads$h2 <- as.numeric(x$comunalidades)
    add(tr("Cargas factoriales", "Loadings"), loads)
    add(tr("Distribucion de items", "Item allocation"), x$distribucion)
    if (ncol(x$loadings)) add(tr("Identificadores de factores", "Factor identifiers"),
      paste(paste0("F", seq_len(ncol(x$loadings)), " = ", colnames(x$loadings)), collapse = "; "))
    add(tr("Correlaciones entre factores", "Factor correlations"), .pm_matrix_table(x$correlacion_factores, "Factor"))
    add(tr("Confiabilidad por dimension", "Reliability by factor"), x$reliability_table)
    if (any(x$reliability_table$n_items >= 3L)) {
      refs <- unique(c(refs, "Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using structural equation modeling: An alternative to coefficient alpha. Psychometrika, 74(1), 155-167. https://doi.org/10.1007/s11336-008-9099-3"))
      add(tr("Nota sobre omega", "Omega note"), tr(
        "Omega se estima con un modelo unidimensional independiente por dimension, como en reliability_auto (ML para items continuos; WLSMV para ordinales, en la escala de respuesta observada). No utiliza directamente las cargas del AFE ni demuestra unidimensionalidad o buen ajuste. La coincidencia requiere los mismos items y participantes.",
        "Omega uses an independent one-factor model per dimension, as in reliability_auto (ML for continuous items; WLSMV for ordinal items, on the observed response scale). It does not directly use EFA loadings or establish unidimensionality or good fit. Agreement requires the same items and participants."))
    }
    add(tr("Items eliminados", "Removed items"), x$items_eliminados)
    add(tr("Observaciones sobre items", "Item screening"), x$screening_suggestions)
    add(tr("Observaciones metodologicas", "Methodological observations"),
        setdiff(x$sugerencias, x$screening_suggestions$Razon))
  } else if (inherits(x, "cfa_auto")) {
    title <- tr("Analisis factorial confirmatorio", "Confirmatory factor analysis")
    add(tr("Datos y metodo", "Data and method"), meta(
      N = x$data_summary$n, Items = x$data_summary$p, Estimator = x$estimator,
      Model = x$model_type, Converged = x$converged))
    add(tr("Justificacion", "Rationale"), c(x$estimator_reason, x$fi_variant_note))
    if (length(x$mardia$p_skew) && any(is.finite(c(x$mardia$p_skew, x$mardia$p_kurt)))) {
      p <- c(x$mardia$p_skew, x$mardia$p_kurt)
      alpha <- .pm_null_default(x$effective_settings$alpha_norm, .05)
      add(tr("Normalidad multivariada Mardia", "Mardia multivariate normality"), data.frame(
        Test = c("Skewness", "Kurtosis"),
        Statistic = c(x$mardia$skew_stat, x$mardia$kurt_stat),
        N = .pm_null_default(x$mardia$n, x$data_summary$n), p = p,
        Decision = ifelse(is.na(p), "-", ifelse(p < alpha,
          tr("[Rechazada]", "[Rejected]"), tr("[No rechazada]", "[Not rejected]")))))
    } else {
      why <- .pm_null_default(x$mardia$reason, "unavailable")
      detail <- switch(why,
        mvn_not_installed = tr("El paquete MVN no esta instalado.", "Package MVN is not installed."),
        mvn_error = tr("El calculo fallo; revise los errores del analisis.", "Computation failed; inspect the analysis errors."),
        mvn_output_unrecognized = tr("La salida de MVN no pudo interpretarse.", "MVN output could not be interpreted."),
        tr("Datos o resultados insuficientes para evaluar normalidad.", "Insufficient data or results to evaluate normality."))
      add(tr("Normalidad multivariada Mardia", "Mardia multivariate normality"),
        paste(tr("No evaluable:", "Not evaluable:"), detail))
    }
    if (isTRUE(x$mardia$ordinal_codes)) add(tr("Nota sobre Mardia", "Mardia note"), tr(
      "Mardia se calcula sobre los codigos observados. No prueba la normalidad de las respuestas latentes ordinales ni modifica el estimador seleccionado.",
      "Mardia is computed on observed numeric codes. It does not test normality of ordinal latent responses or change the selected estimator."))
    add(tr("Indices de ajuste", "Fit indices"), pick(x$fit_indices, c("Index", "Value", "Variant", "Evaluation")))
    add(tr("Cargas factoriales estandarizadas (IC 95%)", "Standardized loadings (95% CI)"), .pm_cfa_context(pick(x$loadings_table,
      c("group", "level", "Factor", "Item", "Loading", "SE", "p_value", "CI_lower", "CI_upper"))))
    if (length(x$r_squared)) {
      r2 <- unlist(x$r_squared)
      add("R2", data.frame(Item = names(r2), R2 = as.numeric(r2)))
    }
    add(tr("Correlaciones entre factores", "Factor correlations"), pick(x$factor_correlation_table,
      c("group", "Factor1", "Factor2", "Correlation", "SE", "p_value", "Alert")))
    rel <- .cfa_reliability_table(x)
    if (nrow(rel)) {
      rel$group <- .pm_group_label(rel$group)
      rel <- .pm_cfa_context(rel)
      is_ave <- rel$metric == "ave"
      is_structure <- rel$metric %in% c("ecv", "puc", "h_general", "h")
      add(tr("Confiabilidad", "Reliability"), .pm_coefficient_evaluation(rel[!is_ave & !is_structure, , drop = FALSE], "value", "metric", language))
      add(tr("Varianza media extraida", "Average variance extracted"), .pm_coefficient_evaluation(rel[is_ave, , drop = FALSE], "value", "metric", language))
      add(tr("Indices estructurales bifactor", "Bifactor structural indices"), rel[is_structure, , drop = FALSE])
    }
    suggestions <- .cfa_respecifications_to_df(x)
    if (nrow(suggestions) && identical(destination, "report")) {
      add(tr("Sugerencias para revisar", "Suggestions for review"), suggestions)
    } else if (nrow(suggestions)) {
      notes <- suggestions$Note
      suggestions$Note <- NULL
      suggestions <- data.frame(No = seq_len(nrow(suggestions)), suggestions, check.names=FALSE)
      add(tr("Sugerencias para revisar", "Suggestions for review"), suggestions)
      add(tr("Notas", "Notes"), paste0(seq_along(notes), ". ", notes))
    }
    vc <- x$diagnostics$vcov
    if (is.list(vc) && is.finite(vc$minimum_eigenvalue) && vc$minimum_eigenvalue <= vc$tolerance) {
      add(tr("Diagnostico de la matriz de parametros", "Parameter covariance diagnostic"), c(
        sprintf(tr("Autovalor minimo de vcov: %.4e; tolerancia numerica: %.4e.",
                   "Minimum vcov eigenvalue: %.4e; numerical tolerance: %.4e."), vc$minimum_eigenvalue, vc$tolerance),
        tr("Un valor cercano a cero no demuestra por si solo falta de identificacion. Revise convergencia, varianzas negativas, restricciones redundantes y errores estandar; no ignore las advertencias originales.",
           "A near-zero value alone does not prove underidentification. Review convergence, negative variances, redundant constraints and standard errors; do not ignore the original warnings.")))
    }
  } else if (inherits(x, "psychomatic_reliability")) {
    title <- tr("Confiabilidad", "Reliability")
    add(tr("Datos y metodo", "Data and method"), meta(
      N = x$data_summary$n, Items = x$data_summary$p, Complete = x$data_summary$complete,
      Measurement = x$effective_settings$level, Keying = x$effective_settings$keying))
    coefficients <- .pm_null_default(x$coefficients, x$results$coefficients)
    rows <- pick(coefficients, c("group", "factor", "coefficient", "estimate", "status", "reason"))
    global <- is.na(rows$factor) & is.na(rows$group)
    rows$group <- .pm_group_label(rows$group)
    rows$factor[is.na(rows$factor) | !nzchar(rows$factor)] <- tr("Todos los items seleccionados", "All selected items")
    rows$group[global] <- tr("Muestra de datos", "Data sample")
    ave <- rows$coefficient == "ave"
    rows <- .pm_coefficient_evaluation(rows, "estimate", "coefficient", language)
    labels <- c(alpha_observed = tr("Alfa (Cronbach)", "Alpha (Cronbach)"),
                alpha_ordinal = tr("Alfa ordinal", "Ordinal alpha"),
                spearman_brown = "Spearman-Brown", omega = "Omega",
                composite_reliability = tr("Omega (modelo suministrado)", "Omega (supplied model)"),
                ave = "AVE")
    matched <- match(rows$coefficient, names(labels))
    rows$coefficient[!is.na(matched)] <- unname(labels[matched[!is.na(matched)]])
    rows$reason[is.na(rows$reason)] <- ""
    add(tr("Coeficientes de confiabilidad", "Reliability coefficients"), rows[!ave, , drop = FALSE])
    add(tr("Varianza media extraida", "Average variance extracted"), rows[ave, , drop = FALSE])
    item_stats <- x$item_stats
    if (is.data.frame(item_stats) && nrow(item_stats)) {
      fields <- c("item", "n", "item_rest_r", "alpha_if_deleted",
                  if (identical(x$effective_settings$level, "ordinal")) "alpha_ordinal_if_deleted",
                  "omega_if_deleted")
      item_table <- pick(item_stats, fields)
      names(item_table)[names(item_table) == "n"] <- tr("N casos", "N cases")
      add(tr("Analisis de items: confiabilidad al eliminar", "Item analysis: reliability if deleted"), item_table)
      reasons <- intersect(c("item_rest_reason", "alpha_if_deleted_reason", "alpha_ordinal_if_deleted_reason", "omega_if_deleted_reason"), names(item_stats))
      if (!identical(x$effective_settings$level, "ordinal")) reasons <- setdiff(reasons, "alpha_ordinal_if_deleted_reason")
      failures <- lapply(reasons, function(field) {
        keep <- !is.na(item_stats[[field]]) & nzchar(item_stats[[field]])
        data.frame(Item = item_stats$item[keep], Index = rep(sub("_reason$", "", field), sum(keep)), reason = item_stats[[field]][keep])
      })
      add(tr("Estadisticos de items no evaluables", "Non-evaluable item statistics"), do.call(rbind, failures))
    }
    add(tr("Como se calcula e interpreta", "Computation and interpretation"), if (es) c(
      "Alfa (Cronbach): k/(k-1) * (1 - suma de varianzas / varianza de la suma). Usa covarianzas de los items seleccionados y observaciones disponibles por par.",
      if (any(x$coefficients$coefficient == "alpha_ordinal")) "Alfa ordinal: aplica la formula de alfa a la matriz policorica. Estima consistencia en las respuestas latentes continuas, no en la puntuacion observada.",
      if (any(x$coefficients$coefficient == "spearman_brown")) "Spearman-Brown: 2*r/(1+r) para dos items. La correlacion utilizada se registra en los resultados.",
      if (any(x$coefficients$coefficient == "composite_reliability")) "Confiabilidad compuesta: semTools::compRelSEM(model), por factor y grupo del CFA suministrado. Usa su muestra y parametros ajustados; por defecto obs.var=TRUE, tau.eq=FALSE y ord.scale=TRUE.",
      if (any(x$coefficients$coefficient == "omega")) "Omega: modelo congenerico de un factor (ML para datos continuos; WLSMV para ordinales), sin imponer igualdad de cargas. En ordinales estima confiabilidad de la puntuacion observada mediante la correccion de Green y Yang. No demuestra unidimensionalidad ni garantiza buen ajuste.",
      if (is.data.frame(item_stats)) "Correlacion item-test corregida: correlacion del item con la suma de los demas, usando casos completos de la escala. Alfa usa pares disponibles; omega usa casos completos de cada ajuste. Al eliminar un item se reajusta omega y su N puede cambiar. Una mejora del coeficiente no justifica por si sola eliminarlo.",
      if (any(ave)) "AVE: semTools::AVE(model). Es varianza media extraida; no es un coeficiente de confiabilidad.",
      "Los alfas del conjunto seleccionado y los coeficientes por factor del modelo pueden referirse a puntuaciones diferentes. No interpretar un alfa global como evidencia de unidimensionalidad."
    ) else c(
      "Alpha (Cronbach): k/(k-1) * (1 - sum of item variances / variance of the sum). Uses the raw covariance matrix with pairwise available observations.",
      if (any(x$coefficients$coefficient == "alpha_ordinal")) "Ordinal alpha applies the alpha formula to the polychoric matrix. Its estimand is the continuous latent-response scale, not the observed score.",
      if (any(x$coefficients$coefficient == "spearman_brown")) "Spearman-Brown: 2*r/(1+r) for two items; the correlation type is retained in the results.",
      if (any(x$coefficients$coefficient == "composite_reliability")) "Composite reliability: semTools::compRelSEM(model), per CFA factor and group, using its fitted sample and parameters. Defaults: obs.var=TRUE, tau.eq=FALSE, ord.scale=TRUE.",
      if (any(x$coefficients$coefficient == "omega")) "Omega: single-factor congeneric model (ML for continuous data; WLSMV for ordinal data), without equal-loading constraints. Ordinal omega targets observed-score reliability using Green and Yang's correction. It does not establish unidimensionality or good model fit.",
      if (is.data.frame(item_stats)) "Corrected item-total correlation excludes the item from the sum and uses full-scale complete cases. Alpha uses pairwise observations; omega uses complete cases per fit. Item deletion refits omega and may change N. A coefficient increase alone does not justify deleting an item.",
      if (any(ave)) "AVE: semTools::AVE(model). Average variance extracted is not a reliability coefficient.",
      "Selected-item alpha and model factor coefficients may refer to different scores. Overall alpha does not establish unidimensionality."
    ))
    refs <- unique(c(refs, "Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. Psychometrika, 16, 297-334. https://doi.org/10.1007/BF02310555",
              if (any(x$coefficients$coefficient %in% c("composite_reliability", "omega"))) .pm_semtools_reference()))
  } else if (inherits(x, "factorial_invariance_auto")) {
    ordinal_workflow <- inherits(x, "ordinal_invariance_auto")
    title <- if (ordinal_workflow) tr("Invarianza ordinal", "Ordinal measurement invariance") else
      tr("Invariancia factorial", "Measurement invariance")
    context <- meta(
      Measurement = x$measurement, Estimator = x$estimator,
      Missing = x$effective_settings$missing, Criteria = x$criteria$preset)
    if (identical(x$workflow, "factorial")) {
      counts <- sort(unique(x$item_info$category_counts))
      counts <- counts[is.finite(counts) & counts > 0]
      label <- tr("Opciones de respuesta observadas", "Response options (observed)")
      value <- if (!length(counts)) tr("No disponible", "Not available") else
        paste0(paste(counts, collapse = ", "), if (length(counts) > 1L)
          tr(" (segun el item)", " (varies by item)") else "")
      context$Information[context$Information == "Measurement"] <- label
      context$Value[context$Information == label] <- value
    }
    add(tr("Datos y metodo", "Data and method"), context)
    if (!is.null(x$effective_settings$fit_variant)) {
      add(tr("Variante de ajuste", "Fit-index variant"), paste0(
        "CFI / RMSEA: ", x$effective_settings$fit_variant,
        tr(" (base con estimadores sin correccion y variante predeterminada). Chi-cuadrado, gl y p: scaled con estimadores corregidos; SRMR: base. Los indices solicitados no disponibles no se sustituyen.",
          " (base for uncorrected estimators with the default variant). Chi-square, df and p: scaled for corrected estimators; SRMR: base. Unavailable requested indices are not substituted.")))
    }
    if (identical(x$workflow, "factorial")) {
      add(tr("Metodo e interpretacion", "Method and interpretation"), tr(
        "Se reproducen llamadas nativas de lavaan con cuatro etapas: configural, cargas iguales (metrica), cargas e interceptos iguales (escalar), y cargas, interceptos y residuos iguales (estricta). WLSMV/ULSMV usan respuestas ordenadas de forma predeterminada e identificacion por marcador y parametrizacion delta. ML/MLR/ULS analizan puntuaciones numericas. En datos ordinales, estas etiquetas no garantizan comparaciones independientes; ordinal_invariance_auto() implementa la secuencia identificada de umbrales y cargas.",
        "Reproduces native lavaan calls with four stages: configural, equal loadings (metric), equal loadings and intercepts (scalar), and equal loadings, intercepts and residuals (strict). WLSMV/ULSMV default to ordered responses, marker identification and delta parameterization. ML/MLR/ULS analyze numeric scores. With ordinal data these labels do not guarantee independent comparisons; ordinal_invariance_auto() implements the identified threshold/loading sequence."))
    }
    add(tr("Muestra por grupo", "Sample by group"), data.frame(
      Group = names(x$sample_info$group_sizes), N = as.integer(x$sample_info$group_sizes)))
    fit <- pick(x$fit_table, c("level", "chi2", "df", "p", "cfi", "rmsea", "srmr"))
    names(fit) <- c("Stage", "Chi-square", "df", "p", "CFI", "RMSEA", "SRMR")
    if (ordinal_workflow) {
      labels <- c(configural = tr("Configural", "Configural"),
        thresholds = tr("Igualdad de umbrales", "Equal thresholds"),
        thresholds_loadings = tr("Igualdad de umbrales y cargas", "Equal thresholds and loadings"))
      fit$Stage <- unname(labels[x$fit_table$level])
      fit <- cbind(fit, pick(x$fit_table, c("delta_cfi", "delta_rmsea", "delta_srmr", "decision")))
    }
    add(tr("Ajuste por etapa", "Fit by stage"), fit)
    changes <- pick(x$fit_table, c("level", "comparison_from", "delta_cfi", "delta_rmsea", "delta_srmr", "decision"))
    if (!ordinal_workflow) add(tr("Cambios y decisiones", "Changes and decisions"), changes)
    failed <- x$fit_table[x$fit_table$status != "ok", , drop = FALSE]
    add(tr("Etapas no evaluables", "Non-evaluable stages"), pick(failed, c("level", "status")))
    add(tr("Motivos", "Reasons"), unique(pick(failed, "reason")))
    config <- x$fit_table[x$fit_table$level == "configural", , drop = FALSE]
    if (nrow(config)) {
      values <- as.numeric(config[1, c("cfi", "rmsea", "srmr")])
      limits <- unlist(x$criteria[c("configural_cfi", "configural_rmsea", "configural_srmr")])
      add(tr("Evaluacion configural", "Configural evaluation"), data.frame(
        Index = c("CFI", "RMSEA", "SRMR"), Value = values,
        Criterion = paste(c(">=", "<=", "<="), limits),
        Evaluation = ifelse(!is.finite(values), tr("No evaluable", "Not evaluable"),
          ifelse(c(values[1] >= limits[1], values[2] <= limits[2], values[3] <= limits[3]),
            tr("Cumple", "Met"), tr("No cumple", "Not met")))))
    }
    add(tr("Limites de cambio", "Change cutoffs"), data.frame(
      Index = if (es) c("Disminucion CFI", "Aumento RMSEA", "Aumento SRMR: metrica", "Aumento SRMR: otras") else
        c("CFI decrease", "RMSEA increase", "SRMR increase: metric", "SRMR increase: other"),
      Value = unlist(x$criteria[c("cfi", "rmsea", "srmr_loading", "srmr_other")]),
      row.names = NULL))
    add(tr("Regla conjunta", "Joint decision rule"), tr(
      "No se apoya invariancia cuando la caida de CFI supera su limite Y el aumento de RMSEA O SRMR supera el suyo. Una etapa previa sin apoyo impide interpretar favorablemente las siguientes.",
      "Invariance is not supported when the CFI decrease exceeds its cutoff AND either the RMSEA or SRMR increase exceeds its cutoff. An unsupported earlier stage prevents favorable interpretation of later stages."))
    if (ordinal_workflow) add(tr("Metodo ordinal e interpretacion", "Ordinal method and interpretation"), tr(
      "Se prueban configural, igualdad de umbrales e igualdad de umbrales y cargas, con identificacion Wu y Estabrook (2016), segun Svetina et al. (2020). Los cambios son modelo actual menos modelo previo evaluable: una mejora no es deterioro. Los limites de Chen (2007) son guias descriptivas extrapoladas a datos ordinales. Con items binarios o de tres categorias, la etapa de umbrales puede no ser una prueba independiente; no se declara invariancia por ausencia de grados de libertad adicionales. Esta secuencia no prueba igualdad de interceptos ni de residuos de las respuestas latentes.",
      "Tests configural, equal thresholds, and equal thresholds and loadings with Wu and Estabrook (2016) identification, following Svetina et al. (2020). Changes are current minus preceding evaluable model: improved fit is not deterioration. Chen (2007) cutoffs are descriptive guidelines extrapolated to ordinal data. For binary or three-category items, thresholds may not be an independent test; no invariance is declared from the absence of additional degrees of freedom. This sequence does not test equal latent-response intercepts or residuals."))
  } else if (inherits(x, "psychomatic_dimensionality")) {
    title <- tr("Sensibilidad dimensional", "Dimensionality sensitivity")
    add(tr("Configuraciones", "Configurations"), pick(x$configurations,
      c("requested_method", "quantile", "seed", "n_iter", "n_factors", "status", "reason")))
    add(tr("Acuerdo", "Agreement"), x$agreement)
  } else if (inherits(x, "psychomatic_invariance_partial")) {
    title <- tr("Invariancia parcial", "Partial invariance")
    add(tr("Datos y metodo", "Data and method"), meta(
      Stage = x$stage_effective, Status = x$conclusion_status, Measurement = x$measurement,
      Estimator = x$effective_settings$estimator_original,
      Workflow = x$effective_settings$workflow,
      Fit_variant = x$effective_settings$fit_variant))
    add(tr("Indices seleccionados", "Selected indices"), tr(
      "CFI y RMSEA usan la variante solicitada; chi-cuadrado, gl y p usan scaled cuando corresponde y SRMR usa base. Para estimadores sin correccion, scaled utiliza indices base. Los indices robust no disponibles no se sustituyen.",
      "CFI and RMSEA use the requested variant; chi-square, df and p use scaled when applicable and SRMR uses base. For uncorrected estimators, scaled uses base indices. Unavailable robust indices are not substituted."))
    add(tr("Parametros liberados", "Released parameters"), x$tables$freed_constraints)
    add(tr("Justificacion", "Rationale"), x$tables$rationale)
    add(tr("Anclajes", "Anchors"), x$tables$anchors)
    add(tr("Comparaciones", "Comparisons"), x$tables$comparisons)
    add(tr("Ajuste", "Fit"), pick(x$tables$fit, c("role", "model", "level", "chi2", "df", "p", "cfi", "rmsea", "srmr", "status")))
    add(tr("Conclusion", "Conclusion"), x$conclusion_text)
  } else if (inherits(x, "psychomatic_dif")) {
    title <- tr("Funcionamiento diferencial de items", "Differential item functioning")
    add(tr("Resumen", "Summary"), x$summary)
    add(tr("Resultados por item", "Item results"), pick(x$results, c("item", "n", "p_uniform", "p_nonuniform", "effect_size", "status", "reason")))
    add(tr("Interpretacion", "Interpretation"), pick(x$interpretation, c("item", "dif_type", "decision", "severity", "recommendation")))
  } else if (inherits(x, "psychomatic_irt")) {
    title <- tr("Teoria de respuesta al item", "Item response theory")
    add(tr("Resumen del modelo", "Model summary"), x$model_summary)
    add(tr("Parametros de items", "Item parameters"), x$item_parameters)
    add(tr("Interpretacion", "Interpretation"), x$interpretation)
  } else if (inherits(x, "psychomatic_model_comparison")) {
    title <- tr("Comparacion de modelos", "Model comparison")
    add(tr("Datos y metodo", "Data and method"), meta(N=x$data_summary$n,
      Estimator=x$estimator,
      Selection=if (identical(x$effective_settings$estimator_selection, "automatic")) tr("Automatica", "Automatic") else tr("Usuario", "User")))
    add(tr("Justificacion", "Rationale"), x$effective_settings$estimator_reason)
    fit <- pick(x$fit, c("model", "cfi", "tli", "rmsea", "srmr", "aic", "bic", "status"))
    corrected <- any(x$estimator %in% c("MLR", "MLM", "MLMV", "MLMVS", "WLSMV", "WLSM", "ULSMV", "ULSM"))
    variant <- .pm_null_default(x$effective_settings$fit_variant, "scaled")
    if (variant == "scaled" && !corrected) variant <- "base"
    for (metric in c("cfi", "tli", "rmsea")) {
      if (variant != "base") fit[[metric]] <- x$fit[[paste0(metric, "_", variant)]]
      names(fit)[names(fit) == metric] <- paste0(toupper(metric), " (", variant, ")")
    }
    names(fit)[names(fit) == "srmr"] <- "SRMR (base)"
    add(tr("Ajuste", "Fit"), fit)
    variant_columns <- paste0(toupper(c("cfi", "tli", "rmsea")), " (", variant, ")")
    if (variant != "base" && any(vapply(fit[variant_columns], function(v) any(!is.finite(v)), logical(1)))) {
      add(tr("Indices no disponibles", "Unavailable indices"), tr(
        "Los indices solicitados sin estimacion finita se muestran como '-'; no se sustituyen por otra variante.",
        "Requested indices without a finite estimate are shown as '-'; no other variant is substituted."))
    }
    add(tr("Recomendacion", "Recommendation"), if (is.list(x$recommendation)) x$recommendation$message else x$recommendation)
  } else if (inherits(x, "inv_align_auto_result")) {
    title <- tr("Alineamiento", "Alignment")
    summary <- x$summary
    if (is.data.frame(summary)) {
      summary$parameter <- ifelse(summary$parameter == "loadings", tr("Cargas", "Loadings"),
        ifelse(summary$parameter == "intercepts", tr("Interceptos", "Intercepts"), summary$parameter))
      names(summary) <- vapply(names(summary), function(n) switch(n,
        parameter=tr("Parametro", "Parameter"), percent_noninvariant=tr("No invariantes (%)", "Noninvariant (%)"),
        cutoff=tr("Limite (%)", "Cutoff (%)"), invariant=tr("Compatible", "Compatible"), n), character(1))
      add(tr("Resumen", "Summary"), summary)
    }
    for (kind in c("loadings", "intercepts")) {
      tab <- x$item_results[[kind]]$summary
      compact <- pick(tab, c("item", "joint_parameter", "noninvariant_groups", "percent_noninvariant", "noninvariant"))
      if (is.data.frame(compact)) {
        if ("noninvariant_groups" %in% names(compact)) compact$noninvariant_groups <- as.integer(compact$noninvariant_groups)
        labels <- c(item="Item", joint_parameter=tr("Estimacion conjunta", "Joint estimate"),
          noninvariant_groups=tr("Grupos senalados", "Flagged N"),
          percent_noninvariant=tr("No invariantes (%)", "Noninvariant (%)"), noninvariant="noninvariant")
        names(compact) <- unname(labels[names(compact)])
      }
      add(if (kind == "loadings") tr("Cargas", "Loadings") else tr("Interceptos", "Intercepts"), compact)
      if (is.data.frame(tab) && "flagged_groups" %in% names(tab)) {
        flagged <- tab[!is.na(tab$flagged_groups) & nzchar(tab$flagged_groups), c("item", "flagged_groups"), drop=FALSE]
        names(flagged) <- c("Item", tr("Grupos", "Groups"))
        add(if (kind == "loadings") tr("Grupos senalados: cargas", "Flagged groups: loadings") else
          tr("Grupos senalados: interceptos", "Flagged groups: intercepts"), flagged)
      }
    }
    groups <- x$group_parameters
    if (!is.null(groups) && NROW(groups)) {
      groups <- as.data.frame(groups)
      labels <- rownames(groups)
      if (is.null(labels) || identical(labels, as.character(seq_len(nrow(groups))))) labels <- x$dimensions$group_names
      if (length(labels) != nrow(groups)) labels <- rep(NA_character_, nrow(groups))
      groups <- data.frame(Group=labels, groups, row.names=NULL, check.names=FALSE)
      names(groups)[names(groups) == "alpha0"] <- tr("Media latente", "Latent mean")
      names(groups)[names(groups) == "psi0"] <- tr("DE latente", "Latent SD")
      add(tr("Parametros por grupo", "Group parameters"), groups)
    }
  } else if (inherits(x, "psychomatic_item_screen")) {
    title <- tr("Revision de items", "Item screening")
    tab <- x
    class(tab) <- "data.frame"
    for (key in c("missing_pct", "floor_pct", "ceiling_pct")) tab[[key]] <- 100 * tab[[key]]
    blocks <- list(
      c("item", "n", "missing", "missing_pct", "mean", "sd", "min", "max"),
      c("item", "skewness", "kurtosis", "n_categories", "floor_pct", "ceiling_pct"),
      c("item", "item_rest_correlation", "zero_variance", "floor_flag", "ceiling_flag", "probable_reverse"))
    titles <- if (es) c("Descriptivos", "Distribucion y extremos", "Revision de items y alertas") else
      c("Descriptives", "Distribution and endpoints", "Item checks and flags")
    labels <- if (es) c(item="Item", n="N", missing="Perdidos", missing_pct="Perdidos (%)",
      mean="Media", sd="DE", min="Min", max="Max", skewness="Asimetria", kurtosis="Curtosis",
      n_categories="Categorias", floor_pct="Suelo (%)", ceiling_pct="Techo (%)",
      item_rest_correlation="r item-resto", zero_variance="Sin varianza", floor_flag="Suelo", ceiling_flag="Techo",
      probable_reverse="Revisar clave") else c(item="Item", n="N", missing="Missing", missing_pct="Missing (%)",
      mean="Mean", sd="SD", min="Min", max="Max", skewness="Skewness", kurtosis="Kurtosis",
      n_categories="Categories", floor_pct="Floor (%)", ceiling_pct="Ceiling (%)",
      item_rest_correlation="Item-rest r", zero_variance="No variance", floor_flag="Floor", ceiling_flag="Ceiling",
      probable_reverse="Review key")
    for (i in seq_along(blocks)) {
      value <- tab[blocks[[i]]]
      for (key in names(value)[vapply(value, is.logical, logical(1))]) value[[key]] <-
        ifelse(value[[key]], tr("Si", "Yes"), tr("No", "No"))
      names(value) <- unname(labels[names(value)])
      add(titles[i], value)
    }
  } else if (is.data.frame(x)) {
    title <- tr("Descriptivos", "Descriptives")
    add(title, x)
    add(tr("Porcentajes de respuesta", "Response percentages"), attr(x, "response_percentages"))
  } else if (is.matrix(x)) {
    method <- attr(x, "effective_type", exact = TRUE)
    title <- switch(.pm_null_default(method, ""),
      pearson = tr("Correlaciones de Pearson", "Pearson correlations"),
      poly = tr("Correlaciones policoricas", "Polychoric correlations"),
      tetra = tr("Correlaciones tetracoricas", "Tetrachoric correlations"),
      tr("Correlaciones", "Correlations"))
    add(title, .pm_matrix_table(x))
    add(tr("Justificacion", "Rationale"), attr(x, "rationale", exact = TRUE))
    refs <- attr(x, "references", exact = TRUE)
  } else {
    # Compatibility for optional utility objects and explicitly requested audits.
    sections <- .pm_object_tables(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    conclusion <- x$conclusion
    if (is.list(conclusion)) conclusion <- conclusion$reason
    if (inherits(x, "cfa_auto")) conclusion <- .pm_cfa_conclusion(x, language)
    if (inherits(x, "factorial_invariance_auto")) conclusion <- .ifa_conclusion_text(x$fit_table, language)
    add(tr("Conclusion", "Conclusion"), conclusion)
    for (kind in c("errors", "warnings")) {
      cond <- x[[kind]]
      if (is.data.frame(cond) && "message" %in% names(cond)) {
        messages <- unique(cond$message)
        if (es) {
          messages <- sub("^([0-9]+) cases were removed by the item missing-data policy\\.$",
            "Se excluyeron \\1 casos por la politica de datos perdidos en los items.", messages)
          messages <- sub("^([0-9]+) cases with missing group membership were removed\\.$",
            "Se excluyeron \\1 casos sin informacion de grupo.", messages)
        }
        add(if (kind == "errors") tr("Errores del analisis", "Analysis errors") else tr("Advertencias", "Warnings"), messages)
      }
    }
  }
  add(tr("Referencias", "References"), .pm_reference_text(refs))
  sections <- lapply(sections, .pm_label_table, language = language)
  list(title = title, sections = sections, language = language)
}

.pm_label_table <- function(table, language) {
  delta <- c(delta_cfi = "\u0394CFI", delta_rmsea = "\u0394RMSEA", delta_srmr = "\u0394SRMR")
  pos <- match(names(table), names(delta))
  names(table)[!is.na(pos)] <- unname(delta[pos[!is.na(pos)]])
  if (!identical(language, "esp")) {
    reasons <- c(
      configural_fit_failed_at_least_one_descriptive_gate = "At least one configural fit criterion was not met; see Configural evaluation.",
      a_less_restrictive_stage_was_not_supported = "A less restrictive stage was not supported; this comparison is not interpreted.",
      failed_configural_fit = "Insufficient configural fit")
    for (key in intersect(c("reason", "status"), names(table))) {
      idx <- match(table[[key]], names(reasons))
      table[[key]][!is.na(idx)] <- unname(reasons[idx[!is.na(idx)]])
    }
    for (key in intersect(c("metric", "Metric"), names(table))) {
      replacements <- c(alpha_observed = "Alpha (Cronbach)", alpha_ordinal = "Ordinal alpha",
        composite_reliability = "Omega (composite reliability)", ave = "AVE")
      idx <- match(table[[key]], names(replacements))
      table[[key]][!is.na(idx)] <- unname(replacements[idx[!is.na(idx)]])
    }
    labels <- c(item_rest_r = "Item-rest r", alpha_if_deleted = "Alpha if deleted",
      alpha_ordinal_if_deleted = "Ordinal alpha if deleted", omega_if_deleted = "Omega if deleted",
      p_value = "p", CI_lower = "Lower CI", CI_upper = "Upper CI")
    idx <- match(names(table), names(labels))
    names(table)[!is.na(idx)] <- unname(labels[idx[!is.na(idx)]])
    return(table)
  }
  if (identical(names(table), "Text")) return(table)
  labels <- c(Information = "Informacion", Value = "Valor", value = "Valor",
    Group = "Grupo", Evaluation = "Evaluacion", Loading = "Carga", SE = "EE", Alert = "Alerta",
    p_value = "p", CI_lower = "IC inferior", CI_upper = "IC superior",
    Estimator = "Estimador", Factors = "Factores", Iterations = "Simulaciones",
    Quantile = "Cuantil", Seed = "Semilla", Correlation = "Correlacion",
    Rotation = "Rotacion", Missing = "Datos perdidos", Excluded = "Casos excluidos",
    Model = "Modelo", Converged = "Convergencia", Complete = "Casos completos",
    Measurement = "Nivel de medicion", Keying = "Claves inversas", Status = "Estado",
    Criteria = "Criterios", Criterion = "Criterio", Selection = "Seleccion", SameFactor = "Mismo factor", No = "Nro.",
    "Chi-square" = "Chi-cuadrado",
    Test = "Prueba", Statistic = "Estadistico", df = "gl", Sample = "Muestra",
    Component = "Componente", metric = "Indice", Metric = "Indice", Index = "Indice",
    Variant = "Variante", factor = "Factor", group = "Grupo", coefficient = "Coeficiente",
    FactorLevel = "Nivel factorial", TotalPercent = "Varianza total (%)", CumulativePercent = "Acumulada (%)",
    item_rest_r = "r item-resto", alpha_if_deleted = "Alfa sin item",
    alpha_ordinal_if_deleted = "Alfa ordinal sin item", omega_if_deleted = "Omega sin item",
    estimate = "Estimacion", status = "Estado", reason = "Motivo", n_items = "N items",
    items = "Items", Item = "Item", item = "Item", Stage = "Etapa", level = "Etapa",
    comparison_from = "Comparado con", delta_cfi = "Cambio CFI", delta_rmsea = "Cambio RMSEA",
    delta_srmr = "Cambio SRMR", decision = "Decision", chi2 = "Chi-cuadrado",
    chi2_variant = "Variante chi2", cfi_variant = "Variante CFI", rmsea_variant = "Variante RMSEA",
    srmr_variant = "Variante SRMR", requested_method = "Correlacion", quantile = "Cuantil",
    seed = "Semilla", n_iter = "Simulaciones", n_factors = "Factores", n = "N",
    n_configurations = "Configuraciones", n_evaluable = "Evaluables", n_skipped = "Omitidas",
    n_distinct_solutions = "Soluciones", exact_agreement = "Acuerdo total",
    observed_solutions = "Factores retenidos", Type = "Tipo", Syntax = "Sintaxis", Note = "Nota",
    p_uniform = "p uniforme", p_nonuniform = "p no uniforme", effect_size = "Tamano del efecto",
    dif_type = "Tipo DIF", severity = "Magnitud", recommendation = "Interpretacion",
    model = "Modelo", cfi = "CFI", tli = "TLI", rmsea = "RMSEA", srmr = "SRMR",
    aic = "AIC", bic = "BIC", alpha = "Alfa", omega = "Omega", spearman_brown = "Spearman-Brown")
  states <- c(ok = "Evaluable", supported = "Compatible", not_supported = "No compatible",
    not_evaluable = "No evaluable", configural = "Configural", thresholds = "Umbrales",
    metric = "Metrica", scalar = "Escalar", strict = "Estricta", skipped = "Omitida")
  for (name in intersect(c("status", "decision", "level", "Stage", "comparison_from"), names(table))) {
    pos <- match(table[[name]], names(states))
    table[[name]][!is.na(pos)] <- unname(states[pos[!is.na(pos)]])
  }
  reasons <- c(
    configural_fit_failed_at_least_one_descriptive_gate = "El ajuste configural no satisface al menos un criterio descriptivo.",
    a_less_restrictive_stage_was_not_supported = "No se interpreto esta comparacion porque una etapa anterior no obtuvo apoyo.",
    failed_configural_fit = "Ajuste configural insuficiente",
    cross_loading = "Carga cruzada", residual_covariance = "Covarianza residual",
    insufficient_items = "Numero insuficiente de items",
    semtools_not_available = "semTools no esta instalado",
    inadmissible_model = "Modelo no admisible",
    insufficient_complete_cases = "Casos completos insuficientes",
    model_not_converged = "El modelo no convergio", model_fit_error = "Fallo el ajuste del modelo",
    constant_item = "Item constante", non_finite_input = "Datos no finitos",
    supplied_model_requires_explicit_refit = "Requiere especificar el modelo reducido",
    insufficient_or_constant_scores = "Datos insuficientes o puntuaciones constantes",
    polychoric_not_evaluable = "Matriz policorica no evaluable",
    alpha_if_deleted = "Alfa sin item", alpha_ordinal_if_deleted = "Alfa ordinal sin item",
    omega_if_deleted = "Omega sin item", item_rest = "Correlacion item-resto",
    zero_variance = "Varianza nula", mixed_measurement_levels = "Niveles de medicion mixtos",
    non_finite_coefficient = "Coeficiente no finito")
  for (name in intersect(c("reason", "status", "Type", "Index"), names(table))) {
    idx <- match(table[[name]], names(reasons))
    table[[name]][!is.na(idx)] <- unname(reasons[idx[!is.na(idx)]])
  }
  metric_labels <- c(composite_reliability = "Omega (confiabilidad compuesta)", ave = "AVE",
    alpha_observed = "Alfa (Cronbach)", alpha_ordinal = "Alfa ordinal",
    omega_total = "Omega total", omega_h = "Omega jerarquico", omega_hs = "Omega especifico",
    ecv = "ECV", puc = "PUC", h_general = "H general", h = "H")
  for (name in intersect(c("metric", "Metric"), names(table))) {
    idx <- match(table[[name]], names(metric_labels))
    table[[name]][!is.na(idx)] <- unname(metric_labels[idx[!is.na(idx)]])
  }
  if (identical(names(table), c("Information", "Value"))) {
    values <- c(none="Ninguna", listwise="Casos completos", pairwise="Pares disponibles",
      continuous="Continuo", mixed="Mixto", categorical_descriptive_extrapolation="Orientaciones descriptivas para items ordinales",
      chen_informed_guidelines="Orientaciones basadas en Chen (2007)", custom="Personalizados")
    pos <- match(table$Value, names(values))
    table$Value[!is.na(pos)] <- unname(values[pos[!is.na(pos)]])
    idx <- match(table$Information, names(labels))
    table$Information[!is.na(idx)] <- unname(labels[idx[!is.na(idx)]])
  }
  idx <- match(names(table), names(labels))
  names(table)[!is.na(idx)] <- unname(labels[idx[!is.na(idx)]])
  table
}

.pm_reader_table <- function(table) {
  out <- table
  for (name in names(out)) {
    if (is.list(out[[name]])) out[[name]] <- vapply(out[[name]], .pm_audit_display, character(1))
    if (is.numeric(out[[name]])) out[[name]][!is.finite(out[[name]])] <- NA_real_
  }
  out
}

.pm_display_table <- function(table, digits = 3) {
  out <- .pm_reader_table(table)
  for (name in names(out)) {
    if (is.numeric(out[[name]])) {
      value <- out[[name]]
      count <- is.integer(value) || tolower(name) %in% c("no", "nro.", "n casos", "n cases", "df", "gl", "n", "n items", "n_items", "n_factors", "n_iter", "seed", "items", "factores", "iteraciones", "simulaciones", "semilla", "configuraciones", "evaluables", "omitidas", "soluciones", "requested", "valid", "failed", "solicitadas", "v\u00e1lidas", "fallidas", "valid_permutations", "permutaciones v\u00e1lidas")
      out[[name]] <- ifelse(is.na(value), "-", formatC(value, format = "f", digits = if (count) 0 else digits))
      if (tolower(name) %in% c("p", "pvalue", "p.value", "p_valor", "p_value", "p_raw", "p_bh", "p_uniform", "p_nonuniform")) {
        small <- is.finite(value) & value < 10^-digits & value >= 0
        out[[name]][small] <- paste0("<", formatC(10^-digits, format = "f", digits = digits))
      }
    } else {
      out[[name]] <- as.character(out[[name]])
      out[[name]][is.na(out[[name]])] <- "-"
    }
  }
  metric_col <- intersect(c("Metric", "metric", "Index", "Indice"), names(out))
  value_col <- intersect(c("value", "Value", "Valor"), names(out))
  if (length(metric_col) && length(value_col) && is.numeric(table[[value_col[1L]]])) {
    metric <- tolower(out[[metric_col[1L]]])
    values <- table[[value_col[1L]]]
    count <- metric %in% c("df", "gl") & is.finite(values)
    out[[value_col[1L]]][count] <- formatC(values[count], format = "f", digits = 0)
    small_p <- grepl("^pvalue$|^p \\(|p-close", metric) & is.finite(values) & values < 10^-digits
    out[[value_col[1L]]][small_p] <- paste0("<", formatC(10^-digits, format = "f", digits = digits))
  }
  out
}

.pm_console_report <- function(x, digits = NULL, language = NULL, detail = "report") {
  if (is.null(digits)) digits <- .pm_null_default(attr(x, "display_digits"), 3L)
  digits <- .pm_validate_scalar_number(digits, "digits", 0, 15, integer = TRUE)
  report <- if (inherits(x, "psychomatic_network")) .pm_network_sections(x, language, detail = detail) else .pm_analysis_sections(x, language)
  width <- max(60L, min(110L, getOption("width", 80L)))
  cat("\n", strrep("=", width), "\n", toupper(report$title), "\n", strrep("=", width), "\n", sep = "")
  for (name in names(report$sections)) {
    table <- report$sections[[name]]
    if (!is.data.frame(table) || !nrow(table)) next
    if (!identical(name, report$title)) cat("\n", toupper(name), "\n", sep = "") else cat("\n")
    if (identical(names(table), "Text")) {
      for (value in table$Text) {
        text <- if (.pm_reference_section(name)) .pm_reference_console(value, width=width) else paste(strwrap(value, width=width),collapse="\n")
        cat(text, "\n\n", sep="")
      }
    } else if (identical(names(table), c("Information", "Value")) || identical(names(table), c("Informacion", "Valor"))) {
      for (i in seq_len(nrow(table))) cat(paste(strwrap(paste0(table[i, 1], ": ", table[i, 2]), width = width), collapse = "\n"), "\n", sep = "")
    } else {
      table <- .pm_display_table(table, digits)
      # Long explanatory cells get their own wrapped paragraphs instead of wide console rows.
      long <- vapply(table, function(v) max(nchar(v), na.rm = TRUE) > 45L, logical(1))
      if (any(!long)) print(table[!long], row.names = FALSE, right = FALSE, quote = FALSE)
      if (any(long)) for (i in seq_len(nrow(table))) {
        for (col in names(table)[long]) {
          text <- paste0(if (nrow(table) > 1) paste0(i, ". ") else "", col, ": ", table[[col]][i])
          cat(paste(strwrap(text, width = width), collapse = "\n"), "\n", sep = "")
        }
      }
    }
  }
  invisible(x)
}

#' @export
print.psychomatic_dif <- function(x, digits = NULL, ...) .pm_console_report(x, digits)
#' @export
print.psychomatic_irt <- function(x, digits = NULL, ...) .pm_console_report(x, digits)
#' @export
print.psychomatic_model_comparison <- function(x, digits = NULL, ...) .pm_console_report(x, digits)

.pm_write_analysis_excel <- function(reports, path, digits = 3, overwrite = FALSE) {
  path <- .pm_safe_output_path(path, ".xlsx", overwrite)
  wb <- openxlsx::createWorkbook()
  es <- identical(reports[[1]]$language, "esp")
  index <- if (es) "Contenido" else "Contents"
  openxlsx::addWorksheet(wb, index, gridLines = FALSE)
  header <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#243B53", textDecoration = "bold", wrapText = TRUE, valign = "center")
  body <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, wrapText = TRUE, valign = "top")
  title <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold")
  section <- openxlsx::createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold", fgFill = "#EAF0F5")
  number <- openxlsx::createStyle(numFmt = if (digits) paste0("0.", strrep("0", digits)) else "0")
  used <- index
  mapping <- list()
  image_files <- character()
  image_sheets <- integer()
  formatted_references <- character()
  on.exit(unlink(image_files), add = TRUE)
  contents <- data.frame(Analysis = character(), Sheet = character())
  for (i in seq_along(reports)) {
    report <- reports[[i]]
    sheet <- .pm_excel_sheet_name(report$title, used)
    used <- c(used, sheet)
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
    contents <- rbind(contents, data.frame(Analysis = report$title, Sheet = sheet))
    openxlsx::writeData(wb, sheet, report$title, startRow = 2, startCol = 2)
    openxlsx::addStyle(wb, sheet, title, rows = 2, cols = 2)
    row <- 4L
    max_cols <- max(vapply(report$sections, ncol, integer(1)), 4L)
    # Numeric tables keep modest widths; prose spans the same content area.
    openxlsx::setColWidths(wb, sheet, cols = 1, widths = 3)
    openxlsx::setColWidths(wb, sheet, cols = 2:(max_cols + 1L), widths = c(26, rep(18, max_cols - 1L)))
    for (name in names(report$sections)) {
      table <- report$sections[[name]]
      if (!nrow(table)) next
      openxlsx::writeData(wb, sheet, name, startRow = row, startCol = 2)
      openxlsx::mergeCells(wb, sheet, cols = 2:(max_cols + 1L), rows = row)
      openxlsx::addStyle(wb, sheet, section, rows = row, cols = 2:(max_cols + 1L), gridExpand = TRUE)
      row <- row + 1L
      if (identical(names(table), "Text")) {
        if (.pm_reference_section(name)) formatted_references <- c(formatted_references,table$Text)
        for (value in table$Text) {
          openxlsx::writeData(wb, sheet, value, startRow = row, startCol = 2)
          openxlsx::mergeCells(wb, sheet, cols = 2:(max_cols + 1L), rows = row)
          openxlsx::addStyle(wb, sheet, body, rows = row, cols = 2)
          openxlsx::setRowHeights(wb, sheet, rows = row, heights = 16 * max(2, ceiling(nchar(value) / (max_cols * 13))))
          row <- row + 1L
        }
      } else {
        table <- .pm_reader_table(table)
        openxlsx::writeData(wb, sheet, table, startRow = row, startCol = 2, headerStyle = header)
        openxlsx::setRowHeights(wb, sheet, rows = row, heights = 32)
        data_rows <- row + seq_len(nrow(table))
        openxlsx::addStyle(wb, sheet, body, rows = data_rows, cols = 2:(ncol(table) + 1L), gridExpand = TRUE)
        nums <- which(vapply(table, is.numeric, logical(1)))
        if (length(nums)) openxlsx::addStyle(wb, sheet, number, rows = data_rows, cols = nums + 1L, gridExpand = TRUE, stack = TRUE)
        counts <- which(vapply(table, is.integer, logical(1)) | tolower(names(table)) %in% c("df", "gl", "n", "n items", "n_items", "n_factors", "n_iter", "seed", "factores", "simulaciones", "semilla"))
        if (length(counts)) openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = "0"), rows = data_rows, cols = counts + 1L, gridExpand = TRUE, stack = TRUE)
        p_style <- openxlsx::createStyle(numFmt = paste0('[<', 10^-digits, ']"<',
          formatC(10^-digits, format = "f", digits = digits), '";',
          if (digits) paste0("0.", strrep("0", digits)) else "0"))
        p_cols <- which(tolower(names(table)) %in% c("p", "pvalue", "p.value", "p_value", "p_raw", "p_bh", "p uniforme", "p no uniforme"))
        if (length(p_cols)) openxlsx::addStyle(wb, sheet, p_style, rows = data_rows, cols = p_cols + 1L, gridExpand = TRUE, stack = TRUE)
        metrics <- intersect(c("Metric", "metric", "Index", "Indice"), names(table))
        values <- intersect(c("Value", "value", "Valor"), names(table))
        if (length(metrics) && length(values)) {
          metric <- tolower(table[[metrics[1]]])
          value_col <- match(values[1], names(table)) + 1L
          df_rows <- which(metric %in% c("df", "gl"))
          if (length(df_rows)) openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = "0"), rows = data_rows[df_rows], cols = value_col, stack = TRUE)
          p_rows <- which(grepl("^pvalue$|^p \\(|p-close", metric))
          if (length(p_rows)) openxlsx::addStyle(wb, sheet, p_style, rows = data_rows[p_rows], cols = value_col, stack = TRUE)
        }
        for (r in seq_len(nrow(table))) {
          lengths <- vapply(table[r, , drop = FALSE], function(v) nchar(as.character(v)), integer(1))
          openxlsx::setRowHeights(wb, sheet, rows = data_rows[r], heights = max(24, 15 * max(ceiling(lengths / c(24, rep(16, ncol(table)-1))), na.rm = TRUE)))
        }
        row <- row + nrow(table) + 1L
      }
      mapping[[length(mapping) + 1L]] <- data.frame(object_index = i, display_label = report$title, section = name,
        requested_sheet = report$title, effective_sheet = sheet, stringsAsFactors = FALSE)
      row <- row + 2L
    }
    for (figure in report$figures) {
      png_path <- tempfile(fileext = ".png")
      image_files <- c(image_files, png_path)
      figure$draw(png_path)
      openxlsx::writeData(wb, sheet, figure$title, startRow = row, startCol = 2)
      openxlsx::insertImage(wb, sheet, png_path, startRow = row + 1L, startCol = 2,
        width = 10, height = 7, units = "in", dpi = 300)
      row <- row + 38L
      openxlsx::writeData(wb, sheet, figure$caption, startRow = row, startCol = 2)
      row <- row + 4L
      image_sheets <- unique(c(image_sheets, i + 1L))
    }
    openxlsx::pageSetup(wb, sheet, orientation = if (max_cols > 7L) "landscape" else "portrait", fitToWidth = 1, fitToHeight = 0)
  }
  openxlsx::writeData(wb, index, "PsychoMatic", startRow = 2, startCol = 2)
  openxlsx::addStyle(wb, index, title, rows = 2, cols = 2)
  openxlsx::writeData(wb, index, contents, startRow = 4, startCol = 2, headerStyle = header)
  openxlsx::setColWidths(wb, index, cols = 1:3, widths = c(3, 46, 34))
  # These reports contain no drawings/comments. openxlsx creates placeholder
  # relationships for them even when the corresponding ZIP parts do not exist.
  for (i in seq_along(wb$worksheets)) {
    # Reader reports do not contain comments/VML. Preserve actual figure parts,
    # but remove the separate phantom VML relationship also on image sheets.
    wb$worksheets_rels[[i]] <- wb$worksheets_rels[[i]][
      !grepl('/vmlDrawing"', wb$worksheets_rels[[i]])]
    if (i %in% image_sheets) next
    wb$worksheets[[i]]$drawing <- character()
    wb$worksheets_rels[[i]] <- wb$worksheets_rels[[i]][
      !grepl('/(drawing|vmlDrawing)"', wb$worksheets_rels[[i]])]
  }
  .pm_reference_excel(wb, formatted_references)
  openxlsx::saveWorkbook(wb, path, overwrite = overwrite)
  attr(path, "sheet_map") <- do.call(rbind, mapping)
  invisible(path)
}

.pm_write_analysis_word <- function(reports, path, digits = 3, overwrite = FALSE) {
  path <- .pm_safe_output_path(path, ".docx", overwrite)
  doc <- officer::read_docx()
  image_files <- character()
  on.exit(unlink(image_files), add = TRUE)
  doc <- officer::docx_set_paragraph_style(doc, style_id = "PsychoMaticTitle", style_name = "Title",
    fp_p = officer::fp_par(padding.bottom = 10, keep_with_next = TRUE),
    fp_t = officer::fp_text_lite(font.family = "Arial", font.size = 20, bold = TRUE))
  doc <- officer::docx_set_paragraph_style(doc, style_id = "PsychoMaticSection", style_name = "PsychoMatic Section",
    fp_p = officer::fp_par(padding.top = 12, padding.bottom = 6, keep_with_next = TRUE),
    fp_t = officer::fp_text_lite(font.family = "Arial", font.size = 12, bold = TRUE))
  doc <- officer::docx_set_paragraph_style(doc, style_id = "PsychoMaticBody", style_name = "PsychoMatic Body",
    fp_p = officer::fp_par(padding.bottom = 6, line_spacing = 1.1),
    fp_t = officer::fp_text_lite(font.family = "Arial", font.size = 10))
  for (i in seq_along(reports)) {
    if (i > 1L) doc <- officer::body_add_break(doc)
    report <- reports[[i]]
    doc <- officer::body_add_par(doc, report$title, style = "Title")
    for (name in names(report$sections)) {
      table <- report$sections[[name]]
      doc <- officer::body_add_par(doc, name, style = "PsychoMatic Section")
      if (identical(names(table), "Text")) {
        for (value in table$Text) doc <- if (.pm_reference_section(name)) .pm_reference_word(doc,value) else officer::body_add_par(doc,value,style="PsychoMatic Body")
      } else {
        displayed <- .pm_display_table(table, digits)
        ft <- flextable::flextable(displayed)
        ft <- flextable::theme_booktabs(ft)
        ft <- flextable::font(ft, fontname = "Arial", part = "all")
        ft <- flextable::fontsize(ft, size = 10, part = "all")
        ft <- flextable::autofit(ft)
        weights <- vapply(seq_along(displayed), function(j) {
          sqrt(min(45, max(8, nchar(names(displayed)[j]), nchar(displayed[[j]]))))
        }, numeric(1))
        ft <- flextable::width(ft, width = 6.3 * weights / sum(weights))
        ft <- flextable::set_table_properties(ft, layout = "fixed", width = 1)
        doc <- flextable::body_add_flextable(doc, ft)
      }
    }
    for (figure in report$figures) {
      png_path <- tempfile(fileext = ".png")
      image_files <- c(image_files, png_path)
      figure$draw(png_path)
      doc <- officer::body_add_par(doc, figure$title, style = "PsychoMatic Section")
      doc <- officer::body_add_img(doc, png_path, width = 6.3, height = 4.41)
      doc <- officer::body_add_par(doc, figure$caption, style = "PsychoMatic Body")
    }
  }
  print(doc, target = path)
  invisible(path)
}
