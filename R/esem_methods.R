.pm_esem_manuscript <- function(x, language, digits = 3) {
  es <- identical(.pm_normalize_language(language), "esp")
  tr <- function(a, b) if (es) a else b
  f <- function(v) if (length(v) != 1L || !is.finite(v)) {
    tr("no disponible", "unavailable")
  } else formatC(v, format = "f", digits = digits)
  s <- x$effective_settings; d <- x$data_summary
  if (!identical(x$conclusion$status, "evaluable")) {
    return(paste(tr("El modelo ESEM no produjo una soluci\u00f3n interpretable.",
                    "The ESEM model could not produce an interpretable solution; its results are not interpretable."),
      tr("Motivo registrado:", "Recorded reason:"), x$conclusion$reason,
      tr("Se requiere revisar los diagn\u00f3sticos antes de elaborar una interpretaci\u00f3n sustantiva.",
         "Diagnostics require review before a substantive interpretation can be drafted.")))
  }
  fi <- x$results$fit_indices
  v <- function(key) fi$value[match(key, fi$metric)]
  text <- sprintf(tr("Se estim\u00f3 un modelo ESEM de %s factores con %s \u00edtems (N = %s), mediante %s y rotaci\u00f3n oblicua %s.",
                     "An ESEM model with %s factors and %s items (N = %s) was estimated using %s and oblique %s rotation."),
                  s$n_factors, d$n_items, d$n, s$estimator, s$rotation)
  text <- paste(text, if (s$rotation == "target") {
    tr("El patr\u00f3n te\u00f3rico orient\u00f3 la rotaci\u00f3n target; las cargas secundarias se estimaron libremente y la rotaci\u00f3n busc\u00f3 aproximarlas a cero, sin fijarlas a ese valor.",
       "The theoretical pattern guided target rotation; secondary loadings were freely estimated and targeted near zero rather than fixed at zero.")
  } else sprintf(tr("La soluci\u00f3n es exploratoria: geomin favoreci\u00f3 un patr\u00f3n de cargas simple (epsilon = %s); las denominaciones F1, F2, etc. requieren interpretaci\u00f3n te\u00f3rica posterior.",
                    "The solution is exploratory: geomin favored a simple loading pattern (epsilon = %s); labels F1, F2, etc. require subsequent theoretical interpretation."),
                   format(s$rotation_args$geomin_epsilon, scientific = FALSE)))
  p <- v("pvalue")
  p_text <- if (length(p) == 1L && is.finite(p) && p < .001) "p < .001" else paste0("p = ", f(p))
  text <- paste(text, sprintf(tr("Los \u00edndices de ajuste fueron: \u03c7\u00b2(%s) = %s, %s; CFI = %s; TLI = %s; RMSEA = %s (IC del 90%% [%s, %s]); SRMR = %s.",
    "Fit indices were: \u03c7\u00b2(%s) = %s, %s; CFI = %s; TLI = %s; RMSEA = %s (90%% CI [%s, %s]); SRMR = %s."),
    if (length(v("df")) == 1L && is.finite(v("df"))) format(v("df"), trim = TRUE) else f(v("df")),
    f(v("chisq")), p_text, f(v("cfi")), f(v("tli")), f(v("rmsea")),
    f(v("rmsea.ci.lower")), f(v("rmsea.ci.upper")), f(v("srmr"))),
    sprintf(tr("Se solicit\u00f3 la variante %s para CFI, TLI y RMSEA; \u03c7\u00b2 utiliza la variante %s y SRMR la base.",
               "The %s variant was requested for CFI, TLI and RMSEA; \u03c7\u00b2 uses the %s variant and SRMR the base variant."),
            s$fit_variant, fi$variant[fi$metric == "chisq"]))
  if (isTRUE(x$diagnostics$saturated)) text <- paste(text,
    tr("El modelo no tiene grados de libertad positivos; sus \u00edndices globales no permiten valorar la adecuaci\u00f3n de la estructura.",
       "The model has no positive degrees of freedom; its global indices cannot assess the adequacy of the structure."))
  id <- x$results$item_diagnostics
  primary <- if (s$rotation == "target") id$expected_loading else id$largest_loading
  cross <- if (s$rotation == "target") {
    ld <- x$results$loadings
    abs(ld$standardized[ld$factor != id$expected_factor[match(ld$item, id$item)]])
  } else abs(id$second_loading)
  text <- paste(text, sprintf(tr("Las cargas estandarizadas %s oscilaron, en valor absoluto, entre %s y %s; la mayor carga secundaria absoluta fue %s.",
    "Standardized loadings %s ranged in absolute value from %s to %s; the largest absolute secondary loading was %s."),
    tr(if (s$rotation == "target") "en los factores previstos" else "dominantes por \u00edtem",
       if (s$rotation == "target") "on the expected factors" else "with the largest magnitude per item"),
    f(min(abs(primary))), f(max(abs(primary))), f(max(cross))))
  if (s$rotation == "target") text <- paste(text,
    sprintf(tr("En %s de %s \u00edtems, la carga de mayor magnitud correspondi\u00f3 al factor previsto.",
               "For %s of %s items, the largest loading corresponded to the expected factor."),
            sum(id$expected_is_largest), nrow(id)))
  fc <- x$results$factor_correlations$correlation
  if (length(fc)) text <- paste(text,
    sprintf(tr("Las correlaciones entre factores variaron entre %s y %s.",
               "Factor correlations ranged from %s to %s."), f(min(fc)), f(max(fc))))
  cmp <- x$results$cfa_comparison
  if (!is.null(cmp)) text <- paste(text, if (identical(cmp$status, "evaluable")) {
    cv <- function(key) cmp$table$delta[match(key, cmp$table$metric)]
    sprintf(tr("La comparaci\u00f3n descriptiva con el CFA de referencia, sobre las mismas observaciones, arroj\u00f3 \u0394CFI (ESEM \u2212 CFA) = %s y \u0394RMSEA = %s; estas diferencias no constituyen una prueba de selecci\u00f3n entre modelos.",
               "The descriptive comparison with the reference CFA on the same observations yielded \u0394CFI (ESEM \u2212 CFA) = %s and \u0394RMSEA = %s; these differences do not constitute a model selection test."),
            f(cv("cfi")), f(cv("rmsea")))
  } else paste(tr("La comparaci\u00f3n CFA no fue evaluable:", "The CFA comparison was not evaluable:"), cmp$reason))
  text
}

.pm_esem_sections <- function(x, language = NULL, destination = "console") {
  language <- .pm_normalize_language(language, x$effective_settings$language)
  es <- language == "esp"
  tr <- function(a, b) if (es) a else b
  sections <- list()
  add <- function(title, value) {
    if (is.data.frame(value) && nrow(value)) sections[[title]] <<- value
  }
  add(tr("Datos y decisiones", "Data and decisions"), data.frame(
    Information = c("N", tr("\u00cdtems", "Items"), tr("Factores", "Factors"),
      tr("Estimador solicitado / nativo", "Requested / native estimator"),
      tr("Rotaci\u00f3n", "Rotation"), tr("Datos faltantes", "Missing data"),
      tr("Estado", "Status")),
    Value = c(x$data_summary$n, x$data_summary$n_items, x$data_summary$n_factors,
      paste(x$effective_settings$estimator, x$effective_settings$native_estimator, sep = " / "),
      x$effective_settings$rotation, x$effective_settings$missing, x$conclusion$status)))
  settings <- x$effective_settings
  controls <- data.frame(Information = c(tr("Inicios aleatorios", "Random starts"),
    tr("Semilla", "Seed"), "GPA tolerance", tr("Iteraciones m\u00e1ximas", "Maximum iterations")),
    Value = as.character(c(settings$rotation_args$rstarts,
                           if (is.null(settings$seed)) "NULL" else settings$seed,
                           settings$rotation_args$gpa_tol, settings$rotation_args$max_iter)))
  if (settings$rotation == "geomin") controls <- rbind(controls,
    data.frame(Information = tr("Epsilon geomin", "Geomin epsilon"),
               Value = format(settings$rotation_args$geomin_epsilon, scientific = FALSE)))
  first <- names(sections)[1L]
  sections[[first]] <- rbind(sections[[first]], controls)
  rename <- function(tab, labels) {
    if (!nrow(tab)) return(tab)
    names(tab) <- if (es) labels else names(tab)
    tab
  }
  add(tr("\u00cdndices de ajuste", "Fit indices"), rename(x$results$fit_indices,
    c("\u00cdndice", "Valor", "Variante", "Base", "Escalada", "Robusta", "Motivo")))
  ld <- x$results$loadings
  if (nrow(ld)) {
    add(tr("Cargas de patr\u00f3n estandarizadas", "Standardized pattern loadings"),
      rename(ld[, c("factor", "item", "standardized", "se", "p_value", "ci_lower", "ci_upper")],
        c("Factor", "\u00cdtem", "Carga", "EE", "p", "IC 95% inf.", "IC 95% sup.")))
    add(tr("Cargas de patr\u00f3n originales", "Unstandardized pattern loadings"),
      rename(ld[, c("factor", "item", "estimate", "raw_se", "raw_ci_lower", "raw_ci_upper")],
        c("Factor", "\u00cdtem", "Carga", "EE", "IC 95% inf.", "IC 95% sup.")))
  }
  add(tr("Diagn\u00f3sticos por \u00edtem", "Item diagnostics"), rename(x$results$item_diagnostics,
    c("\u00cdtem", "Factor dominante", "Carga dominante", "Segundo factor", "Segunda carga", "Diferencia absoluta",
      "Factor previsto", "Carga prevista", "Previsto dominante")))
  add(tr("Correlaciones entre factores", "Factor correlations"), rename(x$results$factor_correlations,
    c("Factor 1", "Factor 2", "Correlaci\u00f3n", "EE", "IC 95% inferior", "IC 95% superior")))
  add(tr("Residuos de covarianza", "Covariance residuals"), x$results$residuals)
  if (!identical(destination, "console")) {
    add(tr("Preparaci\u00f3n", "Preparation"), rename(x$results$preparation,
      c("\u00cdtem", "Nivel detectado", "Tratamiento", "Invertido", "L\u00edmites usados", "Codificaci\u00f3n de entrada")))
    add(tr("Mapa de \u00edtems", "Item map"), rename(x$results$item_map, c("\u00cdtem", "Alias")))
    add(tr("Mapa de factores", "Factor map"), rename(x$results$factor_map[, c("factor", "label")],
      c("ID", "Etiqueta")))
  }
  if (!is.null(x$results$target)) {
    target <- x$results$target
    rownames(target) <- x$results$item_map$item
    colnames(target) <- x$results$factor_map$label
    add(tr("Patr\u00f3n target", "Target pattern"), .pm_matrix_table(target))
    add(tr("Lectura del patr\u00f3n target", "Reading the target pattern"), data.frame(Text =
      tr("0: carga secundaria orientada hacia cero; 1: carga primaria sin objetivo. Todas se estiman libremente; estos valores no son restricciones ni estimaciones.",
         "0: secondary loading targeted near zero; 1: primary loading left untargeted. All loadings are freely estimated; these values are neither fixed constraints nor estimates.")))
  }
  cmp <- x$results$cfa_comparison
  if (!is.null(cmp)) {
    if (identical(cmp$status, "evaluable")) add(tr("Comparaci\u00f3n CFA descriptiva", "Descriptive CFA comparison"),
      rename(cmp$table, c("\u00cdndice", "Variante", "ESEM", "CFA", "Diferencia", "Motivo"))) else
        add(tr("Comparaci\u00f3n CFA", "CFA comparison"), data.frame(Text = cmp$reason))
  }
  if (nrow(x$warnings)) add(tr("Advertencias", "Warnings"), data.frame(Text = x$warnings$message))
  if (nrow(x$errors)) add(tr("Errores", "Errors"), data.frame(Text = x$errors$message))
  if (!identical(destination, "console")) {
    add(tr("Versiones", "Versions"), data.frame(Information = names(x$versions),
      Value = unname(unlist(x$versions, use.names = FALSE))))
  }
  add(tr("Referencias", "References"), data.frame(Text = x$references$reference))
  add(tr("Interpretaci\u00f3n", "Interpretation"), data.frame(Text = .pm_esem_manuscript(
    x, language, .pm_null_default(attr(x, "display_digits"), 3L))))
  list(title = tr("Ecuaciones estructurales exploratorias (ESEM)", "Exploratory structural equation modeling (ESEM)"),
       sections = sections)
}

#' @export
print.esem_auto <- function(x, digits = NULL, ...) .pm_console_report(x, digits)

#' @export
summary.esem_auto <- function(object, ...) {
  list(settings = object$effective_settings, data = object$data_summary,
       loadings = object$results$loadings, fit_indices = object$results$fit_indices,
       factor_correlations = object$results$factor_correlations,
       item_diagnostics = object$results$item_diagnostics, diagnostics = object$diagnostics,
       conclusion = object$conclusion, manuscript_text = object$manuscript_text)
}
