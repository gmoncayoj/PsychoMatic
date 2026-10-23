# One matrix estimates one association type; never combine incompatible scales.
.cormat_select <- function(data, type, ordinal_levels, language) {
  es <- identical(language, "esp")
  tr <- function(es_text, en_text) if (es) es_text else en_text
  if (type != "auto") {
    return(list(type = type,
      data = .cormat_prepare_data(data, type, ordinal_levels), profile = NULL,
      rationale = tr("Metodo especificado por el usuario; no se aplico la seleccion automatica.",
                     "User-specified method; automatic selection was not applied.")))
  }
  mapped <- .pm_apply_ordinal_levels(data, seq_along(data), names(data), ordinal_levels)$data
  # Logical binary values have a known order; unordered labels still require mapping.
  mapped[] <- lapply(mapped, function(x) if (is.logical(x)) as.integer(x) else x)
  prepared <- .cormat_prepare_data(mapped, "poly")
  .cormat_validate_observed_variation(prepared)
  profile <- do.call(rbind, lapply(seq_along(prepared), function(i) {
    x <- stats::na.omit(prepared[[i]])
    k <- length(unique(x))
    centered <- x - mean(x)
    m2 <- mean(centered^2)
    skew <- mean(centered^3) / m2^1.5
    excess <- mean(centered^4) / m2^2 - 3
    # Ordered factors declare ordinality. Numeric integer scales up to 10
    # observed categories are an explicit operational inference, not proof.
    ordinal <- is.ordered(mapped[[i]]) || (k <= 10L && all(x == trunc(x)))
    data.frame(item = names(prepared)[i], categories = k, ordinal = ordinal,
      skewness = skew, excess_kurtosis = excess,
      endpoint_mass = max(mean(x == min(x)), mean(x == max(x))),
      sparse_category = ordinal && any(table(x) < 5L))
  }))
  if (any(profile$ordinal & profile$categories <= 5L) &&
      (!all(profile$ordinal) || any(profile$categories >= 6L))) {
    stop("Mixed ordinal and continuous items: specify `type` explicitly or analyze homogeneous subsets.", call. = FALSE)
  }
  shape_flag <- any(abs(profile$skewness) > 1 | abs(profile$excess_kurtosis) > 2 |
                      profile$endpoint_mass >= .5)
  type <- if (all(profile$ordinal & profile$categories == 2L)) "tetra" else if (
    all(profile$ordinal) && (any(profile$categories <= 5L) || shape_flag)) "poly" else "pearson"
  reason <- if (type == "tetra") tr("Todos los items son dicotomicos.", "All items are dichotomous.") else if (type == "poly") {
    if (any(profile$categories <= 5L)) tr("Indicadores ordinales con 2-5 categorias presentes.",
      "Ordinal indicators include items with 2-5 categories.") else tr(
      "Indicadores ordinales con asimetria, curtosis o concentracion extrema; se conserva el tratamiento ordinal.",
      "Ordinal indicators show skewness, kurtosis or endpoint concentration; ordinal treatment is retained.")
  } else if (all(profile$ordinal)) tr(
    "Seis o mas categorias sin alertas de forma; se aproximan los codigos a puntuaciones continuas.",
    "Six or more categories without shape flags; numeric codes are approximated as continuous scores.") else tr(
    "Puntuaciones numericas tratadas como continuas; Pearson estima asociacion lineal, sin exigir normalidad para calcular el coeficiente.",
    "Numeric scores treated as continuous; Pearson estimates linear association without requiring normality to compute the coefficient.")
  reason <- paste(reason, sprintf(tr(
    "Max |asimetria|=%.2f, max |curtosis excedente|=%.2f, max proporcion extrema=%.2f. Limites operativos: 1, 2 y 0.50; no son pruebas de normalidad.",
    "Max |skewness|=%.2f, max |excess kurtosis|=%.2f, max endpoint proportion=%.2f. Operational flags: 1, 2 and 0.50; these are not normality tests."),
    max(abs(profile$skewness)), max(abs(profile$excess_kurtosis)), max(profile$endpoint_mass)))
  if (type %in% c("tetra", "poly")) reason <- paste(reason, tr(
    "Se asumen respuestas latentes continuas con normalidad bivariada; la forma observada no demuestra ese supuesto.",
    "Continuous latent responses with bivariate normality are assumed; observed shape does not establish this assumption."))
  if (any(profile$sparse_category)) reason <- paste(reason, tr(
    "Hay categorias con menos de 5 respuestas; interprete con cautela la estabilidad de la estimacion.",
    "Some categories have fewer than 5 responses; interpret estimation stability cautiously."))
  if (type == "pearson" && shape_flag) reason <- paste(reason, tr(
    "Hay alertas de forma: revise valores extremos y linealidad; no se convierten datos continuos en ordinales automaticamente.",
    "Shape flags are present: inspect outliers and linearity; continuous data are not automatically converted to ordinal data."))
  list(type = type, data = prepared, profile = profile, rationale = reason)
}

.cormat_references <- function(type) {
  c(if (type == "pearson")
    "Pearson, K. (1895). Note on regression and inheritance in the case of two parents. Proceedings of the Royal Society of London, 58, 240-242. https://doi.org/10.1098/rspl.1895.0041" else
    "Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation coefficient. Psychometrika, 44, 443-460. https://doi.org/10.1007/BF02296207",
    "Rhemtulla, M., Brosseau-Liard, P. E., & Savalei, V. (2012). When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions. Psychological Methods, 17(3), 354-373. https://doi.org/10.1037/a0029315")
}
