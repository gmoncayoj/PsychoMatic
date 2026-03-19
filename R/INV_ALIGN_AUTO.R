# Automatic alignment invariance -----------------------------------------------
#
# Main usage:
#   res <- inv_align_auto(data = my_items, group = my_group, language = "esp")
#   print(res)
#
# Notes:
# - If `data` and `group` are supplied, the function estimates a one-factor
#   configural model with `sirt::invariance_alignment_cfa_config()`.
# - If `lambda` and `nu` are supplied, they are used directly in
#   `sirt::invariance.alignment()`.
# - Requested defaults:
#     align.scale = c(lambda = 0.40, nu = 0.20)
#     align.pow   = c(lambda = 0.25, nu = 0.25)
# - A scale is flagged as noninvariant when more than 25% of lambda or nu
#   parameters are classified as noninvariant.

normalize_alignment_language <- function(language = "esp") {
  if (missing(language) || is.null(language) || length(language) == 0L) {
    return("esp")
  }

  language <- tolower(as.character(language[[1L]]))

  if (language %in% c("esp", "es", "spa", "spanish", "espanol")) {
    return("esp")
  }

  if (language %in% c("eng", "en", "english")) {
    return("eng")
  }

  stop("`language` must be 'esp' or 'eng'.", call. = FALSE)
}

normalize_alignment_pair <- function(x, arg_name, default = NULL) {
  if (is.null(x)) {
    x <- default
  }

  if (!is.numeric(x) || length(x) != 2L || anyNA(x)) {
    stop(
      sprintf("`%s` must be a numeric vector of length 2.", arg_name),
      call. = FALSE
    )
  }

  x <- as.numeric(x)
  names(x) <- c("lambda", "nu")
  x
}

coerce_numeric_matrix <- function(x, arg_name) {
  x <- as.matrix(x)

  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", arg_name), call. = FALSE)
  }

  storage.mode(x) <- "double"

  if (length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a 2D matrix-like object.", arg_name), call. = FALSE)
  }

  x
}

validate_parameter_matrices <- function(lambda, nu) {
  lambda <- coerce_numeric_matrix(lambda, "lambda")
  nu <- coerce_numeric_matrix(nu, "nu")

  if (!identical(dim(lambda), dim(nu))) {
    stop("`lambda` and `nu` must have the same dimensions.", call. = FALSE)
  }

  lambda_rows <- rownames(lambda)
  nu_rows <- rownames(nu)
  lambda_cols <- colnames(lambda)
  nu_cols <- colnames(nu)

  if (!is.null(lambda_rows) && !is.null(nu_rows) && !identical(lambda_rows, nu_rows)) {
    stop("`lambda` and `nu` must use the same row names.", call. = FALSE)
  }

  if (!is.null(lambda_cols) && !is.null(nu_cols) && !identical(lambda_cols, nu_cols)) {
    stop("`lambda` and `nu` must use the same column names.", call. = FALSE)
  }

  if (is.null(lambda_rows) && is.null(nu_rows)) {
    lambda_rows <- nu_rows <- paste0("Group", seq_len(nrow(lambda)))
  } else if (is.null(lambda_rows)) {
    lambda_rows <- nu_rows
  } else if (is.null(nu_rows)) {
    nu_rows <- lambda_rows
  }

  if (is.null(lambda_cols) && is.null(nu_cols)) {
    lambda_cols <- nu_cols <- paste0("Item", seq_len(ncol(lambda)))
  } else if (is.null(lambda_cols)) {
    lambda_cols <- nu_cols
  } else if (is.null(nu_cols)) {
    nu_cols <- lambda_cols
  }

  rownames(lambda) <- rownames(nu) <- lambda_rows
  colnames(lambda) <- colnames(nu) <- lambda_cols

  list(lambda = lambda, nu = nu)
}

normalize_alignment_weights <- function(wgt, lambda) {
  if (is.null(wgt)) {
    return(NULL)
  }

  if (is.numeric(wgt) && is.null(dim(wgt)) && length(wgt) == nrow(lambda)) {
    wgt <- matrix(wgt, nrow = nrow(lambda), ncol = ncol(lambda))
    rownames(wgt) <- rownames(lambda)
    colnames(wgt) <- colnames(lambda)
    return(wgt)
  }

  wgt <- coerce_numeric_matrix(wgt, "wgt")

  if (!identical(dim(wgt), dim(lambda))) {
    stop(
      "`wgt` must either be a numeric vector with one value per group or a matrix with the same dimensions as `lambda` and `nu`.",
      call. = FALSE
    )
  }

  if (is.null(rownames(wgt))) {
    rownames(wgt) <- rownames(lambda)
  }

  if (is.null(colnames(wgt))) {
    colnames(wgt) <- colnames(lambda)
  }

  wgt
}

require_alignment_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required but not installed. Install it before running `inv_align_auto()`.",
        package_name
      ),
      call. = FALSE
    )
  }
}

format_alignment_number <- function(x, digits = 3) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

r2_proximity <- function(value) {
  if (is.na(value)) {
    return("unknown")
  }

  distance_to_one <- abs(1 - value)
  distance_to_zero <- abs(value)

  if (distance_to_one < distance_to_zero) {
    return("closer_to_1")
  }

  if (distance_to_one > distance_to_zero) {
    return("closer_to_0")
  }

  "midpoint"
}

r2_sentence <- function(value, language) {
  flag <- r2_proximity(value)

  if (language == "esp") {
    if (flag == "closer_to_1") {
      return("mas cercano a 1, lo que sugiere mayor invarianza")
    }

    if (flag == "closer_to_0") {
      return("mas cercano a 0, lo que sugiere menor invarianza")
    }

    return("equidistante entre 0 y 1, por lo que requiere cautela interpretativa")
  }

  if (flag == "closer_to_1") {
    return("closer to 1, suggesting greater invariance")
  }

  if (flag == "closer_to_0") {
    return("closer to 0, suggesting lower invariance")
  }

  "equidistant between 0 and 1, so the interpretation should be cautious"
}

cutoff_sentence <- function(parameter_invariant, cutoff, language) {
  if (language == "esp") {
    if (isTRUE(parameter_invariant)) {
      return(
        sprintf(
          "Como el porcentaje no supera %s%%, este conjunto de parametros se considera invariante.",
          format_alignment_number(cutoff, digits = 1)
        )
      )
    }

    return(
      sprintf(
        "Como el porcentaje supera %s%%, este conjunto de parametros se considera no invariante.",
        format_alignment_number(cutoff, digits = 1)
      )
    )
  }

  if (isTRUE(parameter_invariant)) {
    return(
      sprintf(
        "Because the percentage does not exceed %s%%, this parameter set is treated as invariant.",
        format_alignment_number(cutoff, digits = 1)
      )
    )
  }

  sprintf(
    "Because the percentage exceeds %s%%, this parameter set is treated as noninvariant.",
    format_alignment_number(cutoff, digits = 1)
  )
}

language_pack_alignment <- function(language) {
  if (language == "esp") {
    return(
      list(
        title = "Reporte de invarianza por alineamiento",
        input_mode = "Modo de entrada",
        input_matrices = "matrices lambda/nu preestimadas",
        input_data = "datos + grupo con CFA configural estimado con sirt",
        groups = "Numero de grupos",
        items = "Numero de items",
        model = "Modelo configural",
        settings = "Configuracion",
        parameter_line = "%s: R2 = %s (%s). Parametros no invariantes = %s%%. %s",
        loadings = "Cargas factoriales",
        intercepts = "Interceptos",
        item_results_loadings = "Resultados por item: cargas factoriales",
        item_results_intercepts = "Resultados por item: interceptos",
        item_summary = "Resumen por item",
        item_matrix = "Matriz parm_dif por grupo e item (0 = invariante dentro de la tolerancia)",
        item_col_item = "item",
        item_col_joint = "parametro_conjunto",
        item_col_distinct = "parametros_distintos",
        item_col_noninv_groups = "grupos_no_invariantes",
        item_col_pct = "porcentaje_no_invariante",
        item_col_flagged = "grupos_identificados",
        item_col_invariant = "item_invariante",
        none_label = "ninguno",
        yes_label = "si",
        no_label = "no",
        references = "Referencias APA 7",
        overall_invariant = "Conclusion: la escala cumple el criterio porcentual de invarianza por alineamiento.",
        overall_noninvariant = "Conclusion: la escala no cumple el criterio porcentual de invarianza por alineamiento.",
        note = "Nota: valores de R2 cercanos a 1 indican mayor invarianza, mientras que valores cercanos a 0 indican menor invarianza."
      )
    )
  }

  list(
    title = "Alignment invariance report",
    input_mode = "Input mode",
    input_matrices = "pre-estimated lambda/nu matrices",
    input_data = "data + group with a configural CFA estimated in sirt",
    groups = "Number of groups",
    items = "Number of items",
    model = "Configural model",
    settings = "Settings",
    parameter_line = "%s: R2 = %s (%s). Noninvariant parameters = %s%%. %s",
    loadings = "Factor loadings",
    intercepts = "Intercepts",
    item_results_loadings = "Item-level results: factor loadings",
    item_results_intercepts = "Item-level results: intercepts",
    item_summary = "Item summary",
    item_matrix = "parm_dif matrix by group and item (0 = invariant within the selected tolerance)",
    item_col_item = "item",
    item_col_joint = "joint_parameter",
    item_col_distinct = "distinct_parameters",
    item_col_noninv_groups = "noninvariant_groups",
    item_col_pct = "percent_noninvariant",
    item_col_flagged = "flagged_groups",
    item_col_invariant = "item_invariant",
    none_label = "none",
    yes_label = "yes",
    no_label = "no",
    references = "APA 7 references",
    overall_invariant = "Conclusion: the scale meets the percentage-based alignment invariance criterion.",
    overall_noninvariant = "Conclusion: the scale does not meet the percentage-based alignment invariance criterion.",
    note = "Note: R2 values close to 1 indicate greater invariance, whereas values close to 0 indicate lower invariance."
  )
}

build_alignment_summary <- function(alignment, constraints, noninvariance_cutoff) {
  es <- alignment$es.invariance

  if (is.null(es) || is.null(dim(es)) || !"R2" %in% rownames(es)) {
    stop("The alignment object does not contain an `es.invariance` R2 summary.", call. = FALSE)
  }

  loading_col <- if ("loadings" %in% colnames(es)) "loadings" else colnames(es)[1L]
  intercept_col <- if ("intercepts" %in% colnames(es)) "intercepts" else colnames(es)[2L]

  loadings_r2 <- as.numeric(es["R2", loading_col])
  intercepts_r2 <- as.numeric(es["R2", intercept_col])

  loadings_pct <- as.numeric(constraints$lambda_list$prop_noninvariance)
  intercepts_pct <- as.numeric(constraints$nu_list$prop_noninvariance)

  data.frame(
    parameter = c("loadings", "intercepts"),
    R2 = c(loadings_r2, intercepts_r2),
    percent_noninvariant = c(loadings_pct, intercepts_pct),
    cutoff = rep(noninvariance_cutoff, 2L),
    invariant = c(
      loadings_pct <= noninvariance_cutoff,
      intercepts_pct <= noninvariance_cutoff
    ),
    stringsAsFactors = FALSE
  )
}

build_alignment_item_results <- function(parameter_list) {
  parm_dif <- coerce_numeric_matrix(parameter_list$parm_dif, "parameter_list$parm_dif")
  parm_est <- coerce_numeric_matrix(parameter_list$parm_est, "parameter_list$parm_est")
  parm_joint <- as.numeric(parameter_list$parm_joint)

  group_names <- rownames(parm_dif)
  item_names <- colnames(parm_dif)

  if (is.null(group_names)) {
    group_names <- paste0("Group", seq_len(nrow(parm_dif)))
  }

  if (is.null(item_names)) {
    item_names <- paste0("Item", seq_len(ncol(parm_dif)))
  }

  rownames(parm_dif) <- rownames(parm_est) <- group_names
  colnames(parm_dif) <- colnames(parm_est) <- item_names

  if (length(parm_joint) != length(item_names)) {
    stop("The parameter summary does not match the number of items.", call. = FALSE)
  }

  names(parm_joint) <- item_names

  eps <- sqrt(.Machine$double.eps)
  noninvariant_matrix <- !is.na(parm_dif) & abs(parm_dif) > eps
  noninvariant_groups <- colSums(noninvariant_matrix, na.rm = TRUE)
  total_groups <- colSums(!is.na(parm_dif), na.rm = TRUE)
  percent_noninvariant <- ifelse(
    total_groups > 0,
    100 * noninvariant_groups / total_groups,
    NA_real_
  )

  flagged_groups <- vapply(
    seq_len(ncol(noninvariant_matrix)),
    function(j) {
      groups <- group_names[noninvariant_matrix[, j]]
      if (length(groups) == 0L) "" else paste(groups, collapse = ", ")
    },
    character(1L)
  )
  names(flagged_groups) <- item_names

  distinct_parameters <- as.numeric(parameter_list$N_parm_items)
  if (length(distinct_parameters) == length(item_names)) {
    names(distinct_parameters) <- item_names
  } else {
    distinct_parameters <- rep(NA_real_, length(item_names))
    names(distinct_parameters) <- item_names
  }

  summary_table <- data.frame(
    item = item_names,
    joint_parameter = as.numeric(parm_joint[item_names]),
    distinct_parameters = distinct_parameters[item_names],
    noninvariant_groups = as.numeric(noninvariant_groups[item_names]),
    percent_noninvariant = as.numeric(percent_noninvariant[item_names]),
    flagged_groups = flagged_groups[item_names],
    invariant = as.numeric(noninvariant_groups[item_names]) == 0,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_table,
    parm_dif = parm_dif,
    parm_est = parm_est,
    parm_joint = parm_joint,
    noninvariant_matrix = noninvariant_matrix
  )
}

format_alignment_table_for_report <- function(x, digits = 3) {
  numeric_cols <- vapply(x, is.numeric, logical(1L))
  if (any(numeric_cols)) {
    x[numeric_cols] <- lapply(x[numeric_cols], round, digits = digits)
  }
  out <- capture.output(print(x, row.names = FALSE))
  paste(out, collapse = "\n")
}

format_alignment_matrix_for_report <- function(x, digits = 3) {
  x <- round(as.matrix(x), digits = digits)
  x[!is.na(x) & abs(x) < 10^(-digits)] <- 0
  out <- capture.output(print(x))
  paste(out, collapse = "\n")
}

alignment_apa_references <- function() {
  c(
    "1. Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis alignment. Structural Equation Modeling: A Multidisciplinary Journal, 21(4), 495-508. https://doi.org/10.1080/10705511.2014.919210",
    "2. Marsh, H. W., Guo, J., Parker, P. D., Nagengast, B., Asparouhov, T., Muthen, B., & Dicke, T. (2018). What to do when scalar invariance fails: The extended alignment method for multi-group factor analysis comparison of latent means across many groups. Psychological Methods, 23(3), 524-545. https://doi.org/10.1037/met0000113",
    "3. Robitzsch, A., & Ludtke, O. (2023). Why full, partial, or approximate measurement invariance are not a prerequisite for meaningful and valid group comparisons. Structural Equation Modeling: A Multidisciplinary Journal, 30(6), 859-870. https://doi.org/10.1080/10705511.2023.2191292",
    "4. Fischer, R., & Karl, J. A. (2019). A primer to (cross-cultural) multi-group invariance testing possibilities in R. Frontiers in Psychology, 10, Article 1507. https://doi.org/10.3389/fpsyg.2019.01507",
    "5. Muthen, B., & Asparouhov, T. (2014). IRT studies of many groups: The alignment method. Frontiers in Psychology, 5, Article 978. https://doi.org/10.3389/fpsyg.2014.00978"
  )
}

prepare_item_summary_for_report <- function(item_summary, language) {
  text <- language_pack_alignment(language)

  flagged_groups <- ifelse(
    nzchar(item_summary$flagged_groups),
    item_summary$flagged_groups,
    text$none_label
  )

  invariant_label <- ifelse(item_summary$invariant, text$yes_label, text$no_label)

  out <- data.frame(
    item_summary$item,
    item_summary$joint_parameter,
    item_summary$distinct_parameters,
    item_summary$noninvariant_groups,
    item_summary$percent_noninvariant,
    flagged_groups,
    invariant_label,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  names(out) <- c(
    text$item_col_item,
    text$item_col_joint,
    text$item_col_distinct,
    text$item_col_noninv_groups,
    text$item_col_pct,
    text$item_col_flagged,
    text$item_col_invariant
  )

  out
}

build_alignment_report <- function(summary_table, item_results, language, settings, dimensions) {
  text <- language_pack_alignment(language)

  loadings_row <- summary_table[summary_table$parameter == "loadings", , drop = FALSE]
  intercepts_row <- summary_table[summary_table$parameter == "intercepts", , drop = FALSE]

  loadings_line <- sprintf(
    text$parameter_line,
    text$loadings,
    format_alignment_number(loadings_row$R2, settings$digits),
    r2_sentence(loadings_row$R2, language),
    format_alignment_number(loadings_row$percent_noninvariant, settings$digits),
    cutoff_sentence(loadings_row$invariant, loadings_row$cutoff, language)
  )

  intercepts_line <- sprintf(
    text$parameter_line,
    text$intercepts,
    format_alignment_number(intercepts_row$R2, settings$digits),
    r2_sentence(intercepts_row$R2, language),
    format_alignment_number(intercepts_row$percent_noninvariant, settings$digits),
    cutoff_sentence(intercepts_row$invariant, intercepts_row$cutoff, language)
  )

  settings_line <- paste0(
    text$settings, ": ",
    "align.scale = (lambda = ", format_alignment_number(settings$align.scale["lambda"], settings$digits),
    ", nu = ", format_alignment_number(settings$align.scale["nu"], settings$digits), "); ",
    "align.pow = (lambda = ", format_alignment_number(settings$align.pow["lambda"], settings$digits),
    ", nu = ", format_alignment_number(settings$align.pow["nu"], settings$digits), "); ",
    "parm_tol = (lambda = ", format_alignment_number(settings$parm_tol["lambda"], settings$digits),
    ", nu = ", format_alignment_number(settings$parm_tol["nu"], settings$digits), "); ",
    "cutoff = ", format_alignment_number(settings$noninvariance_cutoff, 1), "%."
  )

  input_text <- if (identical(settings$input_mode, "matrices")) {
    text$input_matrices
  } else {
    text$input_data
  }

  header_lines <- c(
    text$title,
    paste0(text$input_mode, ": ", input_text),
    paste0(text$groups, ": ", dimensions$groups),
    paste0(text$items, ": ", dimensions$items)
  )

  if (!is.null(settings$model)) {
    header_lines <- c(header_lines, paste0(text$model, ": ", settings$model))
  }

  overall_line <- if (all(summary_table$invariant)) {
    text$overall_invariant
  } else {
    text$overall_noninvariant
  }

  loadings_item_table <- format_alignment_table_for_report(
    prepare_item_summary_for_report(item_results$loadings$summary, language = language),
    digits = settings$digits
  )
  intercepts_item_table <- format_alignment_table_for_report(
    prepare_item_summary_for_report(item_results$intercepts$summary, language = language),
    digits = settings$digits
  )

  loadings_item_matrix <- format_alignment_matrix_for_report(
    item_results$loadings$parm_dif,
    digits = settings$digits
  )
  intercepts_item_matrix <- format_alignment_matrix_for_report(
    item_results$intercepts$parm_dif,
    digits = settings$digits
  )

  references_block <- paste(alignment_apa_references(), collapse = "\n")

  paste(
    c(
      header_lines,
      settings_line,
      "",
      loadings_line,
      intercepts_line,
      "",
      overall_line,
      text$note,
      "",
      text$item_results_loadings,
      text$item_summary,
      loadings_item_table,
      "",
      text$item_matrix,
      loadings_item_matrix,
      "",
      text$item_results_intercepts,
      text$item_summary,
      intercepts_item_table,
      "",
      text$item_matrix,
      intercepts_item_matrix,
      "",
      text$references,
      references_block
    ),
    collapse = "\n"
  )
}

#' Invarianza por alineamiento automatizada
#'
#' Ejecuta un analisis de alineamiento a partir de matrices `lambda` y `nu`,
#' o directamente desde datos y un vector de grupos.
#'
#' @param lambda,nu Matrices de parametros para cargas e interceptos.
#' @param data Datos de entrada cuando no se proporcionan `lambda` y `nu`.
#' @param group Vector de grupos asociado a `data`.
#' @param sampling_weights Pesos muestrales opcionales por fila.
#' @param wgt Pesos opcionales para el procedimiento de alineamiento.
#' @param language Idioma del reporte: `"esp"` o `"eng"`.
#' @param model Tipo de modelo para el ajuste configural.
#' @param align.scale,align.pow,parm_tol Parametros de ajuste del alineamiento.
#' @param noninvariance_cutoff Porcentaje a partir del cual se marca no invarianza.
#' @param center,optimizer,fixed,meth,eps,vcov Parametros avanzados del algoritmo.
#' @param digits Numero de decimales del reporte.
#' @param config_args Lista de argumentos adicionales para la etapa configural.
#' @param alignment_args Lista de argumentos adicionales para la etapa de alineamiento.
#'
#' @return Objeto de clase `inv_align_auto_result`.
#'
#' @examples
#' # res <- inv_align_auto(data = mis_items, group = mi_grupo)
#' # summary(res)
#'
#' @export
inv_align_auto <- function(
  lambda = NULL,
  nu = NULL,
  data = NULL,
  group = NULL,
  sampling_weights = NULL,
  wgt = NULL,
  language = c("esp", "eng"),
  model = c("2PM", "1PM"),
  align.scale = c(lambda = 0.40, nu = 0.20),
  align.pow = c(lambda = 0.25, nu = 0.25),
  parm_tol = NULL,
  noninvariance_cutoff = 25,
  center = FALSE,
  optimizer = "optim",
  fixed = NULL,
  meth = 1,
  eps = 1e-03,
  vcov = NULL,
  digits = 3,
  config_args = list(),
  alignment_args = list()
) {
  language <- normalize_alignment_language(language)
  model <- match.arg(model)
  require_alignment_package("sirt")

  align.scale <- normalize_alignment_pair(align.scale, "align.scale")
  align.pow <- normalize_alignment_pair(align.pow, "align.pow")

  if (is.null(parm_tol)) {
    parm_tol <- align.scale
  }
  parm_tol <- normalize_alignment_pair(parm_tol, "parm_tol")

  if (!is.numeric(noninvariance_cutoff) ||
      length(noninvariance_cutoff) != 1L ||
      is.na(noninvariance_cutoff) ||
      noninvariance_cutoff < 0 ||
      noninvariance_cutoff > 100) {
    stop("`noninvariance_cutoff` must be a single numeric value between 0 and 100.", call. = FALSE)
  }

  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0) {
    stop("`digits` must be a single non-negative numeric value.", call. = FALSE)
  }
  digits <- as.integer(digits)

  if (!is.list(config_args)) {
    stop("`config_args` must be a list.", call. = FALSE)
  }

  if (!is.list(alignment_args)) {
    stop("`alignment_args` must be a list.", call. = FALSE)
  }

  input_mode <- if (!is.null(lambda) || !is.null(nu)) "matrices" else "data"
  configural <- NULL

  if (identical(input_mode, "matrices")) {
    if (is.null(lambda) || is.null(nu)) {
      stop("Provide both `lambda` and `nu`, or provide both `data` and `group`.", call. = FALSE)
    }

    if (!is.null(data) || !is.null(group) || !is.null(sampling_weights)) {
      warning(
        "`data`, `group`, and `sampling_weights` are ignored because `lambda` and `nu` were provided.",
        call. = FALSE
      )
    }

    matrices <- validate_parameter_matrices(lambda, nu)
    lambda <- matrices$lambda
    nu <- matrices$nu
  } else {
    if (is.null(data) || is.null(group)) {
      stop("Provide both `lambda` and `nu`, or provide both `data` and `group`.", call. = FALSE)
    }

    data <- as.data.frame(data)

    if (!all(vapply(data, is.numeric, logical(1L)))) {
      stop("`data` must contain only numeric item columns.", call. = FALSE)
    }

    if (nrow(data) != length(group)) {
      stop("`data` and `group` must contain the same number of rows.", call. = FALSE)
    }

    if (anyNA(group)) {
      stop("`group` cannot contain missing values.", call. = FALSE)
    }

    if (ncol(data) < 2L) {
      stop("`data` must contain at least two item columns.", call. = FALSE)
    }

    if (length(unique(group)) < 2L) {
      stop("`group` must contain at least two groups.", call. = FALSE)
    }

    if (!is.null(sampling_weights)) {
      if (!is.numeric(sampling_weights) || length(sampling_weights) != nrow(data)) {
        stop(
          "`sampling_weights` must be a numeric vector with one value per row in `data`.",
          call. = FALSE
        )
      }
    }

    config_defaults <- list(
      dat = data,
      group = group,
      weights = sampling_weights,
      model = model,
      verbose = FALSE
    )
    config_call <- utils::modifyList(config_defaults, config_args, keep.null = TRUE)
    configural <- do.call(sirt::invariance_alignment_cfa_config, config_call)

    matrices <- validate_parameter_matrices(configural$lambda, configural$nu)
    lambda <- matrices$lambda
    nu <- matrices$nu

    if (is.null(vcov) && !is.null(configural$vcov)) {
      vcov <- configural$vcov
    }
  }

  wgt <- normalize_alignment_weights(wgt, lambda)

  alignment_defaults <- list(
    lambda = lambda,
    nu = nu,
    wgt = wgt,
    align.scale = align.scale,
    align.pow = align.pow,
    eps = eps,
    center = center,
    optimizer = optimizer,
    fixed = fixed,
    meth = meth,
    vcov = vcov
  )
  alignment_call <- utils::modifyList(alignment_defaults, alignment_args, keep.null = TRUE)
  alignment <- do.call(sirt::invariance.alignment, alignment_call)

  constraints <- sirt::invariance_alignment_constraints(
    model = alignment,
    lambda_parm_tol = parm_tol["lambda"],
    nu_parm_tol = parm_tol["nu"]
  )

  summary_table <- build_alignment_summary(
    alignment = alignment,
    constraints = constraints,
    noninvariance_cutoff = noninvariance_cutoff
  )

  item_results <- list(
    loadings = build_alignment_item_results(constraints$lambda_list),
    intercepts = build_alignment_item_results(constraints$nu_list)
  )

  settings <- list(
    input_mode = input_mode,
    model = if (identical(input_mode, "data")) model else NULL,
    language = language,
    align.scale = align.scale,
    align.pow = align.pow,
    parm_tol = parm_tol,
    noninvariance_cutoff = noninvariance_cutoff,
    digits = digits
  )

  dimensions <- list(
    groups = nrow(lambda),
    items = ncol(lambda)
  )

  report <- build_alignment_report(
    summary_table = summary_table,
    item_results = item_results,
    language = language,
    settings = settings,
    dimensions = dimensions
  )

  result <- list(
    call = match.call(),
    input_mode = input_mode,
    settings = settings,
    dimensions = dimensions,
    configural = configural,
    alignment = alignment,
    constraints = constraints,
    summary = summary_table,
    item_results = item_results,
    group_parameters = alignment$pars,
    aligned_item_parameters = alignment$itempars.aligned,
    noninvariance = list(
      loadings = constraints$lambda_list,
      intercepts = constraints$nu_list
    ),
    references = alignment_apa_references(),
    report = report
  )

  class(result) <- "inv_align_auto_result"
  result
}

#' Imprime un objeto `inv_align_auto_result`
#'
#' @param x Objeto devuelto por `inv_align_auto()`.
#' @param ... Argumentos adicionales no usados.
#'
#' @rdname inv_align_auto
#' @method print inv_align_auto_result
#' @export
print.inv_align_auto_result <- function(x, ...) {
  cat(x$report, sep = "\n")
  invisible(x)
}

#' Resume un objeto `inv_align_auto_result`
#'
#' @param object Objeto devuelto por `inv_align_auto()`.
#' @param ... Argumentos adicionales no usados.
#'
#' @rdname inv_align_auto
#' @method summary inv_align_auto_result
#' @rawNamespace S3method(summary,inv_align_auto_result)
#' @export
summary.inv_align_auto_result <- function(object, ...) {
  object$summary
}
