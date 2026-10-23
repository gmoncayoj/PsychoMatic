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

normalize_alignment_language <- function(language = "eng") {
  if (missing(language) || is.null(language) || length(language) == 0L) {
    return("eng")
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

normalize_alignment_pair <- function(x,
                                     arg_name,
                                     default = NULL,
                                     allow_zero = FALSE) {
  if (is.null(x)) {
    x <- default
  }

  if (!is.numeric(x) || length(x) != 2L || anyNA(x) || any(!is.finite(x))) {
    stop(
      sprintf("`%s` must be a finite numeric vector of length 2.", arg_name),
      call. = FALSE
    )
  }

  component_names <- names(x)
  if (!is.null(component_names)) {
    if (length(component_names) != 2L || anyNA(component_names) ||
        any(!nzchar(component_names)) || anyDuplicated(component_names) ||
        !setequal(component_names, c("lambda", "nu"))) {
      stop(
        sprintf(
          "Named `%s` values must use the unique components 'lambda' and 'nu'.",
          arg_name
        ),
        call. = FALSE
      )
    }
    x <- x[c("lambda", "nu")]
  } else {
    names(x) <- c("lambda", "nu")
  }

  if (isTRUE(allow_zero)) {
    invalid <- x < 0
  } else {
    invalid <- x <= 0
  }
  if (any(invalid)) {
    qualifier <- if (isTRUE(allow_zero)) "non-negative" else "positive"
    stop(sprintf("`%s` values must be %s.", arg_name, qualifier), call. = FALSE)
  }

  stats::setNames(as.numeric(x), c("lambda", "nu"))
}

coerce_numeric_matrix <- function(x, arg_name) {
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a 2D matrix-like object.", arg_name), call. = FALSE)
  }
  x <- as.matrix(x)

  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", arg_name), call. = FALSE)
  }

  storage.mode(x) <- "double"

  if (length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a 2D matrix-like object.", arg_name), call. = FALSE)
  }

  if (any(dim(x) == 0L)) {
    stop(sprintf("`%s` cannot have an empty dimension.", arg_name), call. = FALSE)
  }

  x
}

validate_parameter_matrices <- function(lambda, nu) {
  lambda <- coerce_numeric_matrix(lambda, "lambda")
  nu <- coerce_numeric_matrix(nu, "nu")

  if (!identical(dim(lambda), dim(nu))) {
    stop("`lambda` and `nu` must have the same dimensions.", call. = FALSE)
  }

  if (nrow(lambda) < 2L || ncol(lambda) < 2L) {
    stop(
      "`lambda` and `nu` must describe at least two groups and two items.",
      call. = FALSE
    )
  }

  if (any(!is.finite(lambda)) || any(!is.finite(nu))) {
    stop("`lambda` and `nu` must contain only finite values.", call. = FALSE)
  }

  lambda_rows <- rownames(lambda)
  nu_rows <- rownames(nu)
  lambda_cols <- colnames(lambda)
  nu_cols <- colnames(nu)

  lambda_named <- !is.null(lambda_rows) || !is.null(lambda_cols)
  nu_named <- !is.null(nu_rows) || !is.null(nu_cols)
  if (xor(lambda_named, nu_named)) {
    stop(
      "A named parameter matrix cannot be paired with an unnamed matrix.",
      call. = FALSE
    )
  }

  if (!lambda_named && !nu_named) {
    group_names <- paste0("G", seq_len(nrow(lambda)))
    item_names <- paste0("item", seq_len(ncol(lambda)))
    rownames(lambda) <- rownames(nu) <- group_names
    colnames(lambda) <- colnames(nu) <- item_names
    return(list(
      lambda = lambda,
      nu = nu,
      names_generated = TRUE,
      identity_matching = "generated"
    ))
  }

  complete_identity <- function(rows, columns, matrix_name) {
    if (is.null(rows) || is.null(columns)) {
      stop(
        sprintf("`%s` must provide complete row and column names.", matrix_name),
        call. = FALSE
      )
    }
    if (anyNA(rows) || anyNA(columns) || any(!nzchar(rows)) ||
        any(!nzchar(columns))) {
      stop(
        sprintf("`%s` row and column names must be nonmissing and nonempty.", matrix_name),
        call. = FALSE
      )
    }
    if (anyDuplicated(rows) || anyDuplicated(columns)) {
      stop(
        sprintf("`%s` row and column names must be unique.", matrix_name),
        call. = FALSE
      )
    }
    invisible(NULL)
  }
  complete_identity(lambda_rows, lambda_cols, "lambda")
  complete_identity(nu_rows, nu_cols, "nu")

  if (!setequal(lambda_rows, nu_rows) || !setequal(lambda_cols, nu_cols)) {
    stop(
      "`lambda` and `nu` must identify the same group and item name sets.",
      call. = FALSE
    )
  }

  nu <- nu[lambda_rows, lambda_cols, drop = FALSE]
  list(
    lambda = lambda,
    nu = nu,
    names_generated = FALSE,
    identity_matching = "named"
  )
}

normalize_alignment_weights <- function(wgt, lambda) {
  if (is.null(wgt)) {
    return(list(value = NULL, matching = "none"))
  }

  if (is.numeric(wgt) && is.null(dim(wgt))) {
    if (length(wgt) != nrow(lambda)) {
      stop("Vector `wgt` must contain exactly one value per group.", call. = FALSE)
    }
    if (anyNA(wgt) || any(!is.finite(wgt)) || any(wgt < 0) || all(wgt == 0)) {
      stop("`wgt` must be finite, non-negative, and not all zero.", call. = FALSE)
    }

    weight_names <- names(wgt)
    matching <- "positional"
    if (!is.null(weight_names)) {
      if (anyNA(weight_names) || any(!nzchar(weight_names)) ||
          anyDuplicated(weight_names) ||
          !setequal(weight_names, rownames(lambda))) {
        stop("Named vector `wgt` must uniquely identify every group.", call. = FALSE)
      }
      wgt <- wgt[rownames(lambda)]
      matching <- "named"
    }
    matrix_value <- matrix(
      rep(as.numeric(wgt), ncol(lambda)),
      nrow = nrow(lambda),
      dimnames = dimnames(lambda)
    )
    return(list(value = matrix_value, matching = matching))
  }

  wgt <- coerce_numeric_matrix(wgt, "wgt")

  if (!identical(dim(wgt), dim(lambda))) {
    stop(
      "`wgt` must either be a numeric vector with one value per group or a matrix with the same dimensions as `lambda` and `nu`.",
      call. = FALSE
    )
  }

  if (anyNA(wgt) || any(!is.finite(wgt)) || any(wgt < 0) || all(wgt == 0)) {
    stop("`wgt` must be finite, non-negative, and not all zero.", call. = FALSE)
  }

  has_rows <- !is.null(rownames(wgt))
  has_columns <- !is.null(colnames(wgt))
  if (xor(has_rows, has_columns)) {
    stop("Matrix `wgt` must be wholly named or wholly unnamed.", call. = FALSE)
  }

  if (!has_rows) {
    rownames(wgt) <- rownames(lambda)
    colnames(wgt) <- colnames(lambda)
    return(list(value = wgt, matching = "positional"))
  }

  rows <- rownames(wgt)
  columns <- colnames(wgt)
  if (anyNA(rows) || anyNA(columns) || any(!nzchar(rows)) ||
      any(!nzchar(columns)) || anyDuplicated(rows) || anyDuplicated(columns) ||
      !setequal(rows, rownames(lambda)) ||
      !setequal(columns, colnames(lambda))) {
    stop(
      "Named matrix `wgt` must uniquely identify the same groups and items as `lambda`.",
      call. = FALSE
    )
  }

  list(
    value = wgt[rownames(lambda), colnames(lambda), drop = FALSE],
    matching = "named"
  )
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

.pm_alignment_validate_named_list <- function(x, name, protected) {
  if (!is.list(x)) {
    stop(sprintf("`%s` must be a list.", name), call. = FALSE)
  }
  if (length(x) == 0L) return(x)

  fields <- names(x)
  if (is.null(fields) || anyNA(fields) || any(!nzchar(fields))) {
    stop(sprintf("Every element of `%s` must be named.", name), call. = FALSE)
  }
  if (anyDuplicated(fields)) {
    stop(sprintf("`%s` contains duplicate names.", name), call. = FALSE)
  }
  protected_found <- intersect(fields, protected)
  if (length(protected_found)) {
    stop(
      sprintf(
        "`%s` cannot override protected argument(s): %s.",
        name,
        paste(protected_found, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  x
}

.pm_alignment_validate_scalar_settings <- function(center,
                                                   optimizer,
                                                   fixed,
                                                   meth,
                                                   eps,
                                                   vcov,
                                                   digits,
                                                   noninvariance_cutoff) {
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("`center` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(optimizer) || length(optimizer) != 1L ||
      is.na(optimizer) || !optimizer %in% c("optim", "nlminb")) {
    stop("`optimizer` must be 'optim' or 'nlminb'.", call. = FALSE)
  }
  if (!is.null(fixed) &&
      (!is.logical(fixed) || length(fixed) != 1L || is.na(fixed))) {
    stop("`fixed` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  if (!is.numeric(meth) || length(meth) != 1L || is.na(meth) ||
      !is.finite(meth) || !meth %in% c(0, 0.5, 1, 2, 3, 4)) {
    stop("`meth` must be one of 0, 0.5, 1, 2, 3, or 4.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) ||
      !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite number.", call. = FALSE)
  }
  if (!is.null(vcov)) {
    if (is.null(dim(vcov)) || length(dim(vcov)) != 2L ||
        !is.numeric(vcov) || any(dim(vcov) == 0L) ||
        anyNA(vcov) || any(!is.finite(vcov))) {
      stop("`vcov` must be NULL or a nonempty finite numeric matrix.", call. = FALSE)
    }
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits != trunc(digits)) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(noninvariance_cutoff) ||
      length(noninvariance_cutoff) != 1L ||
      is.na(noninvariance_cutoff) || !is.finite(noninvariance_cutoff) ||
      noninvariance_cutoff < 0 || noninvariance_cutoff > 100) {
    stop(
      "`noninvariance_cutoff` must be a single finite numeric value between 0 and 100.",
      call. = FALSE
    )
  }

  list(
    center = center,
    optimizer = optimizer,
    fixed = fixed,
    meth = meth,
    eps = eps,
    vcov = vcov,
    digits = as.integer(digits),
    noninvariance_cutoff = noninvariance_cutoff
  )
}

.pm_alignment_missing_settings <- function(config_args) {
  requested <- if ("missing" %in% names(config_args)) {
    config_args[["missing"]]
  } else {
    NULL
  }
  normalized <- if (is.character(requested) && length(requested) == 1L &&
                    !is.na(requested) && nzchar(trimws(requested))) {
    tolower(trimws(requested))
  } else {
    NULL
  }
  effective <- if (is.null(normalized) || identical(normalized, "default")) {
    "listwise"
  } else if (normalized %in% c("ml", "fiml", "direct")) {
    "fiml"
  } else {
    normalized
  }
  list(engine_request = requested, effective = effective)
}

.pm_alignment_group_counts <- function(group, group_names, rows = NULL) {
  if (is.null(rows)) rows <- rep(TRUE, length(group))
  counts <- table(factor(
    as.character(group)[rows],
    levels = group_names
  ))
  stats::setNames(as.integer(counts), group_names)
}

.pm_alignment_expected_vcov_names <- function(lambda) {
  unlist(lapply(seq_len(nrow(lambda)), function(group_index) {
    c(
      paste0(colnames(lambda), "_lam_Gr", group_index),
      paste0(colnames(lambda), "_nu_Gr", group_index)
    )
  }), use.names = FALSE)
}

.pm_alignment_validate_vcov <- function(vcov,
                                        lambda,
                                        tolerance = sqrt(.Machine$double.eps)) {
  if (is.null(vcov)) {
    return(list(
      value = NULL,
      matching = "none",
      minimum_eigenvalue = NULL
    ))
  }

  expected_size <- as.integer(2L * nrow(lambda) * ncol(lambda))
  if (is.null(dim(vcov)) || length(dim(vcov)) != 2L ||
      !is.numeric(vcov) ||
      !identical(dim(vcov), c(expected_size, expected_size)) ||
      anyNA(vcov) || any(!is.finite(vcov))) {
    stop(
      sprintf(
        "`vcov` must be a finite numeric %d x %d matrix for these parameters.",
        expected_size,
        expected_size
      ),
      call. = FALSE
    )
  }

  matrix_scale <- max(1, max(abs(vcov)))
  if (max(abs(vcov - t(vcov))) > tolerance * matrix_scale) {
    stop("`vcov` must be symmetric within numerical tolerance.", call. = FALSE)
  }

  row_ids <- rownames(vcov)
  column_ids <- colnames(vcov)
  if (xor(is.null(row_ids), is.null(column_ids))) {
    stop("`vcov` must have both row and column names, or neither.", call. = FALSE)
  }
  matching <- "positional"
  if (!is.null(row_ids)) {
    expected_ids <- .pm_alignment_expected_vcov_names(lambda)
    invalid_ids <- anyNA(row_ids) || anyNA(column_ids) ||
      any(!nzchar(row_ids)) || any(!nzchar(column_ids)) ||
      anyDuplicated(row_ids) || anyDuplicated(column_ids) ||
      !identical(row_ids, column_ids) ||
      !identical(row_ids, expected_ids)
    if (invalid_ids) {
      stop(
        paste0(
          "Named `vcov` rows and columns must be identical, unique, and ",
          "follow the sirt parameter order for `lambda` and `nu`."
        ),
        call. = FALSE
      )
    }
    matching <- "named"
  }

  eigenvalues <- tryCatch(
    eigen((vcov + t(vcov)) / 2, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eigenvalues) || !length(eigenvalues) ||
      anyNA(eigenvalues) || any(!is.finite(eigenvalues))) {
    stop("`vcov` eigenvalues could not be evaluated.", call. = FALSE)
  }
  eigen_scale <- max(1, max(abs(eigenvalues)))
  minimum_eigenvalue <- min(eigenvalues)
  if (minimum_eigenvalue < -tolerance * eigen_scale) {
    stop("`vcov` must be positive semidefinite.", call. = FALSE)
  }

  list(
    value = vcov,
    matching = matching,
    minimum_eigenvalue = minimum_eigenvalue
  )
}

.pm_alignment_effective_fixed <- function(alignment) {
  if (!is.list(alignment)) return(NULL)
  value <- alignment[["fixed"]]
  if (!is.logical(value) || length(value) != 1L || is.na(value)) return(NULL)
  value
}

.pm_alignment_config_engine <- function(args) {
  do.call(sirt::invariance_alignment_cfa_config, args)
}

.pm_alignment_engine <- function(args) {
  do.call(sirt::invariance.alignment, args)
}

.pm_alignment_constraints_engine <- function(args) {
  do.call(sirt::invariance_alignment_constraints, args)
}

.pm_alignment_bound_stdout <- function(output,
                                       max_lines = 50L,
                                       max_characters = 4000L) {
  output <- as.character(output)
  original_lines <- length(output)
  output <- utils::head(output, max_lines)

  kept <- character()
  characters_kept <- 0L
  character_truncated <- FALSE
  for (line in output) {
    available <- max_characters - characters_kept
    if (available <= 0L) {
      character_truncated <- TRUE
      break
    }
    line_characters <- nchar(line, type = "chars")
    if (line_characters > available) {
      kept <- c(kept, substr(line, 1L, available))
      characters_kept <- max_characters
      character_truncated <- TRUE
      break
    }
    kept <- c(kept, line)
    characters_kept <- characters_kept + line_characters
  }

  list(
    output = kept,
    truncated = original_lines > max_lines || character_truncated
  )
}

.pm_alignment_capture_engine <- function(expr) {
  captured <- NULL
  stdout <- utils::capture.output(
    captured <- .pm_capture_conditions(expr),
    type = "output"
  )
  bounded <- .pm_alignment_bound_stdout(stdout)
  captured$stdout <- bounded$output
  captured$stdout_truncated <- bounded$truncated
  captured
}

.pm_alignment_convergence_status <- function(alignment) {
  if (!is.list(alignment)) {
    return(list(
      status = "not_reported",
      reason = "alignment_output_missing"
    ))
  }
  signals <- list(
    alignment$converged,
    alignment$convergence,
    alignment$status
  )
  if (is.list(alignment$res_optim)) {
    signals <- c(
      signals,
      list(
        alignment$res_optim$converged,
        alignment$res_optim$convergence,
        alignment$res_optim$status
      )
    )
  }
  signals <- Filter(Negate(is.null), signals)
  if (!length(signals)) {
    return(list(status = "not_reported", reason = NULL))
  }

  recognized <- FALSE
  for (signal in signals) {
    if (length(signal) != 1L || is.na(signal)) next
    if (is.logical(signal)) {
      recognized <- TRUE
      if (!isTRUE(signal)) {
        return(list(status = "not_converged", reason = "alignment_nonconvergence"))
      }
    } else if (is.numeric(signal)) {
      recognized <- TRUE
      if (!is.finite(signal) || signal != 0) {
        return(list(status = "not_converged", reason = "alignment_nonconvergence"))
      }
    } else if (is.character(signal)) {
      value <- tolower(trimws(signal))
      if (grepl("fail|error|not[ _-]?converg|nonconverg", value)) {
        return(list(status = "not_converged", reason = "alignment_nonconvergence"))
      }
      if (grepl("converg|success|^ok$", value)) recognized <- TRUE
    }
  }
  list(
    status = if (recognized) "converged" else "not_reported",
    reason = NULL
  )
}

.pm_alignment_objective_problem <- function(alignment) {
  if (!is.list(alignment)) return("alignment_objective_not_evaluable")
  candidates <- list(alignment[["fopt"]])
  res_optim <- alignment[["res_optim"]]
  if (is.list(res_optim)) {
    candidates <- c(candidates, list(res_optim[["value"]]))
  } else {
    candidates <- c(candidates, list(NULL))
  }
  available <- !vapply(candidates, is.null, logical(1L))
  if (!any(available)) return("alignment_objective_not_evaluable")
  valid <- vapply(candidates[available], function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value)
  }, logical(1L))
  if (!all(valid)) return("alignment_objective_not_evaluable")
  values <- vapply(candidates[available], as.numeric, numeric(1L))
  if (any(values < 0)) return("alignment_objective_out_of_range")
  NULL
}

.pm_alignment_numeric_problem <- function(x,
                                          name,
                                          expected_dim = NULL,
                                          expected_rows = NULL,
                                          expected_columns = NULL) {
  if (is.null(x)) return(paste0(name, "_missing"))
  converted <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(converted) || !is.numeric(converted) || !length(converted) ||
      anyNA(converted) || any(!is.finite(converted))) {
    return(paste0(name, "_nonfinite_or_empty"))
  }
  if (!is.null(expected_dim) && !identical(dim(converted), expected_dim)) {
    return(paste0(name, "_dimension_mismatch"))
  }
  if (!is.null(expected_rows) && !identical(rownames(converted), expected_rows)) {
    return(paste0(name, "_row_names_mismatch"))
  }
  if (!is.null(expected_columns) &&
      !identical(colnames(converted), expected_columns)) {
    return(paste0(name, "_column_names_mismatch"))
  }
  NULL
}

.pm_alignment_constraint_problem <- function(family,
                                             family_name,
                                             group_item_dim,
                                             group_names,
                                             item_names) {
  if (!is.list(family) || !is.numeric(family$prop_noninvariance) ||
      length(family$prop_noninvariance) != 1L ||
      is.na(family$prop_noninvariance) ||
      !is.finite(family$prop_noninvariance)) {
    return(paste0(family_name, "_summary_not_evaluable"))
  }
  proportion <- as.numeric(family$prop_noninvariance)
  if (proportion < 0 || proportion > 100) {
    return(paste0(family_name, "_prop_noninvariance_out_of_range"))
  }

  for (field in c("parm_dif", "parm_est")) {
    problem <- .pm_alignment_numeric_problem(
      family[[field]],
      paste0(family_name, "_", field),
      group_item_dim,
      group_names,
      item_names
    )
    if (!is.null(problem)) return(problem)
  }
  for (field in c("parm_joint", "N_parm_items")) {
    values <- family[[field]]
    if (!is.numeric(values) || length(values) != length(item_names) ||
        anyNA(values) || any(!is.finite(values)) ||
        !identical(names(values), item_names)) {
      return(paste0(family_name, "_", field, "_not_evaluable"))
    }
  }

  unique_by_group <- family$N_unique_parm_groups
  total <- family$N_total
  if (!is.numeric(unique_by_group) ||
      length(unique_by_group) != length(group_names) ||
      anyNA(unique_by_group) || any(!is.finite(unique_by_group)) ||
      any(unique_by_group < 0) ||
      !identical(names(unique_by_group), group_names) ||
      !is.numeric(total) || length(total) != 1L || is.na(total) ||
      !is.finite(total) || total <= 0) {
    return(paste0(family_name, "_noninvariance_counts_not_evaluable"))
  }
  miss_items <- family$miss_items
  miss_problem <- .pm_alignment_numeric_problem(
    miss_items,
    paste0(family_name, "_miss_items"),
    group_item_dim,
    group_names,
    item_names
  )
  if (!is.null(miss_problem)) {
    return(paste0(family_name, "_noninvariance_counts_not_evaluable"))
  }
  if (any(!(miss_items %in% c(0, 1)))) {
    return(paste0(family_name, "_miss_items_not_binary"))
  }
  flags <- abs(family$parm_dif) > 1e-5
  expected_unique_by_group <- rowSums(flags, na.rm = TRUE)
  expected_n_parm_items <- as.numeric(
    colSums(abs(family$parm_dif) <= 1e-5, na.rm = TRUE) > 0
  ) + colSums(flags, na.rm = TRUE)
  if (!isTRUE(all.equal(
        as.numeric(unique_by_group),
        as.numeric(expected_unique_by_group),
        tolerance = sqrt(.Machine$double.eps),
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        as.numeric(family$N_parm_items),
        as.numeric(expected_n_parm_items),
        tolerance = sqrt(.Machine$double.eps),
        check.attributes = FALSE
      )) ||
      abs(total - sum(miss_items)) > sqrt(.Machine$double.eps)) {
    return(paste0(family_name, "_noninvariance_counts_inconsistent"))
  }
  expected_proportion <- 100 * sum(unique_by_group) / total
  tolerance <- sqrt(.Machine$double.eps) * max(
    1,
    abs(expected_proportion),
    abs(proportion)
  )
  if (abs(proportion - expected_proportion) > tolerance) {
    return(paste0(family_name, "_prop_noninvariance_inconsistent"))
  }

  if (!is.null(family$N_parm_all)) {
    n_all <- family$N_parm_all
    if (!is.numeric(n_all) || length(n_all) != 1L || is.na(n_all) ||
        !is.finite(n_all) ||
        abs(n_all - sum(family$N_parm_items)) > sqrt(.Machine$double.eps)) {
      return(paste0(family_name, "_parameter_counts_inconsistent"))
    }
  }
  NULL
}

.pm_alignment_status <- function(alignment, constraints, lambda, nu) {
  if (!is.list(alignment)) {
    return(list(
      status = "not_evaluable",
      reason = "alignment_output_missing",
      convergence_status = "not_reported"
    ))
  }

  convergence <- .pm_alignment_convergence_status(alignment)
  if (!is.null(convergence$reason)) {
    return(list(
      status = "not_evaluable",
      reason = convergence$reason,
      convergence_status = convergence$status
    ))
  }

  objective_problem <- .pm_alignment_objective_problem(alignment)
  if (!is.null(objective_problem)) {
    return(list(
      status = "not_evaluable",
      reason = objective_problem,
      convergence_status = convergence$status
    ))
  }

  group_names <- rownames(lambda)
  item_names <- colnames(lambda)
  group_item_dim <- dim(lambda)
  problems <- c(
    .pm_alignment_numeric_problem(
      alignment$pars,
      "group_parameters",
      expected_rows = group_names
    ),
    .pm_alignment_numeric_problem(
      alignment$itempars.aligned,
      "aligned_item_parameters",
      expected_rows = item_names
    ),
    .pm_alignment_numeric_problem(
      alignment$lambda,
      "alignment_lambda",
      group_item_dim,
      group_names,
      item_names
    ),
    .pm_alignment_numeric_problem(
      alignment$nu,
      "alignment_nu",
      group_item_dim,
      group_names,
      item_names
    ),
    .pm_alignment_numeric_problem(
      alignment$lambda.aligned,
      "aligned_lambda",
      group_item_dim,
      group_names,
      item_names
    ),
    .pm_alignment_numeric_problem(
      alignment$nu.aligned,
      "aligned_nu",
      group_item_dim,
      group_names,
      item_names
    ),
    .pm_alignment_numeric_problem(
      alignment$wgt,
      "alignment_weights",
      group_item_dim,
      group_names,
      item_names
    )
  )
  problems <- problems[!vapply(problems, is.null, logical(1L))]
  if (length(problems)) {
    return(list(
      status = "not_evaluable",
      reason = as.character(problems[[1L]]),
      convergence_status = convergence$status
    ))
  }

  es <- alignment$es.invariance
  if (is.null(es) || is.null(dim(es)) ||
      !"R2" %in% rownames(es) ||
      !all(c("loadings", "intercepts") %in% colnames(es)) ||
      any(!is.finite(as.numeric(es["R2", c("loadings", "intercepts")]))) ) {
    return(list(
      status = "not_evaluable",
      reason = "alignment_r2_not_evaluable",
      convergence_status = convergence$status
    ))
  }

  if (!is.list(constraints)) {
    return(list(
      status = "not_evaluable",
      reason = "constraint_output_missing",
      convergence_status = convergence$status
    ))
  }
  for (family_name in c("lambda_list", "nu_list")) {
    family <- constraints[[family_name]]
    problem <- .pm_alignment_constraint_problem(
      family,
      family_name,
      group_item_dim,
      group_names,
      item_names
    )
    if (!is.null(problem)) {
      return(list(
        status = "not_evaluable",
        reason = problem,
        convergence_status = convergence$status
      ))
    }
  }

  list(
    status = "interpretable",
    reason = "alignment_outputs_valid",
    convergence_status = convergence$status
  )
}

.pm_alignment_empty_summary <- function(cutoff, reason) {
  data.frame(
    parameter = c("loadings", "intercepts"),
    R2 = c(NA_real_, NA_real_),
    percent_noninvariant = c(NA_real_, NA_real_),
    cutoff = rep(cutoff, 2L),
    invariant = c(NA, NA),
    status = rep("not_evaluable", 2L),
    reason = rep(reason, 2L),
    stringsAsFactors = FALSE
  )
}

.pm_alignment_not_evaluable_text <- function(reason, language) {
  if (identical(language, "esp")) {
    return(sprintf(
      paste0(
        "No fue posible evaluar la invarianza por alineamiento (%s). ",
        "No interprete porcentajes de invarianza ni valores R2."
      ),
      reason
    ))
  }
  sprintf(
    paste0(
      "Alignment invariance could not be evaluated (%s). ",
      "Do not interpret invariance percentages or R2 values."
    ),
    reason
  )
}

.pm_alignment_collect_conditions <- function(captures, field) {
  unlist(lapply(captures, `[[`, field), recursive = FALSE, use.names = FALSE)
}

.pm_alignment_build_result <- function(call,
                                       input_mode,
                                       requested_settings,
                                       effective_settings,
                                       data_summary,
                                       dimensions,
                                       configural,
                                       alignment,
                                       constraints,
                                       engine_calls,
                                       captures,
                                       status_info,
                                       language,
                                       settings,
                                       input_warnings = list()) {
  interpretable <- identical(status_info$status, "interpretable")
  item_results <- NULL
  if (interpretable) {
    summary_table <- build_alignment_summary(
      alignment = alignment,
      constraints = constraints,
      noninvariance_cutoff = settings$noninvariance_cutoff
    )
    item_results <- list(
      loadings = build_alignment_item_results(constraints$lambda_list),
      intercepts = build_alignment_item_results(constraints$nu_list)
    )
    report <- build_alignment_report(
      summary_table = summary_table,
      item_results = item_results,
      language = language,
      settings = settings,
      dimensions = dimensions
    )
    text <- language_pack_alignment(language)
    conclusion <- if (all(summary_table$invariant)) {
      text$overall_invariant
    } else {
      text$overall_noninvariant
    }
  } else {
    summary_table <- .pm_alignment_empty_summary(
      settings$noninvariance_cutoff,
      status_info$reason
    )
    conclusion <- .pm_alignment_not_evaluable_text(status_info$reason, language)
    report <- paste(
      language_pack_alignment(language)$title,
      conclusion,
      sep = "\n"
    )
  }

  captured_warnings <- c(
    input_warnings,
    .pm_alignment_collect_conditions(captures, "warnings")
  )
  captured_errors <- Filter(
    Negate(is.null),
    lapply(captures, `[[`, "error")
  )
  captured_messages <- .pm_alignment_collect_conditions(captures, "messages")
  engine_names <- c("configural", "alignment", "constraints")
  engine_stdout <- stats::setNames(lapply(engine_names, function(engine_name) {
    capture <- captures[[engine_name]]
    if (is.null(capture$stdout)) character() else capture$stdout
  }), engine_names)
  engine_stdout_truncated <- stats::setNames(
    lapply(engine_names, function(engine_name) {
      isTRUE(captures[[engine_name]]$stdout_truncated)
    }),
    engine_names
  )

  decision_log <- rbind(
    .pm_decision(
      "alignment_identity",
      requested = "lambda_nu",
      effective = effective_settings$identity_matching,
      reason = if (isTRUE(dimensions$names_generated)) {
        "stable_names_generated"
      } else {
        "matched_by_identifiers"
      }
    ),
    .pm_decision(
      "alignment_weights",
      requested = if (is.null(requested_settings$wgt)) "none" else "supplied",
      effective = effective_settings$weight_matching,
      reason = "validated_before_engine"
    ),
    .pm_decision(
      "alignment_evaluability",
      requested = "alignment",
      effective = status_info$status,
      reason = status_info$reason,
      severity = if (interpretable) "info" else "warning"
    )
  )

  noninvariance <- list(
    loadings = if (is.list(constraints)) constraints$lambda_list else NULL,
    intercepts = if (is.list(constraints)) constraints$nu_list else NULL
  )
  group_parameters <- if (is.list(alignment)) alignment$pars else NULL
  aligned_item_parameters <- if (is.list(alignment)) {
    alignment$itempars.aligned
  } else {
    NULL
  }
  results <- list(
    summary = summary_table,
    item_results = item_results,
    group_parameters = group_parameters,
    aligned_item_parameters = aligned_item_parameters,
    alignment = alignment,
    constraints = constraints,
    noninvariance = noninvariance
  )

  result <- list(
    call = call,
    input_mode = input_mode,
    settings = settings,
    dimensions = dimensions,
    configural = configural,
    alignment = alignment,
    constraints = constraints,
    summary = summary_table,
    item_results = item_results,
    group_parameters = group_parameters,
    aligned_item_parameters = aligned_item_parameters,
    noninvariance = noninvariance,
    references = alignment_apa_references(),
    report = report,
    requested_settings = requested_settings,
    effective_settings = effective_settings,
    data_summary = data_summary,
    results = results,
    warnings = .pm_conditions_to_table(captured_warnings, "warnings"),
    errors = .pm_conditions_to_table(captured_errors, "errors"),
    diagnostics = list(
      status = status_info$status,
      reason = status_info$reason,
      convergence_status = status_info$convergence_status,
      messages = .pm_conditions_to_table(captured_messages, "messages"),
      engine_stdout = engine_stdout,
      engine_stdout_truncated = engine_stdout_truncated
    ),
    decision_log = decision_log,
    conclusion = conclusion,
    conclusion_status = status_info$status,
    engine_calls = engine_calls
  )
  class(result) <- c("inv_align_auto_result", "psychomatic_result")
  result
}

#' Automated alignment invariance
#'
#' Runs alignment optimization from loading and intercept matrices (`lambda` and
#' `nu`) or from item-level data with a grouping vector.
#'
#' @param lambda,nu Parameter matrices for loadings and intercepts.
#' @inheritParams efa_auto
#' @param data Item-level data used when `lambda` and `nu` are not supplied.
#' @param group Grouping vector associated with `data`.
#' @param sampling_weights Optional sampling weights, one per row.
#' @param wgt Optional weights for the alignment procedure.
#' @param language Report language: `"esp"` or `"eng"`.
#' @param model Model type for the configural alignment step.
#' @param align.scale,align.pow,parm_tol Alignment tuning parameters.
#' @param noninvariance_cutoff Percentage used to flag non-invariance.
#' @param center,optimizer,fixed,meth,eps,vcov Advanced algorithm parameters.
#' @param digits Number of decimal places used in the report.
#' @param config_args Additional arguments for the configural step.
#' @param alignment_args Additional arguments for the alignment step.
#'
#' @return Object of class `inv_align_auto_result`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' inv_align_auto(data = MY_DATA[, 28:39], group = MY_DATA$Country,
#'                language = "eng")
#' }
#'
#' @export
inv_align_auto <- function(
  lambda = NULL,
  nu = NULL,
  data = NULL,
  group = NULL,
  sampling_weights = NULL,
  wgt = NULL,
  language = c("eng", "esp"),
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
  alignment_args = list(),
  report = FALSE, file_name = NULL, format = "excel", overwrite = FALSE, verbose = TRUE
) {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  finish_alignment <- function(...) {
    .pm_finish_analysis(.pm_alignment_build_result(...), output, language)
  }
  requested_call <- match.call()
  requested_settings <- list(
    input_mode = if (!is.null(lambda) || !is.null(nu)) "matrices" else "data",
    lambda = lambda,
    nu = nu,
    data = data,
    group = group,
    sampling_weights = sampling_weights,
    wgt = wgt,
    language = language,
    model = model,
    align.scale = align.scale,
    align.pow = align.pow,
    parm_tol = parm_tol,
    noninvariance_cutoff = noninvariance_cutoff,
    center = center,
    optimizer = optimizer,
    fixed = fixed,
    meth = meth,
    eps = eps,
    vcov = vcov,
    digits = digits,
    config_args = config_args,
    alignment_args = alignment_args
  )
  language <- normalize_alignment_language(language)
  model <- match.arg(model)

  config_args <- .pm_alignment_validate_named_list(
    config_args, "config_args", c("dat", "group", "weights", "model")
  )
  alignment_args <- .pm_alignment_validate_named_list(
    alignment_args,
    "alignment_args",
    c(
      "lambda", "nu", "wgt", "align.scale", "align.pow", "eps",
      "center", "optimizer", "fixed", "meth", "vcov"
    )
  )
  missing_settings <- if (identical(requested_settings$input_mode, "data")) {
    .pm_alignment_missing_settings(config_args)
  } else {
    list(engine_request = NULL, effective = "not_applicable")
  }
  requested_settings$missing_policy <- missing_settings$engine_request
  align.scale <- normalize_alignment_pair(
    align.scale, "align.scale", allow_zero = FALSE
  )
  align.pow <- normalize_alignment_pair(
    align.pow, "align.pow", allow_zero = TRUE
  )

  if (is.null(parm_tol)) {
    parm_tol <- align.scale
  }
  parm_tol <- normalize_alignment_pair(parm_tol, "parm_tol", allow_zero = TRUE)
  scalar_settings <- .pm_alignment_validate_scalar_settings(
    center, optimizer, fixed, meth, eps, vcov, digits, noninvariance_cutoff
  )
  center <- scalar_settings$center
  optimizer <- scalar_settings$optimizer
  fixed <- scalar_settings$fixed
  meth <- scalar_settings$meth
  eps <- scalar_settings$eps
  vcov <- scalar_settings$vcov
  digits <- scalar_settings$digits
  noninvariance_cutoff <- scalar_settings$noninvariance_cutoff

  require_alignment_package("sirt")

  input_mode <- requested_settings$input_mode
  configural <- NULL
  effective_settings_base <- list(
    input_mode = input_mode,
    model = if (identical(input_mode, "data")) model else NULL,
    language = language,
    align.scale = align.scale,
    align.pow = align.pow,
    parm_tol = parm_tol,
    noninvariance_cutoff = noninvariance_cutoff,
    center = center,
    optimizer = optimizer,
    fixed = NULL,
    fixed_engine_request = fixed,
    fixed_effective = NULL,
    meth = meth,
    eps = eps,
    vcov = NULL,
    vcov_engine_request = vcov,
    vcov_matching = if (is.null(vcov)) "none" else "pending",
    digits = digits,
    sampling_weights = if (identical(input_mode, "data")) sampling_weights else NULL,
    wgt = NULL,
    weights = list(
      sampling_weights = if (identical(input_mode, "data")) sampling_weights else NULL,
      alignment = NULL
    ),
    config_args = if (identical(input_mode, "data")) config_args else NULL,
    alignment_args = alignment_args,
    configural_status = if (identical(input_mode, "data")) "pending" else "not_applicable",
    sampling_weights_status = if (identical(input_mode, "data")) {
      if (is.null(sampling_weights)) "not_supplied" else "pending"
    } else {
      "not_applicable"
    },
    config_args_status = if (identical(input_mode, "data")) {
      if (length(config_args)) "pending" else "not_supplied"
    } else {
      "not_applicable"
    },
    missing_engine_request = if (identical(input_mode, "data")) {
      missing_settings$engine_request
    } else {
      NULL
    },
    missing_policy = if (identical(input_mode, "data")) NULL else "not_applicable"
  )
  engine_calls <- list(configural = NULL, alignment = NULL, constraints = NULL)
  captures <- list()
  input_warnings <- list()

  if (identical(input_mode, "matrices")) {
    if (is.null(lambda) || is.null(nu)) {
      stop("Provide both `lambda` and `nu`, or provide both `data` and `group`.", call. = FALSE)
    }

    if (!is.null(data) || !is.null(group) || !is.null(sampling_weights)) {
      input_warnings <- list(simpleWarning(paste0(
        "`data`, `group`, and `sampling_weights` were ignored because ",
        "`lambda` and `nu` were provided."
      )))
    }

    matrices <- validate_parameter_matrices(lambda, nu)
    lambda <- matrices$lambda
    nu <- matrices$nu
    group_names <- rownames(lambda)
    item_names <- colnames(lambda)
    names_generated <- matrices$names_generated
    identity_matching <- matrices$identity_matching
    data_summary <- list(
      input_mode = input_mode,
      rows = NA_integer_,
      groups = nrow(lambda),
      items = ncol(lambda),
      group_names = group_names,
      item_names = item_names,
      names_generated = names_generated
    )
  } else {
    if (is.null(data) || is.null(group)) {
      stop("Provide both `lambda` and `nu`, or provide both `data` and `group`.", call. = FALSE)
    }

    if (is.null(dim(data)) || length(dim(data)) != 2L) {
      stop("`data` must be a two-dimensional data frame or matrix.", call. = FALSE)
    }
    data <- as.data.frame(data, check.names = FALSE)
    if (nrow(data) == 0L || ncol(data) < 2L) {
      stop("`data` must contain rows and at least two item columns.", call. = FALSE)
    }
    if (!all(vapply(data, is.numeric, logical(1L)))) {
      stop("`data` must contain only numeric item columns.", call. = FALSE)
    }
    if (any(vapply(
      data,
      function(column) any(is.nan(column) | is.infinite(column)),
      logical(1L)
    ))) {
      stop("`data` cannot contain NaN or infinite item values.", call. = FALSE)
    }
    if (nrow(data) != length(group)) {
      stop("`data` and `group` must contain the same number of rows.", call. = FALSE)
    }
    if (anyNA(group)) {
      stop("`group` cannot contain missing values.", call. = FALSE)
    }
    group_names <- unique(as.character(group))
    if (length(group_names) < 2L || any(!nzchar(group_names))) {
      stop("`group` must contain at least two nonempty groups.", call. = FALSE)
    }
    item_names <- names(data)
    if (is.null(item_names) || anyNA(item_names) || any(!nzchar(item_names)) ||
        anyDuplicated(item_names)) {
      stop("`data` item names must be unique, nonmissing, and nonempty.", call. = FALSE)
    }
    if (!is.null(sampling_weights) &&
        (!is.numeric(sampling_weights) ||
         length(sampling_weights) != nrow(data) ||
         anyNA(sampling_weights) || any(!is.finite(sampling_weights)) ||
         any(sampling_weights < 0) || all(sampling_weights == 0))) {
      stop(
        paste0(
          "`sampling_weights` must contain one finite, non-negative value ",
          "per row and cannot be all zero."
        ),
        call. = FALSE
      )
    }

    group_values <- as.character(group)
    entirely_missing <- vapply(
      group_names,
      function(group_name) {
        vapply(
          data[group_values == group_name, , drop = FALSE],
          function(column) all(is.na(column)),
          logical(1L)
        )
      },
      logical(ncol(data))
    )
    missing_pairs <- which(entirely_missing, arr.ind = TRUE)
    if (nrow(missing_pairs)) {
      pair_labels <- paste0(
        "group '",
        colnames(entirely_missing)[missing_pairs[, "col"]],
        "' / item '",
        rownames(entirely_missing)[missing_pairs[, "row"]],
        "'"
      )
      stop(
        paste0(
          "Every group-item cell must contain at least one observed value; ",
          "entirely missing pair(s): ",
          paste(pair_labels, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    missing_by_item <- vapply(
      data,
      function(column) as.integer(sum(is.na(column))),
      integer(1L)
    )
    missing_values <- as.integer(sum(missing_by_item))
    complete_case <- stats::complete.cases(data)
    missing_by_row <- as.integer(rowSums(is.na(data)))
    missing_rows <- as.integer(sum(!complete_case))
    complete_rows <- as.integer(nrow(data) - missing_rows)
    n_input_by_group <- .pm_alignment_group_counts(group, group_names)
    n_complete_by_group <- .pm_alignment_group_counts(
      group,
      group_names,
      complete_case
    )
    data_summary <- list(
      input_mode = input_mode,
      rows = nrow(data),
      groups = length(group_names),
      items = ncol(data),
      group_names = group_names,
      item_names = item_names,
      names_generated = FALSE,
      missing_values = missing_values,
      missing_by_item = missing_by_item,
      missing_by_row = missing_by_row,
      missing_rows = missing_rows,
      complete_rows = complete_rows,
      n_input_by_group = n_input_by_group,
      n_complete_by_group = n_complete_by_group,
      missing_policy_requested = missing_settings$engine_request,
      missing_policy_effective = NULL,
      missing_policy = NULL,
      configural_n_by_group = NULL,
      n_effective_by_group = NULL
    )
    effective_settings_base$configural_status <- "requested"
    effective_settings_base$sampling_weights_status <- if (is.null(sampling_weights)) {
      "not_supplied"
    } else {
      "applied"
    }
    effective_settings_base$config_args_status <- if (length(config_args)) {
      "applied"
    } else {
      "not_supplied"
    }
    names_generated <- FALSE
    identity_matching <- "data_identifiers"
    if (!is.null(wgt)) {
      identity_template <- matrix(
        0,
        nrow = length(group_names),
        ncol = length(item_names),
        dimnames = list(group_names, item_names)
      )
      normalize_alignment_weights(wgt, identity_template)
    }
    provisional_dimensions <- list(
      groups = length(group_names),
      items = length(item_names),
      group_names = group_names,
      item_names = item_names,
      names_generated = FALSE
    )
    provisional_settings <- list(
      input_mode = input_mode,
      model = model,
      language = language,
      align.scale = align.scale,
      align.pow = align.pow,
      parm_tol = parm_tol,
      noninvariance_cutoff = noninvariance_cutoff,
      digits = digits
    )

    config_call <- utils::modifyList(
      list(
        dat = data,
        group = group,
        weights = sampling_weights,
        model = model,
        verbose = FALSE
      ),
      config_args,
      keep.null = TRUE
    )
    engine_calls$configural <- list(
      function_name = "sirt::invariance_alignment_cfa_config",
      arguments = config_call
    )
    config_capture <- .pm_alignment_capture_engine(
      .pm_alignment_config_engine(config_call)
    )
    captures <- c(captures, list(configural = config_capture))
    configural <- config_capture$value

    if (!is.null(config_capture$error)) {
      effective_settings <- utils::modifyList(effective_settings_base, list(
        lambda = NULL,
        nu = NULL,
        wgt = NULL,
        configural_status = "failed",
        identity_matching = identity_matching,
        weight_matching = if (is.null(sampling_weights)) "none" else "rowwise"
      ), keep.null = TRUE)
      return(finish_alignment(
        requested_call, input_mode, requested_settings, effective_settings,
        data_summary, provisional_dimensions, NULL, NULL, NULL,
        engine_calls, captures,
        list(
          status = "not_evaluable",
          reason = "configural_engine_error",
          convergence_status = "not_reported"
        ),
        language, provisional_settings
      ))
    }

    data_summary$missing_policy_effective <- missing_settings$effective
    data_summary$missing_policy <- missing_settings$effective
    effective_settings_base$missing_policy <- missing_settings$effective

    config_convergence <- .pm_alignment_convergence_status(configural)
    if (!is.null(config_convergence$reason)) {
      effective_settings <- utils::modifyList(effective_settings_base, list(
        lambda = NULL,
        nu = NULL,
        wgt = NULL,
        configural_status = "nonconverged",
        identity_matching = identity_matching,
        weight_matching = if (is.null(sampling_weights)) "none" else "rowwise"
      ), keep.null = TRUE)
      return(finish_alignment(
        requested_call, input_mode, requested_settings, effective_settings,
        data_summary, provisional_dimensions, configural, NULL, NULL,
        engine_calls, captures,
        list(
          status = "not_evaluable",
          reason = if (identical(
            config_convergence$reason,
            "alignment_output_missing"
          )) {
            "configural_output_invalid"
          } else {
            "configural_nonconvergence"
          },
          convergence_status = config_convergence$status
        ),
        language, provisional_settings
      ))
    }

    matrix_capture <- .pm_capture_conditions({
      validated_matrices <- validate_parameter_matrices(
        configural$lambda,
        configural$nu
      )
      configural_group_names <- rownames(validated_matrices$lambda)
      configural_n_names <- names(configural$N)
      invalid_configural_n_names <- !is.null(configural_n_names) &&
        (anyNA(configural_n_names) || any(!nzchar(configural_n_names)) ||
         anyDuplicated(configural_n_names) ||
         !setequal(configural_n_names, configural_group_names))
      if (!setequal(
            names(data_summary$n_input_by_group),
            configural_group_names
          ) || invalid_configural_n_names) {
        stop(
          "Configural output must identify the same groups as `group`.",
          call. = FALSE
        )
      }
      validated_matrices
    })
    captures <- c(captures, list(configural_matrices = matrix_capture))
    if (!is.null(matrix_capture$error)) {
      effective_settings <- utils::modifyList(effective_settings_base, list(
        lambda = NULL,
        nu = NULL,
        wgt = NULL,
        configural_status = "invalid_output",
        identity_matching = identity_matching,
        weight_matching = if (is.null(sampling_weights)) "none" else "rowwise"
      ), keep.null = TRUE)
      return(finish_alignment(
        requested_call, input_mode, requested_settings, effective_settings,
        data_summary, provisional_dimensions, configural, NULL, NULL,
        engine_calls, captures,
        list(
          status = "not_evaluable",
          reason = "configural_parameter_matrices_invalid",
          convergence_status = config_convergence$status
        ),
        language, provisional_settings
      ))
    }

    matrices <- matrix_capture$value
    lambda <- matrices$lambda
    nu <- matrices$nu
    group_names <- rownames(lambda)
    item_names <- colnames(lambda)
    names_generated <- matrices$names_generated
    identity_matching <- matrices$identity_matching
    data_summary$groups <- nrow(lambda)
    data_summary$items <- ncol(lambda)
    data_summary$group_names <- group_names
    data_summary$item_names <- item_names
    data_summary$names_generated <- names_generated
    data_summary$n_input_by_group <-
      data_summary$n_input_by_group[group_names]
    data_summary$n_complete_by_group <-
      data_summary$n_complete_by_group[group_names]
    configural_n <- configural$N
    if (is.numeric(configural_n) && length(configural_n) == length(group_names)) {
      if (!is.null(names(configural_n))) {
        configural_n <- configural_n[group_names]
      } else if (is.null(names(configural_n))) {
        names(configural_n) <- group_names
      }
    }
    data_summary$configural_n_by_group <- configural_n
    data_summary$n_effective_by_group <- configural_n
    effective_settings_base$configural_status <- "completed"

    if (is.null(vcov) && !is.null(configural$vcov)) {
      vcov <- configural$vcov
    }
  }

  weights <- normalize_alignment_weights(wgt, lambda)
  wgt <- weights$value
  dimensions <- list(
    groups = nrow(lambda),
    items = ncol(lambda),
    group_names = rownames(lambda),
    item_names = colnames(lambda),
    names_generated = names_generated
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

  if (!is.null(vcov)) {
    vcov_capture <- .pm_capture_conditions(
      .pm_alignment_validate_vcov(vcov, lambda)
    )
    captures <- c(captures, list(vcov = vcov_capture))
    if (!is.null(vcov_capture$error)) {
      effective_settings <- utils::modifyList(effective_settings_base, list(
        lambda = lambda,
        nu = nu,
        wgt = wgt,
        configural_status = if (identical(input_mode, "data")) {
          "invalid_output"
        } else {
          "not_applicable"
        },
        vcov = NULL,
        vcov_matching = "invalid",
        identity_matching = identity_matching,
        weight_matching = weights$matching
      ), keep.null = TRUE)
      return(finish_alignment(
        requested_call, input_mode, requested_settings, effective_settings,
        data_summary, dimensions, configural, NULL, NULL, engine_calls,
        captures,
        list(
          status = "not_evaluable",
          reason = "alignment_vcov_invalid",
          convergence_status = "not_reported"
        ),
        language, settings, input_warnings
      ))
    }
    vcov_validation <- vcov_capture$value
    vcov <- vcov_validation$value
    vcov_matching <- vcov_validation$matching
    vcov_minimum_eigenvalue <- vcov_validation$minimum_eigenvalue
  } else {
    vcov_matching <- "none"
    vcov_minimum_eigenvalue <- NULL
  }

  effective_settings <- list(
    input_mode = input_mode,
    model = if (identical(input_mode, "data")) model else NULL,
    language = language,
    align.scale = align.scale,
    align.pow = align.pow,
    parm_tol = parm_tol,
    noninvariance_cutoff = noninvariance_cutoff,
    center = center,
    optimizer = optimizer,
    fixed = NULL,
    fixed_engine_request = fixed,
    fixed_effective = NULL,
    meth = meth,
    eps = eps,
    vcov = vcov,
    vcov_engine_request = vcov,
    vcov_matching = vcov_matching,
    vcov_minimum_eigenvalue = vcov_minimum_eigenvalue,
    digits = digits,
    sampling_weights = if (identical(input_mode, "data")) sampling_weights else NULL,
    wgt = wgt,
    weights = list(
      sampling_weights = if (identical(input_mode, "data")) sampling_weights else NULL,
      alignment = wgt
    ),
    config_args = if (identical(input_mode, "data")) config_args else NULL,
    alignment_args = alignment_args,
    lambda = lambda,
    nu = nu,
    identity_matching = identity_matching,
    weight_matching = weights$matching,
    configural_status = if (identical(input_mode, "data")) "completed" else "not_applicable",
    sampling_weights_status = if (identical(input_mode, "data")) {
      if (is.null(sampling_weights)) "not_supplied" else "applied"
    } else {
      "not_applicable"
    },
    config_args_status = if (identical(input_mode, "data")) {
      if (length(config_args)) "applied" else "not_supplied"
    } else {
      "not_applicable"
    },
    missing_engine_request = if (identical(input_mode, "data")) {
      missing_settings$engine_request
    } else {
      NULL
    },
    missing_policy = if (identical(input_mode, "data")) {
      data_summary$missing_policy_effective
    } else {
      "not_applicable"
    }
  )
  data_summary$weight_matching <- weights$matching

  alignment_call <- c(
    list(
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
    ),
    alignment_args
  )
  engine_calls$alignment <- list(
    function_name = "sirt::invariance.alignment",
    arguments = alignment_call
  )
  alignment_capture <- .pm_alignment_capture_engine(
    .pm_alignment_engine(alignment_call)
  )
  captures <- c(captures, list(alignment = alignment_capture))
  alignment <- alignment_capture$value
  if (!is.null(alignment_capture$error)) {
    return(finish_alignment(
      requested_call, input_mode, requested_settings, effective_settings,
      data_summary, dimensions, configural, NULL, NULL, engine_calls,
      captures,
      list(
        status = "not_evaluable",
        reason = "alignment_engine_error",
        convergence_status = "not_reported"
      ),
      language, settings, input_warnings
    ))
  }

  fixed_effective <- .pm_alignment_effective_fixed(alignment)
  effective_settings$fixed <- fixed_effective
  effective_settings$fixed_effective <- fixed_effective

  convergence <- .pm_alignment_convergence_status(alignment)
  pre_problems <- character()
  if (!is.null(convergence$reason)) {
    pre_problems <- convergence$reason
  } else {
    objective_problem <- .pm_alignment_objective_problem(alignment)
    if (!is.null(objective_problem)) {
      pre_problems <- objective_problem
    }
  }
  if (!length(pre_problems)) {
    for (spec in list(
      list(alignment$pars, "group_parameters", NULL, rownames(lambda), NULL),
      list(alignment$itempars.aligned, "aligned_item_parameters", NULL, colnames(lambda), NULL),
      list(alignment$lambda, "alignment_lambda", dim(lambda), rownames(lambda), colnames(lambda)),
      list(alignment$nu, "alignment_nu", dim(lambda), rownames(lambda), colnames(lambda)),
      list(alignment$lambda.aligned, "aligned_lambda", dim(lambda), rownames(lambda), colnames(lambda)),
      list(alignment$nu.aligned, "aligned_nu", dim(lambda), rownames(lambda), colnames(lambda)),
      list(alignment$wgt, "alignment_weights", dim(lambda), rownames(lambda), colnames(lambda))
    )) {
      problem <- .pm_alignment_numeric_problem(
        spec[[1L]], spec[[2L]], spec[[3L]], spec[[4L]], spec[[5L]]
      )
      if (!is.null(problem)) {
        pre_problems <- problem
        break
      }
    }
    es <- alignment$es.invariance
    if (!length(pre_problems) &&
        (is.null(es) || is.null(dim(es)) || !"R2" %in% rownames(es) ||
         !all(c("loadings", "intercepts") %in% colnames(es)) ||
         any(!is.finite(as.numeric(es["R2", c("loadings", "intercepts")]))))) {
      pre_problems <- "alignment_r2_not_evaluable"
    }
  }
  if (length(pre_problems)) {
    return(finish_alignment(
      requested_call, input_mode, requested_settings, effective_settings,
      data_summary, dimensions, configural, alignment, NULL, engine_calls,
      captures,
      list(
        status = "not_evaluable",
        reason = pre_problems[[1L]],
        convergence_status = convergence$status
      ),
      language, settings, input_warnings
    ))
  }

  constraint_call <- list(
    model = alignment,
    lambda_parm_tol = unname(parm_tol["lambda"]),
    nu_parm_tol = unname(parm_tol["nu"])
  )
  engine_calls$constraints <- list(
    function_name = "sirt::invariance_alignment_constraints",
    arguments = constraint_call
  )
  constraint_capture <- .pm_alignment_capture_engine(
    .pm_alignment_constraints_engine(constraint_call)
  )
  captures <- c(captures, list(constraints = constraint_capture))
  constraints <- constraint_capture$value
  if (!is.null(constraint_capture$error)) {
    return(finish_alignment(
      requested_call, input_mode, requested_settings, effective_settings,
      data_summary, dimensions, configural, alignment, NULL, engine_calls,
      captures,
      list(
        status = "not_evaluable",
        reason = "constraint_engine_error",
        convergence_status = convergence$status
      ),
      language, settings, input_warnings
    ))
  }

  status_info <- .pm_alignment_status(alignment, constraints, lambda, nu)
  finish_alignment(
    requested_call, input_mode, requested_settings, effective_settings,
    data_summary, dimensions, configural, alignment, constraints, engine_calls,
    captures, status_info, language, settings, input_warnings
  )
}

#' Print an `inv_align_auto_result` object
#'
#' @param x Object returned by `inv_align_auto()`.
#' @param ... Additional unused arguments.
#'
#' @rdname print.inv_align_auto_result
#' @method print inv_align_auto_result
#' @keywords internal
#' @export
print.inv_align_auto_result <- function(x, ...) {
  report <- x$report
  if (.pm_reference_ansi()) for (reference in x$references)
    report <- gsub(reference, .pm_reference_console(reference), report, fixed=TRUE)
  cat(report, sep = "\n")
  invisible(x)
}

#' Summarize an `inv_align_auto_result` object
#'
#' @param object Object returned by `inv_align_auto()`.
#' @param ... Additional unused arguments.
#'
#' @rdname summary.inv_align_auto_result
#' @method summary inv_align_auto_result
#' @keywords internal
#' @rawNamespace S3method(summary,inv_align_auto_result)
#' @export
summary.inv_align_auto_result <- function(object, ...) {
  object$summary
}
