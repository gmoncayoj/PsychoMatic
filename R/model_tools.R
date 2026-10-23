#' Fit basic IRT models automatically
#'
#' Fits a one-dimensional Rasch, 2PL, or graded response model using `mirt`.
#' The dependency is optional and is loaded only when this function is called.
#'
#' @inheritParams efa_auto
#' @param data Data frame or matrix containing item responses.
#' @param model IRT model: `"auto"`, `"rasch"`, `"2pl"`, or `"grm"`.
#' @param factors Number of latent dimensions passed to `mirt::mirt()`.
#' @param itemtype Optional `mirt` item type. If `NULL`, it is selected from
#'   `model`.
#' @param technical Optional technical arguments passed to `mirt::mirt()`.
#' @param language Interpretation language: `"eng"` or `"esp"`.
#' @param digits Number of decimal places used in interpretation tables.
#'
#' @return A list of class `psychomatic_irt` with fitted model, item
#'   parameters, item-level interpretation, a model-level summary, and
#'   supporting APA 7 references.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' irt_auto(MY_DATA[, 3:14])
#' }
#'
#' @export
irt_auto <- function(data,
                     model = c("auto", "rasch", "2pl", "grm"),
                     factors = 1,
                     itemtype = NULL,
                     technical = list(),
                     language = c("eng", "esp"),
                     digits = 3, report = FALSE, file_name = NULL,
                     format = "excel", overwrite = FALSE, verbose = TRUE) {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  call <- match.call()
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package `mirt` is required for `irt_auto()`.", call. = FALSE)
  }
  requested_model <- match.arg(model)
  language <- match.arg(language)
  factors <- .pm_validate_scalar_number(
    factors,
    name = "factors",
    lower = 1,
    integer = TRUE
  )
  digits <- .pm_validate_scalar_number(
    digits,
    name = "digits",
    lower = 0,
    integer = TRUE
  )
  if (!is.list(technical)) {
    stop("`technical` must be a list.", call. = FALSE)
  }

  prepared <- .pm_prepare_irt_data(data)
  resolved <- .pm_resolve_irt_model(
    requested_model,
    itemtype,
    prepared$n_categories,
    prepared$item_names
  )

  fit_capture <- .pm_capture_conditions(
    mirt::mirt(
      data = prepared$data,
      model = as.integer(factors),
      itemtype = unname(resolved$itemtype),
      verbose = FALSE,
      technical = technical
    )
  )
  fit <- fit_capture$value
  fit_error <- fit_capture$error

  empty_capture <- list(value = NULL, warnings = list(), messages = list(), error = NULL)
  convergence_capture <- empty_capture
  raw_capture <- empty_capture
  legacy_capture <- empty_capture
  raw_coefficients <- NULL
  legacy_coefficients <- NULL
  raw_item_parameters <- data.frame()
  item_parameters <- data.frame()
  convergence <- FALSE
  initial_fit_reason <- if (is.null(fit_error)) {
    "model fitting failed"
  } else {
    conditionMessage(fit_error)
  }
  assessment <- list(admissible = FALSE, reason = initial_fit_reason)
  fit_statistics <- list(
    status = "not_available",
    reason = initial_fit_reason,
    warnings = list(),
    error = NULL
  )

  if (!is.null(fit)) {
    convergence_capture <- .pm_capture_conditions(mirt::extract.mirt(fit, "converged"))
    convergence <- isTRUE(convergence_capture$value)
    raw_capture <- .pm_capture_conditions(
      mirt::coef(fit, IRTpars = FALSE, simplify = TRUE)
    )
    raw_coefficients <- raw_capture$value

    if (!is.null(raw_coefficients)) {
      raw_item_parameters <- .pm_mirt_item_parameters(raw_coefficients)
      raw_item_parameters <- .pm_add_discrimination_magnitude(raw_item_parameters)
    }

    if (as.integer(factors) == 1L) {
      legacy_capture <- .pm_capture_conditions(
        mirt::coef(fit, IRTpars = TRUE, simplify = TRUE)
      )
      legacy_coefficients <- legacy_capture$value
      if (!is.null(legacy_coefficients)) {
        legacy_items <- .pm_mirt_item_parameters(legacy_coefficients)
        item_parameters <- .pm_merge_item_parameter_aliases(
          raw_item_parameters,
          legacy_items
        )
      } else {
        item_parameters <- raw_item_parameters
      }
    } else {
      item_parameters <- raw_item_parameters
    }

    assessment <- .pm_assess_irt_solution(
      converged = convergence,
      item_parameters = raw_item_parameters,
      error = if (!is.null(convergence_capture$error)) {
        convergence_capture$error
      } else {
        raw_capture$error
      },
      itemtype = resolved$itemtype,
      n_categories = prepared$n_categories
    )
    fit_statistics <- .pm_mirt_fit_statistics(fit, convergence)
  }

  interpretation <- .pm_irt_interpretation(
    item_parameters,
    model = resolved$model,
    itemtype = resolved$display_itemtype,
    language = language,
    digits = digits
  )
  model_summary <- .pm_irt_model_summary(
    interpretation,
    model = resolved$model,
    itemtype = resolved$display_itemtype,
    n = nrow(prepared$data),
    language = language,
    digits = digits,
    raw_discrimination = if (nrow(item_parameters)) {
      vapply(
        seq_len(nrow(item_parameters)),
        function(i) .pm_item_discrimination(item_parameters[i, , drop = FALSE]),
        numeric(1)
      )
    } else {
      numeric()
    }
  )
  if (!isTRUE(assessment$admissible) && nrow(model_summary)) {
    model_summary$recommendation <- if (language == "esp") {
      paste("Modelo no evaluable:", assessment$reason)
    } else {
      paste("Model not evaluable:", assessment$reason)
    }
  }

  conclusion <- if (isTRUE(assessment$admissible)) {
    list(status = "evaluable", reason = NA_character_)
  } else {
    .pm_not_evaluable(assessment$reason)
  }
  result <- .pm_result(
    class = "psychomatic_irt",
    call = call,
    requested_settings = list(
      model = requested_model,
      factors = factors,
      itemtype = itemtype,
      technical = technical,
      language = language,
      digits = digits
    ),
    effective_settings = list(
      model = resolved$model,
      factors = as.integer(factors),
      itemtype = resolved$itemtype,
      language = language
    ),
    data_summary = list(
      n = nrow(prepared$data),
      n_items = ncol(prepared$data),
      categories = prepared$n_categories,
      category_map = prepared$category_map
    ),
    results = list(
      raw_coefficients = raw_coefficients,
      raw_item_parameters = raw_item_parameters,
      item_parameters = item_parameters,
      fit_statistics = fit_statistics
    ),
    warnings = c(
      fit_capture$warnings,
      convergence_capture$warnings,
      raw_capture$warnings,
      legacy_capture$warnings,
      fit_statistics$warnings
    ),
    errors = Filter(
      Negate(is.null),
      list(
        fit_capture$error,
        convergence_capture$error,
        raw_capture$error,
        legacy_capture$error,
        fit_statistics$error
      )
    ),
    diagnostics = list(
      converged = convergence,
      admissible = assessment$admissible,
      reason = assessment$reason,
      category_map = prepared$category_map
    ),
    decision_log = .pm_decision(
      code = "irt_model",
      requested = requested_model,
      effective = resolved$model,
      reason = resolved$reason
    ),
    conclusion = conclusion
  )
  legacy <- list(
    fit = fit,
    coefficients = if (!is.null(legacy_coefficients)) legacy_coefficients else raw_coefficients,
    raw_coefficients = raw_coefficients,
    coefficients_raw = raw_coefficients,
    coefficients_irt = legacy_coefficients,
    raw_item_parameters = raw_item_parameters,
    item_parameters = item_parameters,
    interpretation = interpretation,
    model_summary = model_summary,
    fit_statistics = fit_statistics,
    model = resolved$model,
    itemtype = resolved$display_itemtype,
    effective_itemtype = resolved$itemtype,
    factors = as.integer(factors),
    language = language,
    n = nrow(prepared$data),
    n_items = ncol(prepared$data),
    categories = prepared$n_categories,
    category_map = prepared$category_map,
    converged = convergence,
    admissible = assessment$admissible,
    reason = assessment$reason,
    references = .pm_irt_references(),
    call = call
  )
  for (field in names(legacy)) result[field] <- list(legacy[[field]])
  .pm_finish_analysis(result, output, language)
}

#' Screen differential item functioning by group
#'
#' Performs a lightweight DIF screening using logistic regression. By default,
#' polytomous items are reported as not evaluable because dichotomization can
#' obscure category-specific DIF. Median dichotomization is available only as
#' an explicit compatibility option and should be followed by a dedicated IRT
#' or ordinal DIF analysis.
#'
#' @inheritParams efa_auto
#' @param data Data frame or matrix containing item responses. If `group` is a
#'   column name in `data`, that column is used as the grouping variable.
#' @param group Grouping vector or the name of the grouping column in `data`.
#' @param items Optional item names. Defaults to all numeric item columns except
#'   the grouping column.
#' @param alpha Significance level used for DIF flags.
#' @param language Interpretation language: `"eng"` or `"esp"`.
#' @param digits Number of decimal places used in interpretation tables.
#' @param polytomous Handling of polytomous focal items. The default,
#'   `"not_evaluable"`, does not dichotomize them. Use `"median"` only to opt
#'   into a compatibility screen that is labeled and accompanied by a
#'   methodological warning.
#'
#' @return A list of class `psychomatic_dif` with screening results,
#'   item-level interpretation (`dif_type`, decision, severity, and
#'   recommendation), an overall DIF summary, and supporting APA 7 references.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' dif_auto(MY_DATA[, 3:14], group = MY_DATA$Country)
#' }
#'
#' @export
dif_auto <- function(data,
                     group,
                     items = NULL,
                     alpha = 0.01,
                     language = c("eng", "esp"),
                     digits = 3,
                     polytomous = c("not_evaluable", "median"),
                     report = FALSE, file_name = NULL, format = "excel",
                     overwrite = FALSE, verbose = TRUE) {
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  call <- match.call()
  language <- match.arg(language)
  polytomous <- match.arg(polytomous)
  alpha <- .pm_validate_scalar_number(alpha, "alpha", lower = 0, upper = 1)
  if (alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be strictly between 0 and 1.", call. = FALSE)
  }
  digits <- .pm_validate_scalar_number(digits, "digits", lower = 0, integer = TRUE)
  prepared <- .pm_prepare_dif_data(data, group, items)

  analyses <- lapply(
    prepared$focal_items,
    function(item) {
      .pm_analyse_dif_item(
        item = item,
        normalized_items = prepared$normalized_items,
        group = prepared$group,
        alpha = alpha,
        polytomous = polytomous,
        original_indices = seq_len(nrow(prepared$data)),
        min_group_n = 10L
      )
    }
  )
  rows <- lapply(analyses, `[[`, "row")

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  interpretation <- .pm_dif_interpretation(
    results,
    alpha = alpha,
    language = language,
    digits = digits
  )
  summary <- .pm_dif_summary(
    interpretation,
    alpha = alpha,
    group_levels = levels(prepared$group),
    language = language,
    digits = digits
  )
  captured_warnings <- unlist(lapply(analyses, `[[`, "warnings"), recursive = FALSE)
  if (any(vapply(
    analyses,
    function(x) identical(as.character(x$row$method[[1L]]), "median_dichotomized"),
    logical(1)
  ))) {
    captured_warnings[[length(captured_warnings) + 1L]] <- structure(
      list(
        message = paste(
          "Polytomous items were median-dichotomized for compatibility;",
          "confirm any DIF finding with a dedicated ordinal or IRT method."
        ),
        call = NULL,
        code = "dif_median_dichotomization"
      ),
      class = c("psychomatic_warning", "warning", "condition")
    )
  }
  captured_errors <- unlist(lapply(analyses, `[[`, "errors"), recursive = FALSE)
  conclusion <- if (summary$not_evaluable_items > 0L) {
    .pm_not_evaluable("One or more requested items were not evaluable; inspect item-level reasons.")
  } else {
    list(status = "evaluable", reason = NA_character_)
  }
  out <- .pm_result(
    class = "psychomatic_dif",
    call = call,
    requested_settings = list(
      group = group,
      items = items,
      alpha = alpha,
      language = language,
      digits = digits,
      polytomous = polytomous
    ),
    effective_settings = list(
      group_levels = levels(prepared$group),
      items = prepared$focal_items,
      alpha = alpha,
      polytomous = polytomous,
      min_group_n = 10L,
      language = language
    ),
    data_summary = list(
      n = nrow(prepared$data),
      n_items = length(prepared$focal_items),
      group_counts = table(prepared$group, useNA = "no")
    ),
    results = results,
    warnings = captured_warnings,
    errors = captured_errors,
    diagnostics = list(
      items = lapply(analyses, `[[`, "diagnostics"),
      category_map = prepared$category_map
    ),
    decision_log = .pm_decision(
      code = "dif_polytomous",
      requested = polytomous,
      effective = polytomous,
      reason = if (polytomous == "median") {
        "Explicit compatibility median dichotomization was requested."
      } else {
        "Polytomous items require a dedicated ordinal or IRT DIF method."
      }
    ),
    conclusion = conclusion
  )
  legacy <- list(
    results = results,
    interpretation = interpretation,
    summary = summary,
    alpha = alpha,
    group_levels = levels(prepared$group),
    language = language,
    polytomous = polytomous,
    references = .pm_dif_references(),
    call = call
  )
  for (field in names(legacy)) out[field] <- list(legacy[[field]])
  .pm_finish_analysis(out, output, language)
}

.pm_prepare_irt_data <- function(data) {
  if (inherits(data, "psychomatic_spec")) {
    item_names <- data$items
    data <- data$data[data$item_index]
    names(data) <- item_names
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame, matrix, or `psychomatic_spec`.", call. = FALSE)
  }
  data <- as.data.frame(data, check.names = FALSE, optional = TRUE)
  if (ncol(data) == 0L) stop("`data` must contain at least one item.", call. = FALSE)
  item_names <- names(data)
  if (is.null(item_names) || anyNA(item_names) || any(!nzchar(item_names))) {
    stop("IRT item names must be non-empty.", call. = FALSE)
  }
  if (anyDuplicated(item_names)) stop("IRT item names must not be duplicate.", call. = FALSE)

  normalized <- vector("list", length(data))
  maps <- vector("list", length(data))
  n_categories <- integer(length(data))
  names(normalized) <- names(maps) <- names(n_categories) <- item_names
  for (i in seq_along(data)) {
    item <- .pm_normalize_ordered_item(data[[i]], item_names[i], context = "IRT")
    normalized[[i]] <- item$values
    maps[[i]] <- item$map
    n_categories[i] <- item$n_categories
  }
  normalized <- as.data.frame(normalized, check.names = FALSE, optional = TRUE)
  names(normalized) <- item_names
  list(
    data = normalized,
    item_names = item_names,
    n_categories = n_categories,
    category_map = maps
  )
}

.pm_normalize_ordered_item <- function(x, item, context = "analysis") {
  if (is.logical(x) || is.character(x) || (is.factor(x) && !is.ordered(x))) {
    stop(
      sprintf("`%s` has an unordered or unsupported encoding for %s.", item, context),
      call. = FALSE
    )
  }
  if (!is.numeric(x) && !is.ordered(x)) {
    stop(sprintf("`%s` has an unsupported item type for %s.", item, context), call. = FALSE)
  }
  if (is.ordered(x)) {
    declared <- levels(x)
    observed_labels <- unique(as.character(stats::na.omit(x)))
    if (length(observed_labels) == 0L) {
      stop(sprintf("`%s` is entirely missing.", item), call. = FALSE)
    }
    if (length(observed_labels) < 2L) {
      stop(sprintf("`%s` is constant or has fewer than two observed categories.", item), call. = FALSE)
    }
    used <- declared %in% observed_labels
    if (any(!used)) {
      stop(
        sprintf("`%s` has empty ordered categories; category gaps cannot be normalized safely.", item),
        call. = FALSE
      )
    }
    requested <- declared
    values <- match(as.character(x), declared) - 1L
    values[is.na(x)] <- NA_integer_
  } else {
    if (any(is.nan(x) | (!is.na(x) & !is.finite(x)))) {
      stop(sprintf("`%s` must contain only finite or missing values.", item), call. = FALSE)
    }
    observed <- sort(unique(x[!is.na(x)]))
    if (!length(observed)) stop(sprintf("`%s` is entirely missing.", item), call. = FALSE)
    if (length(observed) < 2L) stop(sprintf("`%s` is constant.", item), call. = FALSE)
    if (any(observed != trunc(observed))) {
      stop(sprintf("`%s` must use integer response categories.", item), call. = FALSE)
    }
    requested <- observed
    values <- match(x, observed) - 1L
    values[is.na(x)] <- NA_integer_
  }
  map <- data.frame(
    requested = I(as.list(requested)),
    effective = seq_along(requested) - 1L,
    stringsAsFactors = FALSE
  )
  list(values = as.integer(values), map = map, n_categories = length(requested))
}

.pm_auto_irt_model <- function(n_categories) {
  if (!is.numeric(n_categories) || !length(n_categories) || anyNA(n_categories) ||
      any(n_categories < 2L)) {
    stop("`n_categories` must contain category counts of at least two.", call. = FALSE)
  }
  item_names <- names(n_categories)
  if (is.null(item_names)) item_names <- paste0("item", seq_along(n_categories))
  if (all(n_categories == 2L)) {
    return(list(model = "2pl", itemtype = stats::setNames(rep("2PL", length(n_categories)), item_names)))
  }
  if (all(n_categories > 2L)) {
    return(list(model = "grm", itemtype = stats::setNames(rep("graded", length(n_categories)), item_names)))
  }
  list(
    model = "mixed",
    itemtype = stats::setNames(ifelse(n_categories == 2L, "2PL", "graded"), item_names)
  )
}

.pm_normalize_itemtype <- function(itemtype, item_names) {
  if (!is.character(itemtype) || !length(itemtype) || anyNA(itemtype) ||
      any(!nzchar(trimws(itemtype)))) {
    stop("`itemtype` must contain supported non-empty mirt item types.", call. = FALSE)
  }
  if (length(itemtype) != length(item_names)) {
    stop("`itemtype` must have exactly one value per item; it is never recycled.", call. = FALSE)
  }
  if (!is.null(names(itemtype))) {
    if (anyNA(names(itemtype)) || any(!nzchar(names(itemtype))) || anyDuplicated(names(itemtype)) ||
        !setequal(names(itemtype), item_names)) {
      stop("A named `itemtype` vector must have unique names matching every item name.", call. = FALSE)
    }
    itemtype <- itemtype[match(item_names, names(itemtype))]
  } else {
    names(itemtype) <- item_names
  }
  normalized <- tolower(trimws(itemtype))
  normalized[normalized == "rasch"] <- "Rasch"
  normalized[normalized == "2pl"] <- "2PL"
  normalized[normalized %in% c("graded", "grm")] <- "graded"
  supported <- normalized %in% c("Rasch", "2PL", "graded")
  if (!all(supported)) {
    stop("`itemtype` supports only Rasch, 2PL, and graded in `irt_auto()`.", call. = FALSE)
  }
  stats::setNames(normalized, item_names)
}

.pm_resolve_irt_model <- function(model, itemtype, n_categories, item_names) {
  automatic <- .pm_auto_irt_model(n_categories)
  if (is.null(itemtype)) {
    if (model == "auto") {
      return(c(automatic, list(
        display_itemtype = if (length(unique(automatic$itemtype)) == 1L) unname(automatic$itemtype[1L]) else automatic$itemtype,
        reason = "Automatic item types were selected from the observed category counts."
      )))
    }
    types <- switch(
      model,
      rasch = stats::setNames(rep("Rasch", length(item_names)), item_names),
      `2pl` = stats::setNames(rep("2PL", length(item_names)), item_names),
      grm = stats::setNames(rep("graded", length(item_names)), item_names)
    )
  } else {
    types <- .pm_normalize_itemtype(itemtype, item_names)
  }

  binary <- n_categories == 2L
  incompatible <- (types %in% c("Rasch", "2PL") & !binary) |
    (types == "graded" & binary)
  if (any(incompatible)) {
    stop(
      "The requested `model`/`itemtype` is incompatible with the observed item categories.",
      call. = FALSE
    )
  }
  if (model != "auto") {
    allowed <- switch(model, rasch = "Rasch", `2pl` = "2PL", grm = "graded")
    if (any(types != allowed)) {
      stop("The requested `model` and `itemtype` are incompatible.", call. = FALSE)
    }
    effective_model <- model
  } else if (length(unique(types)) == 1L) {
    effective_model <- switch(types[1L], Rasch = "rasch", `2PL` = "2pl", graded = "grm")
  } else {
    effective_model <- "mixed"
  }
  list(
    model = effective_model,
    itemtype = types,
    display_itemtype = if (length(unique(types)) == 1L) unname(types[1L]) else types,
    reason = if (model == "auto") {
      "User-supplied item types were validated against each item's categories."
    } else {
      "The explicit model and item types were validated against each item's categories."
    }
  )
}

.pm_discrimination_components <- function(row) {
  if (is.data.frame(row)) row <- unlist(row[1L, , drop = TRUE], use.names = TRUE)
  slope_names <- grep("^a[0-9]+$", names(row), value = TRUE)
  dimensions <- as.numeric(row[slope_names])
  names(dimensions) <- slope_names
  finite <- is.finite(dimensions)
  magnitude <- if (!length(dimensions) || !all(finite)) NA_real_ else sqrt(sum(dimensions^2))
  list(dimensions = dimensions, magnitude = magnitude)
}

.pm_add_discrimination_magnitude <- function(item_parameters) {
  if (!is.data.frame(item_parameters) || !nrow(item_parameters)) return(item_parameters)
  item_parameters$discrimination_magnitude <- vapply(
    seq_len(nrow(item_parameters)),
    function(i) .pm_discrimination_components(item_parameters[i, , drop = FALSE])$magnitude,
    numeric(1)
  )
  item_parameters
}

.pm_merge_item_parameter_aliases <- function(raw, legacy) {
  if (!nrow(raw)) return(legacy)
  if (!nrow(legacy)) return(raw)
  if (!identical(raw$item, legacy$item)) {
    legacy <- legacy[match(raw$item, legacy$item), , drop = FALSE]
  }
  alias_names <- setdiff(names(legacy), names(raw))
  for (name in alias_names) raw[[name]] <- legacy[[name]]
  raw
}

.pm_assess_irt_solution <- function(converged,
                                    item_parameters,
                                    error = NULL,
                                    itemtype = NULL,
                                    n_categories = NULL) {
  if (!is.null(error)) return(list(admissible = FALSE, reason = conditionMessage(error)))
  if (!isTRUE(converged)) return(list(admissible = FALSE, reason = "The mirt optimizer did not converge."))
  if (!is.data.frame(item_parameters) || !nrow(item_parameters)) {
    return(list(admissible = FALSE, reason = "No estimable item parameters were returned."))
  }
  slope_names <- grep("^a[0-9]+$", names(item_parameters), value = TRUE)
  intercept_names <- grep("^d[0-9]*$", names(item_parameters), value = TRUE)
  if (!length(slope_names) || !length(intercept_names)) {
    return(list(admissible = FALSE, reason = "No authoritative slope/intercept parameters were returned."))
  }

  item_names <- if ("item" %in% names(item_parameters)) {
    as.character(item_parameters$item)
  } else {
    paste0("item", seq_len(nrow(item_parameters)))
  }
  if (is.null(itemtype)) {
    slopes <- as.matrix(item_parameters[slope_names])
    storage.mode(slopes) <- "double"
    intercepts <- as.matrix(item_parameters[intercept_names])
    storage.mode(intercepts) <- "double"
    invalid <- any(!is.finite(slopes)) || any(apply(intercepts, 1L, function(x) {
      finite <- is.finite(x)
      !any(finite) || any(is.infinite(x))
    }))
    if (invalid) {
      return(list(admissible = FALSE, reason = "The fitted solution contains nonfinite item parameters."))
    }
    return(list(admissible = TRUE, reason = NA_character_))
  }

  if (length(itemtype) != nrow(item_parameters)) {
    return(list(admissible = FALSE, reason = "Effective item types did not align with fitted item parameters."))
  }
  if (!is.null(names(itemtype))) itemtype <- itemtype[match(item_names, names(itemtype))]
  if (is.null(n_categories) || length(n_categories) != nrow(item_parameters)) {
    return(list(admissible = FALSE, reason = "Effective category counts did not align with fitted item parameters."))
  }
  if (!is.null(names(n_categories))) {
    n_categories <- n_categories[match(item_names, names(n_categories))]
  }

  invalid_item <- vapply(seq_len(nrow(item_parameters)), function(i) {
    row <- item_parameters[i, , drop = FALSE]
    slopes <- as.numeric(row[1L, slope_names, drop = TRUE])
    if (!length(slopes) || any(!is.finite(slopes))) return(TRUE)
    type <- unname(itemtype[i])
    category_count <- as.integer(n_categories[i])
    required <- if (type %in% c("Rasch", "2PL")) {
      c("d", intersect(c("g", "u"), names(row)))
    } else if (identical(type, "graded") && is.finite(category_count) && category_count >= 3L) {
      paste0("d", seq_len(category_count - 1L))
    } else {
      character()
    }
    if (!length(required) || !all(required %in% names(row))) return(TRUE)
    values <- as.numeric(row[1L, required, drop = TRUE])
    any(!is.finite(values))
  }, logical(1))
  if (any(invalid_item)) {
    return(list(admissible = FALSE, reason = "The fitted solution contains nonfinite item parameters."))
  }
  list(admissible = TRUE, reason = NA_character_)
}

.pm_mirt_fit_statistics <- function(fit, converged) {
  if (!isTRUE(converged)) {
    return(list(
      status = "not_evaluable",
      reason = "M2 was not computed because the model did not converge.",
      warnings = list(),
      error = NULL
    ))
  }
  captured <- .pm_capture_conditions(mirt::M2(fit))
  if (!is.null(captured$error) || is.null(captured$value)) {
    return(list(
      status = "not_available",
      reason = if (is.null(captured$error)) "M2 was unavailable." else conditionMessage(captured$error),
      warnings = captured$warnings,
      error = captured$error
    ))
  }
  values <- captured$value
  numeric_values <- unlist(values[vapply(values, is.numeric, logical(1))], use.names = TRUE)
  if (!length(numeric_values) || any(!is.finite(numeric_values))) {
    return(list(
      status = "not_available",
      reason = "M2 returned nonfinite diagnostics.",
      warnings = captured$warnings,
      error = NULL
    ))
  }
  list(
    status = "available",
    M2 = values,
    reason = NA_character_,
    warnings = captured$warnings,
    error = NULL
  )
}

#' Compare alternative CFA models
#'
#' Fits a named set of lavaan CFA models and returns a compact comparison table.
#' The recommendation is descriptive only and favors converged models with
#' stronger comparative fit indices. Information criteria are reported when
#' available, but they are not required and are not used for WLSMV-style
#' categorical estimators.
#'
#' @inheritParams efa_auto
#' @param data Data frame containing observed variables.
#' @param models Named list or named character vector of lavaan model strings.
#' @param ordered Optional ordered variables passed to `lavaan::cfa()`. If
#'   `NULL` and a categorical estimator such as WLSMV is used, observed items
#'   referenced in `models` are treated as ordered.
#' @param estimator Estimator passed to `lavaan::cfa()`. `NULL` (default) or
#'   `"auto"` uses the same category-count, sample-size and Mardia criteria as
#'   `cfa_auto()`, applied to the observed indicators across the candidate models.
#' @param fit_variant Fit-index variant for both reporting and comparison:
#'   `"scaled"` (default), `"robust"`, or `"base"`. SRMR remains base.
#'   Uncorrected estimators use base indices with the default scaled setting.
#'   Otherwise unavailable requested indices remain missing, without substitution.
#' @param language Recommendation language: `"eng"` or `"esp"`.
#'
#' @return A list of class `psychomatic_model_comparison`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' models <- list(one_factor = "F =~ R1 + R2 + R3 + R4 + R5 + R6",
#'                two_factors = "F1 =~ R1 + R2 + R3\nF2 =~ R4 + R5 + R6")
#' compare_models_auto(MY_DATA, models = models)
#' compare_models_auto(MY_DATA, models = models, fit_variant = "robust")
#' }
#'
#' @export
compare_models_auto <- function(data,
                                models,
                                ordered = NULL,
                                estimator = NULL,
                                language = c("eng", "esp"),
                                digits = 3, report = FALSE, file_name = NULL,
                                format = "excel", overwrite = FALSE, verbose = TRUE,
                                fit_variant = c("scaled", "robust", "base")) {
  fit_variant <- match.arg(fit_variant)
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  call <- match.call()
  requested_ordered <- ordered
  requested_estimator <- estimator
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package `lavaan` is required for `compare_models_auto()`.", call. = FALSE)
  }
  language <- match.arg(language)
  if (is.character(models)) {
    models <- as.list(models)
  }
  if (!is.list(models) || length(models) == 0L || is.null(names(models)) ||
      anyNA(names(models)) || any(!nzchar(trimws(names(models))))) {
    stop("`models` must be a named list or named character vector.", call. = FALSE)
  }
  invalid_syntax_values <- !vapply(
    models,
    function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x)),
    logical(1)
  )
  if (any(invalid_syntax_values)) {
    stop("Each entry in `models` must be one non-empty lavaan syntax string.", call. = FALSE)
  }
  automatic_estimator <- is.null(estimator) || identical(tolower(estimator), "auto")
  if (!automatic_estimator && (!is.character(estimator) || length(estimator) != 1L || is.na(estimator) ||
      !nzchar(trimws(estimator)))) {
    stop("`estimator` must be one non-empty character value.", call. = FALSE)
  }
  if (!automatic_estimator) estimator <- toupper(estimator)
  if (inherits(data, "psychomatic_spec")) data <- data$data
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame, matrix, or `psychomatic_spec`.", call. = FALSE)
  }
  data <- as.data.frame(data, check.names = FALSE, optional = TRUE)
  if (is.null(names(data)) || anyNA(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("`data` must have unique, non-empty column names.", call. = FALSE)
  }
  estimator_reason <- if (language == "esp") "Estimador especificado por el usuario" else "Estimator specified by the user"
  if (automatic_estimator) {
    observed <- .pm_model_observed_variables(models, data)
    if (!length(observed)) stop("No observed indicators found for automatic estimator selection.", call. = FALSE)
    items <- data[, observed, drop = FALSE]
    categories <- vapply(items, function(v) length(unique(stats::na.omit(v))), integer(1))
    complete <- sum(stats::complete.cases(items))
    diagnostic <- list(status = "unavailable")
    if (max(categories) > 7 || (max(categories) > 5 && min(categories) > 5)) {
      numeric_items <- as.data.frame(lapply(items, function(v) if (is.ordered(v)) as.integer(v) else v))
      diagnostic <- .pm_mvn_mardia(numeric_items[stats::complete.cases(numeric_items), , drop=FALSE], alpha = .05)
    }
    selection <- .pm_select_cfa_estimator(max(categories), min(categories), complete,
      diagnostic, ordered = ordered, language = language)
    estimator <- selection$estimator
    estimator_reason <- selection$reason
    if (is.null(ordered) && isTRUE(selection$ordered)) ordered <- observed
    if (length(ordered) && !estimator %in% c("WLSMV", "ULSMV")) {
      estimator <- "WLSMV"
      estimator_reason <- if (language == "esp") "WLSMV: indicadores ordinales especificados por el usuario." else "WLSMV: ordered indicators specified by the user."
    }
  }
  if (is.null(ordered) && estimator %in% c("WLSMV", "ULSMV")) {
    ordered <- .pm_model_observed_variables(models, data)
  }
  if (!is.null(ordered)) {
    if (!is.character(ordered) || anyNA(ordered) || any(!nzchar(ordered)) || anyDuplicated(ordered)) {
      stop("`ordered` must contain unique, non-empty column names.", call. = FALSE)
    }
    missing_ordered <- setdiff(ordered, names(data))
    if (length(missing_ordered)) {
      stop("Ordered variables not found in `data`: ", paste(missing_ordered, collapse = ", "), call. = FALSE)
    }
  }

  model_ids <- sprintf("model_%03d", seq_along(models))
  display_names <- names(models)
  fits <- stats::setNames(vector("list", length(models)), model_ids)
  diagnostics <- stats::setNames(vector("list", length(models)), model_ids)
  condition_warnings <- list()
  condition_errors <- list()
  rows <- lapply(seq_along(models), function(i) {
    model_id <- model_ids[i]
    model_name <- display_names[i]
    model_syntax <- models[[i]]
    parsed <- .pm_capture_conditions(lavaan::lavaanify(model_syntax))
    condition_warnings <<- c(condition_warnings, parsed$warnings)
    if (!is.null(parsed$error)) {
      condition_errors[[length(condition_errors) + 1L]] <<- parsed$error
      diagnostics[[model_id]] <<- list(
        model_id = model_id,
        model = model_name,
        syntax = model_syntax,
        observed_variables = character(),
        case_indices = integer(),
        status = "not_evaluable",
        reason = conditionMessage(parsed$error)
      )
      return(.pm_empty_model_fit_row(
        model_id,
        model_name,
        estimator,
        conditionMessage(parsed$error)
      ))
    }
    fit_args <- list(model = model_syntax, data = data, estimator = estimator)
    if (!is.null(ordered)) fit_args$ordered <- ordered
    fit_capture <- .pm_capture_conditions(do.call(lavaan::cfa, fit_args))
    condition_warnings <<- c(condition_warnings, fit_capture$warnings)
    fit <- fit_capture$value
    fits[model_id] <<- list(fit)
    if (!is.null(fit_capture$error)) {
      condition_errors[[length(condition_errors) + 1L]] <<- fit_capture$error
      diagnostics[[model_id]] <<- list(
        model_id = model_id,
        model = model_name,
        syntax = model_syntax,
        observed_variables = .pm_observed_from_parameter_table(parsed$value, names(data)),
        case_indices = integer(),
        status = "not_evaluable",
        reason = conditionMessage(fit_capture$error)
      )
      return(.pm_empty_model_fit_row(
        model_id,
        model_name,
        estimator,
        conditionMessage(fit_capture$error)
      ))
    }

    record <- .pm_extract_lavaan_model_record(
      fit,
      model_id = model_id,
      model_name = model_name,
      requested_estimator = estimator,
      syntax = model_syntax
    )
    condition_warnings <<- c(condition_warnings, record$warnings)
    condition_errors <<- c(condition_errors, record$errors)
    diagnostics[[model_id]] <<- record$diagnostics
    record$row
  })

  fit_table <- do.call(rbind, rows)
  rownames(fit_table) <- NULL
  comparison <- .pm_apply_model_comparability(fit_table, diagnostics)
  fit_table <- comparison$fit
  if (length(comparison$indices)) {
    variants <- .pm_select_common_fit_variants(fit_table[comparison$indices, , drop = FALSE], fit_variant)
    for (metric in names(variants$values)) {
      fit_table[[paste0("comparison_", metric)]][comparison$indices] <- variants$values[[metric]]
      fit_table[[paste0("comparison_", metric, "_variant")]][comparison$indices] <- variants$variants[[metric]]
    }
    fit_table$fit_score[comparison$indices] <- mapply(
      .pm_fit_score,
      fit_table$comparison_cfi[comparison$indices],
      fit_table$comparison_tli[comparison$indices],
      fit_table$comparison_rmsea[comparison$indices],
      fit_table$comparison_srmr[comparison$indices]
    )
  }
  recommendation <- .pm_model_recommendation(fit_table, language)
  effective_model_settings <- .pm_summarize_model_effective_settings(
    diagnostics,
    model_ids
  )

  conclusion <- if (identical(recommendation$basis, "none")) {
    .pm_not_evaluable(recommendation$message)
  } else {
    list(status = "evaluable", reason = NA_character_)
  }
  out <- .pm_result(
    class = "psychomatic_model_comparison",
    call = call,
    requested_settings = list(
      models = display_names,
      fit_variant = fit_variant,
      ordered = requested_ordered,
      estimator = requested_estimator,
      language = language
    ),
    effective_settings = list(
      model_ids = model_ids,
      fit_variant = fit_variant,
      ordered = effective_model_settings$ordered,
      estimator = effective_model_settings$estimator_original,
      estimator_selection = if (automatic_estimator) "automatic" else "user",
      estimator_reason = estimator_reason,
      estimator_original = effective_model_settings$estimator_original,
      estimator_computational = effective_model_settings$estimator_computational,
      meanstructure = effective_model_settings$meanstructure,
      estimand = effective_model_settings$estimand,
      by_model = effective_model_settings$by_model,
      comparison_reason = comparison$reason
    ),
    data_summary = list(n = nrow(data), p = ncol(data)),
    results = fit_table,
    warnings = condition_warnings,
    errors = condition_errors,
    diagnostics = list(models = diagnostics, comparison = comparison),
    decision_log = .pm_decision(
      code = "model_comparison",
      requested = paste(display_names, collapse = ", "),
      effective = if (length(comparison$indices)) {
        paste(model_ids[comparison$indices], collapse = ", ")
      } else {
        "none"
      },
      reason = comparison$reason
    ),
    conclusion = conclusion
  )
  legacy <- list(
    fit = fit_table,
    recommendation = recommendation,
    fitted_models = fits,
    estimator = estimator,
    ordered = ordered,
    call = call
  )
  for (field in names(legacy)) out[field] <- list(legacy[[field]])
  .pm_finish_analysis(out, output, language)
}

.pm_mirt_item_parameters <- function(coefficients) {
  items <- if (is.list(coefficients) && "items" %in% names(coefficients)) {
    coefficients$items
  } else {
    coefficients
  }
  out <- as.data.frame(items)
  out$item <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("item", setdiff(names(out), "item")), drop = FALSE]
  out
}

.pm_irt_interpretation <- function(item_parameters,
                                   model,
                                   itemtype,
                                   language = "eng",
                                   digits = 3) {
  if (!is.data.frame(item_parameters) || nrow(item_parameters) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(item_parameters)), function(i) {
    row <- item_parameters[i, , drop = FALSE]
    discrimination <- .pm_item_discrimination(row)
    thresholds <- .pm_item_thresholds(row)
    difficulty_center <- if (length(thresholds) > 0L) {
      mean(thresholds, na.rm = TRUE)
    } else {
      .pm_item_difficulty(row)
    }
    threshold_min <- if (length(thresholds) > 0L) min(thresholds, na.rm = TRUE) else difficulty_center
    threshold_max <- if (length(thresholds) > 0L) max(thresholds, na.rm = TRUE) else difficulty_center
    threshold_ordered <- if (length(thresholds) <= 1L) TRUE else all(diff(thresholds) > 0)

    data.frame(
      item = as.character(row$item),
      discrimination = round(discrimination, digits),
      discrimination_level = .pm_discrimination_level(discrimination, language),
      difficulty_center = round(difficulty_center, digits),
      threshold_min = round(threshold_min, digits),
      threshold_max = round(threshold_max, digits),
      threshold_ordered = threshold_ordered,
      recommendation = .pm_irt_recommendation(
        discrimination = discrimination,
        difficulty_center = difficulty_center,
        threshold_min = threshold_min,
        threshold_max = threshold_max,
        threshold_ordered = threshold_ordered,
        language = language
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.pm_irt_model_summary <- function(interpretation,
                                  model,
                                  itemtype,
                                  n,
                                  language = "eng",
                                  digits = 3,
                                  raw_discrimination = NULL) {
  itemtype_label <- if (!length(itemtype)) {
    NA_character_
  } else {
    paste(unique(unname(as.character(itemtype))), collapse = ", ")
  }
  if (!is.data.frame(interpretation) || nrow(interpretation) == 0L) {
    return(data.frame(
      model = model,
      itemtype = itemtype_label,
      n = n,
      n_items = 0L,
      items_to_review = 0L,
      mean_discrimination = NA_real_,
      recommendation = if (language == "esp") {
        "No hay parametros interpretables disponibles."
      } else {
        "No interpretable item parameters were available."
      },
      stringsAsFactors = FALSE
    ))
  }

  review_flags <- !grepl(
    if (language == "esp") "^Sin alerta principal" else "^No major flag",
    interpretation$recommendation
  )
  items_to_review <- sum(review_flags, na.rm = TRUE)
  discrimination_evidence <- if (is.null(raw_discrimination)) {
    interpretation$discrimination
  } else {
    if (!is.numeric(raw_discrimination) || length(raw_discrimination) != nrow(interpretation)) {
      stop("`raw_discrimination` must have one numeric value per interpreted item.", call. = FALSE)
    }
    raw_discrimination
  }
  mean_discrimination <- mean(discrimination_evidence, na.rm = TRUE)
  recommendation <- .pm_irt_overall_recommendation(
    items_to_review = items_to_review,
    n_items = nrow(interpretation),
    mean_discrimination = mean_discrimination,
    language = language
  )

  data.frame(
    model = model,
    itemtype = itemtype_label,
    n = n,
    n_items = nrow(interpretation),
    items_to_review = items_to_review,
    mean_discrimination = round(mean_discrimination, digits),
    recommendation = recommendation,
    stringsAsFactors = FALSE
  )
}

.pm_irt_overall_recommendation <- function(items_to_review,
                                           n_items,
                                           mean_discrimination,
                                           language) {
  review_pct <- if (n_items > 0L) items_to_review / n_items else 1
  if (language == "esp") {
    if (items_to_review == 0L && is.finite(mean_discrimination) && mean_discrimination >= 0.65) {
      return("Sin alertas principales; revisar contenido e informacion del test antes de concluir.")
    }
    if (review_pct <= 0.25) {
      return("Pocos items requieren revision; priorice contenido, discriminacion y umbrales.")
    }
    return("Varios items requieren revision; use estos resultados como cribado antes de decisiones definitivas.")
  }
  if (items_to_review == 0L && is.finite(mean_discrimination) && mean_discrimination >= 0.65) {
    return("No major item flags; review item content and test information before concluding.")
  }
  if (review_pct <= 0.25) {
    return("Few items require review; prioritize content, discrimination, and thresholds.")
  }
  "Several items require review; treat these results as screening evidence before final decisions."
}

.pm_item_discrimination <- function(row) {
  if ("discrimination_magnitude" %in% names(row)) {
    magnitude <- as.numeric(row[1L, "discrimination_magnitude", drop = TRUE])
    if (length(magnitude) == 1L && is.finite(magnitude)) return(magnitude)
  }
  candidates <- c(
    intersect("a", names(row)),
    grep("^a[0-9]+$", names(row), value = TRUE)
  )
  values <- as.numeric(row[1, candidates, drop = TRUE])
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  sqrt(sum(values^2))
}

.pm_item_difficulty <- function(row) {
  candidates <- intersect(c("b", "b1"), names(row))
  values <- as.numeric(row[1, candidates, drop = TRUE])
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  values[1L]
}

.pm_item_thresholds <- function(row) {
  threshold_names <- grep("^b[0-9]+$", names(row), value = TRUE)
  values <- as.numeric(row[1, threshold_names, drop = TRUE])
  values[is.finite(values)]
}

.pm_discrimination_level <- function(discrimination, language) {
  if (!is.finite(discrimination)) {
    return(if (language == "esp") "No disponible" else "Not available")
  }
  if (discrimination < 0.65) {
    if (language == "esp") "Baja" else "Low"
  } else if (discrimination < 1.35) {
    if (language == "esp") "Moderada" else "Moderate"
  } else if (discrimination < 1.70) {
    if (language == "esp") "Alta" else "High"
  } else {
    if (language == "esp") "Muy alta" else "Very high"
  }
}

.pm_irt_recommendation <- function(discrimination,
                                   difficulty_center,
                                   threshold_min,
                                   threshold_max,
                                   threshold_ordered,
                                   language) {
  issues <- character(0)
  if (!is.finite(discrimination) || discrimination < 0.65) {
    issues <- c(issues, if (language == "esp") "discriminacion baja" else "low discrimination")
  }
  if (is.finite(difficulty_center) && abs(difficulty_center) > 3) {
    issues <- c(issues, if (language == "esp") "dificultad extrema" else "extreme difficulty")
  }
  if (is.finite(threshold_min) && is.finite(threshold_max) &&
      (threshold_min < -4 || threshold_max > 4)) {
    issues <- c(issues, if (language == "esp") "umbrales extremos" else "extreme thresholds")
  }
  if (!isTRUE(threshold_ordered)) {
    issues <- c(issues, if (language == "esp") "umbrales desordenados" else "disordered thresholds")
  }
  if (length(issues) == 0L) {
    if (language == "esp") {
      "Sin alerta principal; revisar junto con contenido y precision del test."
    } else {
      "No major flag; review together with item content and test information."
    }
  } else if (language == "esp") {
    paste("Revisar item por", paste(unique(issues), collapse = ", "), ".")
  } else {
    paste("Review item for", paste(unique(issues), collapse = ", "), ".")
  }
}

.pm_prepare_dif_data <- function(data, group, items) {
  specification <- if (inherits(data, "psychomatic_spec")) data else NULL
  if (!is.null(specification)) data <- specification$data
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame, matrix, or `psychomatic_spec`.", call. = FALSE)
  }
  data <- as.data.frame(data, check.names = FALSE, optional = TRUE)
  if (!ncol(data) || is.null(names(data)) || anyNA(names(data)) ||
      any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("`data` must have unique, non-empty column names.", call. = FALSE)
  }

  group_name <- NULL
  if (is.character(group) && length(group) == 1L && !is.na(group) && group %in% names(data)) {
    group_name <- group
    group_vector <- data[[group_name]]
  } else {
    group_vector <- group
  }
  if (length(group_vector) != nrow(data)) {
    stop("`group` must have the same length as the number of rows in `data`.", call. = FALSE)
  }
  if (is.data.frame(group_vector) || is.matrix(group_vector) || is.list(group_vector)) {
    stop("`group` must be an atomic grouping vector or one column name.", call. = FALSE)
  }
  group_factor <- droplevels(factor(group_vector))
  if (nlevels(group_factor) < 2L) stop("`group` must contain at least two groups.", call. = FALSE)

  candidate_names <- setdiff(names(data), group_name)
  eligible <- vapply(
    data[candidate_names],
    function(x) is.numeric(x) || is.ordered(x),
    logical(1)
  )
  pool_names <- candidate_names[eligible]
  if (is.null(items)) {
    if (!is.null(specification)) {
      focal_items <- intersect(specification$items, pool_names)
    } else {
      focal_items <- pool_names
    }
  } else {
    if (!is.character(items) || !length(items) || anyNA(items) ||
        any(!nzchar(items)) || anyDuplicated(items)) {
      stop("`items` must contain unique, non-empty item names.", call. = FALSE)
    }
    missing_items <- setdiff(items, names(data))
    if (length(missing_items)) {
      stop("Items not found in `data`: ", paste(missing_items, collapse = ", "), call. = FALSE)
    }
    unsupported <- items[!items %in% pool_names]
    if (length(unsupported)) {
      stop("Requested DIF items must be numeric or ordered: ", paste(unsupported, collapse = ", "), call. = FALSE)
    }
    focal_items <- items
  }
  if (!length(focal_items)) stop("No numeric or ordered items were available for DIF screening.", call. = FALSE)

  # The matching score is defined by the requested/specification item set.
  # Numeric metadata outside that set must never enter a focal item's score.
  matching_names <- focal_items
  normalized <- vector("list", length(matching_names))
  category_map <- vector("list", length(matching_names))
  names(normalized) <- names(category_map) <- matching_names
  for (i in seq_along(matching_names)) {
    name <- matching_names[i]
    captured <- tryCatch(
      .pm_normalize_ordered_item(data[[name]], name, context = "DIF"),
      error = identity
    )
    if (inherits(captured, "error")) {
      normalized[[i]] <- list(values = rep(NA_integer_, nrow(data)), n_categories = 0L,
                              reason = conditionMessage(captured), valid = FALSE)
      category_map[[i]] <- data.frame()
    } else {
      normalized[[i]] <- c(captured, list(reason = NA_character_, valid = TRUE))
      category_map[[i]] <- captured$map
    }
  }
  list(
    data = data,
    group = group_factor,
    group_name = group_name,
    focal_items = focal_items,
    normalized_items = normalized,
    category_map = category_map
  )
}

.pm_dif_outcome <- function(item, polytomous) {
  if (!isTRUE(item$valid)) {
    return(list(status = "not_evaluable", reason = item$reason, y = rep(NA_integer_, length(item$values)),
                method = "not_evaluable", note = item$reason))
  }
  if (item$n_categories == 2L) {
    return(list(
      status = "evaluable",
      reason = NA_character_,
      y = as.integer(item$values == 1L),
      method = "binary_logistic",
      note = ""
    ))
  }
  if (polytomous == "not_evaluable") {
    reason <- "Polytomous or ordinal items require a dedicated ordinal or IRT DIF method."
    return(list(status = "not_evaluable", reason = reason, y = rep(NA_integer_, length(item$values)),
                method = "not_evaluable", note = reason))
  }
  cutoff <- stats::median(item$values, na.rm = TRUE)
  list(
    status = "evaluable",
    reason = NA_character_,
    y = as.integer(item$values > cutoff),
    method = "median_dichotomized",
    note = "Methodological warning: the polytomous item was median-dichotomized; confirm with a dedicated ordinal or IRT DIF method."
  )
}

.pm_empty_dif_row <- function(item, method, note, reason) {
  out <- data.frame(
    item = item,
    status = "not_evaluable",
    reason = reason,
    p_uniform = NA_real_,
    p_nonuniform = NA_real_,
    lr_uniform = NA_real_,
    df_uniform = NA_real_,
    lr_nonuniform = NA_real_,
    df_nonuniform = NA_real_,
    nagelkerke_base = NA_real_,
    nagelkerke_uniform = NA_real_,
    nagelkerke_nonuniform = NA_real_,
    delta_r2_uniform = NA_real_,
    delta_r2_nonuniform = NA_real_,
    effect_size = NA_real_,
    flagged = NA,
    method = method,
    note = note,
    n_used = 0L,
    base_converged = FALSE,
    uniform_converged = FALSE,
    nonuniform_converged = FALSE,
    stringsAsFactors = FALSE
  )
  out$original_indices <- I(list(integer()))
  out$base_indices <- I(list(integer()))
  out$uniform_indices <- I(list(integer()))
  out$nonuniform_indices <- I(list(integer()))
  out$group_counts <- I(list(integer()))
  out$coefficients <- I(list(stats::setNames(numeric(), character())))
  out$odds_ratios <- I(list(stats::setNames(numeric(), character())))
  out
}

.pm_nagelkerke_r2 <- function(model) {
  ll_model <- as.numeric(stats::logLik(model))
  response <- model$y
  weights <- model$prior.weights
  if (is.null(response) || is.null(weights) || !length(response) || length(response) != length(weights)) {
    return(NA_real_)
  }
  probability <- stats::weighted.mean(response, weights)
  if (!is.finite(probability) || probability <= 0 || probability >= 1) return(NA_real_)
  ll_null <- sum(weights * stats::dbinom(response, size = 1, prob = probability, log = TRUE))
  n <- stats::nobs(model)
  denominator <- 1 - exp(2 * ll_null / n)
  if (!is.finite(ll_model) || !is.finite(ll_null) || !is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  (1 - exp((2 / n) * (ll_null - ll_model))) / denominator
}

.pm_glm_solution_reason <- function(models, warnings) {
  if (!all(vapply(models, function(x) inherits(x, "glm"), logical(1)))) return("A logistic DIF model failed.")
  if (!all(vapply(models, function(x) isTRUE(x$converged), logical(1)))) {
    return("At least one logistic DIF model did not converge.")
  }
  coefficients <- unlist(lapply(models, stats::coef), use.names = FALSE)
  likelihoods <- vapply(models, function(x) as.numeric(stats::logLik(x)), numeric(1))
  if (any(!is.finite(coefficients)) || any(!is.finite(likelihoods))) {
    return("The logistic DIF solution contains a singular or nonfinite coefficient or likelihood.")
  }
  warning_text <- paste(vapply(warnings, conditionMessage, character(1)), collapse = " ")
  fitted <- unlist(lapply(models, stats::fitted), use.names = FALSE)
  separated <- any(fitted < 1e-8 | fitted > 1 - 1e-8) || any(abs(coefficients) > 20) ||
    grepl("0 or 1|did not converge|separation", warning_text, ignore.case = TRUE)
  if (separated) return("Complete or quasi separation produced degenerate fitted probabilities.")
  NA_character_
}

.pm_fit_logistic_dif <- function(analysis_data) {
  captures <- list(
    base = .pm_capture_conditions(stats::glm(y ~ score, data = analysis_data, family = stats::binomial())),
    uniform = .pm_capture_conditions(stats::glm(y ~ score + group, data = analysis_data, family = stats::binomial())),
    nonuniform = .pm_capture_conditions(stats::glm(y ~ score * group, data = analysis_data, family = stats::binomial()))
  )
  errors <- Filter(Negate(is.null), lapply(captures, `[[`, "error"))
  warnings <- unlist(lapply(captures, `[[`, "warnings"), recursive = FALSE)
  models <- lapply(captures, `[[`, "value")
  if (length(errors)) {
    return(list(status = "not_evaluable", reason = conditionMessage(errors[[1L]]),
                models = models, warnings = warnings, errors = errors))
  }
  reason <- .pm_glm_solution_reason(models, warnings)
  if (!is.na(reason)) {
    return(list(status = "not_evaluable", reason = reason, models = models,
                warnings = warnings, errors = errors))
  }
  lr_uniform <- .pm_capture_conditions(stats::anova(models$base, models$uniform, test = "Chisq"))
  lr_nonuniform <- .pm_capture_conditions(stats::anova(models$uniform, models$nonuniform, test = "Chisq"))
  warnings <- c(warnings, lr_uniform$warnings, lr_nonuniform$warnings)
  errors <- c(errors, Filter(Negate(is.null), list(lr_uniform$error, lr_nonuniform$error)))
  if (length(errors)) {
    return(list(status = "not_evaluable", reason = conditionMessage(errors[[1L]]), models = models,
                warnings = warnings, errors = errors))
  }
  uniform_table <- lr_uniform$value
  nonuniform_table <- lr_nonuniform$value
  lr_u <- as.numeric(uniform_table$Deviance[2L])
  df_u <- as.numeric(uniform_table$Df[2L])
  p_u <- as.numeric(uniform_table$`Pr(>Chi)`[2L])
  lr_n <- as.numeric(nonuniform_table$Deviance[2L])
  df_n <- as.numeric(nonuniform_table$Df[2L])
  p_n <- as.numeric(nonuniform_table$`Pr(>Chi)`[2L])
  if (any(!is.finite(c(lr_u, df_u, p_u, lr_n, df_n, p_n))) || df_u <= 0 || df_n <= 0) {
    return(list(status = "not_evaluable", reason = "The LR comparison was nonfinite or not correctly nested.",
                models = models, warnings = warnings, errors = errors))
  }
  list(
    status = "evaluable",
    reason = NA_character_,
    models = models,
    warnings = warnings,
    errors = errors,
    lr_uniform = lr_u,
    df_uniform = df_u,
    p_uniform = p_u,
    lr_nonuniform = lr_n,
    df_nonuniform = df_n,
    p_nonuniform = p_n
  )
}

.pm_analyse_dif_item <- function(item,
                                 normalized_items,
                                 group,
                                 alpha,
                                 polytomous,
                                 original_indices,
                                 min_group_n) {
  focal <- normalized_items[[item]]
  outcome <- .pm_dif_outcome(focal, polytomous)
  if (outcome$status != "evaluable") {
    return(list(
      row = .pm_empty_dif_row(item, outcome$method, outcome$note, outcome$reason),
      warnings = list(), errors = list(), diagnostics = list(status = outcome$status, reason = outcome$reason)
    ))
  }
  rest_names <- setdiff(names(normalized_items), item)
  usable_rest <- rest_names[vapply(
    normalized_items[rest_names],
    function(x) isTRUE(x$valid) && x$n_categories >= 2L,
    logical(1)
  )]
  if (!length(usable_rest)) {
    reason <- "No usable rest item was available for a focal-excluded matching score."
    return(list(row = .pm_empty_dif_row(item, outcome$method, outcome$note, reason),
                warnings = list(), errors = list(), diagnostics = list(status = "not_evaluable", reason = reason)))
  }
  rest_data <- as.data.frame(lapply(normalized_items[usable_rest], `[[`, "values"))
  score <- rowMeans(rest_data, na.rm = TRUE)
  score[!is.finite(score)] <- NA_real_
  complete <- stats::complete.cases(outcome$y, score, group)
  indices <- original_indices[complete]
  analysis_data <- data.frame(
    y = outcome$y[complete],
    score = score[complete],
    group = droplevels(group[complete])
  )
  group_counts <- table(analysis_data$group, useNA = "no")
  reason <- NULL
  if (nlevels(analysis_data$group) < 2L) {
    reason <- "At least one group was empty after item-specific filtering."
  } else if (any(group_counts < min_group_n)) {
    reason <- sprintf("At least one operational group has fewer than %d usable observations.", min_group_n)
  } else if (length(unique(analysis_data$y)) < 2L) {
    reason <- "The focal item is constant after item-specific filtering."
  } else if (length(unique(analysis_data$score)) < 2L) {
    reason <- "The focal-excluded rest matching score is degenerate."
  }
  if (!is.null(reason)) {
    row <- .pm_empty_dif_row(item, outcome$method, outcome$note, reason)
    row$original_indices[[1L]] <- indices
    row$base_indices[[1L]] <- indices
    row$uniform_indices[[1L]] <- indices
    row$nonuniform_indices[[1L]] <- indices
    row$group_counts[[1L]] <- group_counts
    row$n_used <- nrow(analysis_data)
    return(list(row = row, warnings = list(), errors = list(),
                diagnostics = list(status = "not_evaluable", reason = reason, rest_items = usable_rest,
                                   score = score, analysis_data = analysis_data)))
  }

  fitted <- .pm_fit_logistic_dif(analysis_data)
  if (fitted$status != "evaluable") {
    row <- .pm_empty_dif_row(item, outcome$method, outcome$note, fitted$reason)
    row$original_indices[[1L]] <- indices
    row$base_indices[[1L]] <- indices
    row$uniform_indices[[1L]] <- indices
    row$nonuniform_indices[[1L]] <- indices
    row$group_counts[[1L]] <- group_counts
    row$n_used <- nrow(analysis_data)
    row$base_converged <- isTRUE(fitted$models$base$converged)
    row$uniform_converged <- isTRUE(fitted$models$uniform$converged)
    row$nonuniform_converged <- isTRUE(fitted$models$nonuniform$converged)
    return(list(row = row, warnings = fitted$warnings, errors = fitted$errors,
                diagnostics = list(status = fitted$status, reason = fitted$reason,
                                   rest_items = usable_rest, score = score, analysis_data = analysis_data,
                                   models = fitted$models)))
  }
  r2 <- c(
    base = .pm_nagelkerke_r2(fitted$models$base),
    uniform = .pm_nagelkerke_r2(fitted$models$uniform),
    nonuniform = .pm_nagelkerke_r2(fitted$models$nonuniform)
  )
  delta_u <- r2["uniform"] - r2["base"]
  delta_n <- r2["nonuniform"] - r2["uniform"]
  coefficient <- stats::coef(fitted$models$nonuniform)
  row <- data.frame(
    item = item,
    status = "evaluable",
    reason = NA_character_,
    p_uniform = fitted$p_uniform,
    p_nonuniform = fitted$p_nonuniform,
    lr_uniform = fitted$lr_uniform,
    df_uniform = fitted$df_uniform,
    lr_nonuniform = fitted$lr_nonuniform,
    df_nonuniform = fitted$df_nonuniform,
    nagelkerke_base = unname(r2["base"]),
    nagelkerke_uniform = unname(r2["uniform"]),
    nagelkerke_nonuniform = unname(r2["nonuniform"]),
    delta_r2_uniform = unname(delta_u),
    delta_r2_nonuniform = unname(delta_n),
    effect_size = if (any(is.finite(c(delta_u, delta_n)))) {
      max(c(delta_u, delta_n), na.rm = TRUE)
    } else {
      NA_real_
    },
    flagged = isTRUE(fitted$p_uniform < alpha || fitted$p_nonuniform < alpha),
    method = outcome$method,
    note = outcome$note,
    n_used = nrow(analysis_data),
    base_converged = isTRUE(fitted$models$base$converged),
    uniform_converged = isTRUE(fitted$models$uniform$converged),
    nonuniform_converged = isTRUE(fitted$models$nonuniform$converged),
    stringsAsFactors = FALSE
  )
  row$original_indices <- I(list(indices))
  row$base_indices <- I(list(indices))
  row$uniform_indices <- I(list(indices))
  row$nonuniform_indices <- I(list(indices))
  row$group_counts <- I(list(group_counts))
  row$coefficients <- I(list(coefficient))
  row$odds_ratios <- I(list(exp(coefficient)))
  list(
    row = row,
    warnings = fitted$warnings,
    errors = fitted$errors,
    diagnostics = list(
      status = "evaluable", reason = NA_character_, rest_items = usable_rest,
      score = score, analysis_data = analysis_data, models = fitted$models
    )
  )
}

.pm_dif_interpretation <- function(results,
                                   alpha,
                                   language = "eng",
                                   digits = 3) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    return(data.frame())
  }
  rows <- lapply(seq_len(nrow(results)), function(i) {
    row <- results[i, , drop = FALSE]
    p_uniform <- as.numeric(row$p_uniform)
    p_nonuniform <- as.numeric(row$p_nonuniform)
    uniform_flag <- isTRUE(is.finite(p_uniform) && p_uniform < alpha)
    nonuniform_flag <- isTRUE(is.finite(p_nonuniform) && p_nonuniform < alpha)
    assessable <- identical(as.character(row$status), "evaluable") &&
      is.finite(p_uniform) && is.finite(p_nonuniform)
    dif_type <- .pm_dif_type(uniform_flag, nonuniform_flag, assessable)
    severity <- .pm_dif_severity(dif_type)
    data.frame(
      item = row$item,
      dif_type = dif_type,
      p_uniform = round(p_uniform, digits),
      p_nonuniform = round(p_nonuniform, digits),
      decision = .pm_dif_decision(dif_type, language),
      severity = severity,
      recommendation = .pm_dif_recommendation(dif_type, row$method, row$note, language),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.pm_dif_summary <- function(interpretation,
                            alpha,
                            group_levels,
                            language = "eng",
                            digits = 3) {
  n_items <- nrow(interpretation)
  uniform_flags <- sum(interpretation$dif_type %in% c("uniform", "uniform_and_nonuniform"), na.rm = TRUE)
  nonuniform_flags <- sum(interpretation$dif_type %in% c("nonuniform", "uniform_and_nonuniform"), na.rm = TRUE)
  flagged_items <- sum(interpretation$severity == "review", na.rm = TRUE)
  not_evaluable_items <- sum(interpretation$severity == "not_evaluable", na.rm = TRUE)
  flagged_pct <- if (n_items > 0L) flagged_items / n_items else NA_real_
  data.frame(
    n_items = n_items,
    flagged_items = flagged_items,
    uniform_flags = uniform_flags,
    nonuniform_flags = nonuniform_flags,
    not_evaluable_items = not_evaluable_items,
    flagged_pct = round(flagged_pct, digits),
    alpha = alpha,
    groups = paste(group_levels, collapse = ", "),
    recommendation = .pm_dif_overall_recommendation(
      flagged_items,
      n_items,
      language,
      not_evaluable_items = not_evaluable_items
    ),
    stringsAsFactors = FALSE
  )
}

.pm_dif_type <- function(uniform_flag, nonuniform_flag, assessable) {
  if (!isTRUE(assessable)) return("not_evaluable")
  if (uniform_flag && nonuniform_flag) return("uniform_and_nonuniform")
  if (uniform_flag) return("uniform")
  if (nonuniform_flag) return("nonuniform")
  "none"
}

.pm_dif_severity <- function(dif_type) {
  if (dif_type == "not_evaluable") return("not_evaluable")
  if (dif_type == "none") return("no_flag")
  "review"
}

.pm_dif_decision <- function(dif_type, language) {
  if (language == "esp") {
    return(switch(
      dif_type,
      none = "Sin evidencia estadistica de DIF con este cribado.",
      uniform = "Posible DIF uniforme.",
      nonuniform = "Posible DIF no uniforme.",
      uniform_and_nonuniform = "Posible DIF uniforme y no uniforme.",
      not_evaluable = "No evaluable con este cribado."
    ))
  }
  switch(
    dif_type,
    none = "No statistical DIF evidence in this screening.",
    uniform = "Possible uniform DIF.",
    nonuniform = "Possible nonuniform DIF.",
    uniform_and_nonuniform = "Possible uniform and nonuniform DIF.",
    not_evaluable = "Not assessable with this screening."
  )
}

.pm_dif_recommendation <- function(dif_type, method, note, language) {
  extra_note <- if (!is.na(note) && nzchar(note)) paste0(" ", note) else ""
  if (language == "esp") {
    return(switch(
      dif_type,
      none = paste0("Conservar salvo evidencia teorica o clinica adicional.", extra_note),
      uniform = paste0("Revisar contenido del item y diferencias sistematicas entre grupos.", extra_note),
      nonuniform = paste0("Revisar interaccion item-grupo y considerar analisis DIF ordinal/IRT.", extra_note),
      uniform_and_nonuniform = paste0("Priorizar revision del item y confirmar con un metodo DIF dedicado.", extra_note),
      not_evaluable = paste0("Aumentar variabilidad o tamano muestral antes de interpretar.", extra_note)
    ))
  }
  switch(
    dif_type,
    none = paste0("Retain unless theory or substantive evidence suggests review.", extra_note),
    uniform = paste0("Review item content and systematic group differences.", extra_note),
    nonuniform = paste0("Review item-by-group interaction and consider ordinal or IRT DIF follow-up.", extra_note),
    uniform_and_nonuniform = paste0("Prioritize item review and confirm with a dedicated DIF method.", extra_note),
    not_evaluable = paste0("Increase variation or sample size before interpreting.", extra_note)
  )
}

.pm_dif_overall_recommendation <- function(flagged_items,
                                           n_items,
                                           language,
                                           not_evaluable_items = 0L) {
  flagged_pct <- if (n_items > 0L) flagged_items / n_items else 1
  if (language == "esp") {
    if (not_evaluable_items > 0L) {
      return(sprintf(
        "%d item(s) no fueron evaluables; no interprete la ausencia de alertas como evidencia de invariancia.",
        not_evaluable_items
      ))
    }
    if (flagged_items == 0L) {
      return("No se detectaron items con DIF en este cribado; interpretar junto con teoria y evidencia externa.")
    }
    if (flagged_pct <= 0.25) {
      return("Se detectaron pocos items con posible DIF; revisar esos items antes de decidir eliminacion o ajuste.")
    }
    return("Se detectaron varios items con posible DIF; confirmar con analisis DIF dedicado antes de decisiones finales.")
  }
  if (not_evaluable_items > 0L) {
    return(sprintf(
      "%d item(s) were not evaluable; absence of flags must not be interpreted as equivalent item functioning.",
      not_evaluable_items
    ))
  }
  if (flagged_items == 0L) {
    return("No DIF-flagged items in this screening; interpret together with theory and external evidence.")
  }
  if (flagged_pct <= 0.25) {
    return("Few items show possible DIF; review those items before deletion or adjustment decisions.")
  }
  "Several items show possible DIF; confirm with a dedicated DIF analysis before final decisions."
}

.pm_irt_references <- function() {
  data.frame(
    method = c(
      "Ordered categories and Rasch-family response models",
      "Nominal and polytomous item response models",
      "Marginal maximum likelihood and EM estimation",
      "Generalized partial credit and polytomous IRT",
      "Rasch, 2PL, graded response, and multidimensional IRT implementation"
    ),
    reference = c(
      "Andrich, D. (1978). A rating formulation for ordered response categories. Psychometrika, 43, 561-573. https://doi.org/10.1007/BF02293814",
      "Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or more nominal categories. Psychometrika, 37(1), 29-51. https://doi.org/10.1007/BF02291411",
      "Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm. Psychometrika, 46, 443-459. https://doi.org/10.1007/BF02293801",
      "Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. Applied Psychological Measurement, 16(2), 159-176. https://doi.org/10.1177/014662169201600206",
      "Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment. Journal of Statistical Software, 48(6), 1-29. https://doi.org/10.18637/jss.v048.i06"
    ),
    stringsAsFactors = FALSE
  )
}

.pm_dif_references <- function() {
  data.frame(
    method = c(
      "Contingency-table and logit item-bias models",
      "Logistic regression DIF screening",
      "Logistic regression DIF effect-size evaluation",
      "IRT-based DIF effect-size reporting",
      "Ordinal logistic regression and IRT DIF implementation"
    ),
    reference = c(
      "Mellenbergh, G. J. (1982). Contingency table models for assessing item bias. Journal of Educational Statistics, 7(2), 105-118. https://doi.org/10.3102/10769986007002105",
      "Swaminathan, H., & Rogers, H. J. (1990). Detecting differential item functioning using logistic regression procedures. Journal of Educational Measurement, 27(4), 361-370. https://doi.org/10.1111/j.1745-3984.1990.tb00754.x",
      "Jodoin, M. G., & Gierl, M. J. (2001). Evaluating Type I error and power rates using an effect size measure with the logistic regression procedure for DIF detection. Applied Measurement in Education, 14(4), 329-349. https://doi.org/10.1207/S15324818AME1404_2",
      "Steinberg, L., & Thissen, D. (2006). Using effect sizes for research reporting: Examples using item response theory to analyze differential item functioning. Psychological Methods, 11(4), 402-415. https://doi.org/10.1037/1082-989X.11.4.402",
      "Choi, S. W., Gibbons, L. E., & Crane, P. K. (2011). lordif: An R package for detecting differential item functioning using iterative hybrid ordinal logistic regression/item response theory and Monte Carlo simulations. Journal of Statistical Software, 39(8), 1-30. https://doi.org/10.18637/jss.v039.i08"
    ),
    stringsAsFactors = FALSE
  )
}

.pm_empty_model_fit_row <- function(model_id, model, estimator, reason) {
  data.frame(
    model_id = model_id,
    model = model,
    status = "not_evaluable",
    reason = reason,
    converged = FALSE,
    post_check = FALSE,
    admissible = FALSE,
    comparable = FALSE,
    comparison_reason = reason,
    estimator = estimator,
    estimator_effective = NA_character_,
    estimator_original = NA_character_,
    estimand = NA_character_,
    meanstructure = NA,
    ordered_variables = NA_character_,
    missing = NA_character_,
    cfi = NA_real_,
    tli = NA_real_,
    rmsea = NA_real_,
    srmr = NA_real_,
    cfi_scaled = NA_real_,
    tli_scaled = NA_real_,
    rmsea_scaled = NA_real_,
    cfi_robust = NA_real_,
    tli_robust = NA_real_,
    rmsea_robust = NA_real_,
    comparison_cfi = NA_real_,
    comparison_tli = NA_real_,
    comparison_rmsea = NA_real_,
    comparison_srmr = NA_real_,
    comparison_cfi_variant = NA_character_,
    comparison_tli_variant = NA_character_,
    comparison_rmsea_variant = NA_character_,
    comparison_srmr_variant = NA_character_,
    fit_score = NA_real_,
    aic = NA_real_,
    bic = NA_real_,
    error = reason,
    stringsAsFactors = FALSE
  )
}

.pm_observed_from_parameter_table <- function(parameter_table, data_names) {
  if (is.null(parameter_table) || !is.data.frame(parameter_table)) return(character())
  mentioned <- unique(c(
    parameter_table$lhs[parameter_table$lhs %in% data_names],
    parameter_table$rhs[parameter_table$rhs %in% data_names]
  ))
  data_names[data_names %in% mentioned]
}

.pm_lavaan_option <- function(options, name, default = NA) {
  value <- options[[name]]
  if (is.null(value) || !length(value)) default else value[[1L]]
}

.pm_extract_lavaan_model_record <- function(fit,
                                            model_id,
                                            model_name,
                                            requested_estimator,
                                            syntax) {
  measure_names <- c(
    "cfi", "tli", "rmsea", "srmr", "aic", "bic",
    "cfi.scaled", "tli.scaled", "rmsea.scaled",
    "cfi.robust", "tli.robust", "rmsea.robust"
  )
  measure_capture <- .pm_capture_conditions(lavaan::fitMeasures(fit, measure_names))
  measures <- measure_capture$value
  if (is.null(measures)) {
    measures <- stats::setNames(rep(NA_real_, length(measure_names)), measure_names)
  }
  measures <- stats::setNames(
    vapply(measure_names, function(name) {
      value <- as.numeric(measures[name])
      if (length(value) == 1L && is.finite(value)) value else NA_real_
    }, numeric(1)),
    measure_names
  )
  converged_capture <- .pm_capture_conditions(lavaan::lavInspect(fit, "converged"))
  post_check_capture <- .pm_capture_conditions(lavaan::lavInspect(fit, "post.check"))
  observed_capture <- .pm_capture_conditions(lavaan::lavNames(fit, "ov"))
  ordered_capture <- .pm_capture_conditions(lavaan::lavNames(fit, "ov.ord"))
  case_capture <- .pm_capture_conditions(lavaan::lavInspect(fit, "case.idx"))
  options_capture <- .pm_capture_conditions(lavaan::lavInspect(fit, "options"))
  estimates_capture <- .pm_capture_conditions(lavaan::parameterEstimates(fit))
  inspection_captures <- list(
    converged_capture,
    post_check_capture,
    observed_capture,
    ordered_capture,
    case_capture,
    options_capture,
    estimates_capture
  )
  inspection_errors <- Filter(
    Negate(is.null),
    lapply(inspection_captures, `[[`, "error")
  )
  converged <- isTRUE(converged_capture$value)
  post_check <- isTRUE(post_check_capture$value)
  observed <- if (is.null(observed_capture$value)) character() else observed_capture$value
  ordered <- if (is.null(ordered_capture$value)) character() else ordered_capture$value
  case_idx <- if (is.null(case_capture$value)) integer() else case_capture$value
  case_indices <- sort(unique(as.integer(unlist(case_idx, use.names = FALSE))))
  options <- if (is.null(options_capture$value)) list() else options_capture$value
  estimator_effective <- toupper(as.character(.pm_lavaan_option(options, "estimator", requested_estimator)))
  estimator_original <- toupper(as.character(.pm_lavaan_option(options, "estimator.orig", requested_estimator)))
  missing <- tolower(as.character(.pm_lavaan_option(options, "missing", "listwise")))
  parameterization <- tolower(as.character(.pm_lavaan_option(options, "parameterization", "delta")))
  std_lv <- isTRUE(.pm_lavaan_option(options, "std.lv", FALSE))
  meanstructure <- isTRUE(.pm_lavaan_option(options, "meanstructure", FALSE))
  estimand <- if (length(ordered)) {
    if (meanstructure) "thresholds_polychoric_correlations" else "polychoric_correlations"
  } else if (meanstructure) {
    "means_covariances"
  } else {
    "covariances"
  }
  estimates <- if (is.null(estimates_capture$value)) data.frame() else estimates_capture$value
  estimable <- if (!nrow(estimates)) FALSE else {
    free_rows <- if ("free" %in% names(estimates)) estimates$free > 0 else rep(TRUE, nrow(estimates))
    all(is.finite(estimates$est[free_rows]))
  }
  admissible <- !length(inspection_errors) && converged && post_check && estimable
  reason <- if (length(inspection_errors)) {
    paste("Could not inspect the fitted lavaan model:", conditionMessage(inspection_errors[[1L]]))
  } else if (!converged) {
    "The lavaan optimizer did not converge."
  } else if (!post_check) {
    "lavaan post-check rejected the solution (for example, a Heywood case)."
  } else if (!estimable) {
    "The solution contains nonfinite estimable parameters."
  } else {
    NA_character_
  }
  diagnostics <- list(
    model_id = model_id,
    model = model_name,
    syntax = syntax,
    observed_variables = observed,
    case_indices = case_indices,
    ordered_variables = ordered,
    estimator_effective = estimator_effective,
    estimator_original = estimator_original,
    estimand = estimand,
    meanstructure = meanstructure,
    missing = missing,
    parameterization = parameterization,
    std_lv = std_lv,
    convergence = converged,
    post_check = post_check,
    admissible = admissible,
    reason = reason
  )
  row <- .pm_empty_model_fit_row(model_id, model_name, requested_estimator, reason)
  row$status <- if (admissible) "evaluable" else "not_evaluable"
  row$reason <- reason
  row$converged <- converged
  row$post_check <- post_check
  row$admissible <- admissible
  row$estimator_effective <- estimator_effective
  row$estimator_original <- estimator_original
  row$estimand <- estimand
  row$meanstructure <- meanstructure
  row$ordered_variables <- paste(ordered, collapse = ",")
  row$missing <- missing
  row$cfi <- measures["cfi"]
  row$tli <- measures["tli"]
  row$rmsea <- measures["rmsea"]
  row$srmr <- measures["srmr"]
  row$cfi_scaled <- measures["cfi.scaled"]
  row$tli_scaled <- measures["tli.scaled"]
  row$rmsea_scaled <- measures["rmsea.scaled"]
  row$cfi_robust <- measures["cfi.robust"]
  row$tli_robust <- measures["tli.robust"]
  row$rmsea_robust <- measures["rmsea.robust"]
  row$aic <- .pm_information_criterion(measures["aic"], estimator_original)
  row$bic <- .pm_information_criterion(measures["bic"], estimator_original)
  row$error <- if (admissible) NA_character_ else reason
  warnings <- c(
    measure_capture$warnings,
    unlist(lapply(inspection_captures, `[[`, "warnings"), recursive = FALSE)
  )
  errors <- c(
    Filter(Negate(is.null), list(measure_capture$error)),
    inspection_errors
  )
  list(row = row, diagnostics = diagnostics, warnings = warnings, errors = errors)
}

.pm_summarize_model_effective_settings <- function(diagnostics, model_ids) {
  by_model <- stats::setNames(lapply(model_ids, function(model_id) {
    record <- diagnostics[[model_id]]
    list(
      ordered = if (is.null(record$ordered_variables)) NULL else record$ordered_variables,
      estimator_original = if (is.null(record$estimator_original)) {
        NA_character_
      } else {
        as.character(record$estimator_original)
      },
      estimator_computational = if (is.null(record$estimator_effective)) {
        NA_character_
      } else {
        as.character(record$estimator_effective)
      },
      meanstructure = if (is.null(record$meanstructure)) NA else isTRUE(record$meanstructure),
      estimand = if (is.null(record$estimand)) NA_character_ else as.character(record$estimand)
    )
  }), model_ids)

  available <- by_model[vapply(
    by_model,
    function(x) !is.na(x$estimator_original) && !is.na(x$estimator_computational),
    logical(1)
  )]
  common_value <- function(field, unavailable) {
    if (!length(available)) return(unavailable)
    values <- lapply(available, `[[`, field)
    if (!all(vapply(values[-1L], identical, logical(1), values[[1L]]))) {
      return(unavailable)
    }
    values[[1L]]
  }

  list(
    ordered = common_value("ordered", NULL),
    estimator_original = common_value("estimator_original", NA_character_),
    estimator_computational = common_value("estimator_computational", NA_character_),
    meanstructure = common_value("meanstructure", NA),
    estimand = common_value("estimand", NA_character_),
    by_model = by_model
  )
}

.pm_model_metadata_equal <- function(a, b) {
  set_fields <- c("observed_variables", "ordered_variables")
  same_sets <- all(vapply(
    set_fields,
    function(field) identical(sort(a[[field]]), sort(b[[field]])),
    logical(1)
  ))
  scalar_fields <- c(
    "case_indices", "estimator_effective", "estimator_original", "estimand",
    "meanstructure", "missing", "parameterization", "std_lv"
  )
  same_scalars <- all(vapply(
    scalar_fields,
    function(field) identical(a[[field]], b[[field]]),
    logical(1)
  ))
  same_sets && same_scalars
}

.pm_apply_model_comparability <- function(fit_table, diagnostics) {
  finalize_information_criteria <- function(table, comparable_indices = integer()) {
    all_indices <- seq_len(nrow(table))
    comparable_indices <- intersect(as.integer(comparable_indices), all_indices)
    excluded <- setdiff(all_indices, comparable_indices)

    table$information_criteria_status <- "not_comparable"
    table$information_criteria_reason <- table$comparison_reason
    fallback_reason <- "Information criteria were suppressed because this candidate is not in a methodologically comparable set."
    missing_reason <- is.na(table$information_criteria_reason) |
      !nzchar(table$information_criteria_reason)
    table$information_criteria_reason[missing_reason] <- fallback_reason

    if (length(excluded)) {
      table$aic[excluded] <- NA_real_
      table$bic[excluded] <- NA_real_
    }
    if (length(comparable_indices)) {
      available <- is.finite(table$aic[comparable_indices]) |
        is.finite(table$bic[comparable_indices])
      table$information_criteria_status[comparable_indices] <- ifelse(
        available,
        "comparable",
        "unavailable"
      )
      table$information_criteria_reason[comparable_indices] <- ifelse(
        available,
        NA_character_,
        "Information criteria are unavailable for the effective estimator or fitted model."
      )
    }
    table
  }

  fit_table$comparable <- FALSE
  fit_table$comparison_reason <- ifelse(
    fit_table$admissible,
    "Fewer than two admissible models were available.",
    fit_table$reason
  )
  admissible <- which(fit_table$admissible)
  if (length(admissible) < 2L) {
    fit_table <- finalize_information_criteria(fit_table)
    return(list(fit = fit_table, indices = integer(),
                reason = "At least two admissible models are required for comparison."))
  }
  groups <- list()
  remaining <- admissible
  while (length(remaining)) {
    reference_index <- remaining[1L]
    reference <- diagnostics[[fit_table$model_id[reference_index]]]
    members <- remaining[vapply(
      remaining,
      function(i) {
        .pm_model_metadata_equal(reference, diagnostics[[fit_table$model_id[i]]])
      },
      logical(1)
    )]
    groups[[length(groups) + 1L]] <- members
    remaining <- setdiff(remaining, members)
  }
  group_sizes <- lengths(groups)
  largest <- max(group_sizes)
  mismatch_reason <- paste(
    "Candidate differs from the selected comparison set in observed variables,",
    "actual rows, estimand, ordered treatment, estimator, or effective settings."
  )
  if (largest < 2L) {
    fit_table$comparison_reason[admissible] <- mismatch_reason
    fit_table <- finalize_information_criteria(fit_table)
    return(list(
      fit = fit_table,
      indices = integer(),
      reason = "No pair of admissible candidates shares the same data, estimand, and effective settings."
    ))
  }
  if (sum(group_sizes == largest) != 1L) {
    reason <- "Multiple equally large compatible comparison sets exist; no set was selected by position."
    fit_table$comparison_reason[admissible] <- reason
    fit_table <- finalize_information_criteria(fit_table)
    return(list(fit = fit_table, indices = integer(), reason = reason))
  }
  selected <- groups[[which.max(group_sizes)]]
  excluded <- setdiff(admissible, selected)
  fit_table$comparable[selected] <- TRUE
  fit_table$comparison_reason[selected] <- NA_character_
  if (length(excluded)) fit_table$comparison_reason[excluded] <- mismatch_reason
  fit_table <- finalize_information_criteria(fit_table, selected)
  list(
    fit = fit_table,
    indices = selected,
    reason = if (length(excluded)) {
      "A compatible subset was compared; incompatible admissible candidates were retained with diagnostics."
    } else {
      "Admissible candidates share observed variables, actual rows, estimand, and effective settings."
    }
  )
}

.pm_select_common_fit_variants <- function(fit_table, fit_variant = "scaled") {
  metrics <- c("cfi", "tli", "rmsea", "srmr")
  variants <- stats::setNames(rep(NA_character_, length(metrics)), metrics)
  values <- stats::setNames(vector("list", length(metrics)), metrics)
  for (metric in metrics) {
    candidates <- if (metric == "srmr") {
      c(base = metric)
    } else {
      corrected <- any(fit_table$estimator %in% c("MLR", "MLM", "MLMV", "MLMVS", "WLSMV", "WLSM", "ULSMV", "ULSM")) ||
        any(is.finite(fit_table[[paste0(metric, "_scaled")]]))
      variant <- if (fit_variant == "scaled" && !corrected) "base" else fit_variant
      stats::setNames(if (variant == "base") metric else paste0(metric, "_", variant), variant)
    }
    selected <- NA_character_
    for (variant in names(candidates)) {
      column <- candidates[[variant]]
      if (column %in% names(fit_table) && all(is.finite(fit_table[[column]]))) {
        selected <- variant
        values[[metric]] <- as.numeric(fit_table[[column]])
        break
      }
    }
    variants[metric] <- selected
    if (is.na(selected)) values[[metric]] <- rep(NA_real_, nrow(fit_table))
  }
  list(variants = variants, values = values)
}

.pm_no_model_recommendation <- function(language, message = NULL) {
  if (is.null(message)) {
    message <- if (language == "esp") {
      "No hay al menos dos modelos admisibles y comparables con indices de ajuste finitos."
    } else {
      "At least two admissible, comparable models with finite fit indices are required."
    }
  }
  list(model = NA_character_, model_id = NA_character_, basis = "none", message = message)
}

.pm_model_recommendation <- function(fit_table, language) {
  required <- c("converged", "admissible", "comparable")
  if (!all(required %in% names(fit_table))) return(.pm_no_model_recommendation(language))
  eligible <- fit_table$converged & fit_table$admissible & fit_table$comparable
  candidates <- fit_table[eligible, , drop = FALSE]
  if (nrow(candidates) < 2L) return(.pm_no_model_recommendation(language))
  ranking_fields <- intersect(
    c("fit_score", "comparison_cfi", "comparison_tli", "comparison_rmsea", "comparison_srmr"),
    names(candidates)
  )
  finite_by_row <- if (length(ranking_fields)) {
    rowSums(vapply(candidates[ranking_fields], is.finite, logical(nrow(candidates)))) > 0L
  } else {
    rep(FALSE, nrow(candidates))
  }
  candidates <- candidates[finite_by_row, , drop = FALSE]
  if (nrow(candidates) < 2L) return(.pm_no_model_recommendation(language))
  ranked <- .pm_rank_fit_models(candidates)
  best <- candidates[ranked[1L], , drop = FALSE]
  second <- candidates[ranked[2L], , drop = FALSE]
  available <- ranking_fields[vapply(
    ranking_fields,
    function(field) is.finite(best[[field]]) && is.finite(second[[field]]),
    logical(1)
  )]
  tied <- length(available) > 0L && all(vapply(
    available,
    function(field) isTRUE(all.equal(best[[field]], second[[field]], tolerance = 1e-7)),
    logical(1)
  ))
  if (tied) {
    message <- if (language == "esp") {
      "Los mejores modelos est\u00e1n empatados en todos los criterios disponibles; no se emite una recomendacion posicional."
    } else {
      "The best models tie on every available criterion; no positional recommendation is issued."
    }
    return(.pm_no_model_recommendation(language, message))
  }
  message <- if (language == "esp") {
    sprintf(
      "El modelo `%s` presenta la combinacion descriptivamente mas favorable de CFI/TLI altos y RMSEA/SRMR bajos; interprete esta recomendacion junto con teoria, cargas y residuos.",
      best$model
    )
  } else {
    sprintf(
      "Model `%s` has the most favorable descriptive combination of higher CFI/TLI and lower RMSEA/SRMR; interpret this recommendation together with theory, loadings, and residuals.",
      best$model
    )
  }
  list(model = best$model, model_id = best$model_id, basis = "fit_indices", message = message)
}

.pm_model_observed_variables <- function(models, data) {
  vars <- unique(unlist(lapply(models, function(model) {
    pt <- .pm_capture_conditions(lavaan::lavaanify(model))$value
    if (is.null(pt)) return(character(0))
    unique(pt$rhs[pt$op == "=~"])
  }), use.names = FALSE))
  intersect(vars, names(data))
}

.pm_first_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else x[1L]
}

.pm_information_criterion <- function(x, estimator) {
  if (estimator %in% c("WLSMV", "DWLS", "ULSMV", "ULS")) return(NA_real_)
  value <- as.numeric(x)
  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

.pm_fit_score <- function(cfi, tli, rmsea, srmr) {
  if (!any(is.finite(c(cfi, tli, rmsea, srmr)))) return(NA_real_)
  score <- 0
  if (is.finite(cfi)) score <- score + ifelse(cfi >= .95, 2, ifelse(cfi >= .90, 1, 0))
  if (is.finite(tli)) score <- score + ifelse(tli >= .95, 2, ifelse(tli >= .90, 1, 0))
  if (is.finite(rmsea)) score <- score + ifelse(rmsea <= .06, 2, ifelse(rmsea <= .08, 1, 0))
  if (is.finite(srmr)) score <- score + ifelse(srmr <= .08, 2, ifelse(srmr <= .10, 1, 0))
  score
}

.pm_rank_fit_models <- function(candidates) {
  order(
    -candidates$fit_score,
    -ifelse(is.finite(candidates$comparison_cfi), candidates$comparison_cfi, -Inf),
    -ifelse(is.finite(candidates$comparison_tli), candidates$comparison_tli, -Inf),
    ifelse(is.finite(candidates$comparison_rmsea), candidates$comparison_rmsea, Inf),
    ifelse(is.finite(candidates$comparison_srmr), candidates$comparison_srmr, Inf)
  )
}
