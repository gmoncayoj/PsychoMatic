# Preparation is separate from fitting so the effective native model is auditable.
.pm_esem_prepare <- function(data, factors, n_factors, items, rotation, ordered,
                             estimator, missing, rotation_args, language) {
  language <- .pm_normalize_language(language)
  if (!is.character(rotation) || length(rotation) != 1L || is.na(rotation) ||
      !rotation %in% c("target", "geomin")) {
    stop("`rotation` must be 'target' or 'geomin'.", call. = FALSE)
  }
  is_spec <- inherits(data, "psychomatic_spec")
  if (is_spec && !is.null(data$group)) {
    stop("ESEM currently supports a single sample; grouped specifications are not supported.", call. = FALSE)
  }
  raw <- if (is_spec) data$data else .pm_as_data_frame(data)
  if (anyDuplicated(names(raw)) || anyNA(names(raw))) {
    stop("ESEM requires unique column names.", call. = FALSE)
  }
  if (is_spec && !is.null(items)) stop("`items` must be NULL for a specification.", call. = FALSE)
  if (!is.null(factors)) {
    if (!is.list(factors) || length(factors) < 2L || is.null(names(factors)) ||
        anyNA(names(factors)) || any(!nzchar(names(factors))) || anyDuplicated(names(factors))) {
      stop("`factors` must be a list of at least two unique named factors.", call. = FALSE)
    }
    if (!is.null(items)) stop("`items` must be NULL with `factors`.", call. = FALSE)
    if (!is.null(n_factors)) stop("`n_factors` must be NULL with `factors`.", call. = FALSE)
    if (any(!vapply(factors, function(x) is.character(x) && length(x) > 0L &&
                   !anyNA(x) && all(nzchar(x)), logical(1)))) {
      stop("Each component of `factors` must contain non-empty item names.", call. = FALSE)
    }
    selected <- unname(unlist(factors, use.names = FALSE))
    if (anyDuplicated(selected)) stop("Each item must appear once in `factors`.", call. = FALSE)
    positions <- .pm_resolve_columns(selected, names(raw), "factors")
    k <- length(factors)
    if (is_spec && !setequal(selected, data$items)) {
      stop("`factors` must cover exactly the specification items.", call. = FALSE)
    }
  } else {
    if (rotation == "target") stop("Target rotation requires `factors`.", call. = FALSE)
    k <- as.integer(.pm_validate_scalar_number(n_factors, "n_factors", 2,
                                              .Machine$integer.max, integer = TRUE))
    if (is_spec) items <- data$items
    if (is.null(items)) stop("Geomin requires explicit `items` for raw data.", call. = FALSE)
    positions <- .pm_resolve_columns(items, names(raw), "items")
    selected <- names(raw)[positions]
  }
  if (k >= length(selected)) stop("`n_factors` must be smaller than the item count.", call. = FALSE)
  if (is_spec) {
    spec <- data
  } else {
    spec <- psychomatic_spec(raw, items = selected, language = language)
  }
  x <- spec$data[, selected, drop = FALSE]
  if (is.null(ordered)) {
    ordinal_names <- selected[spec$measurement_level[selected] == "ordinal"]
  } else if (is.logical(ordered) && length(ordered) == 1L && !is.na(ordered)) {
    ordinal_names <- if (ordered) selected else character()
  } else if (is.character(ordered) && !anyNA(ordered) && !anyDuplicated(ordered) &&
             all(nzchar(ordered)) && all(ordered %in% selected)) {
    ordinal_names <- ordered
  } else stop("`ordered` must be TRUE, FALSE or unique selected item names.", call. = FALSE)
  category_order <- stats::setNames(lapply(selected, function(nm) {
    if (!nm %in% ordinal_names) return(NULL)
    if (is.ordered(x[[nm]])) levels(x[[nm]]) else sort(unique(x[[nm]][!is.na(x[[nm]])]))
  }), selected)
  decisions <- rbind(spec$decision_log, .pm_decision("ordered",
    if (is.null(ordered)) "auto" else paste(ordered, collapse = ", "),
    paste(ordinal_names, collapse = ", "),
    "Ordinal selection used by the native engine; explicit selection overrides detected levels."))
  for (nm in selected) {
    v <- x[[nm]]
    if (is.ordered(v) && !nm %in% ordinal_names) {
      stop("Continuous treatment of ordered factors requires explicit numeric coding.", call. = FALSE)
    }
    if (is.ordered(v)) v <- as.integer(v)
    if (!is.numeric(v) || inherits(v, c("Date", "POSIXt")) || any(!is.finite(v[!is.na(v)]))) {
      stop(sprintf("`%s` requires numeric or ordered responses.", nm), call. = FALSE)
    }
    if (length(unique(v[!is.na(v)])) < 2L) stop("Items need at least two observed values.", call. = FALSE)
    x[[nm]] <- v
  }
  effective_ranges <- list()
  for (nm in intersect(spec$reverse, selected)) {
    limits <- spec$ranges[[nm]]
    if (is.null(limits)) limits <- range(x[[nm]], na.rm = TRUE)
    effective_ranges[[nm]] <- limits
    x[[nm]] <- sum(limits) - x[[nm]]
    decisions <- rbind(decisions, .pm_decision("reverse", nm,
      paste0(nm, ": ", sum(limits), " - score"),
      if (is.null(spec$ranges[[nm]])) "Bounds inferred from observed encoded responses." else
        "Declared bounds applied to encoded responses."))
  }
  if (is.null(estimator)) estimator <- if (length(ordinal_names)) "WLSMV" else "MLR"
  if (!is.character(estimator) || length(estimator) != 1L || is.na(estimator)) {
    stop("Invalid `estimator`.", call. = FALSE)
  }
  estimator <- toupper(estimator)
  allowed <- if (length(ordinal_names)) c("WLSMV", "ULSMV") else c("ML", "MLR")
  if (!estimator %in% allowed) stop("`estimator` is incompatible with the effective item levels.", call. = FALSE)
  allowed_missing <- if (length(ordinal_names)) c("listwise", "pairwise") else c("listwise", "fiml")
  if (!is.character(missing) || length(missing) != 1L || is.na(missing) || !missing %in% allowed_missing) {
    stop("`missing` is incompatible with the estimator and item levels.", call. = FALSE)
  }
  decisions <- rbind(decisions, .pm_decision("missing", spec$missing, missing,
    "ESEM missing-data policy is explicit; no imputation or policy substitution."))
  preparation <- list(source = if (is_spec) "psychomatic_spec" else "raw data",
    measurement_level = spec$measurement_level[selected],
    ordinal_levels = spec$ordinal_levels, category_order = category_order,
    reverse = spec$reverse, ranges_declared = spec$ranges,
    ranges_effective = effective_ranges, missing_spec = spec$missing)
  preparation_table <- data.frame(item = selected,
    level = unname(spec$measurement_level[selected]),
    treatment = ifelse(selected %in% ordinal_names, "ordinal", "continuous"),
    reversed = selected %in% spec$reverse,
    bounds_used = vapply(selected, function(nm) paste(effective_ranges[[nm]], collapse = ", "), character(1)),
    coding = vapply(category_order, function(v) {
      if (!length(v)) "" else paste(paste(seq_along(v), v, sep = " = "), collapse = "; ")
    }, character(1)), row.names = NULL, stringsAsFactors = FALSE)
  if (!is.list(rotation_args) || (length(rotation_args) &&
      (is.null(names(rotation_args)) || anyNA(names(rotation_args)) ||
       any(!nzchar(names(rotation_args))) || anyDuplicated(names(rotation_args)) ||
       any(!names(rotation_args) %in% c("rstarts", "geomin_epsilon", "gpa_tol", "max_iter"))))) {
    stop("Invalid or unsupported `rotation_args`.", call. = FALSE)
  }
  if (rotation == "target" && "geomin_epsilon" %in% names(rotation_args)) {
    stop("`geomin_epsilon` applies only to geomin.", call. = FALSE)
  }
  rot <- list(orthogonal = FALSE, algorithm = "gpa", std_ov = TRUE,
              row_weights = "none", reflect = TRUE, order_lv_by = "index",
              warn = TRUE, verbose = FALSE, rstarts = 30L, gpa_tol = 1e-5, max_iter = 10000L)
  if (rotation == "geomin") rot$geomin_epsilon <- .001
  rot[names(rotation_args)] <- rotation_args
  for (nm in intersect(c("rstarts", "max_iter"), names(rot))) {
    rot[[nm]] <- as.integer(.pm_validate_scalar_number(rot[[nm]], nm, 1,
                                      .Machine$integer.max, integer = TRUE))
  }
  for (nm in intersect(c("geomin_epsilon", "gpa_tol"), names(rot))) {
    rot[[nm]] <- .pm_validate_scalar_number(rot[[nm]], nm, .Machine$double.eps)
  }
  item_map <- data.frame(item = selected, alias = sprintf("pm_i%03d", seq_along(selected)),
                         stringsAsFactors = FALSE)
  factor_map <- data.frame(factor = paste0("F", seq_len(k)),
    label = if (rotation == "target") names(factors) else paste0("F", seq_len(k)),
    reference_label = if (!is.null(factors)) names(factors) else rep(NA_character_, k),
    stringsAsFactors = FALSE)
  names(x) <- item_map$alias
  target <- NULL
  if (rotation == "target") {
    target <- matrix(0, nrow = ncol(x), ncol = k,
                     dimnames = list(item_map$alias, factor_map$factor))
    for (i in seq_len(k)) target[match(factors[[i]], selected), i] <- 1
    rot$target <- target
  }
  syntax <- paste(paste0('efa("measurement")*', factor_map$factor, collapse = " + "),
                  "=~", paste(item_map$alias, collapse = " + "))
  cfa_syntax <- if (!is.null(factors)) paste(vapply(seq_len(k), function(i) {
    paste(factor_map$factor[i], "=~", paste(item_map$alias[match(factors[[i]], selected)], collapse = " + "))
  }, character(1)), collapse = "\n") else NULL
  native_args <- list(model = syntax, data = x, estimator = estimator,
                       std_lv = TRUE, missing = missing,
                       rotation = c(list(rotation), rot))
  if (length(ordinal_names)) native_args$ordered <- item_map$alias[match(ordinal_names, selected)]
  list(analysis_data = x, syntax = syntax, target = target, cfa_syntax = cfa_syntax,
       item_map = item_map, factor_map = factor_map, factors = factors,
       case_indices = seq_len(nrow(x)), native_args = native_args,
       spec = spec, preparation_table = preparation_table,
       settings = list(rotation = rotation, n_factors = k,
       estimator = estimator, ordered = ordinal_names, missing = missing,
       std_lv = TRUE, rotation_args = rot, language = language, preparation = preparation),
       decisions = decisions)
}

.pm_esem_fit <- function(prepared, seed) {
  captured <- .efa_with_seed(seed, .pm_capture_conditions(
    do.call(lavaan::sem, prepared$native_args)))
  list(fit = captured$value, warnings = captured$warnings, error = captured$error,
       native_args = prepared$native_args)
}

