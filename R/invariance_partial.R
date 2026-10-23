.ipi_namespace_available <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

.ipi_scalar_choice <- function(x, name, choices) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !nzchar(trimws(x))) {
    stop(sprintf("`%s` must be one non-empty character value.", name),
         call. = FALSE)
  }
  value <- tolower(trimws(x))
  if (!value %in% choices) {
    stop(
      sprintf("`%s` must be one of: %s.", name, paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }
  value
}

.ipi_parse_expression <- function(expression) {
  if (!is.character(expression) || length(expression) != 1L ||
      is.na(expression) || !nzchar(trimws(expression))) {
    stop("Every `free` expression must be a non-empty character scalar without NA.",
         call. = FALSE)
  }
  expression <- trimws(expression)
  operand_ok <- function(x) {
    nzchar(x) && !grepl("[=~|]", x)
  }
  capture <- function(pattern) {
    match <- regexec(pattern, expression, perl = TRUE)
    value <- regmatches(expression, match)[[1L]]
    if (!length(value)) NULL else value[-1L]
  }

  pieces <- capture("^(.+?)\\s*=~\\s*(.+?)$")
  if (!is.null(pieces)) {
    lhs <- trimws(pieces[1L])
    rhs <- trimws(pieces[2L])
    if (operand_ok(lhs) && operand_ok(rhs)) {
      return(list(
        family = "loadings", lhs = lhs, op = "=~", rhs = rhs,
        canonical = paste(lhs, "=~", rhs)
      ))
    }
  }

  pieces <- capture("^(.+?)\\s*~~\\s*(.+?)$")
  if (!is.null(pieces)) {
    lhs <- trimws(pieces[1L])
    rhs <- trimws(pieces[2L])
    if (!operand_ok(lhs) || !operand_ok(rhs)) {
      stop("Residual `free` expressions contain an unsupported operand.", call. = FALSE)
    }
    if (!identical(lhs, rhs)) {
      stop("Residual `free` expressions must use identical item names on both sides.",
           call. = FALSE)
    }
    return(list(
      family = "residuals", lhs = lhs, op = "~~", rhs = rhs,
      canonical = paste(lhs, "~~", rhs)
    ))
  }

  pieces <- capture("^(.+?)\\s*\\|\\s*t([1-9][0-9]*)$")
  if (!is.null(pieces)) {
    lhs <- trimws(pieces[1L])
    threshold <- paste0("t", pieces[2L])
    if (operand_ok(lhs)) {
      return(list(
        family = "thresholds", lhs = lhs, op = "|", rhs = threshold,
        canonical = paste(lhs, "|", threshold)
      ))
    }
  }

  pieces <- capture("^(.+?)\\s*~\\s*1$")
  if (!is.null(pieces)) {
    lhs <- trimws(pieces[1L])
    if (operand_ok(lhs)) {
      return(list(
        family = "intercepts", lhs = lhs, op = "~1", rhs = "",
        canonical = paste(lhs, "~ 1")
      ))
    }
  }

  stop(
    paste0(
      "Unsupported `free` grammar: '", expression,
      "'. Use factor =~ item, item ~ 1, item | t#, or item ~~ item."
    ),
    call. = FALSE
  )
}

.ipi_parse_free <- function(free) {
  supplied_family <- NULL
  if (is.list(free) && !is.data.frame(free)) {
    family_names <- names(free)
    allowed <- c("loadings", "intercepts", "thresholds", "residuals")
    if (is.null(family_names) || !length(family_names) ||
        anyNA(family_names) || any(!nzchar(family_names)) ||
        anyDuplicated(family_names) || !all(family_names %in% allowed)) {
      stop(
        "Named-list `free` family names must be a unique subset of loadings, intercepts, thresholds, and residuals.",
        call. = FALSE
      )
    }
    lengths_ok <- vapply(
      free,
      function(x) is.character(x) && length(x) > 0L && !anyNA(x) &&
        all(nzchar(trimws(x))),
      logical(1)
    )
    if (!all(lengths_ok)) {
      stop("Every named-list `free` element must contain non-empty expressions without NA.",
           call. = FALSE)
    }
    supplied_family <- rep(family_names, lengths(free))
    expressions <- unlist(free, use.names = FALSE)
  } else {
    if (!is.character(free) || !length(free) || anyNA(free) ||
        any(!nzchar(trimws(free)))) {
      stop("`free` must contain one or more non-empty character expressions without NA.",
           call. = FALSE)
    }
    expressions <- unname(free)
  }

  parsed <- lapply(expressions, .ipi_parse_expression)
  out <- data.frame(
    input_order = seq_along(parsed),
    family = vapply(parsed, `[[`, character(1), "family"),
    lhs = vapply(parsed, `[[`, character(1), "lhs"),
    op = vapply(parsed, `[[`, character(1), "op"),
    rhs = vapply(parsed, `[[`, character(1), "rhs"),
    canonical = vapply(parsed, `[[`, character(1), "canonical"),
    stringsAsFactors = FALSE
  )
  if (!is.null(supplied_family) && any(out$family != supplied_family)) {
    bad <- which(out$family != supplied_family)[1L]
    stop(
      sprintf(
        "`%s` belongs to family '%s', not named-list family '%s'.",
        out$canonical[bad], out$family[bad], supplied_family[bad]
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(out$canonical)) {
    stop("`free` expressions must remain unique after canonicalization; duplicate release found.",
         call. = FALSE)
  }
  out
}

.ipi_normalize_rationale <- function(rationale, parsed) {
  if (!is.data.frame(parsed) || !all(c("canonical", "input_order") %in% names(parsed))) {
    stop("Parsed `free` expressions are required to normalize `rationale`.", call. = FALSE)
  }
  if (!is.character(rationale) || !length(rationale) || anyNA(rationale) ||
      any(!nzchar(trimws(rationale)))) {
    stop("`rationale` is required and must contain non-empty text without NA.",
         call. = FALSE)
  }
  if (length(rationale) == 1L) {
    text <- rep(unname(rationale), nrow(parsed))
  } else {
    rationale_names <- names(rationale)
    if (is.null(rationale_names) || anyNA(rationale_names) ||
        any(!nzchar(trimws(rationale_names))) || anyDuplicated(rationale_names)) {
      stop("A multi-value `rationale` must have unique expression names.", call. = FALSE)
    }
    canonical_names <- vapply(
      rationale_names,
      function(x) .ipi_parse_expression(x)$canonical,
      character(1)
    )
    if (anyDuplicated(canonical_names) ||
        !setequal(canonical_names, parsed$canonical) ||
        length(canonical_names) != nrow(parsed)) {
      stop("Named `rationale` expressions must match the canonical `free` set exactly.",
           call. = FALSE)
    }
    text <- unname(rationale[match(parsed$canonical, canonical_names)])
  }
  data.frame(
    input_order = parsed$input_order,
    canonical = parsed$canonical,
    rationale = text,
    stringsAsFactors = FALSE
  )
}

.ipi_resolve_stage <- function(stage, measurement, category_profile, workflow = "ordinal") {
  stage_requested <- .ipi_scalar_choice(
    stage, "stage", c("metric", "loadings", "thresholds", "scalar", "strict")
  )
  measurement <- .ipi_scalar_choice(
    measurement, "measurement", c("continuous", "ordinal", "mixed")
  )
  if (workflow == "factorial" || measurement == "continuous") {
    allowed <- c("metric", "scalar", "strict")
    if (!stage_requested %in% allowed) {
      stop("For workflow = 'factorial' or continuous data, `stage` must be metric, scalar, or strict. Use workflow = 'ordinal' for identified threshold stages.",
           call. = FALSE)
    }
    stage_key <- stage_requested
    stage_effective <- stage_requested
  } else if (identical(category_profile, "binary")) {
    if (!stage_requested %in% c("scalar", "strict")) {
      stop(
        "Binary indicators permit only the jointly identified scalar or strict stage; separate thresholds/loadings are not supported.",
        call. = FALSE
      )
    }
    stage_key <- stage_requested
    stage_effective <- stage_requested
  } else {
    allowed <- c("thresholds", "loadings", "metric", "scalar", "strict")
    if (!stage_requested %in% allowed) {
      stop("Unsupported ordinal/mixed partial-invariance stage.", call. = FALSE)
    }
    stage_key <- if (stage_requested %in% c("loadings", "metric")) "metric" else stage_requested
    stage_effective <- if (stage_key == "metric") "loadings" else stage_key
  }
  sequence <- .ifa_stage_sequence(if (workflow == "factorial") "continuous" else measurement, category_profile)
  position <- match(stage_key, names(sequence))
  if (is.na(position) || position <= 1L) {
    stop("The requested partial-invariance stage has no supported prerequisite.",
         call. = FALSE)
  }
  prerequisite_key <- names(sequence)[position - 1L]
  list(
    stage_requested = stage_requested,
    stage_effective = stage_effective,
    stage_key = stage_key,
    prerequisite_key = prerequisite_key,
    constraints = sequence[[stage_key]]$constraints,
    prerequisite_constraints = sequence[[prerequisite_key]]$constraints,
    stage_detail = sequence[[stage_key]]$detail,
    prerequisite_detail = sequence[[prerequisite_key]]$detail,
    sequence = sequence
  )
}

.ipi_model_membership <- function(model) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package `lavaan` is required for `invariance_partial()`.", call. = FALSE)
  }
  captured <- .pm_capture_conditions(lavaan::lavaanify(model))
  if (!is.null(captured$error) || is.null(captured$value)) {
    stop(
      paste0("The measurement model could not be parsed: ",
             if (is.null(captured$error)) "unknown syntax error" else conditionMessage(captured$error)),
      call. = FALSE
    )
  }
  table <- captured$value
  loading_rows <- table[table$op == "=~", c("lhs", "rhs"), drop = FALSE]
  if (!nrow(loading_rows)) {
    stop("The model must contain at least one first-order measurement loading.",
         call. = FALSE)
  }
  latent <- unique(as.character(loading_rows$lhs))
  if (any(as.character(loading_rows$rhs) %in% latent)) {
    stop("Higher-order measurement models are not supported by `invariance_partial()` 1.0.0.",
         call. = FALSE)
  }
  membership <- unique(data.frame(
    factor = as.character(loading_rows$lhs),
    item = as.character(loading_rows$rhs),
    stringsAsFactors = FALSE
  ))
  if (!nrow(membership) || any(!nzchar(membership$factor)) ||
      any(!nzchar(membership$item))) {
    stop("First-order factor/item membership could not be established.", call. = FALSE)
  }
  list(membership = membership, parameter_table = table, conditions = captured$warnings)
}

.ipi_validate_structural_releases <- function(parsed, membership) {
  if (!is.data.frame(parsed) ||
      !all(c("family", "lhs", "rhs", "canonical") %in% names(parsed))) {
    stop("Parsed `free` expressions are required for structural validation.",
         call. = FALSE)
  }
  if (!is.data.frame(membership) ||
      !all(c("factor", "item") %in% names(membership)) || !nrow(membership)) {
    stop("First-order factor/item membership is required for structural validation.",
         call. = FALSE)
  }

  present <- vapply(seq_len(nrow(parsed)), function(i) {
    if (identical(parsed$family[i], "loadings")) {
      return(any(
        membership$factor == parsed$lhs[i] & membership$item == parsed$rhs[i]
      ))
    }
    parsed$lhs[i] %in% membership$item
  }, logical(1L))

  if (any(!present)) {
    offending <- parsed$canonical[which(!present)[1L]]
    stop(
      sprintf(
        "Release `%s` does not identify a parameter in the measurement model.",
        offending
      ),
      call. = FALSE
    )
  }
  invisible(parsed)
}

.ipi_anchor_table <- function(membership, parsed, active_constraints) {
  if (!is.data.frame(membership) ||
      !all(c("factor", "item") %in% names(membership)) || !nrow(membership)) {
    stop("A fitted first-order factor/item membership table is required.", call. = FALSE)
  }
  active_constraints <- .ifa_normalize_constraints(active_constraints)
  rows <- lapply(seq_len(nrow(membership)), function(i) {
    factor <- as.character(membership$factor[i])
    item <- as.character(membership$item[i])
    loading_release <- parsed$family == "loadings" &
      parsed$lhs == factor & parsed$rhs == item
    item_release <- parsed$family %in% c("intercepts", "thresholds", "residuals") &
      parsed$lhs == item & parsed$family %in% active_constraints
    affected <- parsed$canonical[loading_release | item_release]
    data.frame(
      factor = factor,
      item = item,
      retained = length(affected) == 0L,
      released_by = paste(affected, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.ipi_validate_anchors <- function(anchors, minimum) {
  counts <- aggregate(retained ~ factor, anchors, function(x) sum(x, na.rm = TRUE))
  names(counts)[2L] <- "retained_count"
  counts$required_count <- as.integer(minimum)
  failing <- counts$retained_count < minimum
  if (any(failing)) {
    details <- vapply(which(failing), function(i) {
      factor <- counts$factor[i]
      affected <- anchors$released_by[anchors$factor == factor & !anchors$retained]
      affected <- unique(affected[nzchar(affected)])
      sprintf(
        "%s retained %d anchor item(s), requires %d; affected releases: %s",
        factor, counts$retained_count[i], minimum,
        if (length(affected)) paste(affected, collapse = ", ") else "none"
      )
    }, character(1))
    stop(paste("Insufficient retained anchors:", paste(details, collapse = " | ")),
         call. = FALSE)
  }
  list(items = anchors, counts = counts)
}

.ipi_prepare_analysis <- function(data,
                                  group_expression,
                                  environment,
                                  model,
                                  ordered,
                                  estimator,
                                  criteria,
                                  measurement_assertion,
                                  workflow = "factorial",
                                  fit_variant = "scaled",
                                  missing = "auto",
                                  min_group_n = 20L) {
  input <- .ifa_extract_input(data)
  resolved_group <- .ifa_resolve_group(
    input$data, group_expression, environment
  )
  indicators <- .ifa_extract_indicators(model)
  item_preparation <- .ifa_prepare_workflow_items(
    resolved_group$data,
    indicators,
    ordered = ordered,
    specification_levels = input$measurement_levels,
    estimator = estimator,
    workflow = workflow
  )
  items <- item_preparation$items
  measurement_requested <- .ipi_scalar_choice(
    measurement_assertion,
    "measurement",
    c("auto", "continuous", "ordinal", "mixed")
  )
  if (measurement_requested != "auto" &&
      !identical(measurement_requested, items$measurement)) {
    stop(
      sprintf(
        paste0(
          "`measurement = \"%s\"` contradicts the prepared indicators, which are '%s'. ",
          "Set `ordered` or a `psychomatic_spec` measurement level explicitly."
        ),
        measurement_requested, items$measurement
      ),
      call. = FALSE
    )
  }

  numeric_wls <- workflow == "factorial" && is.character(estimator) &&
    length(estimator) == 1L && !is.na(estimator) &&
    toupper(trimws(estimator)) %in% c("WLSMV", "ULSMV")
  missing_request <- .ifa_resolve_missing(missing,
    if (numeric_wls) "ordinal" else items$measurement)
  filtered <- .ifa_filter_analysis_rows(
    items$data,
    indicators,
    resolved_group$group_name,
    resolved_group$original_indices,
    missing_request
  )
  items <- .ifa_refresh_analytic_items(items, filtered$data, indicators)
  analysis_data <- items$data
  sample_validation <- .ifa_validate_analysis_data(
    analysis_data,
    indicators,
    resolved_group$group_name,
    min_group_n = min_group_n
  )

  group_sizes <- sample_validation$group_sizes
  ratio <- as.numeric(max(group_sizes) / min(group_sizes))
  normality <- if (items$measurement == "continuous") {
    .ifa_group_normality(
      analysis_data, indicators, resolved_group$group_name
    )
  } else {
    list(
      overall_normal = NA,
      overall_test = list(
        normal = NA, skew_p = NA_real_, kurt_p = NA_real_,
        status = "not_applicable_to_categorical_measurement",
        n_complete = NA_integer_
      ),
      by_group = data.frame(
        group = levels(analysis_data[[resolved_group$group_name]]),
        n_complete = as.integer(group_sizes),
        skew_p = NA_real_, kurtosis_p = NA_real_,
        status = "not_applicable", message = "",
        stringsAsFactors = FALSE
      ),
      conditions = list()
    )
  }
  criteria_effective <- .ifa_resolve_criteria(
    sample_validation$total_n,
    ratio,
    normality$overall_normal,
    items$measurement,
    criteria
  )
  estimator_info <- .ifa_resolve_estimator(
    estimator,
    items$measurement,
    normality$overall_normal,
    ordered_argument = ordered,
    allow_numeric_wls = workflow == "factorial"
  )
  parameterization_request <- if (items$measurement == "continuous") {
    NA_character_
  } else {
    if (workflow == "factorial") "delta" else "theta"
  }
  id_cat <- if (items$measurement == "continuous" || workflow == "factorial") {
    NA_character_
  } else {
    "Wu.Estabrook.2016"
  }
  estimand <- if (items$measurement == "continuous") {
    "means_covariances"
  } else {
    "thresholds_polychoric_correlations"
  }
  blueprint <- list(
    case_indices = filtered$original_indices,
    variables = indicators,
    group_name = resolved_group$group_name,
    group_source = resolved_group$group_source,
    group_coding = levels(analysis_data[[resolved_group$group_name]]),
    estimator = estimator_info$estimator,
    estimator_engine_request = estimator_info$estimator,
    estimator_effective = NA_character_,
    estimator_orig = NA_character_,
    estimand = estimand,
    ordered = items$ordered,
    ordered_engine_request = items$ordered,
    missing_engine_request = missing_request,
    missing = missing_request,
    identification_engine_request = if (workflow == "factorial") "marker" else "std.lv",
    identification = if (workflow == "factorial") "marker" else "std.lv",
    parameterization_engine_request = parameterization_request,
    parameterization = if (workflow == "factorial") "delta" else "theta",
    ID.fac = if (workflow == "factorial") "marker" else "std.lv",
    ID.cat = id_cat,
    measurement = items$measurement,
    category_profile = items$category_profile,
    model = model,
    effective_syntax = list()
  )
  list(
    workflow = workflow,
    fit_variant = fit_variant,
    source_measurement = item_preparation$source_measurement,
    score_coding = item_preparation$score_coding,
    input = input,
    resolved_group = resolved_group,
    indicators = indicators,
    items = items,
    analysis_data = analysis_data,
    sample_validation = sample_validation,
    group_sizes = group_sizes,
    normality = normality,
    criteria = criteria_effective,
    estimator_info = estimator_info,
    missing_request = missing_request,
    blueprint = blueprint,
    dropped_missing_items = filtered$dropped_missing_items
  )
}

.ipi_fit_record <- function(prepared,
                            stage_key,
                            stage_detail,
                            constraints,
                            group_partial = character()) {
  blueprint <- prepared$blueprint
  built <- .ifa_build_syntax(
    model = blueprint$model,
    data = prepared$analysis_data,
    group_name = blueprint$group_name,
    measurement = blueprint$measurement,
    ordered = blueprint$ordered,
    estimator = prepared$estimator_info$estimator,
    constraints = constraints,
    group_partial = group_partial,
    missing = prepared$missing_request,
    parameterization = blueprint$parameterization,
    ID.fac = blueprint$ID.fac,
    ID.cat = if (blueprint$measurement == "continuous") {
      "Wu.Estabrook.2016"
    } else {
      blueprint$ID.cat
    },
    native = prepared$workflow == "factorial"
  )
  record <- .ifa_extract_fit_record(
    built,
    stage = stage_key,
    stage_detail = stage_detail,
    blueprint = blueprint,
    indicators = prepared$indicators,
    fit_variant = prepared$fit_variant
  )
  record$verified_releases <- character()
  record
}

.ipi_sync_blueprint <- function(prepared, record) {
  actual <- .ifa_effective_lavaan_settings(record$fit)
  prepared$blueprint$missing <- actual$missing
  prepared$blueprint$parameterization <- actual$parameterization
  prepared$blueprint$estimator_effective <- actual$estimator
  prepared$blueprint$estimator_orig <- actual$estimator_orig
  if (is.list(record$fit_metadata)) {
    prepared$blueprint$ordered <- record$fit_metadata$ordered
    prepared$blueprint$identification <- record$fit_metadata$identification
    prepared$blueprint$std.lv <- record$fit_metadata$std_lv
    prepared$blueprint$estimand <- record$fit_metadata$estimand
    prepared$blueprint$group_coding <- record$fit_metadata$group_coding
  }
  prepared
}

.ipi_sync_record <- function(record) {
  metadata <- record$fit_metadata
  if (!is.list(metadata)) return(record)
  record$case_indices <- metadata$case_indices
  record$variables <- metadata$variables
  record$estimator_effective <- metadata$estimator_effective
  record$estimator_orig <- metadata$estimator_orig
  record$estimand <- metadata$estimand
  record$ordered <- metadata$ordered
  record$missing <- metadata$missing
  record$identification <- metadata$identification
  record$parameterization <- metadata$parameterization
  record$group_coding <- metadata$group_coding
  record$measurement <- metadata$measurement
  record
}

.ipi_membership_from_record <- function(record, indicators) {
  table <- record$parameter_table
  if (!is.data.frame(table) || !nrow(table)) {
    stop("The fitted configural parameter table is unavailable; anchors cannot be established.",
         call. = FALSE)
  }
  rows <- table$op == "=~" & table$rhs %in% indicators
  if ("group" %in% names(table)) rows <- rows & table$group == min(table$group)
  membership <- unique(data.frame(
    factor = as.character(table$lhs[rows]),
    item = as.character(table$rhs[rows]),
    stringsAsFactors = FALSE
  ))
  if (!nrow(membership)) {
    stop("No first-order factor/item memberships were found in the fitted model.",
         call. = FALSE)
  }
  membership
}

.ipi_parameter_rows <- function(record, expression, group_coding) {
  table <- record$parameter_table
  if (!is.data.frame(table) || !nrow(table)) return(data.frame())
  rhs <- if (expression$op == "~1") "" else expression$rhs
  selected <- table[
    table$lhs == expression$lhs & table$op == expression$op & table$rhs == rhs,
    , drop = FALSE
  ]
  if (!nrow(selected)) return(selected)
  if (!"group" %in% names(selected)) selected$group <- 1L
  selected$group_label <- group_coding[as.integer(selected$group)]
  selected
}

.ipi_globally_tied <- function(rows, group_coding) {
  if (!is.data.frame(rows) || nrow(rows) != length(group_coding) ||
      !setequal(as.character(rows$group_label), as.character(group_coding))) {
    return(FALSE)
  }
  labels <- as.character(rows$label)
  if (anyNA(labels)) return(FALSE)
  all(nzchar(labels)) && length(unique(labels)) == 1L
}

.ipi_validate_complete_families <- function(parsed,
                                            record,
                                            group_coding,
                                            active_constraints) {
  inactive <- !parsed$family %in% active_constraints
  if (any(inactive)) {
    stop(
      sprintf(
        "Release `%s` belongs to an inactive family at the requested stage.",
        parsed$canonical[which(inactive)[1L]]
      ),
      call. = FALSE
    )
  }
  matched_groups <- character(nrow(parsed))
  equality <- logical(nrow(parsed))
  for (i in seq_len(nrow(parsed))) {
    expression <- as.list(parsed[i, , drop = FALSE])
    expression <- lapply(expression, function(x) x[[1L]])
    rows <- .ipi_parameter_rows(record, expression, group_coding)
    counts <- if (nrow(rows)) table(factor(rows$group_label, levels = group_coding)) else integer()
    if (length(counts) != length(group_coding) || any(counts != 1L)) {
      stop(
        sprintf(
          "Release `%s` must match exactly one pertinent parameter row in every analysis group.",
          parsed$canonical[i]
        ),
        call. = FALSE
      )
    }
    equality[i] <- .ipi_globally_tied(rows, group_coding)
    if (!equality[i]) {
      stop(
        sprintf(
          "Release `%s` is not effectively equality-constrained in the complete stage.",
          parsed$canonical[i]
        ),
        call. = FALSE
      )
    }
    matched_groups[i] <- paste(group_coding, collapse = ", ")
  }
  data.frame(
    canonical = parsed$canonical,
    matched_groups = matched_groups,
    equality_in_full = equality,
    stringsAsFactors = FALSE
  )
}

.ipi_verify_releases <- function(full_record,
                                 partial_record,
                                 parsed,
                                 group_coding) {
  if (!nrow(parsed)) {
    return(data.frame(
      canonical = character(), effective_release = logical(),
      release_reason = character(), stringsAsFactors = FALSE
    ))
  }
  effective <- logical(nrow(parsed))
  reasons <- character(nrow(parsed))
  for (i in seq_len(nrow(parsed))) {
    expression <- as.list(parsed[i, , drop = FALSE])
    expression <- lapply(expression, function(x) x[[1L]])
    full_rows <- .ipi_parameter_rows(full_record, expression, group_coding)
    partial_rows <- .ipi_parameter_rows(partial_record, expression, group_coding)
    exact <- nrow(full_rows) == length(group_coding) &&
      nrow(partial_rows) == length(group_coding) &&
      setequal(full_rows$group_label, group_coding) &&
      setequal(partial_rows$group_label, group_coding)
    if (!exact) {
      reasons[i] <- "release_rows_missing_or_duplicated"
      next
    }
    full_rows <- full_rows[match(group_coding, full_rows$group_label), , drop = FALSE]
    partial_rows <- partial_rows[match(group_coding, partial_rows$group_label), , drop = FALSE]
    full_labels <- as.character(full_rows$label)
    partial_labels <- as.character(partial_rows$label)
    full_free <- as.integer(full_rows$free)
    partial_free <- as.integer(partial_rows$free)
    if (anyNA(c(full_labels, partial_labels)) ||
        anyNA(c(full_free, partial_free))) {
      reasons[i] <- "release_metadata_incomplete"
      next
    }
    changed <- any(full_labels != partial_labels) || any(full_free != partial_free)
    still_tied <- .ipi_globally_tied(partial_rows, group_coding)
    pair_independent <- TRUE
    if (length(group_coding) > 1L) {
      pairs <- utils::combn(seq_along(group_coding), 2L)
      pair_independent <- all(apply(pairs, 2L, function(pair) {
        labels <- partial_labels[pair]
        frees <- partial_free[pair]
        !(all(nzchar(labels)) && identical(labels[1L], labels[2L])) &&
          !(all(!nzchar(labels)) && identical(frees[1L], frees[2L]) && frees[1L] > 0L)
      }))
    }
    at_least_one_estimable <- any(partial_free > 0L)
    effective[i] <- isTRUE(changed) && !isTRUE(still_tied) &&
      isTRUE(pair_independent) && isTRUE(at_least_one_estimable)
    reasons[i] <- if (effective[i]) {
      "requested_equality_family_released"
    } else if (!changed) {
      "group_partial_was_ignored"
    } else if (still_tied) {
      "family_remains_globally_tied"
    } else {
      "independent_group_identification_not_established"
    }
  }
  data.frame(
    canonical = parsed$canonical,
    effective_release = effective,
    release_reason = reasons,
    stringsAsFactors = FALSE
  )
}

.ipi_partial_nesting_gate <- function(prerequisite_releases,
                                      partial_releases,
                                      requested,
                                      prerequisite_constraints) {
  prerequisite_releases <- sort(unique(as.character(prerequisite_releases)))
  partial_releases <- sort(unique(as.character(partial_releases)))
  expected_prerequisite <- sort(requested$canonical[
    requested$family %in% prerequisite_constraints
  ])
  expected_partial <- sort(requested$canonical)
  comparable <- identical(prerequisite_releases, expected_prerequisite) &&
    identical(partial_releases, expected_partial)
  list(
    comparable = comparable,
    reason = if (comparable) "group_partial_releases_nested" else "mismatch_group_partial_nesting",
    expected_prerequisite = expected_prerequisite,
    actual_prerequisite = prerequisite_releases,
    expected_partial = expected_partial,
    actual_partial = partial_releases
  )
}

.ipi_metadata_gate <- function(first, second) {
  if (!isTRUE(first$admissible)) {
    return(list(comparable = FALSE, reason = "first_model_not_admissible"))
  }
  if (!isTRUE(second$admissible)) {
    return(list(comparable = FALSE, reason = "second_model_not_admissible"))
  }
  first_metadata <- .ifa_comparison_metadata(first)
  second_metadata <- .ifa_comparison_metadata(second)
  fields <- c(
    "case_indices", "variables", "estimator_effective", "estimator_orig",
    "estimand", "ordered", "missing", "identification", "parameterization",
    "std_lv", "group_coding", "measurement"
  )
  for (field in fields) {
    if (!identical(first_metadata[[field]], second_metadata[[field]])) {
      return(list(comparable = FALSE, reason = paste0("mismatch_", field)))
    }
  }
  list(comparable = TRUE, reason = "same_frozen_settings")
}

.ipi_compare_same_stage <- function(partial, full, requested_releases) {
  metadata <- .ipi_metadata_gate(partial, full)
  if (!metadata$comparable) return(metadata)
  if (!identical(
    .ifa_normalize_constraints(partial$effective_constraints),
    .ifa_normalize_constraints(full$effective_constraints)
  )) {
    return(list(comparable = FALSE, reason = "mismatch_same_stage_constraints"))
  }
  actual_partial <- sort(unique(as.character(partial$verified_releases)))
  expected <- sort(unique(as.character(requested_releases)))
  if (!identical(actual_partial, expected) || length(full$verified_releases)) {
    return(list(comparable = FALSE, reason = "mismatch_same_stage_releases"))
  }
  normalize_group_partial <- function(value) {
    if (is.null(value) || !length(value) || identical(value, "")) {
      return(character())
    }
    if (!is.character(value) || anyNA(value) || any(!nzchar(trimws(value)))) {
      return(NA_character_)
    }
    sort(unique(trimws(value)))
  }
  partial_argument <- normalize_group_partial(partial$group_partial)
  full_argument <- normalize_group_partial(full$group_partial)
  if (anyNA(partial_argument) || anyNA(full_argument) ||
      !identical(partial_argument, expected) || length(full_argument)) {
    return(list(
      comparable = FALSE,
      reason = "mismatch_same_stage_group_partial"
    ))
  }
  indices <- c("cfi", "rmsea", "srmr", "df")
  common <- stats::setNames(lapply(indices, function(index) {
    .ifa_common_variant(partial$measures[[index]], full$measures[[index]],
      intersect(.ifa_record_priority(partial, index), .ifa_record_priority(full, index)))
  }), indices)
  missing <- names(common)[vapply(common, is.null, logical(1))]
  if (length(missing)) {
    return(list(
      comparable = FALSE,
      reason = paste0("no_common_variant_", paste(missing, collapse = "_"))
    ))
  }
  delta_df <- common$df$current - common$df$previous
  if (!is.finite(delta_df) || delta_df <= 0) {
    return(list(
      comparable = FALSE, reason = "full_minus_partial_df_not_positive",
      delta_df = unname(delta_df)
    ))
  }
  list(
    comparable = TRUE,
    reason = "same_stage_partial_to_full_diagnostic",
    delta_df = unname(delta_df),
    improvement_cfi = unname(common$cfi$previous - common$cfi$current),
    improvement_rmsea = unname(common$rmsea$current - common$rmsea$previous),
    improvement_srmr = unname(common$srmr$current - common$srmr$previous),
    variants = vapply(common, `[[`, character(1), "variant")
  )
}

.ipi_configural_gate_record <- function(record, criteria) {
  if (!isTRUE(record$admissible)) {
    reason <- record$admissibility_reason
    if (!is.character(reason) || length(reason) != 1L ||
        is.na(reason) || !nzchar(reason)) {
      reason <- "configural_model_not_admissible"
    }
    return(list(
      passed = FALSE,
      decision = "not_evaluable",
      status = "not_evaluable",
      reason = reason
    ))
  }
  selected <- vapply(
    c("cfi", "rmsea", "srmr"),
    function(index) .ifa_select_variant(record$measures[[index]], .ifa_record_priority(record, index))$value,
    numeric(1)
  )
  .ifa_configural_gate(selected, criteria)
}

.ipi_decide_change <- function(comparison, criteria, stage_key) {
  if (!isTRUE(comparison$comparable)) {
    return(list(
      decision = "not_evaluable", status = "not_evaluable",
      reason = comparison$reason
    ))
  }
  .ifa_change_decision(
    comparison$delta_cfi,
    comparison$delta_rmsea,
    comparison$delta_srmr,
    criteria,
    stage_key
  )
}

.ipi_resolve_conclusion_status <- function(releases_effective,
                                           nesting_comparable,
                                           configural_gate,
                                           full_decision,
                                           partial_comparison,
                                           partial_decision) {
  if (!isTRUE(releases_effective) || !isTRUE(nesting_comparable)) {
    return("not_evaluable")
  }
  if (!isTRUE(configural_gate$passed)) {
    if (identical(configural_gate$decision, "not_supported")) {
      return("configural_not_supported")
    }
    return("not_evaluable")
  }
  if (identical(full_decision$decision, "supported")) {
    return("full_invariance_supported_partial_not_needed")
  }
  if (!isTRUE(partial_comparison$comparable)) {
    return("not_evaluable")
  }
  if (identical(partial_decision$decision, "supported")) {
    return("partial_invariance_supported")
  }
  if (identical(partial_decision$decision, "not_supported")) {
    return("partial_invariance_not_supported")
  }
  "not_evaluable"
}

.ipi_selected_fit_values <- function(record) {
  indices <- c("chi2", "df", "p", "cfi", "rmsea", "srmr")
  selected <- lapply(indices, function(index) {
    .ifa_select_variant(record$measures[[index]], .ifa_record_priority(record, index))
  })
  names(selected) <- indices
  list(
    values = vapply(selected, `[[`, numeric(1), "value"),
    variants = vapply(selected, `[[`, character(1), "variant")
  )
}

.ipi_fit_table <- function(records) {
  rows <- lapply(names(records), function(role) {
    record <- records[[role]]
    selected <- .ipi_selected_fit_values(record)
    data.frame(
      role = role,
      stage = record$stage,
      stage_detail = record$stage_detail,
      constraints = paste(record$effective_constraints, collapse = "+"),
      group_partial = paste(
        if (identical(record$group_partial, "")) character() else record$group_partial,
        collapse = "; "
      ),
      chi2 = selected$values[["chi2"]],
      df = selected$values[["df"]],
      p = selected$values[["p"]],
      cfi = selected$values[["cfi"]],
      rmsea = selected$values[["rmsea"]],
      srmr = selected$values[["srmr"]],
      chi2_variant = selected$variants[["chi2"]],
      df_variant = selected$variants[["df"]],
      p_variant = selected$variants[["p"]],
      cfi_variant = selected$variants[["cfi"]],
      rmsea_variant = selected$variants[["rmsea"]],
      srmr_variant = selected$variants[["srmr"]],
      converged = record$converged,
      post_check = record$post_check,
      admissible = record$admissible,
      reason = record$admissibility_reason,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.ipi_model_envelope <- function(record) {
  fit <- record$fit
  nobs <- if (!is.null(fit) && inherits(fit, "lavaan")) {
    value <- tryCatch(lavaan::lavInspect(fit, "nobs"), error = function(e) NA)
    sum(as.numeric(unlist(value, use.names = FALSE)), na.rm = TRUE)
  } else {
    NA_real_
  }
  metadata <- record$fit_metadata
  list(
    fit = fit,
    fit_record = record,
    nobs = nobs,
    case_idx = if (is.list(metadata)) metadata$case_indices else integer(),
    variables = if (is.list(metadata)) metadata$variables else character(),
    estimator = record$estimator,
    estimator_effective = record$estimator_effective,
    estimator_original = record$estimator_orig,
    ordered = record$ordered,
    missing = record$missing,
    estimand = record$estimand,
    identification = record$identification,
    parameterization = record$parameterization,
    group_coding = record$group_coding
  )
}

.ipi_comparison_row <- function(name, comparison, decision = NULL) {
  get <- function(field, default = NA) {
    value <- comparison[[field]]
    if (is.null(value) || !length(value)) default else value[[1L]]
  }
  variants <- comparison$variants
  data.frame(
    comparison = name,
    comparable = isTRUE(comparison$comparable),
    reason = as.character(get("reason", "not_evaluable")),
    delta_df = as.numeric(get("delta_df", NA_real_)),
    delta_cfi = as.numeric(get("delta_cfi", NA_real_)),
    delta_rmsea = as.numeric(get("delta_rmsea", NA_real_)),
    delta_srmr = as.numeric(get("delta_srmr", NA_real_)),
    improvement_cfi = as.numeric(get("improvement_cfi", NA_real_)),
    improvement_rmsea = as.numeric(get("improvement_rmsea", NA_real_)),
    improvement_srmr = as.numeric(get("improvement_srmr", NA_real_)),
    cfi_variant = if (is.null(variants)) NA_character_ else as.character(variants[["cfi"]]),
    rmsea_variant = if (is.null(variants)) NA_character_ else as.character(variants[["rmsea"]]),
    srmr_variant = if (is.null(variants)) NA_character_ else as.character(variants[["srmr"]]),
    decision = if (is.null(decision)) NA_character_ else as.character(decision$decision),
    decision_reason = if (is.null(decision)) NA_character_ else as.character(decision$reason),
    stringsAsFactors = FALSE
  )
}

.ipi_conclusion_text <- function(status, language) {
  if (language == "esp") {
    return(switch(
      status,
      full_invariance_supported_partial_not_needed = paste(
        "La invariancia completa en la etapa solicitada fue compatible con las heur\u00edsticas de cambio;",
        "la liberaci\u00f3n parcial no era necesaria. Esto no constituye prueba formal de invariancia."
      ),
      partial_invariance_supported = paste(
        "La etapa de invariancia parcial fue compatible con las heur\u00edsticas de cambio usando",
        "las liberaciones justificadas. Esto no prueba invariancia y requiere juicio sustantivo."
      ),
      partial_invariance_not_supported = paste(
        "La etapa de invariancia parcial no fue compatible con las heur\u00edsticas de cambio;",
        "no deben a\u00f1adirse liberaciones autom\u00e1ticas sin justificaci\u00f3n independiente."
      ),
      configural_not_supported = paste(
        "El ajuste configural no satisfizo los criterios descriptivos;",
        "las comparaciones de invariancia parcial no deben interpretarse."
      ),
      "La evidencia necesaria no fue evaluable; revise ajustes, admisibilidad, restricciones y diagn\u00f3sticos."
    ))
  }
  switch(
    status,
    full_invariance_supported_partial_not_needed = paste(
      "Full invariance at the requested stage was compatible with the fit-change heuristics;",
      "the partial release was not needed. This is not formal proof of invariance."
    ),
    partial_invariance_supported = paste(
      "The partial-invariance stage was compatible with the fit-change heuristics using",
      "the justified releases. This is not proof of invariance and requires substantive judgment."
    ),
    partial_invariance_not_supported = paste(
      "The partial-invariance stage was not compatible with the fit-change heuristics;",
      "do not add automated releases without independent justification."
    ),
    configural_not_supported = paste(
      "Configural fit did not meet the descriptive gates; partial-invariance comparisons",
      "must not be interpreted."
    ),
    "Required evidence was not evaluable; inspect fits, admissibility, constraints, and diagnostics."
  )
}

.ipi_configural_fit_failed <- function(record) {
  !is.null(record$error) || is.null(record$fit) || !isTRUE(record$admissible)
}

.ipi_not_attempted_model <- function(prepared, role, reason) {
  blueprint <- prepared$blueprint
  list(
    fit = NULL,
    fit_record = list(
      stage = role,
      status = "not_attempted",
      reason = reason,
      fit = NULL,
      error = NULL
    ),
    nobs = NA_real_,
    case_idx = blueprint$case_indices,
    variables = blueprint$variables,
    estimator = NA_character_,
    estimator_effective = NA_character_,
    estimator_original = NA_character_,
    ordered = blueprint$ordered,
    missing = blueprint$missing,
    estimand = blueprint$estimand,
    identification = blueprint$identification,
    parameterization = blueprint$parameterization,
    group_coding = blueprint$group_coding,
    status = "not_attempted",
    reason = reason
  )
}

.ipi_configural_failure_result <- function(call,
                                           requested_settings,
                                           prepared,
                                           parsed_model_conditions,
                                           parsed_free,
                                           rationale_table,
                                           stage_info,
                                           configural,
                                           report_path,
                                           language_effective) {
  reason <- configural$admissibility_reason
  if (!is.character(reason) || length(reason) != 1L ||
      is.na(reason) || !nzchar(reason)) {
    reason <- "configural_fit_not_evaluable"
  }
  original_error <- configural$error
  message <- if (inherits(original_error, "condition")) {
    conditionMessage(original_error)
  } else {
    paste0("The configural fit was not evaluable: ", reason, ".")
  }
  failure <- if (inherits(original_error, "condition")) {
    original_error
  } else {
    .ifa_condition(message, "configural_fit_not_evaluable", type = "error")
  }
  unavailable <- list(
    comparable = FALSE,
    reason = "configural_fit_not_evaluable"
  )
  unavailable_decision <- list(
    decision = "not_evaluable",
    status = "not_evaluable",
    reason = "configural_fit_not_evaluable"
  )
  configural_gate <- .ipi_configural_gate_record(configural, prepared$criteria)

  roles <- c(
    "prerequisite_full", "prerequisite_matched", "full", "partial"
  )
  models <- list(configural = .ipi_model_envelope(configural))
  for (role in roles) {
    models[[role]] <- .ipi_not_attempted_model(prepared, role, reason)
  }
  models$prerequisite <- models$prerequisite_matched

  freed_constraints <- parsed_free
  freed_constraints$rationale <- rationale_table$rationale[
    match(freed_constraints$canonical, rationale_table$canonical)
  ]
  freed_constraints$matched_groups <- NA_character_
  freed_constraints$equality_in_full <- NA
  freed_constraints$effective_release <- FALSE
  freed_constraints$release_reason <- "not_attempted_configural_failure"
  anchors <- data.frame(
    factor = character(), item = character(), retained = logical(),
    released_by = character(), stringsAsFactors = FALSE
  )
  anchor_counts <- data.frame(
    factor = character(), retained_count = integer(),
    required_count = integer(), stringsAsFactors = FALSE
  )
  fit_table <- .ipi_fit_table(list(configural = configural))
  comparison_table <- do.call(rbind, list(
    .ipi_comparison_row(
      "full_prerequisite_vs_full", unavailable, unavailable_decision
    ),
    .ipi_comparison_row(
      "prerequisite_vs_partial", unavailable, unavailable_decision
    ),
    .ipi_comparison_row(
      "full_vs_partial_diagnostic", unavailable, NULL
    )
  ))
  rownames(comparison_table) <- NULL
  comparisons <- list(
    configural_gate = configural_gate,
    full_prerequisite_vs_full = unavailable,
    prerequisite_vs_partial = unavailable,
    full_vs_partial_diagnostic = unavailable,
    group_partial_nesting = unavailable
  )
  conclusion_status <- "not_evaluable"
  conclusion <- .ipi_conclusion_text(conclusion_status, language_effective)
  fit_conditions <- c(configural$warnings, configural$messages)

  result <- .pm_result(
    class = "psychomatic_invariance_partial",
    call = call,
    requested_settings = requested_settings,
    effective_settings = list(
      language = language_effective,
      workflow = prepared$workflow,
      fit_variant = prepared$fit_variant,
      source_measurement = prepared$source_measurement,
      stage = stage_info$stage_effective,
      stage_key = stage_info$stage_key,
      measurement = prepared$items$measurement,
      category_profile = prepared$items$category_profile,
      estimator_engine_request = prepared$estimator_info$estimator,
      estimator = configural$estimator_effective,
      estimator_original = configural$estimator_orig,
      ordered_engine_request = prepared$blueprint$ordered_engine_request,
      ordered = configural$ordered,
      ordered_effective = configural$ordered,
      missing_engine_request = prepared$missing_request,
      missing = configural$missing,
      parameterization_engine_request = prepared$blueprint$parameterization_engine_request,
      parameterization = configural$parameterization,
      identification = configural$identification,
      std_lv = configural$std_lv,
      ID.fac = configural$ID.fac_effective,
      ID.cat = configural$ID.cat_effective,
      ID.fac_engine_request = configural$ID.fac_engine_request,
      ID.fac_effective = configural$ID.fac_effective,
      ID.cat_engine_request = configural$ID.cat_engine_request,
      ID.cat_effective = configural$ID.cat_effective,
      min_anchor_items = as.integer(requested_settings$min_anchor_items)
    ),
    data_summary = list(
      n = prepared$sample_validation$total_n,
      group_sizes = prepared$group_sizes,
      indicators = prepared$indicators,
      category_counts = prepared$items$category_counts,
      dropped_missing_group = prepared$resolved_group$dropped_missing_group,
      dropped_missing_items = prepared$dropped_missing_items
    ),
    results = list(
      fit_table = fit_table,
      freed_constraints = freed_constraints,
      anchors = list(items = anchors, counts = anchor_counts),
      comparison = comparisons
    ),
    warnings = c(
      parsed_model_conditions, prepared$normality$conditions, fit_conditions
    ),
    errors = list(failure),
    diagnostics = list(
      configural_failure = list(
        code = "configural_fit_not_evaluable",
        message = message,
        admissibility_reason = reason,
        converged = isTRUE(configural$converged),
        post_check = isTRUE(configural$post_check),
        finite_parameters = isTRUE(configural$finite_parameters),
        error = original_error
      ),
      fit = .ifa_diagnostics_long(list(configural = configural)),
      release_verification = data.frame(),
      lower_release_verification = data.frame(),
      comparison_table = comparison_table
    ),
    decision_log = do.call(rbind, list(
      .pm_decision(
        "partial_measurement", requested_settings$measurement,
        prepared$items$measurement,
        "Shared complete-invariance preparation determined the effective measurement type."
      ),
      .pm_decision(
        "partial_stage", requested_settings$stage, stage_info$stage_effective,
        paste("Effective stage:", stage_info$stage_key)
      ),
      .pm_decision(
        "partial_conclusion", "researcher_directed_release",
        conclusion_status,
        "The configural fit failed its admissibility gate; later fits were not attempted."
      )
    )),
    conclusion = list(status = conclusion_status, reason = conclusion)
  )

  references <- .ifa_reference_table()
  settings_for_report <- result$effective_settings
  settings_table <- data.frame(
    setting = names(settings_for_report),
    value = vapply(settings_for_report, function(value) {
      if (is.null(value) || !length(value)) return("")
      if (all(is.na(value))) return(NA_character_)
      paste(as.character(value), collapse = "; ")
    }, character(1L)),
    stringsAsFactors = FALSE
  )
  tables <- list(
    fit = fit_table,
    freed_constraints = freed_constraints,
    rationale = rationale_table,
    anchors = anchors,
    anchor_counts = anchor_counts,
    comparisons = comparison_table,
    criteria = data.frame(
      criterion = names(prepared$criteria),
      value = unlist(prepared$criteria, use.names = FALSE),
      stringsAsFactors = FALSE
    ),
    warnings = .ifa_conditions_for_report(result$warnings),
    errors = .ifa_conditions_for_report(result$errors),
    conclusion = data.frame(
      status = conclusion_status,
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ),
    settings = settings_table,
    references = references
  )
  prepared$blueprint$effective_syntax <- list(
    configural = configural$syntax,
    prerequisite_full = NULL,
    prerequisite_matched = NULL,
    full = NULL,
    partial = NULL
  )
  prepared$blueprint$group_partial <- list(
    configural = character(),
    prerequisite_full = character(),
    prerequisite_matched = character(),
    full = character(),
    partial = character()
  )
  direct <- list(
    measurement = prepared$items$measurement,
    stage = stage_info$stage_key,
    stage_requested = stage_info$stage_requested,
    stage_effective = stage_info$stage_effective,
    prerequisite_stage = stage_info$prerequisite_key,
    analysis_blueprint = prepared$blueprint,
    analysis_data = prepared$analysis_data,
    models = models,
    fitted_models = lapply(
      models[setdiff(names(models), "prerequisite")], `[[`, "fit"
    ),
    freed_constraints = freed_constraints,
    rationale = requested_settings$rationale,
    rationale_table = rationale_table,
    anchors = anchors,
    anchor_counts = anchor_counts,
    comparison = comparisons,
    fit_table = fit_table,
    criteria = prepared$criteria,
    tables = tables,
    references = references,
    conclusion_status = conclusion_status,
    conclusion_text = conclusion,
    report_path = report_path
  )
  for (field in names(direct)) result[[field]] <- direct[[field]]

  if (!is.null(report_path)) {
    written <- .pm_write_excel_report(
      result$tables,
      report_path,
      overwrite = FALSE
    )
    result$report_path <- as.character(written)
  }
  result
}

#' Fit researcher-directed partial measurement invariance
#'
#' Fits a complete target invariance stage and a matched partial stage using
#' explicit, researcher-justified parameter releases. The function automates
#' the lengthy lavaan/semTools syntax while deliberately avoiding automated
#' modification-index searches or iterative freeing.
#'
#' @inheritParams efa_auto
#' @param data A data frame, matrix, or `psychomatic_spec`.
#' @param group Group column name/symbol/index or a vector with one value per row.
#' @param model A first-order lavaan measurement model.
#' @param free Canonical parameter expressions to release, supplied as a
#'   character vector or a family-named list.
#' @param rationale Required substantive rationale, either one global value or
#'   an exactly expression-named character vector.
#' @param stage Target stage: `"metric"`, `"scalar"`, or `"strict"` in the
#'   default factorial workflow. The explicit ordinal workflow additionally
#'   accepts `"thresholds"` and `"loadings"`; binary ordinal data require
#'   joint `"scalar"` or `"strict"` stages.
#' @param measurement Assertion (`"auto"`, `"continuous"`, `"ordinal"`, or
#'   `"mixed"`) checked against prepared item types.
#' @param estimator Estimator passed through the complete-invariance engine.
#' @param ordered Ordered-item assertion accepted by the complete-invariance
#'   engine.
#' @param criteria Fit-change criteria accepted by
#'   `factorial_invariance_auto()`.
#' @param min_anchor_items Minimum number of fully constrained anchor items per
#'   first-order factor; must be an integer of at least two.
#' @param report `FALSE`, `TRUE`, or an Excel output path. Existing files are
#'   protected unless `overwrite = TRUE` explicitly permits replacement.
#' @param language Presentation language, `"eng"` by default; `"esp"` for Spanish.
#'   Explicit `NULL` inherits a specification or the package option.
#' @param fit_variant `"scaled"` (default) or `"robust"` for CFI and RMSEA,
#'   changes and decisions. Unavailable robust values are not substituted.
#'   Corrected chi-square, df and p remain scaled; SRMR remains base.
#' @param workflow `"factorial"` (default) follows the native lavaan settings and
#'   four-stage constraint sequence of `factorial_invariance_auto()`, with
#'   researcher-specified `group.partial` releases. WLSMV/ULSMV default to
#'   ordered indicators, marker identification and delta parameterization.
#'   `"ordinal"` retains the identified semTools theta/std.lv/Wu-Estabrook
#'   procedure, including threshold releases. An ineffective release or a
#'   comparison without additional degrees of freedom cannot support invariance.
#' @param missing Missing-data option accepted by `factorial_invariance_auto()`.
#' @param min_group_n Minimum effective group size; default 20.
#'
#' @return Invisibly, an object of classes
#'   `psychomatic_invariance_partial` and `psychomatic_result`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' model <- "F1 =~ R1 + R2 + R3 + R4"
#' # Continuous indicators; justify the freed intercept before testing.
#' invariance_partial(MY_DATA, group = "Country", model = model,
#'                    free = "R2 ~ 1", estimator = "ML",
#'                    rationale = "A preregistered item-wording hypothesis")
#' }
#'
#' @export
invariance_partial <- function(
  data, group, model, free, rationale,
  stage = "scalar", measurement = "auto", estimator = "auto",
  ordered = NULL, criteria = "auto", min_anchor_items = 2L,
  report = FALSE, language = "eng",
  digits = 3, file_name = NULL, format = "excel", overwrite = FALSE, verbose = TRUE,
  fit_variant = c("scaled", "robust"), workflow = c("factorial", "ordinal"),
  missing = "auto", min_group_n = 20L
) {
  call <- match.call()
  fit_variant <- match.arg(fit_variant, c("scaled", "robust"))
  workflow <- match.arg(workflow, c("factorial", "ordinal"))
  requested_language <- language
  requested_stage <- stage
  requested_measurement <- measurement
  requested_estimator <- estimator
  requested_ordered <- ordered
  requested_criteria <- criteria
  requested_free <- free
  requested_rationale <- rationale

  if (workflow == "ordinal" && !.ipi_namespace_available("semTools")) {
    stop(
      "Package `semTools` is required by `invariance_partial(workflow = 'ordinal')`. Install it with install.packages('semTools').",
      call. = FALSE
    )
  }
  if (!.ipi_namespace_available("lavaan")) {
    stop("Package `lavaan` is required by `invariance_partial()`.", call. = FALSE)
  }
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  report_path <- NULL
  min_anchor_items <- .pm_validate_scalar_number(
    min_anchor_items,
    "min_anchor_items",
    lower = 2,
    integer = TRUE
  )
  parsed_model <- .ipi_model_membership(model)
  parsed_free <- .ipi_parse_free(free)
  rationale_table <- .ipi_normalize_rationale(rationale, parsed_free)
  .ipi_validate_structural_releases(parsed_free, parsed_model$membership)
  prepared <- .ipi_prepare_analysis(
    data = data,
    group_expression = substitute(group),
    environment = parent.frame(),
    model = model,
    ordered = ordered,
    estimator = estimator,
    criteria = criteria,
    measurement_assertion = measurement,
    workflow = workflow,
    fit_variant = fit_variant,
    missing = missing,
    min_group_n = min_group_n
  )
  language_effective <- .pm_normalize_language(
    language, prepared$input$language
  )
  stage_info <- .ipi_resolve_stage(
    stage,
    prepared$items$measurement,
    prepared$items$category_profile,
    workflow = workflow
  )
  inactive <- !parsed_free$family %in% stage_info$constraints
  if (any(inactive)) {
    stop(
      sprintf(
        "Release `%s` belongs to a family inactive at stage '%s'.",
        parsed_free$canonical[which(inactive)[1L]], stage_info$stage_effective
      ),
      call. = FALSE
    )
  }
  requested_settings <- list(
    workflow = workflow,
    fit_variant = fit_variant,
    missing = missing,
    min_group_n = min_group_n,
    group = paste(deparse(substitute(group)), collapse = " "),
    free = requested_free,
    rationale = requested_rationale,
    stage = requested_stage,
    measurement = requested_measurement,
    estimator = requested_estimator,
    ordered = requested_ordered,
    criteria = requested_criteria,
    min_anchor_items = min_anchor_items,
    report = report,
    language = requested_language
  )

  configural <- .ipi_fit_record(
    prepared,
    "configural",
    stage_info$sequence$configural$detail,
    character(),
    character()
  )
  prepared <- .ipi_sync_blueprint(prepared, configural)
  configural <- .ipi_sync_record(configural)
  if (.ipi_configural_fit_failed(configural)) {
    failed <- .ipi_configural_failure_result(
      call = call,
      requested_settings = requested_settings,
      prepared = prepared,
      parsed_model_conditions = parsed_model$conditions,
      parsed_free = parsed_free,
      rationale_table = rationale_table,
      stage_info = stage_info,
      configural = configural,
      report_path = report_path,
      language_effective = language_effective
    )
    return(.pm_finish_analysis(failed, output, language_effective))
  }

  if (stage_info$prerequisite_key == "configural") {
    prerequisite_full <- configural
  } else {
    prerequisite_full <- .ipi_sync_record(.ipi_fit_record(
      prepared,
      stage_info$prerequisite_key,
      stage_info$prerequisite_detail,
      stage_info$prerequisite_constraints,
      character()
    ))
  }
  lower <- parsed_free[parsed_free$family %in% stage_info$prerequisite_constraints, , drop = FALSE]
  if (nrow(lower)) {
    .ipi_validate_complete_families(
      lower, prerequisite_full, prepared$blueprint$group_coding,
      stage_info$prerequisite_constraints
    )
    prerequisite_matched <- .ipi_sync_record(.ipi_fit_record(
      prepared,
      stage_info$prerequisite_key,
      stage_info$prerequisite_detail,
      stage_info$prerequisite_constraints,
      lower$canonical
    ))
    lower_verification <- .ipi_verify_releases(
      prerequisite_full, prerequisite_matched, lower,
      prepared$blueprint$group_coding
    )
    prerequisite_matched$verified_releases <- lower_verification$canonical[
      lower_verification$effective_release
    ]
  } else {
    prerequisite_matched <- prerequisite_full
    lower_verification <- data.frame(
      canonical = character(), effective_release = logical(),
      release_reason = character(), stringsAsFactors = FALSE
    )
  }

  full <- .ipi_sync_record(.ipi_fit_record(
    prepared,
    stage_info$stage_key,
    stage_info$stage_detail,
    stage_info$constraints,
    character()
  ))
  complete_match <- .ipi_validate_complete_families(
    parsed_free,
    full,
    prepared$blueprint$group_coding,
    stage_info$constraints
  )

  membership <- .ipi_membership_from_record(configural, prepared$indicators)
  anchors_raw <- .ipi_anchor_table(
    membership, parsed_free, stage_info$constraints
  )
  anchors <- .ipi_validate_anchors(anchors_raw, as.integer(min_anchor_items))

  partial <- .ipi_sync_record(.ipi_fit_record(
    prepared,
    stage_info$stage_key,
    stage_info$stage_detail,
    stage_info$constraints,
    parsed_free$canonical
  ))
  release_verification <- .ipi_verify_releases(
    full, partial, parsed_free, prepared$blueprint$group_coding
  )
  partial$verified_releases <- release_verification$canonical[
    release_verification$effective_release
  ]

  freed_constraints <- Reduce(
    function(x, y) merge(x, y, by = "canonical", all.x = TRUE, sort = FALSE),
    list(
      parsed_free,
      rationale_table[c("canonical", "rationale")],
      complete_match,
      release_verification
    )
  )
  freed_constraints <- freed_constraints[
    match(parsed_free$canonical, freed_constraints$canonical), , drop = FALSE
  ]
  freed_constraints$input_order <- parsed_free$input_order

  nesting <- .ipi_partial_nesting_gate(
    prerequisite_releases = prerequisite_matched$verified_releases,
    partial_releases = partial$verified_releases,
    requested = parsed_free,
    prerequisite_constraints = stage_info$prerequisite_constraints
  )
  configural_gate <- .ipi_configural_gate_record(configural, prepared$criteria)
  full_comparison <- .ifa_compare_records(prerequisite_full, full)
  full_decision <- .ipi_decide_change(
    full_comparison, prepared$criteria, stage_info$stage_key
  )
  partial_comparison <- if (isTRUE(nesting$comparable)) {
    .ifa_compare_records(prerequisite_matched, partial)
  } else {
    list(comparable = FALSE, reason = nesting$reason)
  }
  partial_decision <- .ipi_decide_change(
    partial_comparison, prepared$criteria, stage_info$stage_key
  )
  diagnostic_comparison <- .ipi_compare_same_stage(
    partial, full, parsed_free$canonical
  )

  releases_effective <- nrow(release_verification) == nrow(parsed_free) &&
    all(release_verification$effective_release) &&
    (!nrow(lower_verification) || all(lower_verification$effective_release))
  conclusion_status <- .ipi_resolve_conclusion_status(
    releases_effective,
    nesting$comparable,
    configural_gate,
    full_decision,
    partial_comparison,
    partial_decision
  )
  conclusion <- .ipi_conclusion_text(conclusion_status, language_effective)

  role_records <- list(
    configural = configural,
    prerequisite_full = prerequisite_full,
    prerequisite_matched = prerequisite_matched,
    full = full,
    partial = partial
  )
  prepared$blueprint$effective_syntax <- lapply(role_records, `[[`, "syntax")
  prepared$blueprint$group_partial <- lapply(role_records, `[[`, "group_partial")
  fit_table <- .ipi_fit_table(role_records)
  comparison_table <- do.call(rbind, list(
    .ipi_comparison_row(
      "full_prerequisite_vs_full", full_comparison, full_decision
    ),
    .ipi_comparison_row(
      "prerequisite_vs_partial", partial_comparison, partial_decision
    ),
    .ipi_comparison_row(
      "full_vs_partial_diagnostic", diagnostic_comparison, NULL
    )
  ))
  rownames(comparison_table) <- NULL

  fit_conditions <- unlist(lapply(role_records, function(record) {
    c(record$warnings, record$messages)
  }), recursive = FALSE)
  fit_errors <- Filter(
    Negate(is.null),
    lapply(role_records, `[[`, "error")
  )
  methodological_conditions <- list()
  if (!releases_effective) {
    methodological_conditions[[length(methodological_conditions) + 1L]] <-
      .ifa_condition(
        "At least one requested equality family was not effectively released by the engine.",
        "partial_release_not_effective",
        type = "error"
      )
  }
  if (!isTRUE(nesting$comparable)) {
    methodological_conditions[[length(methodological_conditions) + 1L]] <-
      .ifa_condition(
        "Matched-prerequisite and target release sets were not nested exactly.",
        "mismatch_group_partial_nesting",
        type = "error"
      )
  }
  if (prepared$resolved_group$dropped_missing_group > 0L) {
    fit_conditions[[length(fit_conditions) + 1L]] <- .ifa_condition(
      sprintf(
        "%d cases with missing group membership were removed.",
        prepared$resolved_group$dropped_missing_group
      ),
      "missing_group_rows_removed"
    )
  }
  if (prepared$dropped_missing_items > 0L) {
    fit_conditions[[length(fit_conditions) + 1L]] <- .ifa_condition(
      sprintf(
        "%d cases were removed by the item missing-data policy.",
        prepared$dropped_missing_items
      ),
      "missing_item_rows_removed"
    )
  }

  models <- lapply(role_records, .ipi_model_envelope)
  models$prerequisite <- models$prerequisite_matched
  comparisons <- list(
    configural_gate = configural_gate,
    full_prerequisite_vs_full = full_comparison,
    prerequisite_vs_partial = partial_comparison,
    full_vs_partial_diagnostic = diagnostic_comparison,
    group_partial_nesting = nesting
  )
  tables <- list(
    fit = fit_table,
    freed_constraints = freed_constraints,
    rationale = rationale_table,
    anchors = anchors$items,
    anchor_counts = anchors$counts,
    comparisons = comparison_table,
    criteria = data.frame(
      criterion = names(prepared$criteria),
      value = unlist(prepared$criteria, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  )

  result <- .pm_result(
    class = "psychomatic_invariance_partial",
    call = call,
    requested_settings = requested_settings,
    effective_settings = list(
      language = language_effective,
      workflow = prepared$workflow,
      fit_variant = prepared$fit_variant,
      source_measurement = prepared$source_measurement,
      stage = stage_info$stage_effective,
      stage_key = stage_info$stage_key,
      measurement = prepared$items$measurement,
      category_profile = prepared$items$category_profile,
      estimator_engine_request = prepared$estimator_info$estimator,
      estimator = configural$estimator_effective,
      estimator_original = configural$estimator_orig,
      ordered_engine_request = prepared$blueprint$ordered_engine_request,
      ordered = configural$ordered,
      ordered_effective = configural$ordered,
      missing_engine_request = prepared$missing_request,
      missing = configural$missing,
      parameterization_engine_request = prepared$blueprint$parameterization_engine_request,
      parameterization = configural$parameterization,
      identification = configural$identification,
      std_lv = configural$std_lv,
      ID.fac = configural$ID.fac_effective,
      ID.cat = configural$ID.cat_effective,
      ID.fac_engine_request = configural$ID.fac_engine_request,
      ID.fac_effective = configural$ID.fac_effective,
      ID.cat_engine_request = configural$ID.cat_engine_request,
      ID.cat_effective = configural$ID.cat_effective,
      min_anchor_items = as.integer(min_anchor_items)
    ),
    data_summary = list(
      n = prepared$sample_validation$total_n,
      group_sizes = prepared$group_sizes,
      indicators = prepared$indicators,
      category_counts = prepared$items$category_counts,
      dropped_missing_group = prepared$resolved_group$dropped_missing_group,
      dropped_missing_items = prepared$dropped_missing_items
    ),
    results = list(
      fit_table = fit_table,
      freed_constraints = freed_constraints,
      anchors = anchors,
      comparison = comparisons
    ),
    warnings = c(parsed_model$conditions, prepared$normality$conditions, fit_conditions),
    errors = c(fit_errors, methodological_conditions),
    diagnostics = list(
      fit = .ifa_diagnostics_long(role_records),
      release_verification = release_verification,
      lower_release_verification = lower_verification,
      comparison_table = comparison_table
    ),
    decision_log = do.call(rbind, list(
      .pm_decision(
        "partial_measurement", requested_measurement,
        prepared$items$measurement,
        "Shared complete-invariance preparation determined the effective measurement type."
      ),
      .pm_decision(
        "partial_stage", requested_stage, stage_info$stage_effective,
        paste("Effective stage:", stage_info$stage_key)
      ),
      .pm_decision(
        "partial_conclusion", "researcher_directed_release", conclusion_status,
        "The conclusion uses configural, admissibility, nesting, and joint fit-change gates."
      )
    )),
    conclusion = list(status = conclusion_status, reason = conclusion)
  )
  references <- .ifa_reference_table()
  settings_for_report <- result$effective_settings[c(
    "language", "stage", "stage_key", "measurement", "category_profile",
    "estimator_engine_request", "estimator", "estimator_original",
    "ordered_engine_request", "ordered_effective",
    "missing_engine_request", "missing",
    "parameterization_engine_request", "parameterization",
    "identification", "std_lv",
    "ID.fac_engine_request", "ID.fac_effective",
    "ID.cat_engine_request", "ID.cat_effective", "min_anchor_items"
  )]
  settings_table <- data.frame(
    setting = names(settings_for_report),
    value = vapply(settings_for_report, function(value) {
      if (is.null(value) || !length(value)) return("")
      if (all(is.na(value))) return(NA_character_)
      paste(as.character(value), collapse = "; ")
    }, character(1L)),
    stringsAsFactors = FALSE
  )
  tables <- c(tables, list(
    warnings = .ifa_conditions_for_report(result$warnings),
    errors = .ifa_conditions_for_report(result$errors),
    conclusion = data.frame(
      status = conclusion_status,
      conclusion = conclusion,
      stringsAsFactors = FALSE
    ),
    settings = settings_table,
    references = references
  ))
  direct <- list(
    measurement = prepared$items$measurement,
    stage = stage_info$stage_key,
    stage_requested = stage_info$stage_requested,
    stage_effective = stage_info$stage_effective,
    prerequisite_stage = stage_info$prerequisite_key,
    analysis_blueprint = prepared$blueprint,
    analysis_data = prepared$analysis_data,
    models = models,
    fitted_models = lapply(models[setdiff(names(models), "prerequisite")], `[[`, "fit"),
    freed_constraints = freed_constraints,
    rationale = requested_rationale,
    rationale_table = rationale_table,
    anchors = anchors$items,
    anchor_counts = anchors$counts,
    comparison = comparisons,
    fit_table = fit_table,
    criteria = prepared$criteria,
    tables = tables,
    references = references,
    conclusion_status = conclusion_status,
    conclusion_text = conclusion,
    report_path = report_path
  )
  for (field in names(direct)) result[[field]] <- direct[[field]]

  .pm_finish_analysis(result, output, language_effective)
}

#' @export
print.psychomatic_invariance_partial <- function(x, digits = NULL, ...) {
  .pm_console_report(x, digits)
}

#' @export
plot.psychomatic_invariance_partial <- function(x, ...) {
  rows <- match(c("full", "partial"), x$fit_table$role)
  plot_data <- x$fit_table[rows, c("role", "cfi", "rmsea", "srmr"), drop = FALSE]
  long <- data.frame(
    model = rep(plot_data$role, times = 3L),
    index = rep(c("CFI", "RMSEA", "SRMR"), each = nrow(plot_data)),
    value = c(plot_data$cfi, plot_data$rmsea, plot_data$srmr),
    stringsAsFactors = FALSE
  )
  finite <- is.finite(long$value)
  if (!any(finite)) {
    graphics::plot.new()
    graphics::title("Partial invariance fit indices unavailable")
    return(invisible(long))
  }
  graphics::matplot(
    t(as.matrix(plot_data[c("cfi", "rmsea", "srmr")])),
    type = "b", pch = c(16, 17), lty = 1,
    xaxt = "n", xlab = "Fit index", ylab = "Value", ...
  )
  graphics::axis(1, at = 1:3, labels = c("CFI", "RMSEA", "SRMR"))
  graphics::legend(
    "topright", legend = plot_data$role,
    pch = c(16, 17), lty = 1, bty = "n"
  )
  invisible(long)
}
