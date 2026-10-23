.pm_measurement_level <- function(x) {
  observed <- x[!is.na(x)]
  if (!length(observed)) return("invalid")
  if (is.ordered(x)) return("ordinal")
  if (is.factor(x) || is.character(x) || is.logical(x)) return("ambiguous")
  if (!is.numeric(x) || inherits(x, c("Date", "POSIXt"))) return("ineligible")
  values <- unique(observed[is.finite(observed)])
  if (length(values) == 1L) return("invalid")
  if (length(values) <= 7L && all(abs(values - round(values)) < 1e-8)) return("ordinal")
  "continuous"
}

.pm_as_data_frame <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }
  out <- as.data.frame(
    data,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    optional = TRUE
  )
  if (ncol(out) == 0L) {
    stop("`data` must contain at least one column.", call. = FALSE)
  }
  if (is.null(names(out)) || length(names(out)) != ncol(out)) {
    names(out) <- paste0("V", seq_len(ncol(out)))
  }
  out
}

.pm_resolve_columns <- function(selection, data_names, name, multiple = TRUE) {
  if (is.null(selection)) return(integer())
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("Internal column-selection name is invalid.", call. = FALSE)
  }

  if (is.character(selection)) {
    if (length(selection) == 0L || anyNA(selection) || any(!nzchar(selection))) {
      stop(sprintf("`%s` must contain non-empty column names.", name), call. = FALSE)
    }
    occurrences <- vapply(
      selection,
      function(value) sum(!is.na(data_names) & data_names == value),
      integer(1)
    )
    if (any(occurrences == 0L)) {
      stop(
        sprintf(
          "Unknown columns in `%s`: %s.",
          name,
          paste(unique(selection[occurrences == 0L]), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (any(occurrences > 1L)) {
      stop(
        sprintf(
          "Duplicated column names in `%s` are ambiguous; use positive indices.",
          name
        ),
        call. = FALSE
      )
    }
    positions <- match(selection, data_names)
  } else if (is.numeric(selection)) {
    if (length(selection) == 0L || anyNA(selection) ||
        any(!is.finite(selection)) || any(selection != floor(selection)) ||
        any(selection < 1L) || any(selection > length(data_names))) {
      stop(
        sprintf("`%s` indices must be positive integers within `data`.", name),
        call. = FALSE
      )
    }
    positions <- as.integer(selection)
  } else {
    stop(sprintf("`%s` must contain column names or positive indices.", name), call. = FALSE)
  }

  if (!multiple && length(positions) != 1L) {
    stop(sprintf("`%s` must select exactly one column.", name), call. = FALSE)
  }
  if (anyDuplicated(positions)) {
    stop(sprintf("`%s` must not select a column more than once.", name), call. = FALSE)
  }
  positions
}

.pm_normalized_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^[:alnum:]]+", "_", x)
  gsub("^_+|_+$", "", x)
}

.pm_identifier_reason <- function(x, name) {
  normalized <- .pm_normalized_name(name)
  by_name <- grepl(
    "(^|_)(id|case|participant|respondent)($|_)",
    normalized,
    perl = TRUE
  )

  by_unique <- FALSE
  by_sequence <- FALSE
  if (is.numeric(x) && !inherits(x, c("Date", "POSIXt"))) {
    finite <- x[is.finite(x)]
    if (length(finite) > 0L) {
      by_unique <- length(unique(finite)) / length(finite) >= 0.95
    }
    if (length(finite) > 1L && all(abs(finite - round(finite)) < 1e-8)) {
      differences <- diff(finite)
      by_sequence <- all(differences > 0) || all(differences < 0)
    }
  }

  reasons <- c(
    if (by_name) "identifier-like name",
    if (by_unique) "at least 95% unique finite values",
    if (by_sequence) "strictly monotone integer sequence"
  )
  if (!length(reasons)) NA_character_ else paste(reasons, collapse = "; ")
}

.pm_level_overrides <- function(level, item_names) {
  if (is.list(level)) {
    level_names <- names(level)
    level <- unlist(level, recursive = FALSE, use.names = TRUE)
    if (!is.null(level_names) && is.null(names(level))) names(level) <- level_names
  }
  if (!is.character(level) || length(level) == 0L || anyNA(level)) {
    stop("`level` must contain 'auto', 'ordinal', or 'continuous'.", call. = FALSE)
  }

  supplied_names <- names(level)
  level <- tolower(trimws(as.character(level)))
  names(level) <- supplied_names
  allowed <- c("auto", "ordinal", "continuous")
  if (any(!level %in% allowed)) {
    stop("`level` must contain only 'auto', 'ordinal', or 'continuous'.", call. = FALSE)
  }

  is_named <- !is.null(supplied_names) && any(nzchar(supplied_names))
  if (!is_named && length(level) == 1L) {
    return(rep(level, length(item_names)))
  }
  if (!is_named && length(level) == length(item_names)) {
    return(unname(level))
  }
  if (is_named && all(nzchar(supplied_names))) {
    if (anyDuplicated(supplied_names)) {
      stop("Named `level` overrides must have unique names.", call. = FALSE)
    }
    unknown <- setdiff(supplied_names, unique(item_names))
    if (length(unknown)) {
      stop(
        sprintf("Unknown item names in `level`: %s.", paste(unknown, collapse = ", ")),
        call. = FALSE
      )
    }
    ambiguous <- supplied_names[
      vapply(supplied_names, function(value) sum(item_names == value), integer(1)) > 1L
    ]
    if (length(ambiguous)) {
      stop("Named `level` overrides cannot target duplicated item names.", call. = FALSE)
    }
    out <- rep("auto", length(item_names))
    out[match(supplied_names, item_names)] <- unname(level)
    return(out)
  }

  stop(
    "`level` must be scalar, align with `items`, or be a fully named override vector.",
    call. = FALSE
  )
}

.pm_apply_ordinal_levels <- function(data, item_positions, item_names, ordinal_levels) {
  if (is.null(ordinal_levels)) {
    return(list(data = data, levels = list(), mapped = rep(FALSE, length(item_names))))
  }
  if (!is.list(ordinal_levels) || is.null(names(ordinal_levels)) ||
      anyNA(names(ordinal_levels)) || any(!nzchar(names(ordinal_levels))) ||
      anyDuplicated(names(ordinal_levels))) {
    stop("`ordinal_levels` must be a uniquely named list.", call. = FALSE)
  }
  unknown <- setdiff(names(ordinal_levels), unique(item_names))
  if (length(unknown)) {
    stop(
      sprintf(
        "Unknown item names in `ordinal_levels`: %s.",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  duplicated_targets <- names(ordinal_levels)[
    vapply(
      names(ordinal_levels),
      function(value) sum(item_names == value),
      integer(1)
    ) > 1L
  ]
  if (length(duplicated_targets)) {
    stop("`ordinal_levels` cannot target duplicated item names.", call. = FALSE)
  }

  normalized <- list()
  mapped <- rep(FALSE, length(item_names))
  for (item in names(ordinal_levels)) {
    declared <- ordinal_levels[[item]]
    if (!is.atomic(declared) || length(declared) == 0L || anyNA(declared)) {
      stop(
        sprintf("Ordinal levels for `%s` must be a non-empty atomic vector.", item),
        call. = FALSE
      )
    }
    declared <- as.character(declared)
    if (anyDuplicated(declared)) {
      stop(sprintf("Ordinal levels for `%s` must be unique.", item), call. = FALSE)
    }

    item_offset <- match(item, item_names)
    data_offset <- item_positions[item_offset]
    observed <- unique(as.character(data[[data_offset]][!is.na(data[[data_offset]])]))
    uncovered <- setdiff(observed, declared)
    if (length(uncovered)) {
      stop(
        sprintf(
          "`ordinal_levels` for `%s` must cover all observed values: %s.",
          item,
          paste(uncovered, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    data[[data_offset]] <- factor(
      as.character(data[[data_offset]]),
      levels = declared,
      ordered = TRUE
    )
    normalized[[item]] <- declared
    mapped[item_offset] <- TRUE
  }
  list(data = data, levels = normalized, mapped = mapped)
}

.pm_normalize_ranges <- function(ranges, item_names) {
  if (is.null(ranges)) return(list())
  if (is.matrix(ranges) || is.data.frame(ranges)) {
    if (ncol(ranges) != 2L || is.null(rownames(ranges)) ||
        any(!nzchar(rownames(ranges)))) {
      stop("Matrix/data-frame `ranges` need two columns and item row names.", call. = FALSE)
    }
    ranges <- lapply(seq_len(nrow(ranges)), function(i) unlist(ranges[i, ], use.names = FALSE))
    names(ranges) <- rownames(ranges)
  }
  if (!is.list(ranges) || is.null(names(ranges)) || anyNA(names(ranges)) ||
      any(!nzchar(names(ranges))) || anyDuplicated(names(ranges))) {
    stop("`ranges` must be a uniquely named list of lower/upper limits.", call. = FALSE)
  }
  unknown <- setdiff(names(ranges), unique(item_names))
  if (length(unknown)) {
    stop(
      sprintf("Unknown item names in `ranges`: %s.", paste(unknown, collapse = ", ")),
      call. = FALSE
    )
  }
  duplicated_targets <- names(ranges)[
    vapply(names(ranges), function(value) sum(item_names == value), integer(1)) > 1L
  ]
  if (length(duplicated_targets)) {
    stop("`ranges` cannot target duplicated item names.", call. = FALSE)
  }

  lapply(names(ranges), function(item) {
    limits <- ranges[[item]]
    if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits)) {
      stop(sprintf("Range for `%s` must contain two finite numbers.", item), call. = FALSE)
    }
    lower <- .pm_validate_scalar_number(limits[1L], paste0("ranges$", item, "[1]"))
    upper <- .pm_validate_scalar_number(limits[2L], paste0("ranges$", item, "[2]"))
    if (lower >= upper) {
      stop(sprintf("Range for `%s` must be strictly increasing.", item), call. = FALSE)
    }
    c(min = lower, max = upper)
  }) |>
    stats::setNames(names(ranges))
}

.pm_measurement_reason <- function(x, effective, mapped, override) {
  if (mapped) return("explicit ordinal levels")
  if (!identical(override, "auto")) return("explicit measurement-level override")
  if (identical(effective, "ordinal") && is.ordered(x)) return("ordered factor")
  if (identical(effective, "ordinal")) return("two to seven finite integer values")
  if (identical(effective, "continuous")) return("more than seven categories or non-integer values")
  if (identical(effective, "invalid")) return("no usable variation")
  if (identical(effective, "ambiguous")) return("unordered categorical encoding")
  "ineligible data class"
}

#' Define a psychometric data specification
#'
#' Selects item columns without silently treating identifiers or unordered
#' categories as measurements. It records item order, measurement levels,
#' ordinal encodings, reverse keys, ranges, grouping, metadata, missing-data
#' policy, language, and the decisions used to build the specification.
#'
#' @details With `level = "auto"`, ordered factors are ordinal. Ordinary
#'   numeric columns with two through seven distinct finite integer values are
#'   ordinal; other varying numeric columns are continuous. Unordered factors,
#'   character vectors, and logical vectors are ambiguous and are never sorted
#'   alphabetically. Explicit `ordinal_levels` take precedence over an explicit
#'   per-item `level`, which in turn takes precedence over automatic
#'   classification.
#'
#'   When `items = NULL`, a candidate is treated as a possible identifier if
#'   its normalized name contains `id`, `case`, `participant`, or `respondent`,
#'   if at least 95 percent of its finite numeric values are unique, or if those
#'   values form a strictly monotone integer sequence. The function then stops
#'   and asks for an explicit `items` or `metadata` declaration.
#'
#' @param data A data frame or matrix containing item responses and optional
#'   metadata/group columns.
#' @param items Item column names or positive integer indices. If `NULL`, all
#'   eligible columns except `metadata` and `group` are considered, and likely
#'   identifiers cause an error rather than being silently selected.
#' @param metadata Optional metadata column names or positive integer indices.
#' @param level Measurement level: `"auto"`, `"ordinal"`, or `"continuous"`.
#'   It may be scalar, aligned with `items`, or a named vector of per-item
#'   overrides.
#' @param ordinal_levels Optional uniquely named list giving the ordered levels
#'   for selected items. Declared levels take precedence over `level` and must
#'   cover every observed non-missing value.
#' @param reverse Optional item names or positive integer indices identifying
#'   reverse-keyed items. This records keys; it does not transform responses.
#' @param ranges Optional uniquely named list of numeric vectors `c(min, max)`,
#'   or a two-column matrix/data frame with item row names.
#' @param group Optional single group-column name or positive integer index.
#' @param missing Character scalar recording the missing-data policy. The
#'   default is `"available"`.
#' @param language Output language, `"eng"` by default; `"esp"` for Spanish.
#'   Explicit `NULL` uses `getOption("PsychoMatic.language", "eng")`.
#'
#' @return An object of classes `psychomatic_spec` and `psychomatic_result`.
#'   In addition to the common result-contract fields, direct fields expose
#'   `data`, the pre-mapping `raw_data`, `items`, `measurement_level`,
#'   `ordinal_levels`, `reverse`, `ranges`, `group`, `metadata`, `missing`,
#'   and `language`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' psychomatic_spec(MY_DATA, items = c("R1", "R2", "R3", "R4"))
#' }
#'
#' @export
psychomatic_spec <- function(data,
                             items = NULL,
                             metadata = NULL,
                             level = "auto",
                             ordinal_levels = NULL,
                             reverse = NULL,
                             ranges = NULL,
                             group = NULL,
                             missing = "available",
                             language = "eng") {
  call <- match.call()
  source_is_matrix <- is.matrix(data)
  data <- .pm_as_data_frame(data)
  raw_data <- data
  data_names <- names(data)
  language_effective <- .pm_normalize_language(language)

  if (!is.character(missing) || length(missing) != 1L || is.na(missing) ||
      !nzchar(trimws(missing))) {
    stop("`missing` must be a non-empty character scalar.", call. = FALSE)
  }
  missing <- trimws(missing)

  metadata_positions <- .pm_resolve_columns(metadata, data_names, "metadata")
  group_position <- .pm_resolve_columns(group, data_names, "group", multiple = FALSE)
  if (length(intersect(metadata_positions, group_position))) {
    stop("`metadata` and `group` must not select the same column.", call. = FALSE)
  }
  excluded <- c(metadata_positions, group_position)

  automatic_items <- is.null(items)
  if (automatic_items) {
    candidate_positions <- setdiff(seq_along(data), excluded)
    candidate_names <- data_names[candidate_positions]
    candidate_identifier_reasons <- vapply(
      seq_along(candidate_positions),
      function(i) {
        .pm_identifier_reason(
          raw_data[[candidate_positions[i]]],
          candidate_names[i]
        )
      },
      character(1)
    )
    possible_candidate_identifiers <- !is.na(candidate_identifier_reasons)
    if (any(possible_candidate_identifiers)) {
      stop(
        sprintf(
          "Columns flagged as possible identifiers: %s. Declare `items` or `metadata` explicitly.",
          paste(candidate_names[possible_candidate_identifiers], collapse = ", ")
        ),
        call. = FALSE
      )
    }
    eligible <- vapply(
      data[candidate_positions],
      function(x) {
        is.ordered(x) ||
          (is.numeric(x) && !inherits(x, c("Date", "POSIXt")))
      },
      logical(1)
    )
    item_positions <- candidate_positions[eligible]
  } else {
    item_positions <- .pm_resolve_columns(items, data_names, "items")
    overlap <- intersect(item_positions, excluded)
    if (length(overlap)) {
      stop("Columns declared as `metadata` or `group` cannot also be items.", call. = FALSE)
    }
  }
  if (!length(item_positions)) {
    stop("No eligible item columns were selected.", call. = FALSE)
  }

  item_names <- data_names[item_positions]
  identifier_reasons <- vapply(
    seq_along(item_positions),
    function(i) .pm_identifier_reason(data[[item_positions[i]]], item_names[i]),
    character(1)
  )
  possible_identifiers <- !is.na(identifier_reasons)

  level_overrides <- .pm_level_overrides(level, item_names)
  ordinal <- .pm_apply_ordinal_levels(
    data,
    item_positions,
    item_names,
    ordinal_levels
  )
  data <- ordinal$data
  normalized_levels <- ordinal$levels

  measurement_level <- character(length(item_positions))
  decision_rows <- vector("list", length(item_positions))
  for (i in seq_along(item_positions)) {
    x <- data[[item_positions[i]]]
    effective <- if (ordinal$mapped[i]) {
      "ordinal"
    } else if (!identical(level_overrides[i], "auto")) {
      level_overrides[i]
    } else {
      .pm_measurement_level(x)
    }
    measurement_level[i] <- effective
    decision_rows[[i]] <- .pm_decision(
      code = "measurement_level",
      requested = if (ordinal$mapped[i]) "ordinal_levels" else level_overrides[i],
      effective = effective,
      reason = .pm_measurement_reason(
        x,
        effective,
        ordinal$mapped[i],
        level_overrides[i]
      ),
      severity = if (effective %in% c("ambiguous", "invalid", "ineligible")) {
        "error"
      } else {
        "info"
      }
    )
  }
  names(measurement_level) <- item_names

  ambiguous <- item_names[measurement_level == "ambiguous"]
  if (length(ambiguous)) {
    stop(
      sprintf(
        "Items with ambiguous measurement encoding: %s. Supply `ordinal_levels` or an explicit `level` override; values were not alphabetically ordered.",
        paste(ambiguous, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unusable <- item_names[measurement_level %in% c("invalid", "ineligible")]
  if (length(unusable)) {
    stop(
      sprintf(
        "Items are invalid or ineligible under automatic classification: %s. Supply an explicit defensible `level` override or revise the data.",
        paste(unusable, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  reverse_positions <- .pm_resolve_columns(reverse, data_names, "reverse")
  unknown_reverse <- setdiff(reverse_positions, item_positions)
  if (length(unknown_reverse)) {
    stop("Every `reverse` key must identify a selected item.", call. = FALSE)
  }
  reverse_names <- data_names[reverse_positions]
  normalized_ranges <- .pm_normalize_ranges(ranges, item_names)
  decisions <- do.call(rbind, decision_rows)
  rownames(decisions) <- NULL

  requested_settings <- list(
    items = items,
    metadata = metadata,
    level = level,
    ordinal_levels = ordinal_levels,
    reverse = reverse,
    ranges = ranges,
    group = group,
    missing = missing,
    language = language
  )
  effective_settings <- list(
    items = item_names,
    metadata = data_names[metadata_positions],
    level = measurement_level,
    ordinal_levels = normalized_levels,
    reverse = reverse_names,
    ranges = normalized_ranges,
    group = if (length(group_position)) data_names[group_position] else NULL,
    missing = missing,
    language = language_effective
  )
  results <- list(
    data = data,
    raw_data = raw_data,
    items = item_names,
    item_index = item_positions,
    measurement_level = measurement_level,
    ordinal_levels = normalized_levels,
    reverse = reverse_names,
    reverse_index = reverse_positions,
    ranges = normalized_ranges,
    group = if (length(group_position)) data_names[group_position] else NULL,
    metadata = data_names[metadata_positions],
    missing = missing,
    language = language_effective,
    source_is_matrix = source_is_matrix
  )

  out <- .pm_result(
    class = "psychomatic_spec",
    call = call,
    requested_settings = requested_settings,
    effective_settings = effective_settings,
    data_summary = list(
      n = nrow(data),
      p = length(item_positions),
      metadata_columns = length(metadata_positions),
      has_group = length(group_position) == 1L
    ),
    results = results,
    warnings = NULL,
    errors = NULL,
    diagnostics = list(
      possible_identifier = data.frame(
        item = item_names[possible_identifiers],
        reason = identifier_reasons[possible_identifiers],
        stringsAsFactors = FALSE
      )
    ),
    decision_log = decisions,
    conclusion = NA_character_
  )
  for (field in names(results)) out[field] <- results[field]
  out
}

#' Print a psychometric data specification
#'
#' @param x A `psychomatic_spec` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#'
#' @method print psychomatic_spec
#' @keywords internal
#' @export
print.psychomatic_spec <- function(x, ...) {
  language <- x$effective_settings$language
  title <- if (identical(language, "esp")) {
    "Especificacion de datos PsychoMatic"
  } else {
    "PsychoMatic data specification"
  }
  cat(title, "\n", sep = "")
  cat("  n = ", x$data_summary$n, ", items = ", x$data_summary$p, "\n", sep = "")
  print(
    data.frame(
      item = x$items,
      measurement_level = unname(x$measurement_level),
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  invisible(x)
}

#' Summarize a psychometric data specification
#'
#' @param object A `psychomatic_spec` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list containing the data summary, item specification table, and
#'   decision log.
#'
#' @method summary psychomatic_spec
#' @keywords internal
#' @export
summary.psychomatic_spec <- function(object, ...) {
  range_min <- rep(NA_real_, length(object$items))
  range_max <- rep(NA_real_, length(object$items))
  if (length(object$ranges)) {
    for (i in seq_along(object$items)) {
      position <- match(object$items[i], names(object$ranges))
      if (!is.na(position)) {
        range_min[i] <- object$ranges[[position]][1L]
        range_max[i] <- object$ranges[[position]][2L]
      }
    }
  }
  list(
    data_summary = object$data_summary,
    items = data.frame(
      item = object$items,
      measurement_level = unname(object$measurement_level),
      reverse = if (!is.null(object$reverse_index)) {
        object$item_index %in% object$reverse_index
      } else {
        object$items %in% object$reverse
      },
      range_min = range_min,
      range_max = range_max,
      stringsAsFactors = FALSE
    ),
    decision_log = object$decision_log
  )
}

.pm_empty_issue_table <- function() {
  data.frame(
    code = character(),
    severity = character(),
    item = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

.pm_add_issue <- function(issues, code, severity, item, message) {
  rbind(
    issues,
    data.frame(
      code = as.character(code),
      severity = as.character(severity),
      item = as.character(item),
      message = as.character(message),
      stringsAsFactors = FALSE
    )
  )
}

.pm_item_label <- function(item) {
  if (length(item) != 1L || is.na(item) || !nzchar(item)) {
    "<unnamed>"
  } else {
    paste0("`", item, "`")
  }
}

.pm_available_values <- function(x) {
  if (is.numeric(x)) is.finite(x) else !is.na(x)
}

.pm_severity_counts <- function(issues) {
  levels <- c("error", "warning", "info")
  counts <- table(factor(issues$severity, levels = levels))
  data.frame(
    severity = levels,
    n = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

#' Check psychometric item data
#'
#' Inspects a `psychomatic_spec`, data frame, or matrix without fitting a model.
#' It returns stable issue codes for missing/non-finite data, low variation,
#' categorical sparsity, unsafe names, range violations, pairwise overlap, and
#' correlation-matrix validity.
#'
#' @details Near-zero variance means that at least 95 percent of usable values
#'   occupy one response value while more than one value is observed. An
#'   observed ordinal category is sparse below five cases; declared factor
#'   levels with zero cases are empty. Pairwise overlap below three usable cases
#'   is insufficient. Numeric square matrices are additionally checked for
#'   symmetry, unit diagonal and correlation bounds, and positive definiteness.
#'
#'   Stable issue codes are `all_missing`, `non_finite`, `zero_variance`,
#'   `near_zero_variance`, `sparse_category`, `empty_category`,
#'   `duplicate_name`, `nonsyntactic_name`, `range_violation`,
#'   `insufficient_overlap`, `non_symmetric_matrix`,
#'   `not_correlation_matrix`, and `not_positive_definite`.
#'
#' @param x A `psychomatic_spec`, data frame, or matrix.
#' @param language Output language. `NULL` inherits the language of a
#'   `psychomatic_spec`; otherwise it uses
#'   `getOption("PsychoMatic.language", "eng")`.
#'
#' @return An object of classes `psychomatic_check` and `psychomatic_result`.
#'   The direct `issues` field is a data frame with `code`, `severity`, `item`,
#'   and `message`; `severity_counts` contains deterministic counts.
#'
#' @examples
#' spec <- psychomatic_spec(
#'   data.frame(q1 = c(1, 1, 1), q2 = c(1, 2, 3)),
#'   items = c("q1", "q2"),
#'   level = "continuous"
#' )
#' check_items(spec)
#'
#' @keywords internal
#' @noRd
check_items <- function(x, language = NULL) {
  call <- match.call()
  is_spec <- inherits(x, "psychomatic_spec")
  if (is_spec) {
    inherited_language <- x$effective_settings$language
    language_effective <- .pm_normalize_language(language, inherited_language)
    data <- x$data
    raw_data <- if (is.null(x$raw_data)) x$data else x$raw_data
    item_positions <- x$item_index
    item_names <- x$items
    measurement_level <- unname(x$measurement_level)
    ranges <- x$ranges
    source_is_matrix <- isTRUE(x$source_is_matrix)
  } else {
    language_effective <- .pm_normalize_language(language)
    source_is_matrix <- is.matrix(x)
    data <- .pm_as_data_frame(x)
    raw_data <- data
    item_positions <- seq_along(data)
    item_names <- names(data)
    measurement_level <- vapply(data, .pm_measurement_level, character(1))
    ranges <- list()
  }

  item_data <- data[item_positions]
  names(item_data) <- item_names
  raw_item_data <- raw_data[item_positions]
  names(raw_item_data) <- item_names
  issues <- .pm_empty_issue_table()
  text_for <- function(esp, eng) {
    if (identical(language_effective, "esp")) esp else eng
  }

  duplicated_names <- item_names[
    duplicated(item_names) | duplicated(item_names, fromLast = TRUE)
  ]
  for (item in unique(duplicated_names)) {
    issues <- .pm_add_issue(
      issues,
      "duplicate_name",
      "error",
      item,
      text_for(
        paste0("El nombre de item ", .pm_item_label(item), " esta duplicado."),
        paste0("Item name ", .pm_item_label(item), " is duplicated.")
      )
    )
  }

  syntactic <- make.names(item_names, unique = FALSE)
  nonsyntactic <- which(is.na(item_names) | !nzchar(item_names) | syntactic != item_names)
  for (i in nonsyntactic) {
    issues <- .pm_add_issue(
      issues,
      "nonsyntactic_name",
      "warning",
      item_names[i],
      text_for(
        paste0("El nombre de item ", .pm_item_label(item_names[i]), " no es sintactico en R."),
        paste0("Item name ", .pm_item_label(item_names[i]), " is not syntactic in R.")
      )
    )
  }

  for (i in seq_along(item_data)) {
    values <- raw_item_data[[i]]
    mapped_values <- item_data[[i]]
    item <- item_names[i]
    label <- .pm_item_label(item)
    all_missing <- length(values) == 0L || all(is.na(values))
    if (all_missing) {
      issues <- .pm_add_issue(
        issues,
        "all_missing",
        "error",
        item,
        text_for(
          paste0("El item ", label, " no tiene valores observados."),
          paste0("Item ", label, " has no observed values.")
        )
      )
    }

    if (is.numeric(values)) {
      non_finite_n <- sum(is.nan(values) | is.infinite(values))
      if (non_finite_n > 0L) {
        issues <- .pm_add_issue(
          issues,
          "non_finite",
          "error",
          item,
          text_for(
            paste0("El item ", label, " contiene ", non_finite_n, " valor(es) no finito(s)."),
            paste0("Item ", label, " contains ", non_finite_n, " non-finite value(s).")
          )
        )
      }
      observed <- values[is.finite(values)]
    } else {
      observed <- values[!is.na(values)]
    }

    observed_unique <- unique(observed)
    if (length(observed) > 0L && length(observed_unique) == 1L) {
      issues <- .pm_add_issue(
        issues,
        "zero_variance",
        "error",
        item,
        text_for(
          paste0("El item ", label, " no presenta variacion."),
          paste0("Item ", label, " has zero variance.")
        )
      )
    } else if (length(observed_unique) > 1L) {
      frequencies <- table(observed, useNA = "no")
      dominant <- max(frequencies) / sum(frequencies)
      if (is.finite(dominant) && dominant >= 0.95) {
        issues <- .pm_add_issue(
          issues,
          "near_zero_variance",
          "warning",
          item,
          text_for(
            paste0("El item ", label, " concentra al menos 95% de respuestas en un valor."),
            paste0("Item ", label, " places at least 95% of responses in one value.")
          )
        )
      }
    }

    is_ordinal <- identical(measurement_level[i], "ordinal") || is.ordered(mapped_values)
    if (is_ordinal) {
      category_counts <- table(mapped_values, useNA = "no")
      empty <- names(category_counts)[category_counts == 0L]
      sparse <- names(category_counts)[category_counts > 0L & category_counts < 5L]
      if (length(empty)) {
        issues <- .pm_add_issue(
          issues,
          "empty_category",
          "warning",
          item,
          text_for(
            paste0("El item ", label, " contiene categorias declaradas sin observaciones: ", paste(empty, collapse = ", "), "."),
            paste0("Item ", label, " has declared categories with no observations: ", paste(empty, collapse = ", "), ".")
          )
        )
      }
      if (length(sparse)) {
        issues <- .pm_add_issue(
          issues,
          "sparse_category",
          "warning",
          item,
          text_for(
            paste0("El item ", label, " contiene categorias con menos de cinco observaciones: ", paste(sparse, collapse = ", "), "."),
            paste0("Item ", label, " has categories with fewer than five observations: ", paste(sparse, collapse = ", "), ".")
          )
        )
      }
    }

    range_position <- match(item, names(ranges))
    if (!is.na(range_position) && is.numeric(values)) {
      limits <- ranges[[range_position]]
      outside <- is.finite(values) & (values < limits[1L] | values > limits[2L])
      if (any(outside)) {
        issues <- .pm_add_issue(
          issues,
          "range_violation",
          "error",
          item,
          text_for(
            paste0("El item ", label, " contiene ", sum(outside), " valor(es) fuera del rango declarado."),
            paste0("Item ", label, " contains ", sum(outside), " value(s) outside its declared range.")
          )
        )
      }
    }
  }

  pairwise_n <- matrix(
    0L,
    nrow = length(item_data),
    ncol = length(item_data),
    dimnames = list(item_names, item_names)
  )
  if (length(item_data)) {
    availability <- vapply(
      raw_item_data,
      .pm_available_values,
      logical(nrow(raw_item_data))
    )
    if (length(item_data) == 1L) availability <- matrix(availability, ncol = 1L)
    pairwise_n[] <- as.integer(crossprod(availability))
  }
  if (length(item_data) >= 2L) {
    for (i in seq_len(length(item_data) - 1L)) {
      for (j in seq.int(i + 1L, length(item_data))) {
        if (pairwise_n[i, j] < 3L) {
          pair <- paste(item_names[c(i, j)], collapse = ", ")
          issues <- .pm_add_issue(
            issues,
            "insufficient_overlap",
            "warning",
            pair,
            text_for(
              paste0("Los items `", pair, "` comparten menos de tres casos utilizables."),
              paste0("Items `", pair, "` share fewer than three usable cases.")
            )
          )
        }
      }
    }
  }

  matrix_diagnostics <- list(evaluated = FALSE)
  if (source_is_matrix) {
    matrix_data <- if (is_spec) as.matrix(raw_item_data) else x
    if (is.numeric(matrix_data) && nrow(matrix_data) == ncol(matrix_data)) {
      matrix_diagnostics$evaluated <- TRUE
      tolerance <- sqrt(.Machine$double.eps)
      finite_matrix <- all(is.finite(matrix_data))
      symmetric <- finite_matrix && isTRUE(all.equal(
        matrix_data,
        t(matrix_data),
        tolerance = tolerance,
        check.attributes = FALSE
      ))
      if (!symmetric) {
        issues <- .pm_add_issue(
          issues,
          "non_symmetric_matrix",
          "error",
          NA_character_,
          text_for(
            "La matriz no es simetrica dentro de la tolerancia numerica.",
            "The matrix is not symmetric within numerical tolerance."
          )
        )
      }
      correlation_matrix <- symmetric &&
        all(abs(diag(matrix_data) - 1) <= tolerance) &&
        all(abs(matrix_data) <= 1 + tolerance)
      if (!correlation_matrix) {
        issues <- .pm_add_issue(
          issues,
          "not_correlation_matrix",
          "error",
          NA_character_,
          text_for(
            "La matriz no satisface el contrato de una matriz de correlaciones.",
            "The matrix does not satisfy the correlation-matrix contract."
          )
        )
      }
      minimum_eigenvalue <- NA_real_
      if (symmetric) {
        minimum_eigenvalue <- min(eigen(matrix_data, symmetric = TRUE, only.values = TRUE)$values)
        if (!is.finite(minimum_eigenvalue) || minimum_eigenvalue <= tolerance) {
          issues <- .pm_add_issue(
            issues,
            "not_positive_definite",
            "error",
            NA_character_,
            text_for(
              "La matriz no es positiva definida.",
              "The matrix is not positive definite."
            )
          )
        }
      }
      matrix_diagnostics$symmetric <- symmetric
      matrix_diagnostics$correlation_matrix <- correlation_matrix
      matrix_diagnostics$minimum_eigenvalue <- minimum_eigenvalue
    }
  }

  rownames(issues) <- NULL
  severity_counts <- .pm_severity_counts(issues)
  error_issues <- issues[issues$severity == "error", , drop = FALSE]
  warning_issues <- issues[issues$severity == "warning", , drop = FALSE]
  requested_language <- if (is.null(language)) "inherited_or_option" else as.character(language)[1L]
  conclusion <- if (nrow(error_issues)) {
    .pm_not_evaluable(
      text_for(
        "La validacion detecto errores que requieren correccion.",
        "Validation found errors that require correction."
      )
    )
  } else {
    NA_character_
  }

  out <- .pm_result(
    class = "psychomatic_check",
    call = call,
    requested_settings = list(language = language),
    effective_settings = list(language = language_effective),
    data_summary = list(
      n = nrow(item_data),
      p = ncol(item_data),
      issues = nrow(issues)
    ),
    results = list(issues = issues, severity_counts = severity_counts),
    warnings = warning_issues,
    errors = error_issues,
    diagnostics = list(
      pairwise_n = pairwise_n,
      matrix = matrix_diagnostics
    ),
    decision_log = .pm_decision(
      "language",
      requested_language,
      language_effective,
      if (is.null(language)) "inherited from specification or package option" else "explicit language"
    ),
    conclusion = conclusion
  )
  out["issues"] <- list(issues)
  out["severity_counts"] <- list(severity_counts)
  out
}

#' Print psychometric item checks
#'
#' @param x A `psychomatic_check` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return `x`, invisibly.
#'
#' @method print psychomatic_check
#' @keywords internal
#' @export
print.psychomatic_check <- function(x, ...) {
  language <- x$effective_settings$language
  title <- if (identical(language, "esp")) {
    "Revision de items PsychoMatic"
  } else {
    "PsychoMatic item check"
  }
  cat(title, "\n", sep = "")
  print(x$severity_counts, row.names = FALSE)
  invisible(x)
}

#' Summarize psychometric item checks
#'
#' @param object A `psychomatic_check` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list containing severity counts and the complete issue table.
#'
#' @method summary psychomatic_check
#' @keywords internal
#' @export
summary.psychomatic_check <- function(object, ...) {
  list(severity_counts = object$severity_counts, issues = object$issues)
}

#' Plot psychometric item issues by severity
#'
#' @param x A `psychomatic_check` object.
#' @param ... Additional graphical arguments passed to [graphics::barplot()].
#'
#' @return The plotted severity-count data frame, invisibly.
#'
#' @method plot psychomatic_check
#' @keywords internal
#' @export
plot.psychomatic_check <- function(x, ...) {
  counts <- summary(x)$severity_counts
  language <- x$effective_settings$language
  y_label <- if (identical(language, "esp")) "Numero de incidencias" else "Number of issues"
  title <- if (identical(language, "esp")) "Incidencias por severidad" else "Issues by severity"
  graphics::barplot(
    counts$n,
    names.arg = counts$severity,
    ylab = y_label,
    main = title,
    ...
  )
  invisible(counts)
}
