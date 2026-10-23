#' Build a structured audit trail from PsychoMatic results
#'
#' Extracts the settings, data provenance, conditions, decisions, diagnostics,
#' conclusions, references, package versions, and report artifacts already
#' stored in one or more PsychoMatic result objects. The function does not
#' refit models or recalculate statistics.
#'
#' @param result A PsychoMatic result object, or a list of such objects.
#' @param language Optional audit language. Use `"eng"` or `"esp"`. When
#'   `NULL`, a single result's effective language is inherited when available.
#'
#' @return A visible object of class `"psychomatic_audit"` containing typed,
#'   tabular audit sections and lightweight result tables.
#'
#' @examples
#' items <- data.frame(i1 = 1:6, i2 = c(1, 2, 2, 3, 4, 4))
#' checked <- check_items(items)
#' audit_report(checked, language = "eng")
#'
#' @keywords internal
#' @noRd
audit_report <- function(result, language = NULL) {
  call <- match.call()
  collected <- .pm_audit_objects(result)
  inherited_language <- .pm_audit_inherited_language(collected$values)
  language_effective <- .pm_normalize_language(language, inherited_language)

  sections <- .pm_build_audit_sections(
    collected$values,
    collected$metadata,
    language = language_effective
  )

  out <- .pm_result(
    class = "psychomatic_audit",
    call = call,
    requested_settings = list(language = language),
    effective_settings = list(
      language = language_effective,
      n_objects = nrow(collected$metadata)
    ),
    data_summary = list(n_objects = nrow(collected$metadata)),
    results = sections$results,
    warnings = sections$warnings,
    errors = sections$errors,
    diagnostics = sections$diagnostics,
    decision_log = sections$decisions,
    conclusion = list(
      status = if (any(sections$conclusions$status == "not_evaluable", na.rm = TRUE)) {
        "not_evaluable"
      } else {
        "evaluable"
      },
      reason = if (any(sections$conclusions$status == "not_evaluable", na.rm = TRUE)) {
        "At least one audited result is not evaluable; inspect `conclusions`."
      } else {
        NA_character_
      }
    )
  )

  out$objects <- sections$objects
  out$settings <- sections$settings
  out$data <- sections$data
  out$decisions <- sections$decisions
  out$conclusions <- sections$conclusions
  out$references <- sections$references
  out$versions <- sections$versions
  out$artifacts <- sections$artifacts
  out
}

#' @keywords internal
#' @export
print.psychomatic_audit <- function(x, ...) {
  language <- .pm_normalize_language(
    x$effective_settings$language,
    inherited = "esp"
  )
  title <- if (language == "esp") "Auditoria PsychoMatic" else "PsychoMatic audit"
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat(if (language == "esp") "Objetos: " else "Objects: ", nrow(x$objects), "\n", sep = "")
  cat(if (language == "esp") "Advertencias: " else "Warnings: ", nrow(x$warnings), "\n", sep = "")
  cat(if (language == "esp") "Errores: " else "Errors: ", nrow(x$errors), "\n", sep = "")
  if (nrow(x$conclusions)) {
    counts <- table(x$conclusions$status, useNA = "ifany")
    cat(if (language == "esp") "Estados: " else "Statuses: ")
    cat(paste(names(counts), as.integer(counts), sep = "=", collapse = ", "), "\n")
  }
  invisible(x)
}

#' @keywords internal
#' @export
summary.psychomatic_audit <- function(object, ...) {
  object[c(
    "objects", "settings", "data", "warnings", "errors", "decisions",
    "diagnostics", "conclusions", "references", "versions", "artifacts",
    "results"
  )]
}

.pm_is_common_result <- function(x) {
  is.list(x) && (
    inherits(x, "psychomatic_result") ||
      all(c(
        "requested_settings", "effective_settings", "data_summary", "results",
        "warnings", "errors", "diagnostics", "decision_log", "conclusion"
      ) %in% names(x))
  )
}

.pm_is_supported_audit_result <- function(x) {
  .pm_is_common_result(x) || inherits(x, c("psychomatic_scores", "psychomatic_item_screen"))
}

.pm_legacy_audit_envelope <- function(x) {
  if (.pm_is_common_result(x)) return(x)
  if (inherits(x, "psychomatic_item_screen")) {
    return(.pm_result(
      class = "psychomatic_item_screen_audit_adapter",
      call = NULL,
      requested_settings = list(),
      effective_settings = list(),
      data_summary = list(n_items = nrow(x)),
      results = list(item_screen = as.data.frame(x)),
      warnings = NULL,
      errors = NULL,
      diagnostics = list(adapter = "psychomatic_item_screen"),
      decision_log = data.frame(),
      conclusion = list(status = "evaluable", reason = NA_character_)
    ))
  }
  if (inherits(x, "psychomatic_scores")) {
    return(.pm_result(
      class = "psychomatic_scores_audit_adapter",
      call = x$call,
      requested_settings = list(
        method = x$method,
        impute = x$impute,
        min_valid = x$min_valid,
        reversed_items = x$reversed_items
      ),
      effective_settings = list(
        method = x$method,
        impute = x$impute,
        scale_keys = names(x$keys)
      ),
      data_summary = list(
        n = if (is.data.frame(x$scores)) nrow(x$scores) else NA_integer_,
        n_scales = if (is.data.frame(x$scores)) ncol(x$scores) else NA_integer_
      ),
      results = list(
        scores = x$scores,
        reliability = x$reliability,
        scale_map = x$scale_map
      ),
      warnings = NULL,
      errors = NULL,
      diagnostics = list(adapter = "psychomatic_scores"),
      decision_log = data.frame(),
      conclusion = list(status = "evaluable", reason = NA_character_)
    ))
  }
  x
}

.pm_audit_objects <- function(result) {
  if (is.data.frame(result) && !inherits(result, "psychomatic_item_screen")) {
    stop(
      "A data.frame is a direct report table, not an auditable PsychoMatic result.",
      call. = FALSE
    )
  }

  single <- .pm_is_supported_audit_result(result)
  values <- if (single) list(result) else result
  if (!is.list(values) || length(values) == 0L) {
    stop("`result` must be a supported PsychoMatic result or a non-empty list of results.", call. = FALSE)
  }

  supplied_names <- if (single) "" else names(values)
  if (is.null(supplied_names)) supplied_names <- rep("", length(values))
  supplied_names[is.na(supplied_names)] <- ""
  unsupported <- which(!vapply(values, .pm_is_supported_audit_result, logical(1)))
  if (length(unsupported)) {
    i <- unsupported[[1L]]
    stop(
      sprintf(
        "Unsupported audit object at index %d (class: %s).",
        i,
        paste(class(values[[i]]), collapse = "/")
      ),
      call. = FALSE
    )
  }

  object_id <- sprintf("result_%03d", seq_along(values))
  generated <- !nzchar(trimws(supplied_names))
  display_label <- supplied_names
  display_label[generated] <- object_id[generated]
  primary_class <- vapply(values, function(x) class(x)[1L], character(1))
  all_classes <- vapply(values, function(x) paste(class(x), collapse = ", "), character(1))
  call_text <- vapply(values, function(x) {
    value <- if (is.list(x)) x[["call"]] else NULL
    if (is.null(value)) return(NA_character_)
    paste(deparse(value), collapse = " ")
  }, character(1L))
  metadata <- data.frame(
    object_index = seq_along(values),
    object_id = object_id,
    display_label = display_label,
    label_generated = generated,
    primary_class = primary_class,
    classes = all_classes,
    call_text = call_text,
    stringsAsFactors = FALSE
  )
  list(
    values = unname(lapply(values, .pm_legacy_audit_envelope)),
    metadata = metadata
  )
}

.pm_audit_inherited_language <- function(objects) {
  languages <- unlist(lapply(objects, function(x) {
    value <- x$effective_settings$language
    if (is.null(value)) value <- x$language
    if (is.null(value)) value <- x$requested_settings$language
    if (is.null(value) || length(value) != 1L || is.na(value)) return(NULL)
    tryCatch(.pm_normalize_language(value), error = function(e) NULL)
  }), use.names = FALSE)
  languages <- unique(languages)
  if (length(languages) == 1L) languages[[1L]] else NULL
}

.pm_audit_value_type <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.environment(x)) return("environment")
  if (inherits(x, "externalptr")) return("externalptr")
  if (inherits(x, "connection")) return("connection")
  paste(class(x), collapse = "/")
}

.pm_audit_unsafe_value <- function(x) {
  is.environment(x) || is.function(x) || isS4(x) ||
    typeof(x) == "externalptr" || inherits(x, "connection") ||
    inherits(x, c("lavaan", "SingleGroupClass", "mirt", "sirt", "psych", "fa"))
}

.pm_audit_contains_unsafe <- function(x, depth = 0L, max_depth = 12L) {
  if (.pm_audit_unsafe_value(x)) return(TRUE)
  attributes_x <- attributes(x)
  if (!is.null(attributes_x) && any(vapply(
    attributes_x,
    function(value) {
      .pm_audit_unsafe_value(value) ||
        (depth < max_depth && is.list(value) && .pm_audit_contains_unsafe(
          value,
          depth = depth + 1L,
          max_depth = max_depth
        ))
    },
    logical(1L)
  ))) {
    return(TRUE)
  }
  if (depth >= max_depth || !is.list(x)) {
    return(FALSE)
  }
  any(vapply(
    x,
    .pm_audit_contains_unsafe,
    logical(1L),
    depth = depth + 1L,
    max_depth = max_depth
  ))
}

.pm_audit_sanitize_value <- function(x, depth = 0L, max_depth = 12L) {
  if (.pm_audit_unsafe_value(x)) {
    return(paste0("<", .pm_audit_value_type(x), " omitted>"))
  }
  if (inherits(x, "call") || is.language(x)) {
    return(paste(deparse(x), collapse = " "))
  }
  if (depth >= max_depth && is.list(x)) {
    return(sprintf("<depth limit: %d elements omitted>", length(x)))
  }
  if (is.data.frame(x)) {
    columns <- lapply(x, function(column) {
      if (is.list(column)) {
        return(I(lapply(
          column,
          .pm_audit_sanitize_value,
          depth = depth + 1L,
          max_depth = max_depth
        )))
      }
      .pm_audit_sanitize_value(
        column,
        depth = depth + 1L,
        max_depth = max_depth
      )
    })
    names(columns) <- names(x)
    return(structure(
      columns,
      class = "data.frame",
      row.names = attr(x, "row.names")
    ))
  }
  if (is.list(x)) {
    out <- lapply(
      x,
      .pm_audit_sanitize_value,
      depth = depth + 1L,
      max_depth = max_depth
    )
    names(out) <- names(x)
    return(out)
  }
  out <- x
  attrs <- attributes(x)
  if (!is.null(attrs)) {
    keep <- !vapply(attrs, .pm_audit_contains_unsafe, logical(1L))
    attributes(out) <- attrs[keep]
  }
  out
}

.pm_audit_sanitize_table <- function(x) {
  if (!is.data.frame(x)) return(x)
  .pm_audit_sanitize_value(x)
}

.pm_audit_display <- function(x, max_values = 20L) {
  if (is.null(x)) return("<NULL>")
  if (.pm_audit_unsafe_value(x)) {
    return(paste0("<", .pm_audit_value_type(x), " omitted>"))
  }
  if (inherits(x, "call") || is.language(x)) return(paste(deparse(x), collapse = " "))
  if (length(x) == 0L) return("<empty>")
  if (is.list(x)) {
    values <- utils::head(x, max_values)
    rendered <- vapply(values, function(value) {
      if (is.null(value)) return("<NULL>")
      if (.pm_audit_unsafe_value(value)) {
        return(paste0("<", .pm_audit_value_type(value), " omitted>"))
      }
      if (is.atomic(value) && length(value) <= max_values) {
        return(paste(as.character(value), collapse = ", "))
      }
      sprintf("<%s: %d values>", .pm_audit_value_type(value), length(value))
    }, character(1))
    suffix <- if (length(x) > max_values) sprintf(" ... [%d entries]", length(x)) else ""
    return(paste0(paste(rendered, collapse = " | "), suffix))
  }
  values <- utils::head(x, max_values)
  rendered <- vapply(values, function(value) {
    if (length(value) == 0L) return("<empty>")
    if (is.numeric(value) && is.nan(value)) return("NaN")
    if (length(value) == 1L && is.na(value)) return("NA")
    if (is.numeric(value) && is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
    as.character(value)
  }, character(1))
  suffix <- if (length(x) > max_values) sprintf(" ... [%d values]", length(x)) else ""
  paste0(paste(rendered, collapse = ", "), suffix)
}

.pm_empty_audit_values <- function() {
  out <- data.frame(
    path = character(),
    value_display = character(),
    value_type = character(),
    truncated = logical(),
    stringsAsFactors = FALSE
  )
  out$value_raw <- I(vector("list", 0L))
  out
}

.pm_audit_flatten <- function(x,
                              max_nodes = 250L,
                              max_depth = 6L,
                              max_values = 20L) {
  rows <- list()
  node_count <- 0L
  limit_reported <- FALSE

  add_row <- function(path, value, truncated = FALSE, keep_raw = TRUE) {
    node_count <<- node_count + 1L
    unsafe <- .pm_audit_contains_unsafe(value)
    raw_value <- if (isTRUE(keep_raw) && !.pm_audit_unsafe_value(value)) {
      .pm_audit_sanitize_value(value)
    } else {
      NULL
    }
    row <- data.frame(
      path = if (nzchar(path)) path else "value",
      value_display = .pm_audit_display(value, max_values = max_values),
      value_type = .pm_audit_value_type(value),
      truncated = isTRUE(truncated) || unsafe,
      stringsAsFactors = FALSE
    )
    row$value_raw <- I(list(raw_value))
    rows[[length(rows) + 1L]] <<- row
  }

  walk <- function(value, path = "", depth = 0L) {
    if (node_count >= max_nodes) {
      if (!limit_reported) {
        add_row(if (nzchar(path)) path else "value", "<node limit reached>", TRUE)
        limit_reported <<- TRUE
      }
      return(invisible(NULL))
    }
    unsafe <- .pm_audit_unsafe_value(value)
    if (unsafe) {
      add_row(path, value, TRUE, keep_raw = FALSE)
      return(invisible(NULL))
    }
    if (depth >= max_depth && (is.list(value) || is.data.frame(value))) {
      add_row(path, sprintf("<depth limit: %d elements>", length(value)), TRUE)
      return(invisible(NULL))
    }
    if (is.data.frame(value)) {
      add_row(path, value, unsafe, keep_raw = TRUE)
      return(invisible(NULL))
    }
    if (is.list(value) && !inherits(value, "condition")) {
      if (length(value) == 0L) {
        if (nzchar(path)) add_row(path, value)
        return(invisible(NULL))
      }
      nms <- names(value)
      if (is.null(nms)) nms <- rep("", length(value))
      child_limit <- min(length(value), max_values)
      for (i in seq_len(child_limit)) {
        raw_label <- nms[[i]]
        if (is.na(raw_label) || !nzchar(raw_label)) {
          label <- paste0("[[", i, "]]" )
        } else {
          label <- .pm_audit_path_component(raw_label)
          same_name <- which(!is.na(nms) & nms == raw_label)
          if (length(same_name) > 1L) {
            label <- paste0(label, "[[", match(i, same_name), "]]" )
          }
        }
        child_path <- if (nzchar(path)) paste(path, label, sep = ".") else label
        walk(value[[i]], child_path, depth + 1L)
      }
      if (length(value) > child_limit) {
        add_row(
          if (nzchar(path)) paste0(path, ".<remaining>") else "<remaining>",
          sprintf("<%d entries omitted>", length(value) - child_limit),
          TRUE
        )
      }
      return(invisible(NULL))
    }
    add_row(
      path,
      value,
      truncated = length(value) > max_values,
      keep_raw = length(value) <= max_values && !is.language(value)
    )
    invisible(NULL)
  }

  walk(x)
  if (!length(rows)) return(.pm_empty_audit_values())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.pm_audit_path_component <- function(x) {
  x <- gsub("\\", "\\\\", as.character(x), fixed = TRUE)
  for (reserved in c(".", "[", "]", "<", ">")) {
    x <- gsub(reserved, paste0("\\", reserved), x, fixed = TRUE)
  }
  x
}

.pm_audit_bind <- function(rows, template) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(template)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.pm_audit_identity <- function(table, metadata, source = NULL) {
  if (!is.data.frame(table)) return(table)
  table <- .pm_audit_sanitize_table(table)
  prefix <- data.frame(
    object_index = rep(metadata$object_index, nrow(table)),
    object_id = rep(metadata$object_id, nrow(table)),
    display_label = rep(metadata$display_label, nrow(table)),
    stringsAsFactors = FALSE
  )
  if (!is.null(source)) prefix$source <- rep(source, nrow(table))
  cbind(prefix, table, stringsAsFactors = FALSE)
}

.pm_audit_settings <- function(object, metadata) {
  requested <- .pm_audit_flatten(object$requested_settings)
  effective <- .pm_audit_flatten(object$effective_settings)
  paths <- union(requested$path, effective$path)
  if (!length(paths)) {
    out <- data.frame(
      object_index = integer(), object_id = character(), display_label = character(),
      path = character(), requested = character(), effective = character(),
      requested_type = character(), effective_type = character(),
      stringsAsFactors = FALSE
    )
    out$requested_raw <- I(vector("list", 0L))
    out$effective_raw <- I(vector("list", 0L))
    return(out)
  }
  lookup <- function(table, path, field, unavailable) {
    hit <- match(path, table$path)
    if (field == "value_raw") {
      return(I(lapply(hit, function(i) if (is.na(i)) NULL else table$value_raw[[i]])))
    }
    value <- rep(unavailable, length(hit))
    available <- !is.na(hit)
    value[available] <- table[[field]][hit[available]]
    value
  }
  table <- data.frame(
    path = paths,
    requested = lookup(requested, paths, "value_display", NA_character_),
    effective = lookup(effective, paths, "value_display", NA_character_),
    requested_type = lookup(requested, paths, "value_type", NA_character_),
    effective_type = lookup(effective, paths, "value_type", NA_character_),
    stringsAsFactors = FALSE
  )
  table$requested_raw <- lookup(requested, paths, "value_raw", NULL)
  table$effective_raw <- lookup(effective, paths, "value_raw", NULL)
  .pm_audit_identity(table, metadata)
}

.pm_audit_conditions <- function(object, metadata, field) {
  table <- .pm_conditions_to_table(object[[field]], field)
  class_text <- if (nrow(table)) {
    vapply(table$class, paste, collapse = "/", character(1))
  } else {
    character()
  }
  if (nrow(table)) {
    call_text <- vapply(table$call, function(value) {
      if (is.null(value)) NA_character_ else paste(deparse(value), collapse = " ")
    }, character(1))
  } else {
    call_text <- character()
  }
  extras <- setdiff(names(table), c("code", "message", "class", "call", "condition"))
  safe <- data.frame(
    code = as.character(table$code),
    message = as.character(table$message),
    class_text = class_text,
    call_text = call_text,
    stringsAsFactors = FALSE
  )
  if (length(extras)) {
    extra_table <- table[extras]
    for (name in names(extra_table)) {
      if (is.list(extra_table[[name]])) {
        extra_table[[name]] <- vapply(
          extra_table[[name]], .pm_audit_display, character(1)
        )
      }
    }
    safe <- cbind(safe, extra_table, stringsAsFactors = FALSE)
  }
  .pm_audit_identity(safe, metadata, source = field)
}

.pm_audit_decisions <- function(object, metadata) {
  value <- object$decision_log
  if (is.null(value)) value <- data.frame()
  if (!is.data.frame(value)) {
    value <- .pm_audit_flatten(value)
  }
  .pm_audit_identity(value, metadata, source = "decision_log")
}

.pm_audit_conclusion <- function(object, metadata, language) {
  value <- object$conclusion
  status <- NULL
  reason <- NA_character_
  conclusion_text <- NA_character_
  if (is.list(value)) {
    if (!is.null(value$status) && length(value$status)) status <- as.character(value$status[[1L]])
    if (!is.null(value$reason) && length(value$reason)) reason <- as.character(value$reason[[1L]])
    if (!is.null(value$text) && length(value$text)) conclusion_text <- as.character(value$text[[1L]])
  } else if (!is.null(value) && length(value)) {
    conclusion_text <- as.character(value[[1L]])
  }
  if (is.null(status) && !is.null(object$conclusion_status) && length(object$conclusion_status)) {
    status <- as.character(object$conclusion_status[[1L]])
  }
  if (is.null(status) && is.list(object$diagnostics) &&
      !is.null(object$diagnostics$status) && length(object$diagnostics$status)) {
    status <- as.character(object$diagnostics$status[[1L]])
  }
  if (is.null(status)) {
    not_evaluable_text <- !is.na(conclusion_text) && grepl(
      "not[_ ]evaluable|no (fue )?evaluable",
      conclusion_text,
      ignore.case = TRUE
    )
    status <- if (isTRUE(not_evaluable_text)) "not_evaluable" else "evaluable"
  }
  if (is.na(reason) && identical(status, "not_evaluable") && !is.na(conclusion_text)) {
    reason <- conclusion_text
  }
  data.frame(
    object_index = metadata$object_index,
    object_id = metadata$object_id,
    display_label = metadata$display_label,
    status = status,
    conclusion_status = status,
    reason = reason,
    conclusion = conclusion_text,
    language = language,
    stringsAsFactors = FALSE
  )
}

.pm_audit_references <- function(object, metadata) {
  value <- object$references
  if (is.null(value) || length(value) == 0L) {
    return(data.frame(
      object_index = integer(), object_id = character(), display_label = character(),
      reference = character(), identifier = character(), stringsAsFactors = FALSE
    ))
  }
  if (is.data.frame(value)) {
    table <- value
  } else {
    table <- data.frame(reference = as.character(unlist(value, use.names = FALSE)), stringsAsFactors = FALSE)
  }
  if (!"reference" %in% names(table)) table$reference <- apply(table, 1L, paste, collapse = " | ")
  if (!"identifier" %in% names(table)) table$identifier <- NA_character_
  .pm_audit_identity(table, metadata, source = "references")
}

.pm_audit_versions <- function(object, metadata) {
  value <- object$versions
  if (is.null(value)) value <- object$package_versions
  if (is.null(value) || length(value) == 0L) {
    return(data.frame(
      object_index = integer(), object_id = character(), display_label = character(),
      package = character(), version = character(), stringsAsFactors = FALSE
    ))
  }
  if (is.data.frame(value)) {
    table <- value
    if (!"package" %in% names(table)) table$package <- rownames(table)
    if (!"version" %in% names(table)) table$version <- NA_character_
  } else {
    nms <- names(value)
    if (is.null(nms)) nms <- rep("unknown", length(value))
    table <- data.frame(
      package = nms,
      version = vapply(value, function(x) paste(as.character(x), collapse = ", "), character(1)),
      stringsAsFactors = FALSE
    )
  }
  .pm_audit_identity(table, metadata, source = "versions")
}

.pm_audit_artifacts <- function(object, metadata) {
  fields <- intersect(
    c("report_path", "report_paths", "export_path", "export_paths", "artifact_path"),
    names(object)
  )
  rows <- lapply(fields, function(field) {
    values <- object[[field]]
    if (is.null(values) || !length(values)) return(NULL)
    data.frame(
      artifact = field,
      path = as.character(unlist(values, use.names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      object_index = integer(), object_id = character(), display_label = character(),
      source = character(), artifact = character(), path = character(),
      stringsAsFactors = FALSE
    ))
  }
  .pm_audit_identity(do.call(rbind, rows), metadata, source = "artifacts")
}

.pm_audit_result_tables <- function(object, object_id) {
  tables <- list()
  seen <- character()
  heavy_names <- c(
    "fit", "fits", "fitted_models", "efa", "alignment", "constraints",
    "configural", "raw_coefficients", "group_parameters",
    "aligned_item_parameters"
  )
  safe_key <- function(path) {
    path <- gsub("[^[:alnum:]_]+", "_", path)
    path <- gsub("_+", "_", path)
    path <- sub("^_+", "", path)
    path <- sub("_+$", "", path)
    if (!nzchar(path)) "table" else path
  }
  add_table <- function(value, semantic_path, output_path = semantic_path) {
    if (semantic_path %in% seen) return(invisible(FALSE))
    table <- .pm_audit_tabular_value(value)
    if (is.null(table)) return(invisible(FALSE))
    key_base <- paste(object_id, safe_key(output_path), sep = "_")
    key <- key_base
    suffix <- 1L
    while (key %in% names(tables)) {
      suffix <- suffix + 1L
      key <- paste0(key_base, "_", suffix)
    }
    tables[[key]] <<- .pm_audit_sanitize_table(table)
    seen <<- c(seen, semantic_path)
    invisible(TRUE)
  }
  walk_tables <- function(value,
                          semantic_path,
                          output_path = semantic_path,
                          depth = 0L,
                          include_atomic = FALSE) {
    leaf_name <- tail(strsplit(semantic_path, "\\.", fixed = FALSE)[[1L]], 1L)
    if (.pm_audit_unsafe_value(value)) {
      return(invisible(NULL))
    }
    if (is.data.frame(value) || is.matrix(value) ||
        (is.array(value) && length(dim(value)) == 2L)) {
      add_table(value, semantic_path, output_path)
      return(invisible(NULL))
    }
    if (tolower(leaf_name) %in% heavy_names) {
      return(invisible(NULL))
    }
    if (is.atomic(value) && is.null(dim(value))) {
      if (isTRUE(include_atomic) && length(value)) {
        add_table(value, semantic_path, output_path)
      }
      return(invisible(NULL))
    }
    if (!is.list(value) || depth >= 6L || !length(value)) {
      return(invisible(NULL))
    }
    labels <- names(value)
    if (is.null(labels)) labels <- rep("", length(value))
    for (i in seq_along(value)) {
      label <- labels[[i]]
      if (is.na(label) || !nzchar(label)) label <- paste0("item_", i)
      child_semantic <- if (nzchar(semantic_path)) {
        paste(semantic_path, label, sep = ".")
      } else {
        label
      }
      child_output <- if (nzchar(output_path)) {
        paste(output_path, label, sep = "_")
      } else {
        label
      }
      walk_tables(
        value[[i]],
        child_semantic,
        child_output,
        depth = depth + 1L,
        include_atomic = include_atomic
      )
    }
    invisible(NULL)
  }

  canonical <- object$tables
  canonical_names <- character()
  if (is.list(canonical)) {
    canonical_names <- names(canonical)
    if (is.null(canonical_names)) canonical_names <- rep("", length(canonical))
    for (i in seq_along(canonical)) {
      label <- canonical_names[[i]]
      if (is.na(label) || !nzchar(label)) label <- paste0("table_", i)
      walk_tables(canonical[[i]], label, label, include_atomic = TRUE)
    }
  }

  results <- object$results
  if (is.data.frame(results)) {
    add_table(results, "results", "results")
    results <- list()
  }
  if (!is.list(results)) results <- list()
  result_names <- names(results)
  if (is.null(result_names)) result_names <- rep("", length(results))
  canonical_alias <- c(
    fit = "fit", fit_table = "fit", fit_indices = "fit",
    loadings = "loadings", tidy_loadings = "loadings",
    reliability = "reliability", reliability_table = "reliability",
    decisions = "decisions", constraint_families = "constraints",
    freed_constraints = "freed_constraints", anchors = "anchors",
    comparison = "comparisons"
  )
  for (i in seq_along(results)) {
    label <- result_names[[i]]
    if (is.na(label) || !nzchar(label)) label <- paste0("results_", i)
    alias <- unname(canonical_alias[label])
    if (length(alias) && !is.na(alias) && alias %in% canonical_names) next
    walk_tables(results[[i]], label, label, include_atomic = FALSE)
  }

  diagnostics <- object$diagnostics
  if (is.list(diagnostics)) {
    diagnostic_names <- names(diagnostics)
    if (is.null(diagnostic_names)) diagnostic_names <- rep("", length(diagnostics))
    for (i in seq_along(diagnostics)) {
      label <- diagnostic_names[[i]]
      if (is.na(label) || !nzchar(label)) label <- paste0("diagnostic_", i)
      walk_tables(
        diagnostics[[i]],
        paste0("diagnostics.", label),
        paste0("diagnostics_", label),
        include_atomic = FALSE
      )
    }
  }

  direct_fields <- c(
    "higher_order_loadings", "modification_indices",
    "respecification_suggestions", "factor_correlation_table"
  )
  for (field in intersect(direct_fields, names(object))) {
    if (field %in% result_names || field %in% canonical_names) next
    walk_tables(object[[field]], field, field, include_atomic = FALSE)
  }

  is_direct_tabular <- vapply(results, function(value) {
    is.data.frame(value) || is.matrix(value) ||
      (is.array(value) && length(dim(value)) == 2L)
  }, logical(1L))
  remaining <- results[!is_direct_tabular]
  if (length(remaining)) {
    remaining_names <- names(remaining)
    if (is.null(remaining_names)) remaining_names <- rep("", length(remaining))
    heavy <- tolower(remaining_names) %in% c(
      "fit", "fits", "fitted_models", "efa", "alignment", "constraints",
      "configural", "raw_coefficients", "group_parameters",
      "aligned_item_parameters", "item_results"
    )
    values <- .pm_audit_flatten(remaining[!heavy])
    if (any(heavy)) {
      omitted <- lapply(which(heavy), function(i) {
        row <- data.frame(
          path = if (nzchar(remaining_names[[i]])) remaining_names[[i]] else paste0("[[", i, "]]"),
          value_display = paste0("<", .pm_audit_value_type(remaining[[i]]), " omitted>"),
          value_type = .pm_audit_value_type(remaining[[i]]),
          truncated = TRUE,
          stringsAsFactors = FALSE
        )
        row$value_raw <- I(list(NULL))
        row
      })
      omitted <- do.call(rbind, omitted)
      values <- rbind(values, omitted)
    }
    if (nrow(values)) tables[[paste0(object_id, "_values")]] <- values
  }
  tables
}

.pm_audit_tabular_value <- function(value) {
  if (.pm_audit_unsafe_value(value)) return(NULL)
  if (is.data.frame(value)) return(.pm_audit_sanitize_table(value))
  if (is.matrix(value) || (is.array(value) && length(dim(value)) == 2L)) {
    matrix_value <- as.matrix(value)
    row_labels <- rownames(matrix_value)
    if (is.null(row_labels)) row_labels <- as.character(seq_len(nrow(matrix_value)))
    table <- data.frame(
      row_index = seq_len(nrow(matrix_value)),
      row_name = row_labels,
      as.data.frame(matrix_value, check.names = FALSE, optional = TRUE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    return(.pm_audit_sanitize_table(table))
  }
  if (is.atomic(value) && is.null(dim(value)) && length(value)) {
    value_names <- names(value)
    if (is.null(value_names)) value_names <- rep("", length(value))
    return(.pm_audit_sanitize_table(data.frame(
      index = seq_along(value),
      name = value_names,
      value = unname(value),
      stringsAsFactors = FALSE
    )))
  }
  NULL
}

.pm_build_audit_sections <- function(objects, metadata, language) {
  settings <- data_rows <- warning_rows <- error_rows <- decision_rows <- list()
  diagnostic_rows <- conclusion_rows <- reference_rows <- version_rows <- artifact_rows <- list()
  result_tables <- list()

  for (i in seq_along(objects)) {
    object <- objects[[i]]
    meta <- metadata[i, , drop = FALSE]
    settings[[i]] <- .pm_audit_settings(object, meta)
    data_rows[[i]] <- .pm_audit_identity(
      .pm_audit_flatten(object$data_summary), meta, source = "data_summary"
    )
    warning_rows[[i]] <- .pm_audit_conditions(object, meta, "warnings")
    error_rows[[i]] <- .pm_audit_conditions(object, meta, "errors")
    decision_rows[[i]] <- .pm_audit_decisions(object, meta)
    diagnostic_rows[[i]] <- .pm_audit_identity(
      .pm_audit_flatten(object$diagnostics), meta, source = "diagnostics"
    )
    conclusion_rows[[i]] <- .pm_audit_conclusion(object, meta, language)
    reference_rows[[i]] <- .pm_audit_references(object, meta)
    version_rows[[i]] <- .pm_audit_versions(object, meta)
    artifact_rows[[i]] <- .pm_audit_artifacts(object, meta)
    result_tables <- c(result_tables, .pm_audit_result_tables(object, meta$object_id))
  }

  empty_values <- .pm_audit_identity(
    .pm_empty_audit_values(), metadata[1L, , drop = FALSE], "data_summary"
  )
  empty_settings <- .pm_audit_settings(
    list(requested_settings = list(), effective_settings = list()),
    metadata[1L, , drop = FALSE]
  )
  empty_conditions <- .pm_audit_conditions(
    list(warnings = NULL), metadata[1L, , drop = FALSE], "warnings"
  )[0L, , drop = FALSE]
  empty_decisions <- .pm_audit_identity(data.frame(
    code = character(), requested = character(), effective = character(),
    reason = character(), severity = character(), stringsAsFactors = FALSE
  ), metadata[1L, , drop = FALSE], "decision_log")
  empty_references <- .pm_audit_references(
    list(references = NULL), metadata[1L, , drop = FALSE]
  )
  empty_versions <- .pm_audit_versions(
    list(versions = NULL), metadata[1L, , drop = FALSE]
  )
  empty_artifacts <- .pm_audit_artifacts(
    list(), metadata[1L, , drop = FALSE]
  )
  list(
    objects = metadata,
    settings = .pm_audit_bind(settings, empty_settings),
    data = .pm_audit_bind(data_rows, empty_values),
    warnings = .pm_audit_bind(warning_rows, empty_conditions),
    errors = .pm_audit_bind(error_rows, transform(empty_conditions, source = character())),
    decisions = .pm_audit_bind(decision_rows, empty_decisions),
    diagnostics = .pm_audit_bind(diagnostic_rows, empty_values),
    conclusions = .pm_audit_bind(conclusion_rows, data.frame()),
    references = .pm_audit_bind(reference_rows, empty_references),
    versions = .pm_audit_bind(version_rows, empty_versions),
    artifacts = .pm_audit_bind(artifact_rows, empty_artifacts),
    results = result_tables
  )
}
