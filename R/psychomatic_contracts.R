.pm_normalize_language <- function(language = NULL, inherited = NULL) {
  value <- language
  if (is.null(value)) value <- inherited
  if (is.null(value)) value <- getOption("PsychoMatic.language", "eng")

  value <- as.character(value)
  if (length(value) == 0L || is.na(value[1L])) {
    stop(
      "`language` must identify Spanish ('esp') or English ('eng').",
      call. = FALSE
    )
  }

  value <- tolower(trimws(value[1L]))
  if (value %in% c("esp", "es", "spa", "spanish", "espanol")) return("esp")
  if (value %in% c("eng", "en", "english", "ingles")) return("eng")
  stop(
    "`language` must identify Spanish ('esp') or English ('eng').",
    call. = FALSE
  )
}

.pm_decision <- function(code, requested, effective, reason, severity = "info") {
  fields <- list(
    code = as.character(code),
    requested = as.character(requested),
    effective = as.character(effective),
    reason = as.character(reason),
    severity = as.character(severity)
  )
  for (field in names(fields)) {
    if (length(fields[[field]]) != 1L) {
      stop(sprintf("`%s` must be a single value.", field), call. = FALSE)
    }
  }

  data.frame(
    code = fields$code,
    requested = fields$requested,
    effective = fields$effective,
    reason = fields$reason,
    severity = fields$severity,
    stringsAsFactors = FALSE
  )
}

.pm_not_evaluable <- function(reason) {
  list(status = "not_evaluable", reason = reason)
}

.pm_validate_scalar_number <- function(x,
                                       name = deparse(substitute(x)),
                                       lower = -Inf,
                                       upper = Inf,
                                       integer = FALSE) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.numeric(lower) || length(lower) != 1L || is.na(lower) ||
      !is.numeric(upper) || length(upper) != 1L || is.na(upper) ||
      lower > upper) {
    stop("Validation bounds must be ordered numeric scalars.", call. = FALSE)
  }
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    stop("`integer` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(
      sprintf("`%s` must be a single finite numeric value.", name),
      call. = FALSE
    )
  }
  if (isTRUE(integer) && x != trunc(x)) {
    stop(sprintf("`%s` must be an integer.", name), call. = FALSE)
  }
  if (x < lower || x > upper) {
    stop(
      sprintf("`%s` must be between %s and %s.", name, lower, upper),
      call. = FALSE
    )
  }

  x
}

.pm_capture_conditions <- function(expr) {
  captured_warnings <- list()
  captured_messages <- list()
  captured_error <- NULL

  value <- tryCatch(
    withCallingHandlers(
      eval.parent(substitute(expr)),
      warning = function(condition) {
        captured_warnings[[length(captured_warnings) + 1L]] <<- condition
        invokeRestart("muffleWarning")
      },
      message = function(condition) {
        captured_messages[[length(captured_messages) + 1L]] <<- condition
        invokeRestart("muffleMessage")
      }
    ),
    error = function(condition) {
      captured_error <<- condition
      NULL
    }
  )

  list(
    value = value,
    warnings = captured_warnings,
    messages = captured_messages,
    error = captured_error
  )
}

.pm_assert_unique_fields <- function(x, name) {
  field_names <- names(x)
  if (is.null(field_names)) return(invisible(NULL))

  named_fields <- field_names[!is.na(field_names) & nzchar(field_names)]
  if (anyDuplicated(named_fields)) {
    stop(sprintf("`%s` contains duplicate field names.", name), call. = FALSE)
  }

  invisible(NULL)
}

.pm_empty_condition_table <- function() {
  out <- data.frame(
    code = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
  out$class <- I(vector("list", 0L))
  out$call <- I(vector("list", 0L))
  out$condition <- I(vector("list", 0L))
  out
}

.pm_complete_condition_table <- function(x, kind) {
  .pm_assert_unique_fields(x, kind)

  if (!"code" %in% names(x)) {
    x$code <- rep(kind, nrow(x))
  }
  if (!"message" %in% names(x)) {
    if (nrow(x) > 0L) {
      stop(sprintf("`%s` must contain a `message` field.", kind), call. = FALSE)
    }
    x$message <- character()
  }
  if (!"class" %in% names(x)) {
    x$class <- I(rep(list(character()), nrow(x)))
  } else if (!is.list(x$class)) {
    x$class <- I(lapply(x$class, as.character))
  }
  if (!"call" %in% names(x)) {
    x$call <- I(rep(list(NULL), nrow(x)))
  } else if (!is.list(x$call)) {
    x$call <- I(as.list(x$call))
  }
  if (!"condition" %in% names(x)) {
    x$condition <- I(rep(list(NULL), nrow(x)))
  } else if (!is.list(x$condition)) {
    x$condition <- I(as.list(x$condition))
  }

  required <- c("code", "message", "class", "call", "condition")
  x[c(required, setdiff(names(x), required))]
}

.pm_conditions_to_table <- function(x, kind) {
  if (is.null(x) || (is.list(x) && !inherits(x, "condition") && length(x) == 0L)) {
    return(.pm_empty_condition_table())
  }
  if (is.data.frame(x)) {
    return(.pm_complete_condition_table(x, kind))
  }
  if (inherits(x, "condition")) {
    x <- list(x)
  }
  if (is.character(x)) {
    out <- data.frame(
      code = rep(kind, length(x)),
      message = x,
      stringsAsFactors = FALSE
    )
    out$class <- I(rep(list(character()), length(x)))
    out$call <- I(rep(list(NULL), length(x)))
    out$condition <- I(rep(list(NULL), length(x)))
    return(out)
  }
  if (!is.list(x) || !all(vapply(x, inherits, logical(1), what = "condition"))) {
    stop(
      sprintf("`%s` must be a data frame, condition, or list of conditions.", kind),
      call. = FALSE
    )
  }

  codes <- vapply(
    x,
    function(condition) {
      code <- condition[["code"]]
      if (is.null(code) || length(code) == 0L || is.na(code[1L])) {
        class(condition)[1L]
      } else {
        as.character(code[1L])
      }
    },
    character(1)
  )
  out <- data.frame(
    code = codes,
    message = vapply(x, conditionMessage, character(1)),
    stringsAsFactors = FALSE
  )
  out$class <- I(lapply(x, class))
  out$call <- I(lapply(x, conditionCall))
  out$condition <- I(unname(x))
  out
}

.pm_result <- function(class,
                       call,
                       requested_settings,
                       effective_settings,
                       data_summary,
                       results,
                       warnings,
                       errors,
                       diagnostics,
                       decision_log,
                       conclusion) {
  requested_class <- as.character(class)
  if (length(requested_class) == 0L || anyNA(requested_class) ||
      any(!nzchar(requested_class))) {
    stop("`class` must contain at least one non-empty class name.", call. = FALSE)
  }

  fields <- list(
    requested_settings = requested_settings,
    effective_settings = effective_settings,
    data_summary = data_summary,
    results = results,
    diagnostics = diagnostics,
    decision_log = decision_log
  )
  for (field_name in names(fields)) {
    .pm_assert_unique_fields(fields[[field_name]], field_name)
  }

  out <- list(
    call = call,
    requested_settings = requested_settings,
    effective_settings = effective_settings,
    data_summary = data_summary,
    results = results,
    warnings = .pm_conditions_to_table(warnings, "warnings"),
    errors = .pm_conditions_to_table(errors, "errors"),
    decision_log = decision_log,
    diagnostics = diagnostics,
    conclusion = conclusion
  )
  class(out) <- unique(c(requested_class, "psychomatic_result"))
  out
}
