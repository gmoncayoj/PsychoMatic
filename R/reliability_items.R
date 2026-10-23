#' Internal reliability calculations
#' @keywords internal
.pm_auto_omega <- function(x, level, keying = "none") {
  n <- sum(stats::complete.cases(x))
  diagnostics <- list(n = n, missing = "listwise", warnings = character(), messages = character())
  finish <- function(value = NA_real_, reason = NA_character_) list(
    row = .pm_reliability_row("omega", value, "observed_score_congeneric",
      if (level == "ordinal") "ordinal_observed_score" else "raw_covariance",
      keying, domain = c(0, 1), reason = reason), diagnostics = diagnostics)
  if (ncol(x) < 3L) return(finish(reason = "insufficient_items"))
  if (level == "mixed") return(finish(reason = "mixed_measurement_levels"))
  if (!.pm_semtools_available()) return(finish(reason = "semtools_not_available"))
  if (any(vapply(x, function(z) any(!is.finite(z) & !is.na(z)), logical(1))))
    return(finish(reason = "non_finite_input"))
  x <- x[stats::complete.cases(x), , drop = FALSE]
  if (n < 3L) return(finish(reason = "insufficient_complete_cases"))
  if (any(vapply(x, function(z) length(unique(z)) < 2L, logical(1))))
    return(finish(reason = "constant_item"))
  names(x) <- paste0("item", seq_along(x))
  fit_capture <- .pm_capture_conditions(lavaan::cfa(
    paste("f =~", paste(names(x), collapse = "+")), data = x, std.lv = TRUE,
    ordered = if (level == "ordinal") names(x) else NULL,
    estimator = if (level == "ordinal") "WLSMV" else "ML", missing = "listwise"))
  diagnostics$warnings <- fit_capture$warnings
  diagnostics$messages <- fit_capture$messages
  diagnostics$error <- fit_capture$error
  if (!is.null(fit_capture$error)) return(finish(reason = "model_fit_error"))
  fit <- fit_capture$value
  diagnostics$fit <- fit
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) return(finish(reason = "model_not_converged"))
  checked <- .pm_capture_conditions(lavaan::lavInspect(fit, "post.check"))
  diagnostics$warnings <- c(diagnostics$warnings, checked$warnings)
  if (!isTRUE(checked$value)) return(finish(reason = "inadmissible_model"))
  cr <- .pm_capture_conditions(semTools::compRelSEM(fit, obs.var = TRUE,
    tau.eq = FALSE, ord.scale = TRUE, return.df = FALSE))
  diagnostics$warnings <- c(diagnostics$warnings, cr$warnings)
  diagnostics$messages <- c(diagnostics$messages, cr$messages)
  diagnostics$error <- cr$error
  if (!is.null(cr$error)) return(finish(reason = "semtools_error"))
  value <- as.numeric(unlist(cr$value))
  if (length(value) != 1L) return(finish(reason = "semtools_output_unrecognized"))
  finish(value)
}

.pm_reliability_items <- function(x, level, covariance, correlation, model = NULL) {
  diagnostics <- list()
  rows <- lapply(seq_along(x), function(i) {
    rest <- x[-i]
    complete <- stats::complete.cases(x)
    r <- if (ncol(rest) && sum(complete) > 1L) suppressWarnings(
      stats::cor(x[[i]][complete], rowSums(rest[complete, , drop = FALSE]))) else NA_real_
    if (!is.finite(r)) r <- NA_real_
    alpha <- if (is.matrix(covariance) && ncol(rest) >= 2L)
      .pm_alpha_from_matrix(covariance[-i, -i, drop = FALSE]) else
        list(value = NA_real_, reason = "insufficient_items")
    ordinal <- if (level == "ordinal" && is.matrix(correlation) && ncol(rest) >= 2L)
      .pm_alpha_from_matrix(correlation[-i, -i, drop = FALSE]) else
        list(value = NA_real_, reason = if (level != "ordinal") "not_ordinal" else
          if (ncol(rest) < 2L) "insufficient_items" else "polychoric_not_evaluable")
    alpha_checked <- .pm_reliability_row("alpha", alpha$value, "observed_score",
      "raw_covariance", "none", domain = c(-Inf, 1), reason = alpha$reason)
    ordinal_checked <- .pm_reliability_row("alpha_ordinal", ordinal$value, "latent_response",
      "polychoric", "none", domain = c(-Inf, 1), reason = ordinal$reason)
    omega <- if (is.null(model)) .pm_auto_omega(rest, level) else list(
      row = .pm_reliability_row("omega", NA_real_, "model_based", "model_implied", "none",
        reason = "supplied_model_requires_explicit_refit"), diagnostics = list(n = sum(stats::complete.cases(rest))))
    diagnostics[[names(x)[i]]] <<- omega$diagnostics
    data.frame(item = names(x)[i], n = sum(complete), item_rest_r = r,
      item_rest_reason = if (is.finite(r)) NA_character_ else "insufficient_or_constant_scores",
      alpha_if_deleted = alpha_checked$estimate, alpha_if_deleted_reason = alpha_checked$reason,
      alpha_ordinal_if_deleted = ordinal_checked$estimate, alpha_ordinal_if_deleted_reason = ordinal_checked$reason,
      omega_if_deleted = omega$row$estimate, omega_if_deleted_reason = omega$row$reason,
      omega_if_deleted_n = omega$diagnostics$n, stringsAsFactors = FALSE)
  })
  list(rows = do.call(rbind, rows), diagnostics = diagnostics)
}

.pm_cfa_alpha_rows <- function(fit) {
  ordinal <- .pm_capture_conditions(semTools::compRelSEM(fit, tau.eq = TRUE,
    obs.var = TRUE, ord.scale = FALSE))
  groups <- lavaan::lavInspect(fit, "group.label")
  ordered <- lavaan::lavNames(fit, "ov.ord")
  raw <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL)
  if (is.matrix(raw)) raw <- list(raw)
  pt <- lavaan::parameterTable(fit)
  observed_names <- lavaan::lavNames(fit, "ov")
  rows <- .pm_semtools_long(NULL, "alpha_observed", groups, "L1")
  if (length(raw)) for (g in seq_along(raw)) {
    if (!is.matrix(raw[[g]]) || !nrow(raw[[g]])) next
    factors <- unique(pt$lhs[pt$op == "=~" & pt$rhs %in% observed_names & pt$group == g])
    for (f in factors) {
      items <- unique(pt$rhs[pt$lhs == f & pt$op == "=~" & pt$rhs %in% observed_names & pt$group == g])
      a <- if (length(items) >= 2L) .pm_alpha_from_matrix(
        stats::cov(raw[[g]][, items, drop = FALSE], use = "pairwise.complete.obs"))$value else NA_real_
      rows <- rbind(rows, data.frame(group = if (length(groups) >= g) groups[g] else "group1",
        level = "L1", factor = f, coefficient = "alpha_observed", value = a))
    }
  }
  if (!length(ordered) && !nrow(rows)) {
    rows <- .pm_semtools_long(ordinal$value, "alpha_observed", groups, "L1")
  }
  if (length(ordered)) {
    factors <- unique(pt$lhs[pt$op == "=~" & pt$rhs %in% ordered])
    factors <- factors[vapply(factors, function(f) {
      indicators <- pt$rhs[pt$op == "=~" & pt$lhs == f & pt$rhs %in% observed_names]
      length(indicators) > 1L && all(indicators %in% ordered)
    }, logical(1))]
    extra <- .pm_semtools_long(ordinal$value, "alpha_ordinal", groups, "L1")
    rows <- rbind(rows, extra[extra$factor %in% factors, , drop = FALSE])
  }
  list(rows = rows, warnings = ordinal$warnings,
    messages = ordinal$messages, errors = ordinal$error)
}
