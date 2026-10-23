.pm_network_wording <- function(raw, keyed, ega, settings) {
  if(settings$corr=="polychoric" && any(vapply(raw,function(x)length(unique(x[!is.na(x)])),integer(1))>7L))
    return(.pm_not_evaluable("The installed riEGA only passes its default seven-category threshold to the initial correlation estimator; riEGA was skipped to avoid mixing correlation methods for items with more than seven categories."))
  original <- .pm_network_ega(raw, settings)
  adjusted <- .pm_network_ega(raw, settings, "riEGA")
  validation <- .pm_network_validate_ri(adjusted)
  # Evaluate both partitions on one ORIGINAL correlation matrix. Comparing the
  # native TEFI of RI residuals to the keyed/raw TEFI would change the data.
  original_fit <- as.numeric(original$TEFI)
  evaluable <- .pm_network_partition_evaluable(adjusted$EGA)
  adjusted_partition_fit <- if (evaluable)
    as.numeric(EGAnet::tefi(original$correlation, adjusted$EGA$wc)$VN.Entropy.Fit) else NA_real_
  signs <- original$correlation
  negative <- data.frame(item = names(raw), negative_links = colSums(signs < 0, na.rm = TRUE),
    declared_reverse = names(raw) %in% settings$reverse)
  list(original = original, adjusted = adjusted, validation = validation, candidates = negative,
    status = if (evaluable) "completed" else "not_evaluable",
    reason = if (evaluable) "" else
      "The riEGA residual network has no estimable communities or unassigned items; its partition TEFI and structural stability are not evaluable. This does not prove wording effects.",
    comparison = data.frame(model = c("raw_EGA", "riEGA"),
      dimensions = c(original$n.dim, adjusted$EGA$n.dim),
      status = c("completed", if (evaluable) "completed" else "not_evaluable"),
      same_partition = c(TRUE, if (evaluable) identical(outer(original$wc, original$wc, "=="),
        outer(adjusted$EGA$wc, adjusted$EGA$wc, "==")) else NA)),
    tefi = data.frame(model = c("raw_EGA", "riEGA_partition", "riEGA_residual"),
      index = "TEFI", value = c(original_fit, adjusted_partition_fit, if(evaluable) as.numeric(adjusted$TEFI) else NA_real_),
      comparison_set = c("raw_retained", "raw_retained", "RI_residual_retained"),
      comparable = c(TRUE, evaluable, FALSE)),
    conclusion = "RI fit and negative associations are diagnostic evidence; they do not identify semantic reverse keys.")
}

.pm_network_hierarchy <- function(data, ega, spec, settings) {
  native <- tryCatch(do.call(EGAnet::hierEGA,c(list(data=as.matrix(data),
    na.data = settings$missing,
    seed = settings$seed, loading.method = "revised", scores = "network", impute = "none", plot.EGA = FALSE),
    .pm_network_native_options(settings,hierarchy=TRUE),
    .pm_network_corr_args(settings$corr,data,settings$ordinal_categories))),error=function(e) {
      if(identical(settings$model,"TMFG") && grepl("invalid 'nrow' value",conditionMessage(e),fixed=TRUE))
        stop("TMFG requires at least four nodes in an estimated network order. hierEGA could not estimate an order with fewer nodes; consider GLASSO or hierarchy=FALSE. No estimator was substituted.",call.=FALSE)
      stop(e)
    })
  groups <- split(seq_len(ncol(data)), ega$wc)
  assignments <- native$dim.variables
  if(!is.null(assignments))names(assignments) <- c("item","lower_dimension","higher_dimension")
  result <- list(native = native, candidates = data.frame(), fits = list(),
    network_structure=.pm_network_hierarchy_table(native),
    network_assignments=assignments,
    network_interpretation=attr(native,"methods")$interpretation,
    recommendation = "indeterminate", reason = "Requires at least three dimensions with three indicators each.",
    tefi = data.frame(model = "hierEGA", index = "genTEFI", value = as.numeric(native$TEFI),
      comparison_set = "hierarchical_entropy", comparable = FALSE))
  if (length(groups) < 3L || any(lengths(groups) < 3L) || anyNA(ega$wc)) return(result)
  # Safe internal names avoid embedding user item labels in lavaan syntax.
  d <- as.data.frame(data); names(d) <- paste0("x", seq_len(ncol(d)))
  factors <- paste0("F", seq_along(groups))
  first <- vapply(seq_along(groups), function(j)
    paste(factors[j], "=~", paste(names(d)[groups[[j]]], collapse = " + ")), character(1))
  syntax <- list(correlated = paste(first, collapse = "\n"),
    second_order = paste(c(first, paste("G =~", paste(factors, collapse = " + "))), collapse = "\n"),
    bifactor = paste(c(paste("G =~", paste(names(d), collapse = " + ")), first), collapse = "\n"))
  encoding <- .pm_network_cfa_encoding(names(data), spec)
  ordered <- encoding$ordered
  estimator <- encoding$estimator
  for (name in names(syntax)) {
    capture <- .pm_capture_conditions(lavaan::cfa(syntax[[name]], data = d, ordered = ordered,
      estimator = estimator, std.lv = TRUE, orthogonal = name == "bifactor",
      missing = "listwise"))
    fit <- capture$value
    row <- data.frame(model = name, estimator = estimator, n = sum(stats::complete.cases(d)),
      df=NA_real_,higher_order_df=if(name=="second_order")length(groups)*(length(groups)-3)/2 else NA_real_,
      converged = FALSE, admissible = FALSE, CFI = NA_real_, TLI = NA_real_,
      RMSEA = NA_real_, SRMR = NA_real_, AIC = NA_real_, BIC = NA_real_,
      reason = if (is.null(capture$error)) "" else conditionMessage(capture$error))
    if (!is.null(fit)) {
      result$fits[[name]] <- fit
      row$converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
      if (row$converged) {
        admissibility <- .pm_cfa_admissibility(fit, lavaan::standardizedSolution(fit), converged = TRUE)
        row$admissible <- isTRUE(admissibility$admissible)
        row$reason <- paste(admissibility$reason_codes, collapse = "; ")
        measures <- lavaan::fitMeasures(fit)
        row$df <- unname(measures["df"])
        for (metric in c("CFI", "TLI", "RMSEA", "SRMR", "AIC", "BIC")) {
          key <- tolower(metric)
          robust <- paste0(key, ".robust")
          if (robust %in% names(measures) && is.finite(measures[robust])) key <- robust
          if (key %in% names(measures)) row[[metric]] <- unname(measures[key])
        }
      }
    }
    result$candidates <- rbind(result$candidates, row)
  }
  tab <- result$candidates
  acceptable <- tab$converged & tab$admissible & is.finite(tab$CFI) & tab$CFI >= .95 &
    is.finite(tab$RMSEA) & tab$RMSEA <= .06 & is.finite(tab$SRMR) & tab$SRMR <= .08
  if (any(acceptable)) {
    eligible <- which(acceptable)
    # BIC only across the same sample and likelihood estimator. A >=10 advantage
    # is a conservative declared heuristic, not proof of a population structure.
    chosen <- eligible[1L]
    if (estimator == "MLR" && all(is.finite(tab$BIC[eligible]))) {
      best <- eligible[which.min(tab$BIC[eligible])]
      if (!(1L %in% eligible) || tab$BIC[1L] - tab$BIC[best] >= 10) chosen <- best
    }
    result$recommendation <- tab$model[chosen]
    result$reason <- paste("Exploratory candidate meets declared admissibility/fit heuristics;",
      if (estimator == "MLR") "BIC advantage >=10 required over an acceptable correlated model." else
        "Ordinal WLSMV: AIC/BIC unavailable; prefer the first acceptable, simpler candidate.")
  } else result$reason <- "No candidate met all declared convergence, admissibility and fit heuristics."
  result
}

.pm_network_hierarchy_table <- function(native) {
  orders <- c("lower_order","higher_order")
  do.call(rbind,lapply(orders,function(order) {
    part <- native[[order]]
    data.frame(order=order,dimensions=if(is.null(part$n.dim))NA_real_ else part$n.dim,
      TEFI=if(is.null(part$TEFI))NA_real_ else as.numeric(part$TEFI),
      separately_estimated=order=="lower_order" || isTRUE(native$lower_order$n.dim>1))
  }))
}

.pm_network_cfa_encoding <- function(items, spec) {
  ordinal <- spec$measurement_level[items] == "ordinal"
  ordered <- if (any(ordinal)) paste0("x", which(ordinal)) else NULL
  list(ordered = ordered, estimator = if (length(ordered)) "WLSMV" else "MLR")
}

.pm_network_validate_ri <- function(adjusted) {
  fit <- adjusted$RI$fit
  if (!inherits(fit,"lavaan") || !isTRUE(lavaan::lavInspect(fit,"converged")))
    stop("Random-intercept model did not converge; riEGA residual interpretation is not evaluable.", call. = FALSE)
  admissibility <- .pm_cfa_admissibility(fit, converged = TRUE)
  parameters <- lavaan::parameterEstimates(fit)
  variance <- parameters$est[parameters$op == "~~" & parameters$lhs == parameters$rhs &
    parameters$lhs %in% lavaan::lavNames(fit,"lv")]
  if (!isTRUE(admissibility$admissible) || !length(variance) ||
      any(!is.finite(variance)) || any(variance < 0))
    stop("Random-intercept model is inadmissible; riEGA residual interpretation is not evaluable.", call. = FALSE)
  list(converged = TRUE, admissible = TRUE, random_intercept_variance = variance)
}
