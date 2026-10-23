test_that("partial targets share native settings and fit values with complete invariance", {
  .pm_skip_extended()
  model <- "F =~ i1+i2+i3+i4+i5+i6"
  for (estimator in c("MLR", "WLSMV", "ULSMV")) {
    d <- simulate_invariance_fixture(180L, ordinal = estimator != "MLR", seed = 7145L)
    for (variant in c("scaled", "robust")) {
      complete <- factorial_invariance_auto(d, "group", model,
        estimator = estimator, fit_variant = variant, verbose = FALSE)
      partial <- invariance_partial(d, "group", model, free = "F =~ i3",
        rationale = "Prior substantive evidence for a group-specific loading.",
        stage = "metric", estimator = estimator, fit_variant = variant, verbose = FALSE)
      expect_identical(partial$effective_settings$workflow, "factorial")
      expect_identical(partial$effective_settings$fit_variant, variant)
      for (stage in c("configural", "metric")) {
        role <- if (stage == "metric") "full" else stage
        expect_equal(lavaan::fitMeasures(partial$models[[role]]$fit),
          lavaan::fitMeasures(complete$fitted_models[[stage]]), tolerance = 1e-8)
      }
      native <- .pm_capture_conditions(lavaan::cfa(model, data = d, group = "group",
        estimator = estimator, ordered = estimator != "MLR", missing = if (estimator == "MLR") "fiml" else "listwise",
        group.equal = "loadings", group.partial = "F =~ i3"))$value
      expect_equal(lavaan::fitMeasures(partial$models$partial$fit),
        lavaan::fitMeasures(native), tolerance = 1e-8)
      expect_true(partial$freed_constraints$effective_release)
      expect_false(lavaan::lavInspect(partial$models$partial$fit, "options")$std.lv)
      fm <- lavaan::fitMeasures(native)
      row <- partial$fit_table[partial$fit_table$role == "partial", ]
      expect_equal(row$cfi, unname(fm[paste0("cfi.", variant)]))
      expect_equal(row$rmsea, unname(fm[paste0("rmsea.", variant)]))
      expect_equal(row$chi2, unname(fm["chisq.scaled"]))
      expect_equal(row$srmr, unname(fm["srmr"]))
      for (language in c("eng", "esp")) {
        text <- paste(unlist(.pm_analysis_sections(partial, language)$sections), collapse = " ")
        expect_match(text, variant, fixed = TRUE)
      }
    }
  }
})

test_that("partial invariance aligns numeric WLSMV and categorical MLR preparation", {
  d <- simulate_invariance_fixture(180L, ordinal = TRUE, seed = 7146L)
  model <- "F =~ i1+i2+i3+i4+i5+i6"
  for (estimator in c("MLR", "WLSMV")) {
    complete <- factorial_invariance_auto(d, "group", model, estimator = estimator,
      ordered = FALSE, verbose = FALSE)
    partial <- invariance_partial(d, "group", model, "i3 ~ 1", "Substantive rationale",
      estimator = estimator, ordered = FALSE, verbose = FALSE)
    expect_identical(partial$measurement, "continuous")
    expect_identical(partial$analysis_data, complete$analysis_data)
    expect_equal(lavaan::fitMeasures(partial$models$full$fit),
      lavaan::fitMeasures(complete$fitted_models$scalar), tolerance = 1e-8)
    expect_true(partial$freed_constraints$effective_release)
  }
})

test_that("native partial targets reject threshold releases and unsupported fit variants", {
  d <- simulate_invariance_fixture(120L, ordinal = TRUE, seed = 7147L)
  model <- "F =~ i1+i2+i3+i4+i5+i6"
  expect_error(invariance_partial(d, "group", model, "i3 | t1", "Prior evidence",
    stage = "metric", estimator = "WLSMV", verbose = FALSE), "inactive")
  expect_error(invariance_partial(d, "group", model, "F =~ i3", "Prior evidence",
    stage = "metric", fit_variant = "invalid", verbose = FALSE), "fit_variant|arg")
  expect_error(invariance_partial(d, "group", model, "F =~ i3", "Prior evidence",
    stage = "loadings", verbose = FALSE), "workflow|stage")
})

test_that("native scalar and strict releases match lavaan and preserve matched prerequisites", {
  d <- simulate_invariance_fixture(180L, ordinal = FALSE, partial = TRUE, seed = 7148L)
  model <- "F =~ i1+i2+i3+i4+i5+i6"
  for (stage in c("scalar", "strict")) {
    free <- c("F =~ i3", "i3 ~ 1")
    equal <- c("loadings", "intercepts")
    if (stage == "strict") {
      free <- c(free, "i4 ~~ i4")
      equal <- c(equal, "residuals")
    }
    out <- invariance_partial(d, "group", model, free, "Prespecified item-content evidence",
      stage = stage, estimator = "MLR", fit_variant = "robust", verbose = FALSE)
    native <- lavaan::cfa(model, data = d, group = "group", estimator = "MLR",
      missing = "fiml", group.equal = equal, group.partial = free)
    expect_equal(lavaan::fitMeasures(out$models$partial$fit),
      lavaan::fitMeasures(native), tolerance = 1e-8)
    expect_true(all(out$freed_constraints$effective_release))
    expect_true(out$comparison$group_partial_nesting$comparable)
    expect_gte(min(out$anchor_counts$retained_count), 2L)
    comparison <- out$comparison$prerequisite_vs_partial
    if (isTRUE(comparison$comparable)) {
      fm_previous <- lavaan::fitMeasures(out$models$prerequisite_matched$fit)
      fm_partial <- lavaan::fitMeasures(native)
      expect_equal(comparison$delta_cfi, unname(fm_partial["cfi.robust"] - fm_previous["cfi.robust"]))
      expect_equal(comparison$delta_rmsea, unname(fm_partial["rmsea.robust"] - fm_previous["rmsea.robust"]))
    }
  }
})

test_that("unavailable robust partial indices remain missing and stop decisions", {
  d <- simulate_invariance_fixture(120L, seed = 7149L)
  extract <- .ifa_extract_measures
  testthat::local_mocked_bindings(.ifa_extract_measures = function(fit) {
    measures <- extract(fit)
    for (index in c("cfi", "rmsea")) measures[[index]]["robust"] <- NA_real_
    measures
  }, .package = "PsychoMatic")
  out <- invariance_partial(d, "group", "F =~ i1+i2+i3+i4+i5+i6", "F =~ i3",
    "Prior substantive evidence", stage = "metric", estimator = "MLR",
    ordered = FALSE, fit_variant = "robust", verbose = FALSE)
  expect_true(all(is.na(out$fit_table$cfi)))
  expect_true(all(is.na(out$fit_table$rmsea)))
  expect_identical(out$conclusion_status, "not_evaluable")
  expect_false(any(out$freed_constraints$effective_release))
  expect_identical(out$effective_settings$fit_variant, "robust")
  expect_true(is.finite(out$models$configural$fit_record$measures$cfi[["scaled"]]))
})
