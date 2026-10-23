test_that("WLSMV defaults reproduce literal native ordered lavaan calls", {
  d <- simulate_invariance_fixture(180L, ordinal = TRUE, seed = 7140L)
  model <- "F =~ i1+i2+i3+i4+i5+i6"
  constraints <- list(configural = character(), metric = "loadings",
    scalar = c("loadings", "intercepts"),
    strict = c("loadings", "intercepts", "residuals"))
  fit <- factorial_invariance_auto(d, "group", model,
    estimator = "WLSMV", verbose = FALSE)
  expect_identical(fit$measurement, "ordinal")
  expect_identical(fit$fit_table$level, names(constraints))
  for (stage in names(constraints)) {
    native <- .pm_capture_conditions(lavaan::cfa(model, data = d, group = "group",
      estimator = "WLSMV", ordered = TRUE, group.equal = constraints[[stage]]))$value
    expect_equal(lavaan::fitMeasures(fit$fitted_models[[stage]]),
      lavaan::fitMeasures(native), tolerance = 1e-8)
    expect_false(lavaan::lavInspect(fit$fitted_models[[stage]], "options")$std.lv)
    expect_identical(fit$fit_records[[stage]]$engine, "lavaan::cfa")
    expect_identical(fit$fit_records[[stage]]$parameterization, "delta")
  }
  expect_true(all(fit$fit_table$cfi_variant == "scaled"))
  expect_false(fit$fit_table$decision[fit$fit_table$level == "strict"] == "supported")
  numeric_data <- d
  numeric_data[paste0("i", 1:6)] <- lapply(d[paste0("i", 1:6)], function(x) as.numeric(as.character(x)))
  numeric_fit <- factorial_invariance_auto(numeric_data, "group", model,
    estimator = "WLSMV", verbose = FALSE)
  expect_identical(numeric_fit$measurement, "ordinal")
  expect_equal(lavaan::fitMeasures(numeric_fit$fitted_models$configural),
    lavaan::fitMeasures(fit$fitted_models$configural))
  continuous <- factorial_invariance_auto(numeric_data, "group", model,
    estimator = "WLSMV", ordered = FALSE, verbose = FALSE)
  expect_identical(continuous$measurement, "continuous")
  expect_false(any(lavaan::parameterTable(continuous$fitted_models$configural)$op == "|"))
})

test_that("requested fit variants govern tables, deltas and reports", {
  d <- simulate_invariance_fixture(180L, ordinal = TRUE, seed = 7141L)
  model <- "F =~ i1+i2+i3+i4+i5+i6"
  for (variant in c("scaled", "robust")) {
    fit <- factorial_invariance_auto(d, "group", model, estimator = "WLSMV",
      fit_variant = variant, verbose = FALSE)
    expect_identical(fit$effective_settings$fit_variant, variant)
    for (i in seq_along(fit$fitted_models)) {
      fm <- lavaan::fitMeasures(fit$fitted_models[[i]])
      expect_equal(fit$fit_table$cfi[i], unname(fm[paste0("cfi.", variant)]))
      expect_equal(fit$fit_table$rmsea[i], unname(fm[paste0("rmsea.", variant)]))
      expect_equal(fit$fit_table$chi2[i], unname(fm["chisq.scaled"]))
      expect_equal(fit$fit_table$srmr[i], unname(fm["srmr"]))
    }
    comparable <- which(fit$fit_table$comparable %in% TRUE)
    for (i in comparable) {
      expect_equal(fit$fit_table$delta_cfi[i], fit$fit_table$cfi[i] - fit$fit_table$cfi[i-1L])
      expect_equal(fit$fit_table$delta_rmsea[i], fit$fit_table$rmsea[i] - fit$fit_table$rmsea[i-1L])
    }
    for (language in c("eng", "esp")) {
      sections <- .pm_analysis_sections(fit, language)$sections
      expect_match(paste(unlist(sections), collapse = " "), variant, fixed = TRUE)
    }
  }
  expect_error(factorial_invariance_auto(d, "group", model,
    fit_variant = "invalid", verbose = FALSE), "fit_variant|arg")
  record <- list(estimator = "DWLS", estimator_orig = "WLSMV", fit_variant = "robust")
  expect_identical(.ifa_record_priority(record, "cfi"), "robust")
  expect_true(is.na(.ifa_select_variant(c(robust = NA, scaled = .95, base = .96),
    .ifa_record_priority(record, "cfi"))$value))
})

test_that("ordinal identification remains distinct and permits robust reporting", {
  d <- simulate_invariance_fixture(180L, ordinal = TRUE, seed = 7142L)
  fit <- ordinal_invariance_auto(d, "group", "F =~ i1+i2+i3+i4+i5+i6",
    fit_variant = "robust", verbose = FALSE)
  expect_identical(fit$fit_table$level, c("configural", "thresholds", "thresholds_loadings"))
  expect_identical(fit$fit_records$configural$engine, "semTools::measEq.syntax")
  expect_identical(fit$fit_records$configural$parameterization, "theta")
  expect_equal(fit$fit_table$cfi[1], unname(lavaan::fitMeasures(fit$fitted_models$configural)["cfi.robust"]))
})
