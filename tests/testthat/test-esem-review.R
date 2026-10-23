test_that('a robust-test warning does not invalidate otherwise finite inference', {
  p <- esem_prepare(rotation_args = list(rstarts = 1))
  fit <- PsychoMatic:::.pm_esem_fit(p, 2026)
  fit$warnings <- c(fit$warnings, list(simpleWarning(paste0(
    'lavaan->lav_test_sb(): could not invert information matrix needed for robust test statistic'))))
  out <- PsychoMatic:::.pm_esem_result(p, fit, quote(probe()), 'scaled', 'eng')
  expect_true(out$diagnostics$finite_se)
  expect_identical(out$conclusion$status, 'evaluable')
  expect_match(paste(out$warnings$message, collapse = ' '), 'robust test statistic')
})

test_that('ESEM preserves declared and inferred preparation provenance', {
  d <- esem_fixture(TRUE)
  d$x1 <- ordered(d$x1, levels = 1:5,
    labels = c('never', 'rarely', 'sometimes', 'often', 'always'))
  s <- psychomatic_spec(d, items = names(d), reverse = c('x1', 'x2'),
    ranges = list(x1 = c(1, 5)),
    ordinal_levels = list(x1 = c('never', 'rarely', 'sometimes', 'often', 'always')))
  out <- esem_auto(s, factors = esem_factors(), rotation_args = list(rstarts = 1), verbose = FALSE)
  prep <- out$effective_settings$preparation
  expect_identical(prep$reverse, s$reverse)
  expect_identical(prep$ordinal_levels, s$ordinal_levels)
  expect_identical(prep$ranges_effective$x1, s$ranges$x1)
  expect_equal(prep$ranges_effective$x2, range(as.integer(s$data$x2), na.rm = TRUE))
  expect_true(any(out$decision_log$code == 'reverse'))
  expect_true(any(out$decision_log$code == 'ordered'))
  expect_true(is.data.frame(out$results$preparation))
})

test_that('reader reports expose target maps versions and ordinal selection', {
  out <- esem_auto(esem_fixture(TRUE), factors = esem_factors(), ordered = c('x1','x2'),
    rotation_args = list(rstarts = 1), verbose = FALSE)
  for (lang in c('eng','esp')) for (destination in c('console','report')) {
    sec <- PsychoMatic:::.pm_analysis_sections(out, lang, destination)$sections
    technical <- if(lang == 'eng') c('Item map','Factor map','Preparation','Versions') else
      c('Mapa de ítems','Mapa de factores','Preparación','Versiones')
    expect_true((if(lang == 'eng') 'Target pattern' else 'Patrón target') %in% names(sec))
    if (destination == 'console') {
      expect_false(any(technical %in% names(sec)))
    } else {
      expect_true(all(technical %in% names(sec)))
      prep <- sec[[if(lang == 'eng') 'Preparation' else 'Preparación']]
      expect_true(all(c('x1','x2') %in% prep[[1]]))
      expect_match(paste(unlist(prep),collapse=' '), 'ordinal')
    }
    expect_identical(tail(names(sec),1), if(lang=='eng') 'Interpretation' else 'Interpretaci\u00f3n')
  }
})

test_that('raw loading intervals match parameterEstimates without relabeling', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(),
    rotation_args = list(rstarts = 1), verbose = FALSE)
  pe <- lavaan::parameterEstimates(out$fit)
  ld <- pe[pe$op == '=~', ]
  expect_equal(out$results$loadings$raw_ci_lower, ld$ci.lower, tolerance = 1e-7)
  expect_equal(out$results$loadings$raw_ci_upper, ld$ci.upper, tolerance = 1e-7)
})

test_that('fit-index conditions are retained without rejecting parameter inference', {
  p <- esem_prepare(rotation_args = list(rstarts = 1))
  fit <- PsychoMatic:::.pm_esem_fit(p, 2026)
  testthat::local_mocked_bindings(fitMeasures = function(...) {
    warning('controlled fit-index warning')
    stop('controlled fit-index failure')
  }, .package = 'lavaan')
  out <- PsychoMatic:::.pm_esem_result(p, fit, quote(probe()), 'robust', 'eng')
  expect_identical(out$conclusion$status, 'evaluable')
  expect_true(all(is.na(out$results$fit_indices$value)))
  expect_match(paste(out$warnings$message,collapse=' '), 'fit-index warning')
  expect_match(out$diagnostics$fit_indices_error, 'fit-index failure')
  expect_match(paste(out$results$fit_indices$reason,collapse=' '), 'fit-index failure')
})
