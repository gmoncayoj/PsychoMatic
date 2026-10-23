test_that('ESEM results retain native standardized inference and manuscript text', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(), verbose = FALSE,
                   rotation_args = list(rstarts = 3), language = 'esp')
  expect_s3_class(out, 'esem_auto')
  expect_s3_class(out, 'psychomatic_result')
  expect_identical(out$fit, out$results$fit)
  native <- lavaan::standardizedSolution(out$fit)
  ld <- native[native$op == '=~', ]
  expect_equal(out$results$loadings$standardized, ld$est.std, tolerance = 1e-7)
  expect_equal(out$results$loadings$ci_lower, ld$ci.lower, tolerance = 1e-7)
  expect_equal(nrow(out$results$loadings), 27)
  expect_identical(out$conclusion$status, 'evaluable')
  expect_match(out$manuscript_text, '400')
  expect_match(out$manuscript_text, 'target')
  expect_match(out$manuscript_text, 'CFI')
  expect_false(grepl('confirm[oó]|causal', out$manuscript_text))
  expect_identical(summary(out)$manuscript_text, out$manuscript_text)
})

test_that('requested unavailable indices remain unavailable without invalidating estimation', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(), estimator = 'ML',
                   fit_variant = 'robust', verbose = FALSE, rotation_args = list(rstarts = 3))
  fi <- out$results$fit_indices
  expect_true(all(is.na(fi$value[fi$metric %in% c('cfi', 'tli', 'rmsea')])))
  expect_true(all(nzchar(fi$reason[fi$metric %in% c('cfi', 'tli', 'rmsea')])))
  expect_identical(out$conclusion$status, 'evaluable')
  expect_match(out$manuscript_text, 'unavailable')
})

test_that('failed estimation cannot produce substantive manuscript interpretation', {
  p <- esem_prepare(rotation_args = list(rstarts = 3))
  failure <- list(fit = NULL, warnings = list(), error = simpleError('native failure'))
  out <- PsychoMatic:::.pm_esem_result(p, failure, quote(esem_auto(d)), 'scaled', 'eng')
  expect_null(out$fit)
  expect_identical(out$conclusion$status, 'not_evaluable')
  expect_match(out$manuscript_text, 'not interpretable|could not')
  expect_false(grepl('supported|adequate fit|confirmed', out$manuscript_text))
})

test_that('geomin does not assign substantive factor names', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(), rotation = 'geomin',
    verbose = FALSE, rotation_args = list(rstarts = 3))
  expect_identical(out$results$factor_map$label, paste0('F', 1:3))
  expect_true(all(is.na(out$results$item_diagnostics$expected_factor)))
  expect_match(out$manuscript_text, 'exploratory')
})
