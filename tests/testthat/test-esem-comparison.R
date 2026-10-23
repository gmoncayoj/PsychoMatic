test_that('ESEM CFA reference has the same data and effective options', {
  d <- esem_fixture(TRUE, TRUE)
  out <- esem_auto(d, factors = esem_factors(), compare_cfa = TRUE,
    missing = 'pairwise', verbose = FALSE, rotation_args = list(rstarts = 3))
  p <- esem_prepare(d, missing = 'pairwise', rotation_args = list(rstarts = 3))
  args <- p$native_args; args$model <- p$cfa_syntax; args$rotation <- NULL
  native <- do.call(lavaan::cfa, args)
  cmp <- out$results$cfa_comparison
  expect_identical(cmp$status, 'evaluable')
  expect_equal(lavaan::coef(cmp$fit), lavaan::coef(native), tolerance = 1e-7)
  expect_equal(lavaan::lavInspect(cmp$fit, 'case.idx'), lavaan::lavInspect(out$fit, 'case.idx'))
  expect_false(any(c('aic', 'bic', 'winner', 'lrt') %in% names(cmp$table)))
  expect_equal(cmp$table$delta, cmp$table$esem - cmp$table$cfa)
  expect_match(out$manuscript_text, 'descriptive')
  expect_error(esem_auto(d, items = 1:9, n_factors = 3, rotation = 'geomin',
    compare_cfa = TRUE, verbose = FALSE), 'factors')
})
test_that('failed, unmatched or unavailable CFA references remain descriptive', {
  p <- esem_prepare(rotation_args = list(rstarts = 1))
  fit <- PsychoMatic:::.pm_esem_fit(p, 2026)$fit
  bad <- p; bad$cfa_syntax <- 'invalid syntax'
  a <- PsychoMatic:::.pm_esem_compare_cfa(bad, fit, 2026, 'scaled')
  expect_identical(a$status, 'not_evaluable')
  expect_equal(nrow(a$table), 0)
  bad <- p; bad$native_args$data <- bad$native_args$data[1:350, ]
  a <- PsychoMatic:::.pm_esem_compare_cfa(bad, fit, 2026, 'scaled')
  expect_identical(a$status, 'not_evaluable')
  expect_equal(nrow(a$table), 0)
  p$native_args$estimator <- 'ML'
  fit <- PsychoMatic:::.pm_esem_fit(p, 2026)$fit
  a <- PsychoMatic:::.pm_esem_compare_cfa(p, fit, 2026, 'robust')
  expect_identical(a$status, 'evaluable')
  expect_true(all(is.na(a$table$delta[a$table$metric %in% c('cfi','tli','rmsea')])))
})
