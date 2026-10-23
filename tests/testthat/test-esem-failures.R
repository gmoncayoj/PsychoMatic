test_that('nonconvergent and Heywood native solutions cannot support a manuscript', {
  .pm_skip_extended()
  d <- esem_fixture()[, 1:6]
  out <- esem_auto(d, items = names(d), n_factors = 3, rotation = 'geomin',
    rotation_args = list(rstarts = 1), seed = 2026, verbose = FALSE)
  expect_false(out$diagnostics$converged)
  expect_false(out$diagnostics$admissible)
  expect_identical(out$conclusion$status, 'not_evaluable')
  expect_match(paste(out$warnings$message, collapse = ' '), 'variances are negative')
  expect_match(out$manuscript_text, 'not interpretable')
})

test_that('a saturated summary cannot treat perfect fit as validation', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(),
    rotation_args = list(rstarts = 1), seed = 2026, verbose = FALSE)
  # Controlled saturation scenario: native saturated probes were inadmissible.
  # This isolates the reporting rule without claiming a valid statistical fit.
  out$diagnostics$saturated <- TRUE
  out$results$fit_indices$value[out$results$fit_indices$metric == 'df'] <- 0
  text <- PsychoMatic:::.pm_esem_manuscript(out, 'eng')
  expect_match(text, 'no positive degrees of freedom')
  expect_match(text, 'cannot assess')
  expect_identical(out$conclusion$status, 'evaluable')
})

test_that('a native underidentified model cannot support interpretation', {
  d <- esem_fixture()[, 1:6]
  out <- esem_auto(d, items = names(d), n_factors = 4, rotation = 'geomin',
    rotation_args = list(rstarts = 1), seed = 2026, verbose = FALSE)
  expect_true(out$diagnostics$converged)
  # Native warning wording depends on BLAS and numerical rounding. Structural
  # identification and the public conclusion must not depend on that wording.
  expect_true(out$diagnostics$underidentified)
  df <- out$results$fit_indices$base[out$results$fit_indices$metric == 'df']
  expect_length(df, 1L)
  expect_lt(df, 0)
  expect_identical(out$conclusion$status, 'not_evaluable')
  expect_match(out$manuscript_text, 'not interpretable')
})

test_that('controlled parameter-inference warnings block interpretation', {
  p <- esem_prepare(rotation_args = list(rstarts = 1))
  fitted <- PsychoMatic:::.pm_esem_fit(p, 2026)
  # Deliberately controlled engine condition, not a statistical reproduction.
  messages <- c(
    'Could not compute standard errors: information matrix could not be inverted',
    'The variance-covariance matrix of the estimated parameters (vcov) does not appear to be positive definite!')
  for (message in messages) {
    probe <- fitted
    probe$warnings <- c(probe$warnings, list(simpleWarning(message)))
    out <- PsychoMatic:::.pm_esem_result(p, probe, quote(probe()), 'scaled', 'eng')
    expect_true(message %in% out$warnings$message)
    expect_true(out$diagnostics$inference_warning)
    expect_identical(out$conclusion$status, 'not_evaluable')
    expect_match(out$manuscript_text, 'not interpretable')
  }
})
