test_that('negative degrees of freedom block interpretation without a native warning', {
  p <- esem_prepare(rotation_args = list(rstarts = 1))
  fitted <- PsychoMatic:::.pm_esem_fit(p, 2026)
  native_indices <- PsychoMatic:::.pm_esem_indices
  result_with_df <- function(df) {
    # Controlled index scenario: isolates structural identification from native warnings.
    testthat::local_mocked_bindings(.pm_esem_indices = function(fit, variant) {
      indices <- native_indices(fit, variant)
      row <- indices$metric == 'df'
      indices$base[row] <- df
      indices$value[row] <- df
      indices
    }, .package = 'PsychoMatic')
    PsychoMatic:::.pm_esem_result(p, fitted, quote(probe()), 'scaled', 'eng')
  }
  negative <- result_with_df(-3)
  expect_true(negative$diagnostics$finite_se)
  expect_identical(negative$conclusion$status, 'not_evaluable')
  expect_true(negative$diagnostics$underidentified)
  expect_false(negative$diagnostics$saturated)
  expect_match(negative$conclusion$reason, 'negative degrees of freedom')
  expect_match(negative$manuscript_text, 'not interpretable')
  saturated <- result_with_df(0)
  expect_identical(saturated$conclusion$status, 'evaluable')
  expect_false(saturated$diagnostics$underidentified)
  expect_true(saturated$diagnostics$saturated)
  expect_match(saturated$manuscript_text, 'no positive degrees of freedom')
})
