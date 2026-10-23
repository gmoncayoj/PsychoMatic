test_that('ESEM default NULL seed follows native rotation RNG without imposing a seed', {
  d <- esem_fixture()
  p <- esem_prepare(d, factors = NULL, n_factors = 3, items = names(d),
    rotation = 'geomin', ordered = FALSE, rotation_args = list(rstarts = 2))
  expect_null(formals(esem_auto)$seed)
  PsychoMatic:::.efa_with_seed(444, {
    before <- .Random.seed
    native <- do.call(lavaan::sem, p$native_args)
    native_after <- .Random.seed
    assign('.Random.seed', before, .GlobalEnv)
    out <- esem_auto(d, n_factors = 3, items = names(d), rotation = 'geomin',
      ordered = FALSE, rotation_args = list(rstarts = 2), verbose = FALSE)
    expect_null(out$effective_settings$seed)
    expect_true('seed' %in% names(out$effective_settings))
    expect_null(out$requested_settings$seed)
    expect_identical(.Random.seed, native_after)
    expect_equal(lavaan::coef(out$fit), lavaan::coef(native), tolerance = 1e-7)
    expect_equal(lavaan::vcov(out$fit), lavaan::vcov(native), tolerance = 1e-7)
  })
  expect_no_error(esem_auto(d, n_factors = 3, items = names(d),
    rotation = 'geomin', ordered = FALSE, seed = NULL,
    rotation_args = list(rstarts = 2), verbose = FALSE))
})

test_that('ESEM explicit seeds remain reproducible and restore caller RNG', {
  d <- esem_fixture()
  PsychoMatic:::.efa_with_seed(81, {
    before <- .Random.seed
    a <- esem_auto(d, factors = esem_factors(), seed = 2026,
      rotation_args = list(rstarts = 2), verbose = FALSE)
    expect_identical(.Random.seed, before)
    b <- esem_auto(d, factors = esem_factors(), seed = 2026,
      rotation_args = list(rstarts = 2), verbose = FALSE)
    expect_identical(.Random.seed, before)
    expect_equal(lavaan::coef(a$fit), lavaan::coef(b$fit), tolerance = 1e-7)
    expect_identical(a$effective_settings$seed, 2026L)
  })
  for (seed in list(-1, NA_real_, c(1, 2), '2026')) {
    expect_error(esem_auto(d, factors = esem_factors(), seed = seed,
      verbose = FALSE), 'seed')
  }
})

test_that('ESEM console is concise and exports retain technical provenance', {
  out <- esem_auto(esem_fixture(TRUE), factors = esem_factors(), seed = 2026,
    rotation_args = list(rstarts = 2), verbose = FALSE)
  for (language in c('eng', 'esp')) {
    technical <- if (language == 'eng') c('Preparation', 'Item map', 'Factor map', 'Versions') else
      c('Preparaci\u00f3n', 'Mapa de \u00edtems', 'Mapa de factores', 'Versiones')
    final <- if (language == 'eng') c('References', 'Interpretation') else
      c('Referencias', 'Interpretaci\u00f3n')
    console <- PsychoMatic:::.pm_analysis_sections(out, language, 'console')$sections
    report <- PsychoMatic:::.pm_analysis_sections(out, language, 'report')$sections
    expect_false(any(technical %in% names(console)))
    expect_true(all(technical %in% names(report)))
    expect_identical(tail(names(console), 2), final)
    expect_identical(tail(names(report), 2), final)
    meta <- console[[1]]
    expect_identical(meta$Value[meta$Information == if (language == 'eng') 'Seed' else 'Semilla'], '2026')
    text <- tail(console, 1)[[1]]$Text
    expect_false(grepl('Este borrador describe|This draft describes', text))
    expect_match(text, 'CFI')
  }
  expect_true(is.data.frame(out$results$preparation))
  expect_true(is.data.frame(out$results$item_map))
  expect_true(length(out$versions) > 0L)
})

test_that('ESEM displays a NULL seed without dropping settings rows', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(), seed = 2026,
    rotation_args = list(rstarts = 2), verbose = FALSE)
  out$effective_settings['seed'] <- list(NULL)
  for (language in c('eng', 'esp')) {
    sec <- PsychoMatic:::.pm_analysis_sections(out, language)$sections
    meta <- sec[[1]]
    expect_identical(meta$Value[meta$Information == if (language == 'eng') 'Seed' else 'Semilla'], 'NULL')
    expect_true(all(c('GPA tolerance', if (language == 'eng') 'Maximum iterations' else
                     'Iteraciones m\u00e1ximas') %in% meta$Information))
  }
})
