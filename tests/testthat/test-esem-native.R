for (ord in c(FALSE, TRUE)) for (rotation in c('target', 'geomin')) {
  for (miss in c('listwise', if (ord) 'pairwise' else 'fiml')) local({
    ordinal <- ord; rot <- rotation; missing_policy <- miss
    test_that(paste('esem matches native', ordinal, rot, missing_policy), {
      # Both rotations and data types run normally; the full missing-data
      # Cartesian product is retained in the explicitly enabled CI suite.
      if (missing_policy != 'listwise') .pm_skip_extended()
      d <- esem_fixture(ordinal, TRUE)
      p <- esem_prepare(d, rotation = rot, missing = missing_policy,
                        rotation_args = list(rstarts = 3))
      a <- PsychoMatic:::.pm_esem_fit(p, 2026)
      b <- PsychoMatic:::.efa_with_seed(2026, do.call(lavaan::sem, p$native_args))
      expect_null(a$error)
      expect_true(lavaan::lavInspect(a$fit, 'converged'))
      expect_equal(lavaan::coef(a$fit), lavaan::coef(b), tolerance = 1e-7)
      expect_equal(lavaan::vcov(a$fit), lavaan::vcov(b), tolerance = 1e-7)
      expect_equal(lavaan::fitted(a$fit), lavaan::fitted(b), tolerance = 1e-7)
      expect_equal(lavaan::fitMeasures(a$fit), lavaan::fitMeasures(b), tolerance = 1e-7)
    })
  })
}

test_that('ESEM restores RNG with no prior seed, a prior seed and failed fitting', {
  p <- esem_prepare(rotation_args = list(rstarts = 3))
  PsychoMatic:::.efa_with_seed(2, {
    before <- .Random.seed
    a <- PsychoMatic:::.pm_esem_fit(p, 2026)
    expect_identical(.Random.seed, before)
    b <- PsychoMatic:::.pm_esem_fit(p, 2026)
    expect_identical(lavaan::coef(a$fit), lavaan::coef(b$fit))
    p$native_args$model <- 'invalid syntax'
    expect_s3_class(PsychoMatic:::.pm_esem_fit(p, 2026)$error, 'error')
    expect_identical(.Random.seed, before)
  })
  old <- if (exists('.Random.seed', .GlobalEnv)) .Random.seed else NULL
  on.exit(if (is.null(old) && exists('.Random.seed', .GlobalEnv)) rm(list = '.Random.seed', envir = .GlobalEnv) else if (!is.null(old))
    assign('.Random.seed', old, .GlobalEnv), add = TRUE)
  if (exists('.Random.seed', .GlobalEnv)) rm(list = '.Random.seed', envir = .GlobalEnv)
  PsychoMatic:::.pm_esem_fit(p, 2026)
  expect_false(exists('.Random.seed', .GlobalEnv))
})

test_that('ESEM fits a mixed continuous and ordinal block', {
  d <- esem_fixture(); d$x1 <- cut(d$x1, c(-Inf, -.5, .5, Inf), labels = FALSE)
  p <- esem_prepare(d, ordered = 'x1', rotation_args = list(rstarts = 3))
  fit <- PsychoMatic:::.pm_esem_fit(p, 2026)$fit
  expect_true(lavaan::lavInspect(fit, 'converged'))
  expect_identical(lavaan::lavNames(fit, 'ov.ord'), 'pm_i001')
})

test_that('rotation changes pattern without choosing by global fit', {
  d <- esem_fixture()
  target <- esem_prepare(d, rotation_args = list(rstarts = 3))
  geomin <- esem_prepare(d, rotation = 'geomin',
    rotation_args = list(rstarts = 3, geomin_epsilon = .01))
  a <- PsychoMatic:::.pm_esem_fit(target, 2026)$fit
  b <- PsychoMatic:::.pm_esem_fit(geomin, 2026)$fit
  expect_equal(lavaan::fitted(a), lavaan::fitted(b), tolerance = 1e-7)
  expect_equal(lavaan::fitMeasures(a), lavaan::fitMeasures(b), tolerance = 1e-7)
  expect_identical(geomin$settings$rotation_args$geomin_epsilon, .01)
  expect_equal(geomin$settings$rotation_args$rstarts, 3)
})
