test_that('ESEM preparation preserves target order and safe names', {
  d <- esem_fixture(); before <- serialize(d, NULL)
  p <- esem_prepare(d, factors = rev(esem_factors()))
  expect_identical(p$item_map$item, unname(unlist(rev(esem_factors()))))
  expect_equal(unname(p$target), kronecker(diag(3), matrix(1, 3, 1)))
  expect_identical(serialize(d, NULL), before)
  names(d)[1] <- 'x 1 + evil'; f <- esem_factors(); f[[1]][1] <- names(d)[1]
  p <- esem_prepare(d, f)
  expect_false(grepl('evil', p$syntax, fixed = TRUE))
  expect_identical(p$item_map$item[1], names(d)[1])
})

test_that('ESEM rejects contradictory or invalid specifications', {
  expect_error(esem_prepare(factors = NULL), 'factors')
  expect_error(esem_prepare(items = 'x1'), 'items')
  expect_error(esem_prepare(n_factors = 3), 'n_factors')
  f <- esem_factors(); f[[2]][1] <- 'x1'
  expect_error(esem_prepare(factors = f), 'once|duplicat')
  f <- esem_factors(); names(f)[2] <- names(f)[1]
  expect_error(esem_prepare(factors = f), 'unique')
  f <- esem_factors(); f[[1]][1] <- 'unknown'
  expect_error(esem_prepare(factors = f), 'Unknown')
  expect_error(esem_prepare(rotation_args = list(geomin_epsilon = .01)), 'geomin')
  expect_error(esem_prepare(rotation_args = list(orthogonal = TRUE)), 'rotation_args')
  expect_error(esem_prepare(rotation_args = list(rstarts = 0)), 'rstarts')
  expect_error(esem_prepare(ordered = TRUE, estimator = 'MLR'), 'estimator')
  expect_error(esem_prepare(ordered = TRUE, missing = 'fiml'), 'missing')
  expect_error(esem_prepare(factors = NULL, rotation = 'geomin', n_factors = 1, items = 1:9), 'n_factors')
})

test_that('ESEM honors data specifications and category meaning', {
  d <- esem_fixture(TRUE)
  s <- psychomatic_spec(d, items = names(d), level = 'ordinal')
  a <- esem_prepare(d); b <- esem_prepare(s)
  expect_identical(a$analysis_data, b$analysis_data)
  expect_identical(a$native_args$ordered, b$native_args$ordered)
  d[[1]] <- factor(d[[1]])
  expect_error(esem_prepare(d), 'ambiguous|nominal|ordered')
  d[[1]] <- ordered(d[[1]])
  expect_error(esem_prepare(d, ordered = FALSE), 'numeric')
  p <- esem_prepare(esem_fixture(), factors = NULL, n_factors = 3,
                    items = 1:9, rotation = 'geomin')
  expect_null(p$target)
  expect_identical(p$factor_map$label, paste0('F', 1:3))
})

test_that('ESEM applies declared reversals without mutating specifications', {
  d <- esem_fixture(TRUE)
  s <- psychomatic_spec(d, items = names(d), reverse = 'x1', ranges = list(x1 = c(1, 5)))
  before <- serialize(s, NULL)
  p <- esem_prepare(s)
  expect_equal(p$analysis_data$pm_i001, 6 - d$x1)
  expect_identical(serialize(s, NULL), before)
  d$group <- rep(c('a', 'b'), each = 200)
  s <- psychomatic_spec(d, items = paste0('x', 1:9), group = 'group')
  expect_error(esem_prepare(s), 'group|single')
})
