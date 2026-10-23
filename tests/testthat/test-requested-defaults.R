test_that("automatic correlations distinguish binary ordinal and symmetric scores", {
  set.seed(920)
  d <- as.data.frame(replicate(3, sample(1:6, 600, TRUE)))
  txt <- capture.output(r <- cormat(d))
  expect_identical(attr(r, "requested_type"), "auto")
  expect_identical(attr(r, "effective_type"), "pearson")
  expect_true(any(grepl("PEARSON CORRELATIONS", txt, fixed = TRUE)))
  expect_true(any(grepl("RATIONALE", txt, fixed = TRUE)))
  expect_gte(length(attr(r, "references")), 2)
  b <- as.data.frame(lapply(d, function(x) as.integer(x > 3)))
  expect_identical(attr(cormat(b, verbose = FALSE), "effective_type"), "tetra")
  o <- as.data.frame(lapply(d, function(x) pmin(x, 5)))
  expect_identical(attr(cormat(o, verbose = FALSE), "effective_type"), "poly")
  skew <- as.data.frame(replicate(3, sample(1:6, 600, TRUE,
    prob = c(.7, .1, .06, .05, .05, .04))))
  expect_identical(attr(cormat(skew, verbose = FALSE), "effective_type"), "poly")
  expect_identical(attr(cormat(skew, type = "pearson", verbose = FALSE), "effective_type"), "pearson")
  mixed <- data.frame(binary = b[[1]], continuous = rnorm(600))
  expect_error(cormat(mixed, verbose = FALSE), "mixed|Mixed")
})

test_that("analysis defaults produce English unless Spanish is explicit", {
  d <- lavaan::HolzingerSwineford1939[paste0("x", 1:4)]
  expect_true(any(grepl("DESCRIPT", capture.output(desc_auto(d)), fixed = TRUE)))
  expect_identical(attr(cormat(d, verbose = FALSE), "display_language"), "eng")
  expect_identical(attr(cormat(d, verbose = FALSE, language = "esp"), "display_language"), "esp")
  x <- reliability_auto(d, verbose = FALSE)
  expect_identical(x$effective_settings$language, "eng")
})

test_that("CFA computes Mardia with manual and ordinal estimators", {
  set.seed(921)
  f <- rnorm(350)
  d <- as.data.frame(replicate(4, as.integer(cut(f + rnorm(350),
    c(-Inf, -1, -.3, .3, 1, Inf)))))
  names(d) <- paste0("i", 1:4)
  x <- cfa_auto(d, "f =~ i1 + i2 + i3 + i4", estimator = "WLSMV", verbose = FALSE)
  expect_identical(x$mardia$status, "ok")
  expect_true(all(is.finite(c(x$mardia$p_skew, x$mardia$p_kurt))))
  expect_identical(x$language, "eng")
  s <- PsychoMatic:::.pm_analysis_sections(x, "esp")$sections
  expect_true("Alfa (Cronbach)" %in% s$Confiabilidad$Indice)
  expect_true(any(grepl("codigos observados", unlist(s), fixed = TRUE)))
  y <- cfa_auto(lavaan::HolzingerSwineford1939, "f =~ x1 + x2 + x3 + x4",
    estimator = "ML", verbose = FALSE)
  expect_identical(y$mardia$status, "ok")
})

test_that("invariance uses signed current minus previous changes and scaled priority", {
  d <- PsychoMatic:::.ifa_deteriorations(c(cfi=.96, rmsea=.04, srmr=.03),
    c(cfi=.94, rmsea=.06, srmr=.07))
  expect_equal(unname(d), c(-.02, .02, .04))
  expect_identical(PsychoMatic:::.ifa_select_variant(c(robust=.80,scaled=.95,base=.97))$variant, "scaled")
  expect_identical(PsychoMatic:::.ifa_change_decision(-.02, .02, .04,
    list(cfi=.01, rmsea=.015, srmr_loading=.03, srmr_other=.01), "metric")$decision,
    "not_supported")
})

test_that("ordinal invariance reports the actual scaled indices and threshold rationale", {
  d <- simulate_invariance_fixture(180L, ordinal = TRUE, seed = 922L)
  x <- ordinal_invariance_auto(d, "group", paste0("f =~ ", paste0("i", 1:6, collapse=" + ")),
    estimator = "WLSMV", verbose = FALSE)
  expect_identical(x$effective_settings$language, "eng")
  expect_gt(length(x$fitted_models), 0L)
  for (stage in names(x$fitted_models)) {
    expect_s4_class(x$fitted_models[[stage]], "lavaan")
    row <- x$fit_table[x$fit_table$level == stage, ]
    expect_equal(row$cfi, unname(lavaan::fitMeasures(x$fitted_models[[stage]], "cfi.scaled")))
    expect_equal(row$rmsea, unname(lavaan::fitMeasures(x$fitted_models[[stage]], "rmsea.scaled")))
  }
  expect_true(all(x$fit_table$cfi_variant[!is.na(x$fit_table$cfi_variant)] == "scaled"))
  sections <- PsychoMatic:::.pm_analysis_sections(x)$sections
  expect_false(any(c("Fit-index variants", "Comparison criteria") %in% names(sections)))
  expect_true("Ordinal method and interpretation" %in% names(sections))
  expect_true("Configural evaluation" %in% names(sections))
  expect_true(any(grepl("Bowen", unlist(sections), fixed=TRUE)))
  record <- x$fit_records[[1]]
  record$measures$cfi[c("scaled")] <- NA_real_
  expect_false(PsychoMatic:::.ifa_admissibility(record)$admissible)
})

test_that("automatic correlation handles ordered labels and continuous skew safely", {
  set.seed(923)
  x <- as.data.frame(replicate(3, ordered(sample(letters[1:6], 600, TRUE), levels=letters[1:6]), simplify=FALSE))
  expect_identical(attr(cormat(x, verbose=FALSE), "effective_type"), "pearson")
  d <- as.data.frame(replicate(3, rexp(600)))
  y <- cormat(d, verbose=FALSE)
  expect_identical(attr(y, "effective_type"), "pearson")
  expect_true(grepl("Shape flags", attr(y, "rationale"), fixed=TRUE))
  expect_error(cormat(data.frame(a=rep(1,10), b=1:10), verbose=FALSE), "two observed")
  expect_error(cormat(data.frame(a=letters[1:6], b=1:6), verbose=FALSE), "mapping|ordered")
})
