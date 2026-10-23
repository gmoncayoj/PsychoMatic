test_that("automatic comparison treats five-category items as ordinal", {
  set.seed(2509)
  f <- rnorm(400)
  d <- as.data.frame(replicate(6, as.integer(cut(f + rnorm(400),
    c(-Inf, -1, -.3, .3, 1, Inf)))))
  names(d) <- paste0("i", 1:6)
  models <- list(one="f =~ i1+i2+i3+i4+i5+i6",
                 two="a =~ i1+i2+i3\nb =~ i4+i5+i6")
  x <- compare_models_auto(d, models, verbose=FALSE)
  expect_identical(x$estimator, "WLSMV")
  expect_identical(x$ordered, names(d))
  explicit <- compare_models_auto(d, models, estimator="ML", verbose=FALSE)
  expect_identical(explicit$estimator, "ML")
  expect_null(explicit$ordered)
  scaled <- compare_models_auto(d, models, estimator="WLSMV", verbose=FALSE)
  direct <- lavaan::cfa(models[[1]], d, estimator="WLSMV", ordered=names(d))
  expect_equal(scaled$fit$comparison_cfi[1], unname(lavaan::fitMeasures(direct, "cfi.scaled")))
  expect_equal(scaled$fit$comparison_rmsea[1], unname(lavaan::fitMeasures(direct, "rmsea.scaled")))
  reader <- .pm_analysis_sections(scaled, "eng")$sections
  expect_true("Data and method" %in% names(reader))
  expect_equal(reader$Fit[["CFI (scaled)"]][1], unname(lavaan::fitMeasures(direct, "cfi.scaled")))
  expect_true(all(is.na(scaled$fit$aic)))
  expect_true(all(is.na(scaled$fit$bic)))
  # A numeric ID outside the models must not alter category detection.
  d$id <- seq_len(nrow(d))
  alias <- compare_models_auto(d, models, estimator="auto", verbose=FALSE)
  expect_identical(alias$estimator, "WLSMV")
  expect_false("id" %in% alias$ordered)
})

test_that("continuous comparison retains the CFA automatic selection policy", {
  d <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  models <- list(one="f =~ x1+x2+x3+x4+x5+x6",
                 two="a =~ x1+x2+x3\nb =~ x4+x5+x6")
  x <- compare_models_auto(d, models, estimator="auto", verbose=FALSE)
  cfa <- cfa_auto(d, model=models[[1]], verbose=FALSE)
  expect_identical(x$estimator, cfa$estimator)
  expect_null(x$ordered)
  expect_identical(x$effective_settings$estimator_selection, "automatic")
  manual <- compare_models_auto(d, models, estimator="MLR", verbose=FALSE)
  expect_identical(manual$effective_settings$estimator_selection, "user")
  expect_equal(manual$fit$comparison_tli, manual$fit$tli_scaled)
})

test_that("missing scaled indices do not silently fall back to base or robust", {
  tab <- data.frame(estimator=c("WLSMV", "WLSMV"), cfi=c(.9,.9),
    cfi_scaled=c(.8,NA), cfi_robust=c(.95,.95), tli=c(.9,.9),
    tli_scaled=c(.8,.8), rmsea=c(.05,.05), rmsea_scaled=c(.06,.06), srmr=c(.05,.05))
  selected <- .pm_select_common_fit_variants(tab)
  expect_true(all(is.na(selected$values$cfi)))
  expect_identical(unname(selected$variants["tli"]), "scaled")
})

test_that("requested index variant controls both reader and recommendation inputs", {
  d <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  models <- list(one="f =~ x1+x2+x3+x4+x5+x6",
                 two="a =~ x1+x2+x3\nb =~ x4+x5+x6")
  for (variant in c("scaled", "robust", "base")) {
    x <- compare_models_auto(d, models, estimator="MLR", fit_variant=variant, verbose=FALSE)
    column <- if (variant == "base") "cfi" else paste0("cfi_", variant)
    expect_equal(x$fit$comparison_cfi, x$fit[[column]])
    reader <- .pm_analysis_sections(x, "eng")$sections$Fit
    expect_equal(reader[[paste0("CFI (", variant, ")")]], x$fit[[column]])
    expect_identical(x$effective_settings$fit_variant, variant)
  }
  expect_error(compare_models_auto(d, models, fit_variant="invalid", verbose=FALSE), "arg")
})
