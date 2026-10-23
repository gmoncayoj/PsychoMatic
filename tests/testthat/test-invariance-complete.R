invariance_model_6 <- paste0("F =~ ", paste0("i", 1:6, collapse = " + "))

make_binary_invariance_fixture <- function(n_per_group = 600L, seed = 333L) {
  d <- simulate_invariance_fixture(n_per_group, ordinal = FALSE, seed = seed)
  d[paste0("i", 1:6)] <- lapply(d[paste0("i", 1:6)], function(x) {
    ordered(as.integer(x > 0), levels = 0:1)
  })
  d[c("group", paste0("i", 1:6))]
}

make_mixed_invariance_fixture <- function(n_per_group = 180L, seed = 6102L) {
  d <- simulate_invariance_fixture(n_per_group, ordinal = FALSE, seed = seed)
  for (item in paste0("i", 1:3)) {
    d[[item]] <- ordered(
      cut(d[[item]], c(-Inf, -.7, 0, .7, Inf), labels = FALSE),
      levels = 1:4
    )
  }
  d[c("group", paste0("i", 1:6))]
}

expect_fit_records_use_effective_metadata <- function(out,
                                                      estimator_effective,
                                                      estimator_orig,
                                                      missing_effective,
                                                      parameterization_effective,
                                                      std_lv_effective,
                                                      identification_effective,
                                                      measurement_effective,
                                                      estimand_effective,
                                                      estimator_engine_request,
                                                      missing_engine_request,
                                                      parameterization_engine_request,
                                                      id_fac_engine_request,
                                                      id_cat_engine_request,
                                                      id_fac_effective,
                                                      id_cat_effective) {
  for (stage in names(out$fit_records)) {
    record <- out$fit_records[[stage]]
    metadata <- record$fit_metadata

    expect_identical(record$case_indices, metadata$case_indices, info = stage)
    expect_identical(record$variables, metadata$variables, info = stage)
    expect_identical(record$estimator, estimator_effective, info = stage)
    expect_identical(record$estimator_effective, estimator_effective, info = stage)
    expect_identical(record$estimator_orig, estimator_orig, info = stage)
    expect_identical(record$estimand, estimand_effective, info = stage)
    expect_identical(record$ordered, metadata$ordered, info = stage)
    expect_identical(record$missing, missing_effective, info = stage)
    expect_identical(record$identification, identification_effective, info = stage)
    expect_identical(record$parameterization, parameterization_effective, info = stage)
    expect_identical(record$std_lv, std_lv_effective, info = stage)
    expect_identical(record$group_coding, metadata$group_coding, info = stage)
    expect_identical(record$measurement, measurement_effective, info = stage)

    expect_identical(record$estimator_engine_request,
                     estimator_engine_request, info = stage)
    expect_identical(record$missing_engine_request,
                     missing_engine_request, info = stage)
    expect_identical(record$parameterization_engine_request,
                     parameterization_engine_request, info = stage)
    expect_identical(record$ID.fac_engine_request,
                     id_fac_engine_request, info = stage)
    expect_identical(record$ID.cat_engine_request,
                     id_cat_engine_request, info = stage)
    expect_identical(record$ID.fac, id_fac_engine_request, info = stage)
    expect_identical(record$ID.cat, id_cat_engine_request, info = stage)
    expect_identical(record$ID.fac_effective,
                     id_fac_effective, info = stage)
    expect_identical(record$ID.cat_effective,
                     id_cat_effective, info = stage)
  }
}

test_that("historical formals and positional calling remain compatible", {
  f <- formals(factorial_invariance_auto)
  expect_identical(
    head(names(f), 11L),
    c(
      "data", "group", "model", "language", "estimator", "ordered",
      "report", "criteria", "missing", "min_group_n", "overwrite"
    )
  )
  expect_identical(deparse(f$language), "c(\"eng\", \"esp\")")
  expect_identical(f$estimator, "auto")
  expect_null(f$ordered)
  expect_identical(f$report, FALSE)
  expect_identical(f$criteria, "auto")
  expect_identical(f$missing, "auto")
  expect_identical(f$min_group_n, 20L)
  expect_identical(f$overwrite, FALSE)

  d <- simulate_invariance_fixture(45L, seed = 6103L)
  d <- d[c("group", paste0("i", 1:6))]
  out <- factorial_invariance_auto(
    d, "group", invariance_model_6, "eng", "ML", FALSE, FALSE
  )
  expect_s3_class(out, "factorial_invariance_auto")
  expect_s3_class(out, "psychomatic_result")
  expect_identical(rownames(out$fit_table), c("configural", "metric", "scalar", "strict"))
  expect_type(out$conclusion, "character")
  expect_invisible(print(out))
})

test_that("continuous stages match direct cumulative lavaan fits", {
  d <- simulate_invariance_fixture(80L, seed = 6113L)
  d <- d[c("group", paste0("i", 1:6))]
  out <- factorial_invariance_auto(
    d, "group", invariance_model_6, estimator = "ML", ordered = FALSE
  )
  constraints <- list(
    configural = character(),
    metric = "loadings",
    scalar = c("loadings", "intercepts"),
    strict = c("loadings", "intercepts", "residuals")
  )
  for (stage in names(constraints)) {
    direct <- lavaan::cfa(
      invariance_model_6, data = d, group = "group", estimator = "ML",
      std.lv = FALSE, meanstructure = TRUE, missing = "fiml",
      group.equal = constraints[[stage]], group.partial = character()
    )
    expect_identical(out$fit_records[[stage]]$effective_constraints,
                     constraints[[stage]], info = stage)
    expect_equal(out$fit_table[stage, "df"],
                 unname(lavaan::fitMeasures(direct, "df")), info = stage)
    expect_equal(out$fit_table[stage, "cfi"],
                 unname(lavaan::fitMeasures(direct, "cfi")), tolerance = 1e-10,
                 info = stage)
  }
  expect_fit_records_use_effective_metadata(
    out,
    estimator_effective = "ML", estimator_orig = "ML",
    missing_effective = "fiml", parameterization_effective = "delta",
    std_lv_effective = FALSE, identification_effective = "marker",
    measurement_effective = "continuous",
    estimand_effective = "means_covariances",
    estimator_engine_request = "ML", missing_engine_request = "fiml",
    parameterization_engine_request = NA_character_,
    id_fac_engine_request = "marker", id_cat_engine_request = NA_character_,
    id_fac_effective = "marker", id_cat_effective = NA_character_
  )
})

test_that("Chen-informed defaults do not infer an unknown noninvariance pattern", {
  favourable <- PsychoMatic:::.ifa_resolve_criteria(
    total_n = 300L, group_ratio = 1.5, normal = TRUE,
    measurement = "continuous", criteria = "auto"
  )
  expect_equal(
    unlist(favourable[c("cfi", "rmsea", "srmr_loading", "srmr_other")]),
    c(cfi = .010, rmsea = .015, srmr_loading = .030, srmr_other = .010)
  )
  expect_identical(favourable$preset, "chen_informed_guidelines")

  cases <- list(
    list(301L, 1.5, TRUE),
    list(300L, 1.5001, TRUE),
    list(300L, 1.5, FALSE),
    list(300L, 1.5, NA)
  )
  for (case in cases) {
    adverse <- PsychoMatic:::.ifa_resolve_criteria(
      total_n = case[[1]], group_ratio = case[[2]], normal = case[[3]],
      measurement = "continuous", criteria = "auto"
    )
    expect_identical(adverse$preset, "chen_informed_guidelines")
    expect_equal(adverse$cfi, .010)
  }

  ordinal <- PsychoMatic:::.ifa_resolve_criteria(
    200L, 1, NA, "ordinal", "auto"
  )
  expect_identical(ordinal$preset, "categorical_descriptive_extrapolation")
  expect_error(
    PsychoMatic:::.ifa_resolve_criteria(
      200L, 1, TRUE, "continuous", c(cfi = .01, rmsea = .015)
    ),
    "fully named|srmr_loading"
  )
  expect_error(
    PsychoMatic:::.ifa_resolve_criteria(
      200L, 1, TRUE, "continuous",
      c(cfi = .01, rmsea = .015, srmr_loading = .03,
        srmr_other = .01, mystery = .2)
    ),
    "unknown"
  )

  custom <- PsychoMatic:::.ifa_resolve_criteria(
    200L, 1, TRUE, "continuous",
    c(cfi = .02, rmsea = .03, srmr_loading = .04, srmr_other = .01,
      configural_cfi = .95, configural_rmsea = .06, configural_srmr = .05)
  )
  expect_identical(custom$preset, "custom")
  expect_equal(custom$configural_cfi, .95)
  expect_false(PsychoMatic:::.ifa_configural_gate(
    c(cfi = .94, rmsea = .05, srmr = .04), custom
  )$passed)
  expect_error(
    PsychoMatic:::.ifa_resolve_criteria(
      200L, 1, TRUE, "continuous",
      c(cfi = -.01, rmsea = .03, srmr_loading = .04, srmr_other = .01)
    ),
    "nonnegative"
  )
  duplicated <- c(cfi = .01, rmsea = .02, srmr_loading = .03, srmr_other = .01)
  names(duplicated)[2] <- "cfi"
  expect_error(
    PsychoMatic:::.ifa_resolve_criteria(200L, 1, TRUE, "continuous", duplicated),
    "unique|fully named"
  )
})

test_that("Mardia boundary is inclusive and any group failure is adverse", {
  ok <- data.frame(
    group = c("g1", "g2"), skew_p = c(.05, .8),
    kurtosis_p = c(.05, .7), status = "ok", stringsAsFactors = FALSE
  )
  expect_true(PsychoMatic:::.ifa_normality_from_groups(ok)$normal)

  low <- ok
  low$skew_p[2] <- .0499
  expect_false(PsychoMatic:::.ifa_normality_from_groups(low)$normal)
  unknown <- ok
  unknown$kurtosis_p[1] <- NA_real_
  expect_true(is.na(PsychoMatic:::.ifa_normality_from_groups(unknown)$normal))
  errored <- ok
  errored$status[1] <- "error"
  expect_true(is.na(PsychoMatic:::.ifa_normality_from_groups(errored)$normal))
})

test_that("a real groupwise Mardia failure selects the adverse path without global fallback", {
  d <- simulate_invariance_fixture(30L, seed = 6117L)
  d <- d[c("group", paste0("i", 1:6))]
  d$i6[d$group == "g1"] <- d$i1[d$group == "g1"]
  out <- factorial_invariance_auto(
    d, "group", invariance_model_6,
    estimator = "auto", ordered = FALSE, min_group_n = 20L
  )
  expect_true(is.na(out$normality$overall_normal))
  expect_identical(out$criteria$preset, "chen_informed_guidelines")
  expect_identical(out$requested_settings$estimator, "auto")
  expect_identical(out$effective_settings$estimator, "MLR")
  failed <- out$normality$by_group[out$normality$by_group$status == "error", , drop = FALSE]
  expect_gt(nrow(failed), 0L)
  expect_true(any(vapply(
    out$warnings$condition,
    inherits,
    logical(1),
    what = "psychomatic_mardia_singular"
  )))
})

test_that("raw deterioration signs and inclusive decision boundaries are exact", {
  deltas <- PsychoMatic:::.ifa_deteriorations(
    previous = c(cfi = .96, rmsea = .04, srmr = .03),
    current = c(cfi = .955, rmsea = .05, srmr = .055)
  )
  expect_equal(deltas, c(delta_cfi = -.005, delta_rmsea = .01, delta_srmr = .025))
  improved <- PsychoMatic:::.ifa_deteriorations(
    c(cfi = .95, rmsea = .05, srmr = .05),
    c(cfi = .96, rmsea = .04, srmr = .03)
  )
  expect_gt(improved[["delta_cfi"]], 0)
  expect_true(all(improved[c("delta_rmsea", "delta_srmr")] < 0))

  cutoffs <- list(cfi = .005, rmsea = .010, srmr_loading = .025, srmr_other = .005)
  expect_identical(
    PsychoMatic:::.ifa_change_decision(-.005, .010, .025, cutoffs, "metric")$decision,
    "supported"
  )
  expect_identical(
    PsychoMatic:::.ifa_change_decision(-.0051, .0101, 0, cutoffs, "metric")$decision,
    "not_supported"
  )
  expect_identical(
    PsychoMatic:::.ifa_change_decision(-.0051, 0, .0251, cutoffs, "metric")$decision,
    "not_supported"
  )
  expect_identical(
    PsychoMatic:::.ifa_change_decision(-.0051, .0101, .0251, cutoffs, "scalar")$decision,
    "not_supported"
  )
  expect_identical(
    PsychoMatic:::.ifa_change_decision(-.0051, 0, .005, cutoffs, "scalar")$decision,
    "supported"
  )
})

test_that("configural fit uses an inclusive three-index gate", {
  expect_true(PsychoMatic:::.ifa_configural_gate(
    c(cfi = .90, rmsea = .08, srmr = .08)
  )$passed)
  bad <- PsychoMatic:::.ifa_configural_gate(c(cfi = .899, rmsea = .08, srmr = .08))
  expect_false(bad$passed)
  expect_identical(bad$status, "failed_configural_fit")
  expect_identical(bad$decision, "not_supported")
  expect_identical(
    PsychoMatic:::.ifa_configural_gate(c(cfi = NA, rmsea = .05, srmr = .04))$status,
    "not_evaluable"
  )
})

test_that("comparability rejects each independently frozen mismatch", {
  base <- list(
    case_indices = 1:20,
    variables = paste0("i", 1:6),
    estimator = "MLR",
    estimand = "means_covariances",
    ordered = character(),
    missing = "fiml",
    identification = "std.lv",
    parameterization = "theta",
    group_coding = c("g1", "g2"),
    measurement = "continuous",
    effective_constraints = "loadings",
    measures = list(
      cfi = c(robust = .95, scaled = .94, base = .93),
      rmsea = c(robust = .04, scaled = .05, base = .06),
      srmr = c(base = .03), df = c(scaled = 20, base = 20)
    ),
    admissible = TRUE
  )
  current <- base
  current$effective_constraints <- c("loadings", "intercepts")
  current$measures$df <- c(scaled = 24, base = 24)
  expect_true(PsychoMatic:::.ifa_compare_records(base, current)$comparable)

  mutations <- list(
    case_indices = 2:21,
    variables = rev(base$variables),
    estimator = "ML",
    estimand = "covariances",
    ordered = "i1",
    missing = "listwise",
    identification = "marker",
    parameterization = "delta",
    group_coding = rev(base$group_coding),
    measurement = "ordinal"
  )
  for (field in names(mutations)) {
    changed <- current
    changed[[field]] <- mutations[[field]]
    cmp <- PsychoMatic:::.ifa_compare_records(base, changed)
    expect_false(cmp$comparable, info = field)
    expect_match(cmp$reason, field, fixed = TRUE)
  }

  no_df <- current
  no_df$measures$df <- c(base = 20)
  expect_false(PsychoMatic:::.ifa_compare_records(base, no_df)$comparable)
  no_common <- current
  no_common$measures$cfi <- c(base = .94)
  base_robust <- base
  base_robust$measures$cfi <- c(robust = .95)
  expect_false(PsychoMatic:::.ifa_compare_records(base_robust, no_common)$comparable)

  non_nested <- current
  non_nested$effective_constraints <- "intercepts"
  expect_match(
    PsychoMatic:::.ifa_compare_records(base, non_nested)$reason,
    "constraints_not_nested"
  )
  bad_previous <- base
  bad_previous$admissible <- FALSE
  expect_match(PsychoMatic:::.ifa_compare_records(bad_previous, current)$reason,
               "previous_admissible")
  bad_current <- current
  bad_current$admissible <- FALSE
  expect_match(PsychoMatic:::.ifa_compare_records(base, bad_current)$reason,
               "current_admissible")

  selected <- PsychoMatic:::.ifa_select_variant(
    c(base = .91, scaled = .92, robust = .93)
  )
  expect_identical(selected$variant, "scaled")
  expect_equal(selected$value, .92)
  selected_scaled <- PsychoMatic:::.ifa_select_variant(
    c(base = .91, scaled = .92, robust = NA_real_)
  )
  expect_identical(selected_scaled$variant, "scaled")

  for (index in c("rmsea", "srmr")) {
    left <- base
    right <- current
    left$measures[[index]] <- c(robust = .04)
    right$measures[[index]] <- c(base = .04)
    cmp <- PsychoMatic:::.ifa_compare_records(left, right)
    expect_false(cmp$comparable, info = index)
    expect_match(cmp$reason, index, fixed = TRUE)
  }
})

test_that("effective lavaan metadata, not copied requests, controls comparability", {
  d <- simulate_invariance_fixture(35L, seed = 6122L)
  indicators <- paste0("i", 1:6)
  d <- d[c("group", indicators)]
  d[2L, indicators] <- NA_real_
  source_indices <- as.integer(seq_len(nrow(d)) * 10L + 7L)

  built <- PsychoMatic:::.ifa_build_syntax(
    model = invariance_model_6,
    data = d,
    group_name = "group",
    measurement = "continuous",
    ordered = character(),
    estimator = "MLR",
    missing = "fiml"
  )
  blueprint <- list(
    case_indices = source_indices,
    variables = rev(indicators),
    estimator = "copied_request_is_not_fit_metadata",
    estimand = "copied_estimand",
    ordered = "i1",
    missing = "listwise",
    identification = "marker",
    parameterization = "theta",
    group_coding = c("wrong_2", "wrong_1"),
    measurement = "ordinal"
  )
  record <- PsychoMatic:::.ifa_extract_fit_record(
    built, "configural", "configural", blueprint, indicators
  )

  expect_true(record$admissible)
  expect_identical(record$fit_metadata$case_indices, source_indices[-2L])
  expect_identical(record$fit_metadata$variables, indicators)
  expect_identical(record$fit_metadata$ordered, character())
  expect_identical(record$fit_metadata$group_coding, c("g1", "g2"))
  expect_identical(record$fit_metadata$estimator_effective, "ML")
  expect_identical(record$fit_metadata$estimator_orig, "MLR")
  expect_identical(record$fit_metadata$missing, "fiml")
  expect_identical(record$fit_metadata$parameterization, "delta")
  expect_identical(record$fit_metadata$std_lv, TRUE)
  expect_identical(record$fit_metadata$identification, "std.lv")
  expect_identical(record$fit_metadata$measurement, "continuous")
  expect_identical(record$fit_metadata$estimand, "means_covariances")

  previous <- record
  previous$effective_constraints <- character()
  current <- record
  current$effective_constraints <- "loadings"
  current$measures$df <- current$measures$df + 5
  current$case_indices <- rev(current$case_indices)
  current$variables <- rev(current$variables)
  current$estimator <- "WLSMV"
  current$estimand <- "thresholds_polychoric_correlations"
  current$ordered <- "i1"
  current$missing <- "listwise"
  current$identification <- "marker"
  current$parameterization <- "theta"
  current$group_coding <- rev(current$group_coding)
  current$measurement <- "ordinal"
  expect_true(PsychoMatic:::.ifa_compare_records(previous, current)$comparable)

  metadata_mutations <- list(
    case_indices = rev(record$fit_metadata$case_indices),
    variables = rev(indicators),
    estimator_effective = "DWLS",
    estimator_orig = "WLSMV",
    estimand = "thresholds_polychoric_correlations",
    ordered = "i1",
    missing = "listwise",
    identification = "marker",
    parameterization = "theta",
    std_lv = FALSE,
    group_coding = rev(record$fit_metadata$group_coding),
    measurement = "ordinal"
  )
  for (field in names(metadata_mutations)) {
    changed <- current
    changed$fit_metadata[[field]] <- metadata_mutations[[field]]
    comparison <- PsychoMatic:::.ifa_compare_records(previous, changed)
    expect_false(comparison$comparable, info = field)
    expect_match(comparison$reason, field, fixed = TRUE)
  }
})

test_that("public continuous settings report lavaan-effective values", {
  d <- simulate_invariance_fixture(35L, seed = 6124L)
  d <- d[c("group", paste0("i", 1:6))]
  out <- factorial_invariance_auto(
    d, "group", invariance_model_6,
    estimator = "ULS", ordered = FALSE, missing = "auto", report = FALSE
  )

  options <- lavaan::lavInspect(out$fitted_models$configural, "options")
  actual_missing <- PsychoMatic:::.ifa_canonical_missing(options$missing)
  actual_parameterization <- as.character(options$parameterization[[1L]])

  expect_identical(actual_missing, "two.stage")
  expect_identical(actual_parameterization, "delta")
  expect_identical(out$requested_settings$missing, "auto")
  expect_identical(out$effective_settings$missing_engine_request, "fiml")
  expect_identical(out$effective_settings$missing, actual_missing)
  expect_identical(
    out$effective_settings$parameterization_engine_request,
    NA_character_
  )
  expect_identical(
    out$effective_settings$parameterization,
    actual_parameterization
  )
  expect_identical(out$analysis_blueprint$missing_engine_request, "fiml")
  expect_identical(out$analysis_blueprint$missing, actual_missing)
  expect_identical(
    out$analysis_blueprint$parameterization_engine_request,
    NA_character_
  )
  expect_identical(
    out$analysis_blueprint$parameterization,
    actual_parameterization
  )
  expect_identical(
    out$tables$settings$value[out$tables$settings$setting == "missing"],
    actual_missing
  )
  expect_identical(
    out$tables$settings$value[
      out$tables$settings$setting == "parameterization"
    ],
    actual_parameterization
  )
  expect_true(all(out$fit_table$comparable[-1L]))
  expect_true(all(vapply(
    out$fit_records,
    function(record) identical(record$fit_metadata$missing, actual_missing) &&
      identical(record$fit_metadata$parameterization, actual_parameterization),
    logical(1)
  )))
  expect_fit_records_use_effective_metadata(
    out,
    estimator_effective = "ULS", estimator_orig = "ULS",
    missing_effective = "two.stage", parameterization_effective = "delta",
    std_lv_effective = FALSE, identification_effective = "marker",
    measurement_effective = "continuous",
    estimand_effective = "means_covariances",
    estimator_engine_request = "ULS", missing_engine_request = "fiml",
    parameterization_engine_request = NA_character_,
    id_fac_engine_request = "marker", id_cat_engine_request = NA_character_,
    id_fac_effective = "marker", id_cat_effective = NA_character_
  )
})

test_that("fit table values remain paired with the common comparison variant", {
  make_record <- function(stage, constraints, cfi, rmsea, srmr, df) {
    list(
      stage = stage,
      stage_detail = stage,
      requested_constraints = constraints,
      effective_constraints = constraints,
      measures = list(
        chi2 = c(base = 20),
        df = df,
        p = c(base = .20),
        cfi = cfi,
        rmsea = rmsea,
        srmr = srmr
      ),
      converged = TRUE,
      post_check = TRUE,
      admissible = TRUE,
      admissibility_reason = "admissible",
      case_indices = 1:40,
      variables = paste0("i", 1:6),
      estimator = "MLR",
      estimand = "means_covariances",
      ordered = character(),
      missing = "fiml",
      identification = "std.lv",
      parameterization = "delta",
      group_coding = c("g1", "g2"),
      measurement = "continuous"
    )
  }
  records <- list(
    configural = make_record(
      "configural", character(),
      c(robust = NA_real_, scaled = .95, base = .90),
      c(robust = NA_real_, scaled = .04, base = .08),
      c(robust = NA_real_, scaled = .03, base = .07),
      c(robust = NA_real_, scaled = 18, base = 17)
    ),
    metric = make_record(
      "metric", "loadings",
      c(robust = .99, scaled = .94, base = .80),
      c(robust = .01, scaled = .05, base = .09),
      c(robust = .01, scaled = .04, base = .08),
      c(robust = 99, scaled = 23, base = 22)
    )
  )
  criteria <- list(
    cfi = .01, rmsea = .015, srmr_loading = .03, srmr_other = .01,
    configural_cfi = .90, configural_rmsea = .08, configural_srmr = .08
  )

  table <- PsychoMatic:::.ifa_fit_table(records, criteria, "continuous")

  expect_identical(table["metric", "cfi_variant"], "scaled")
  expect_identical(table["metric", "rmsea_variant"], "scaled")
  expect_identical(table["metric", "srmr_variant"], "base")
  expect_identical(table["metric", "df_variant"], "scaled")
  expect_equal(table["metric", "cfi"], .94)
  expect_equal(table["metric", "rmsea"], .05)
  expect_equal(table["metric", "srmr"], .08)
  expect_equal(table["metric", "df"], 23)
  expect_equal(table["metric", "delta_cfi"], -.01)
  expect_equal(table["metric", "delta_rmsea"], .01)
  expect_equal(table["metric", "delta_srmr"], .01)
})

test_that("fit checks fail closed for convergence post-check Heywood and nonfinite evidence", {
  good <- list(
    converged = TRUE, post_check = TRUE, finite_parameters = TRUE,
    negative_observed_variance = FALSE, negative_latent_variance = FALSE,
    measures = list(cfi = c(base = .95), rmsea = c(base = .04),
                    srmr = c(base = .03), df = c(base = 20))
  )
  expect_true(PsychoMatic:::.ifa_admissibility(good)$admissible)
  for (field in c(
    "converged", "post_check", "finite_parameters",
    "negative_observed_variance", "negative_latent_variance"
  )) {
    bad <- good
    bad[[field]] <- field %in% c("negative_observed_variance", "negative_latent_variance")
    expect_false(PsychoMatic:::.ifa_admissibility(bad)$admissible, info = field)
  }
  bad_index <- good
  bad_index$measures$cfi <- c(base = NA_real_)
  expect_false(PsychoMatic:::.ifa_admissibility(bad_index)$admissible)
})

test_that("an inadmissible mandatory stage closes the sequential inference chain", {
  make_record <- function(stage, constraints, df, admissible = TRUE) {
    list(
      stage = stage,
      stage_detail = stage,
      requested_constraints = constraints,
      effective_constraints = constraints,
      measures = list(
        chi2 = c(base = 20 + df),
        df = c(base = df),
        p = c(base = .20),
        cfi = c(base = .96),
        rmsea = c(base = .04),
        srmr = c(base = .03)
      ),
      converged = admissible,
      post_check = admissible,
      admissible = admissible,
      admissibility_reason = if (admissible) "admissible" else "converged",
      case_indices = 1:40,
      variables = paste0("i", 1:6),
      estimator = "ML",
      estimand = "means_covariances",
      ordered = character(),
      missing = "fiml",
      identification = "std.lv",
      parameterization = "theta",
      group_coding = c("g1", "g2"),
      measurement = "continuous"
    )
  }
  records <- list(
    configural = make_record("configural", character(), 18),
    metric = make_record("metric", "loadings", 23, admissible = FALSE),
    scalar = make_record("scalar", c("loadings", "intercepts"), 28),
    strict = make_record(
      "strict", c("loadings", "intercepts", "residuals"), 34
    )
  )
  criteria <- list(
    cfi = .01, rmsea = .015, srmr_loading = .03, srmr_other = .01,
    configural_cfi = .90, configural_rmsea = .08, configural_srmr = .08
  )

  table <- PsychoMatic:::.ifa_fit_table(records, criteria, "continuous")

  expect_identical(table["metric", "decision"], "not_evaluable")
  expect_identical(table["scalar", "decision"], "not_evaluable")
  expect_identical(table["strict", "decision"], "not_evaluable")
  expect_match(table["scalar", "reason"], "less_restrictive")
  expect_false(table["scalar", "comparable"])
  expect_false(table["strict", "comparable"])
  expect_identical(table["scalar", "comparison_from"], "metric")
  expect_identical(table["strict", "comparison_from"], "scalar")
  expect_true(all(is.na(table[c("scalar", "strict"), c(
    "delta_df", "delta_cfi", "delta_rmsea", "delta_srmr"
  )])))
})

test_that("group name symbol index and direct vector preserve the same analytic rows", {
  d <- simulate_invariance_fixture(35L, seed = 6104L)
  d <- d[c("group", paste0("i", 1:6))]
  d$group <- factor(d$group, levels = c("g1", "g2", "unused"))
  d$i1[c(2, 39)] <- NA_real_
  group <- rep("collision", nrow(d))

  by_name <- factorial_invariance_auto(
    d, "group", invariance_model_6, estimator = "MLR", ordered = FALSE,
    missing = "listwise", report = FALSE
  )
  by_symbol <- factorial_invariance_auto(
    d, group, invariance_model_6, estimator = "MLR", ordered = FALSE,
    missing = "listwise", report = FALSE
  )
  by_index <- factorial_invariance_auto(
    d, 1L, invariance_model_6, estimator = "MLR", ordered = FALSE,
    missing = "listwise", report = FALSE
  )
  direct <- d$group
  by_vector <- factorial_invariance_auto(
    d, direct, invariance_model_6, estimator = "MLR", ordered = FALSE,
    missing = "listwise", report = FALSE
  )
  expect_identical(by_name$analysis_blueprint$case_indices, by_symbol$analysis_blueprint$case_indices)
  expect_identical(by_name$analysis_blueprint$case_indices, by_index$analysis_blueprint$case_indices)
  expect_identical(by_name$analysis_blueprint$case_indices, by_vector$analysis_blueprint$case_indices)
  expect_identical(
    by_name$analysis_blueprint$case_indices,
    as.integer(setdiff(seq_len(nrow(d)), c(2L, 39L)))
  )
  expect_identical(by_name$sample_info$dropped_missing_items, 2L)
  expect_identical(by_vector$sample_info$dropped_missing_items, 2L)
  expect_identical(by_name$analysis_blueprint$group_coding, c("g1", "g2"))
  expect_false(identical(by_vector$analysis_blueprint$group_name, "group"))
  expect_identical(as.character(by_vector$analysis_data$group), as.character(d$group[-c(2L, 39L)]))
})

test_that("a positive group index disambiguates duplicated column names", {
  d <- simulate_invariance_fixture(25L, seed = 6121L)
  duplicated_group <- data.frame(
    decoy = rep("one_group", nrow(d)),
    selected_group = d$group,
    d[paste0("i", 1:6)],
    check.names = FALSE
  )
  names(duplicated_group)[1:2] <- c("group", "group")

  out <- factorial_invariance_auto(
    duplicated_group, 2L, invariance_model_6,
    estimator = "ML", ordered = FALSE, report = FALSE
  )

  expect_identical(out$analysis_blueprint$group_source, "column_index")
  expect_identical(out$analysis_blueprint$group_coding, c("g1", "g2"))
  expect_equal(sum(names(out$analysis_data) == out$analysis_blueprint$group_name), 1L)
  expect_identical(
    as.character(out$analysis_data[[out$analysis_blueprint$group_name]]),
    as.character(d$group)
  )
})

test_that("missing and invalid groups fail before fitting with counts", {
  d <- simulate_invariance_fixture(25L, seed = 6105L)
  d <- d[c("group", paste0("i", 1:6))]
  d$group[1] <- NA
  out <- factorial_invariance_auto(
    d, "group", invariance_model_6, estimator = "ML", ordered = FALSE,
    missing = "listwise", min_group_n = 20L
  )
  expect_identical(out$sample_info$dropped_missing_group, 1L)
  expect_error(
    factorial_invariance_auto(transform(d, group = "g1"), "group", invariance_model_6),
    "at least two"
  )
  d_small <- d
  d_small$group <- factor(c(rep("g1", 10), rep("g2", nrow(d_small) - 10)))
  expect_error(
    factorial_invariance_auto(d_small, "group", invariance_model_6, min_group_n = 20L),
    "g1=10|min_group_n"
  )
})

test_that("item filtering removes empty group levels and rebuilds analytic coding", {
  d <- simulate_invariance_fixture(30L, seed = 6123L)
  indicators <- paste0("i", 1:6)
  d <- d[c("group", indicators)]
  d$group <- factor(rep(c("g1", "g2", "g3"), each = 20L),
                    levels = c("g1", "g2", "g3", "unused"))
  d[d$group == "g3", indicators] <- NA_real_

  out <- factorial_invariance_auto(
    d, "group", invariance_model_6,
    estimator = "ML", ordered = FALSE, missing = "listwise",
    min_group_n = 20L, report = FALSE
  )

  expect_identical(levels(out$analysis_data$group), c("g1", "g2"))
  expect_identical(out$analysis_blueprint$group_coding, c("g1", "g2"))
  expect_identical(names(out$sample_info$group_sizes), c("g1", "g2"))
  expect_identical(out$sample_info$dropped_missing_items, 20L)

  two_groups <- d[d$group %in% c("g1", "g2"), , drop = FALSE]
  two_groups[two_groups$group == "g2", indicators] <- NA_real_
  expect_error(
    factorial_invariance_auto(
      two_groups, "group", invariance_model_6,
      estimator = "ML", ordered = FALSE, missing = "listwise",
      min_group_n = 2L
    ),
    "analytic sample.*at least two groups"
  )
})

test_that("ambiguous nonfinite constant and group-degenerate indicators are rejected", {
  d <- simulate_invariance_fixture(25L, seed = 6106L)
  d <- d[c("group", paste0("i", 1:6))]
  expect_error(
    factorial_invariance_auto(transform(d, i1 = as.character(round(i1))), "group", invariance_model_6),
    "character|explicit levels"
  )
  expect_error(
    factorial_invariance_auto(transform(d, i1 = factor(round(i1))), "group", invariance_model_6, ordered = TRUE),
    "unordered factor|explicit levels"
  )
  expect_error(
    factorial_invariance_auto(transform(d, i1 = Inf), "group", invariance_model_6),
    "finite|Inf"
  )
  expect_error(
    factorial_invariance_auto(transform(d, i1 = 1), "group", invariance_model_6),
    "constant|two observed"
  )
  d_degenerate <- d
  d_degenerate$i1[d_degenerate$group == "g1"] <- 1
  expect_error(
    factorial_invariance_auto(d_degenerate, "group", invariance_model_6),
    "g1|group"
  )

  all_missing <- d
  all_missing$i1 <- NA_real_
  expect_error(
    factorial_invariance_auto(all_missing, "group", invariance_model_6),
    "two observed|eligible"
  )
  with_nan <- d
  with_nan$i1[1] <- NaN
  expect_error(
    factorial_invariance_auto(with_nan, "group", invariance_model_6),
    "finite|eligible"
  )

  ordinal <- simulate_invariance_fixture(25L, ordinal = TRUE, seed = 6114L)
  ordinal <- ordinal[c("group", paste0("i", 1:6))]
  ordinal$i1[ordinal$group == "g1"] <- ordered(
    rep(1L, sum(ordinal$group == "g1")), levels = 1:5
  )
  expect_error(
    factorial_invariance_auto(ordinal, "group", invariance_model_6, ordered = TRUE),
    "group-degenerate|g1"
  )
})

test_that("ordered FALSE rejects semantic categorical encodings without code coercion", {
  ordered_data <- data.frame(
    item = ordered(rep(c("low", "mid", "high"), 4L),
                   levels = c("low", "mid", "high"))
  )
  factor_data <- data.frame(
    item = factor(rep(c("low", "mid", "high"), 4L),
                  levels = c("low", "mid", "high"))
  )
  logical_data <- data.frame(item = rep(c(FALSE, TRUE), 6L))

  expect_error(
    PsychoMatic:::.ifa_prepare_items(ordered_data, "item", ordered = FALSE),
    "ordered = FALSE|categorical.*continuous|numeric codes"
  )
  expect_error(
    PsychoMatic:::.ifa_prepare_items(factor_data, "item", ordered = FALSE),
    "unordered factor|explicit levels|ordered = FALSE"
  )
  expect_error(
    PsychoMatic:::.ifa_prepare_items(logical_data, "item", ordered = FALSE),
    "ordered = FALSE|categorical.*continuous|numeric codes"
  )

  numeric_codes <- data.frame(item = rep(1:3, 4L))
  numeric_result <- PsychoMatic:::.ifa_prepare_items(
    numeric_codes, "item", ordered = FALSE
  )
  expect_identical(numeric_result$measurement, "continuous")
  expect_identical(numeric_result$ordered, character())
  expect_identical(numeric_result$data$item, as.numeric(numeric_codes$item))
})

test_that("matrix spec and ordered item indices use the same public workflow", {
  d <- simulate_invariance_fixture(25L, seed = 6115L)
  d <- d[c("group", paste0("i", 1:6))]
  matrix_data <- data.matrix(d)
  matrix_out <- factorial_invariance_auto(
    matrix_data, 1L, invariance_model_6,
    estimator = "MLR", ordered = FALSE
  )
  expect_s3_class(matrix_out, "factorial_invariance_auto")
  expect_identical(matrix_out$analysis_blueprint$variables, paste0("i", 1:6))

  spec <- psychomatic_spec(
    d, items = paste0("i", 1:6), group = "group", level = "continuous"
  )
  spec_out <- factorial_invariance_auto(
    spec, "group", invariance_model_6,
    estimator = "MLR", ordered = FALSE
  )
  expect_s3_class(spec_out, "factorial_invariance_auto")
  expect_identical(spec_out$analysis_blueprint$variables, paste0("i", 1:6))

  ordinal <- simulate_invariance_fixture(25L, ordinal = TRUE, seed = 6116L)
  ordinal <- ordinal[c("group", paste0("i", 1:6))]
  prepared <- PsychoMatic:::.ifa_prepare_items(
    ordinal, paste0("i", 1:6), ordered = 2:7
  )
  expect_identical(prepared$ordered, paste0("i", 1:6))
  expect_identical(prepared$measurement, "ordinal")
})

test_that("private categorical estimator and missing guards remain available to partial invariance", {
  continuous <- simulate_invariance_fixture(25L, ordinal = FALSE, seed = 6107L)
  continuous <- continuous[c("group", paste0("i", 1:6))]
  ordinal <- simulate_invariance_fixture(25L, ordinal = TRUE, seed = 6108L)
  ordinal <- ordinal[c("group", paste0("i", 1:6))]

  expect_error(
    legacy_ordinal_identification(continuous, "group", invariance_model_6,
                              estimator = "WLSMV", ordered = FALSE),
    "ordered = FALSE|WLSMV"
  )
  expect_error(
    legacy_ordinal_identification(ordinal, "group", invariance_model_6,
                              estimator = "MLR", ordered = TRUE),
    "ordered indicators|MLR"
  )
  expect_error(
    legacy_ordinal_identification(ordinal, "group", invariance_model_6,
                              estimator = "WLSMV", ordered = TRUE, missing = "fiml"),
    "FIML|listwise"
  )
  expect_error(
    legacy_ordinal_identification(continuous, "group", invariance_model_6,
                              estimator = "ULSMV", ordered = FALSE),
    "ULSMV|continuous"
  )
  expect_error(
    legacy_ordinal_identification(ordinal, "group", invariance_model_6,
                              estimator = "ULS", ordered = TRUE),
    "ULS|ordered indicators"
  )
})

test_that("continuous invariance works without semTools while categorical requires it", {
  absent <- function(package, quietly = TRUE) !identical(package, "semTools")
  expect_invisible(PsychoMatic:::.ifa_validate_packages("continuous", absent))
  expect_error(
    PsychoMatic:::.ifa_validate_packages("ordinal", absent),
    "semTools"
  )
})

test_that("private polytomous primitives retain native semTools identification", {
  skip_if_not_installed("semTools")
  d <- simulate_invariance_fixture(180L, ordinal = TRUE, seed = 6109L)
  d <- d[c("group", paste0("i", 1:6))]
  out <- legacy_ordinal_identification(
    d, "group", invariance_model_6, ordered = TRUE,
    estimator = "auto", report = FALSE
  )

  expect_identical(out$measurement, "ordinal")
  expect_identical(
    out$fit_table$level,
    c("configural", "thresholds", "metric", "scalar", "strict")
  )
  expect_identical(out$effective_settings$parameterization, "theta")
  expect_identical(out$effective_settings$ID.fac, "std.lv")
  expect_identical(out$effective_settings$ID.cat, "Wu.Estabrook.2016")
  expect_identical(out$requested_settings$estimator, "auto")
  expect_identical(out$effective_settings$estimator, "WLSMV")
  expect_identical(out$effective_settings$estimator_engine, "DWLS")
  expect_identical(out$effective_settings$estimator_orig, "WLSMV")
  expect_identical(out$effective_settings$ID.fac_effective, "std.lv")
  expect_identical(out$effective_settings$ID.cat_effective, "wu")
  expect_identical(out$effective_settings$ID.fac_engine_request, "std.lv")
  expect_identical(
    out$effective_settings$ID.cat_engine_request,
    "Wu.Estabrook.2016"
  )
  expect_identical(out$effective_settings$identification, "std.lv")
  expect_identical(out$effective_settings$ordered, paste0("i", 1:6))
  expect_identical(
    out$fit_records$scalar$engine_args$group.equal,
    c("thresholds", "loadings", "intercepts")
  )
  expect_true(isTRUE(out$fit_records$scalar$engine_args$return.fit))
  expect_identical(out$fit_records$scalar$engine_args$group.partial, "")
  expect_identical(out$fit_records$scalar$engine_args$ID.fac, "std.lv")
  expect_identical(out$fit_records$scalar$engine_args$ID.cat, "Wu.Estabrook.2016")

  direct <- semTools::measEq.syntax(
    configural.model = invariance_model_6,
    data = out$analysis_data,
    group = out$analysis_blueprint$group_name,
    ordered = paste0("i", 1:6), estimator = "WLSMV",
    parameterization = "theta", ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016",
    group.equal = c("thresholds", "loadings", "intercepts"),
    group.partial = "", return.fit = TRUE
  )
  expect_equal(
    unname(lavaan::fitMeasures(out$fitted_models$scalar, "df")),
    unname(lavaan::fitMeasures(direct, "df"))
  )
  expect_equal(
    unname(lavaan::fitMeasures(out$fitted_models$scalar, "cfi.scaled")),
    unname(lavaan::fitMeasures(direct, "cfi.scaled")), tolerance = 1e-8
  )
  expect_fit_records_use_effective_metadata(
    out,
    estimator_effective = "DWLS", estimator_orig = "WLSMV",
    missing_effective = "listwise", parameterization_effective = "theta",
    std_lv_effective = TRUE, identification_effective = "std.lv",
    measurement_effective = "ordinal",
    estimand_effective = "thresholds_polychoric_correlations",
    estimator_engine_request = "WLSMV", missing_engine_request = "listwise",
    parameterization_engine_request = "theta",
    id_fac_engine_request = "std.lv",
    id_cat_engine_request = "Wu.Estabrook.2016",
    id_fac_effective = "std.lv", id_cat_effective = "wu"
  )
})

test_that("private binary primitives retain the joint scalar stage for partial invariance", {
  skip_if_not_installed("semTools")
  d <- make_binary_invariance_fixture()
  out <- legacy_ordinal_identification(
    d, "group", invariance_model_6, ordered = TRUE, report = FALSE
  )
  expect_identical(out$measurement, "ordinal")
  expect_identical(out$item_info$category_profile, "binary")
  expect_identical(out$fit_table$level, c("configural", "scalar", "strict"))
  expect_identical(
    out$fit_table$stage_detail[2],
    "binary_thresholds_loadings_intercepts_joint"
  )
  expect_identical(
    out$fit_records$scalar$requested_constraints,
    c("thresholds", "loadings", "intercepts")
  )
  expect_gt(out$fit_table["scalar", "delta_df"], 0)
  expect_false(any(out$fit_table$level %in% c("thresholds", "metric")))

  direct_scalar <- semTools::measEq.syntax(
    configural.model = invariance_model_6,
    data = out$analysis_data,
    group = out$analysis_blueprint$group_name,
    ordered = paste0("i", 1:6), estimator = "WLSMV",
    parameterization = "theta", ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016",
    group.equal = c("thresholds", "loadings", "intercepts"),
    group.partial = "", return.fit = TRUE
  )
  expect_equal(out$fit_table["scalar", "df"],
               unname(lavaan::fitMeasures(direct_scalar, "df")))
})

test_that("ordinal stage selection uses categories in the analytic sample", {
  skip_if_not_installed("semTools")
  d <- make_binary_invariance_fixture(n_per_group = 300L, seed = 6125L)
  rare <- unlist(lapply(split(seq_len(nrow(d)), d$group), head, 8L), use.names = FALSE)
  for (item in paste0("i", 1:6)) {
    values <- as.integer(as.character(d[[item]])) + 1L
    values[rare] <- 3L
    d[[item]] <- ordered(values, levels = 1:3)
  }
  d$i1[rare] <- NA

  before <- PsychoMatic:::.ifa_prepare_items(
    d, paste0("i", 1:6), ordered = TRUE
  )
  expect_identical(before$category_profile, "polytomous")
  expect_true(any(before$category_counts == 3L))

  out <- legacy_ordinal_identification(
    d, "group", invariance_model_6,
    ordered = TRUE, missing = "listwise", report = FALSE
  )

  expect_identical(out$sample_info$dropped_missing_items, 16L)
  expect_identical(out$item_info$category_counts, stats::setNames(
    rep(2L, 6L), paste0("i", 1:6)
  ))
  expect_identical(out$item_info$category_profile, "binary")
  expect_identical(
    vapply(out$analysis_data[paste0("i", 1:6)], nlevels, integer(1)),
    stats::setNames(rep(2L, 6L), paste0("i", 1:6))
  )
  expect_identical(out$fit_table$level, c("configural", "scalar", "strict"))
  expect_false(any(out$fit_table$level %in% c("thresholds", "metric")))
})

test_that("mixed sequence remains categorical and does not collapse to continuous", {
  skip_if_not_installed("semTools")
  d <- make_mixed_invariance_fixture()
  out <- legacy_ordinal_identification(
    d, "group", invariance_model_6,
    ordered = paste0("i", 1:3), report = FALSE
  )
  expect_identical(out$measurement, "mixed")
  expect_identical(out$requested_settings$estimator, "auto")
  expect_identical(out$effective_settings$estimator, "WLSMV")
  expect_identical(out$effective_settings$parameterization, "theta")
  expect_true("thresholds" %in% out$fit_table$level)
  expect_identical(
    out$fit_records$metric$effective_constraints,
    c("thresholds", "loadings")
  )
  expect_identical(out$fit_table["metric", "comparison_from"],
                   if (out$fit_table["thresholds", "status"] == "skipped_not_identified") {
                     "configural"
                   } else {
                     "thresholds"
                   })
  expect_true(out$fit_table$status[out$fit_table$level == "thresholds"] %in%
                c("ok", "skipped_not_identified"))

  binary_mixed <- d
  binary_mixed$i1 <- ordered(as.integer(binary_mixed$i1 > 2), levels = 0:1)
  prepared <- PsychoMatic:::.ifa_prepare_items(
    binary_mixed, paste0("i", 1:6), ordered = "i1"
  )
  expect_identical(prepared$measurement, "mixed")
  expect_identical(prepared$category_profile, "mixed")
})

test_that("constraint families retain one canonical row per group", {
  skip_if_not_installed("semTools")
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 6110L)
  d <- d[c("group", paste0("i", 1:6))]
  out <- legacy_ordinal_identification(d, "group", invariance_model_6, ordered = TRUE)
  families <- out$constraint_families
  target <- families[
    families$stage == "scalar" & families$canonical == "F =~ i2",
    , drop = FALSE
  ]
  expect_identical(sort(unique(target$group_label)), c("g1", "g2"))
  expect_equal(nrow(target), 2L)
  expect_true(all(target$family == "loadings"))
  expect_true(all(target$equality_constrained))
  expect_true(all(nzchar(target$label)))
  expect_length(unique(target$label), 1L)

  for (family in c("thresholds", "loadings", "intercepts")) {
    stage_rows <- families[families$stage == "scalar" & families$family == family,
                           , drop = FALSE]
    expect_gt(nrow(stage_rows), 0L)
    split_rows <- split(stage_rows, stage_rows$canonical)
    expect_true(all(vapply(split_rows, nrow, integer(1)) == 2L), info = family)
    expect_true(all(vapply(
      split_rows,
      function(rows) all(nzchar(rows$label)) && length(unique(rows$label)) == 1L,
      logical(1)
    )), info = family)
  }
  residual_rows <- families[
    families$stage == "strict" & families$family == "residuals",
    , drop = FALSE
  ]
  expect_gt(nrow(residual_rows), 0L)
  expect_true(all(residual_rows$equality_constrained))
})

test_that("result envelope methods and reporting preserve raw results", {
  d <- simulate_invariance_fixture(35L, seed = 6111L)
  d <- d[c("group", paste0("i", 1:6))]
  out <- factorial_invariance_auto(
    d, "group", invariance_model_6,
    estimator = "MLR", ordered = FALSE, report = FALSE
  )
  expect_s3_class(out, "psychomatic_result")
  expect_true(is.data.frame(out$warnings))
  expect_true(is.data.frame(out$errors))
  expect_true(is.character(out$warnings_list))
  expect_true(all(c(
    "requested_settings", "effective_settings", "data_summary", "results",
    "warnings", "errors", "decision_log", "diagnostics", "conclusion"
  ) %in% names(out)))
  expect_true(all(c(
    "requested_constraints", "effective_constraints", "delta_srmr",
    "post_check", "admissible", "comparable", "reason"
  ) %in% names(out$fit_table)))
  legacy_fields <- c(
    "estimator", "estimator_info", "estimator_reason", "indicators",
    "item_info", "sample_info", "normality", "requested_estimator",
    "requested_ordered", "chen_scenario", "fitted_models", "display_table",
    "decisions", "notes", "call", "tables", "report_text", "report_path"
  )
  expect_true(all(legacy_fields %in% names(out)))
  expect_identical(out$notes, out$warnings_list)
  expect_identical(out$estimator, out$effective_settings$estimator)

  integrated_tables <- PsychoMatic:::.pm_object_tables(out)
  expect_true(all(c(
    "fit", "decisions", "diagnostics", "criteria", "constraints",
    "warnings", "errors", "conclusion", "settings", "references"
  ) %in% names(integrated_tables)))
  expect_identical(integrated_tables$constraints, out$constraint_families)
  expect_identical(integrated_tables$conclusion$conclusion, out$conclusion)

  before <- serialize(out, NULL)
  smry <- summary(out)
  expect_true(all(c("fit", "decisions", "diagnostics", "criteria",
                    "conditions", "conclusion") %in% names(smry)))
  printed <- capture.output(value <- print(out))
  expect_identical(value, out)
  expect_true(length(printed) > 0L)
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit({ grDevices::dev.off(); unlink(pdf_file) }, add = TRUE)
  plotting <- plot(out)
  expect_true(is.data.frame(plotting))
  expect_identical(serialize(out, NULL), before)

  skipped <- out
  skipped$fit_table$status[2] <- "skipped_not_identified"
  skipped$fit_table[2, c("delta_df", "delta_cfi", "delta_rmsea", "delta_srmr")] <- NA_real_
  expect_invisible(print(skipped))
  expect_true(is.list(summary(skipped)))
  skipped_pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(skipped_pdf)
  skipped_plot <- plot(skipped)
  grDevices::dev.off()
  unlink(skipped_pdf)
  expect_true(is.data.frame(skipped_plot))

  report_file <- tempfile(fileext = ".xlsx")
  file.create(report_file)
  on.exit(unlink(report_file), add = TRUE)
  expect_error(
    factorial_invariance_auto(
      d, "group", invariance_model_6, estimator = "ML", ordered = FALSE,
      report = report_file, overwrite = FALSE
    ),
    "already exists"
  )

  skip_if_not_installed("openxlsx")
  successful_report <- tempfile(fileext = ".xlsx")
  on.exit(unlink(successful_report), add = TRUE)
  reported <- factorial_invariance_auto(
    d, "group", invariance_model_6, estimator = "ML", ordered = FALSE,
    report = successful_report, overwrite = FALSE
  )
  expect_true(file.exists(successful_report))
  expect_identical(normalizePath(reported$report_path, winslash = "/"),
                   normalizePath(successful_report, winslash = "/"))
})

test_that("verbose FALSE is silent and display rounding does not alter raw values", {
  d <- simulate_invariance_fixture(35L, seed = 6112L)
  d <- d[c("group", paste0("i", 1:6))]
  output <- capture.output(
    out <- factorial_invariance_auto(
      d, "group", invariance_model_6,
      estimator = "MLR", ordered = FALSE, report = FALSE, verbose = FALSE
    )
  )
  expect_length(output, 0L)
  numeric_columns <- names(out$fit_table)[vapply(out$fit_table, is.numeric, logical(1))]
  expect_true(any(vapply(
    numeric_columns,
    function(name) any(
      is.finite(out$fit_table[[name]]) &
        out$fit_table[[name]] != round(out$fit_table[[name]], 3)
    ),
    logical(1)
  )))
})
