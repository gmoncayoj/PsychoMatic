partial_model_6 <- paste0("f =~ ", paste0("i", 1:6, collapse = " + "))

partial_fixture <- function(n_per_group = 90L, ordinal = FALSE, seed = 7301L) {
  d <- simulate_invariance_fixture(
    n_per_group = n_per_group,
    ordinal = ordinal,
    partial = TRUE,
    seed = seed
  )
  d[c("group", paste0("i", 1:6))]
}

test_that("partial-invariance public signature is exact", {
  expect_identical(
    head(names(formals(invariance_partial)), 13L),
    c(
      "data", "group", "model", "free", "rationale", "stage",
      "measurement", "estimator", "ordered", "criteria",
      "min_anchor_items", "report", "language"
    )
  )
  expect_identical(eval(formals(invariance_partial)$stage), "scalar")
  expect_identical(eval(formals(invariance_partial)$measurement), "auto")
  expect_identical(eval(formals(invariance_partial)$estimator), "auto")
  expect_null(eval(formals(invariance_partial)$ordered))
  expect_identical(eval(formals(invariance_partial)$criteria), "auto")
  expect_identical(eval(formals(invariance_partial)$min_anchor_items), 2L)
  expect_false(eval(formals(invariance_partial)$report))
  expect_identical(eval(formals(invariance_partial)$language), "eng")
})

test_that("free grammar canonicalizes exact supported families", {
  parsed <- PsychoMatic:::.ipi_parse_free(c(
    " f =~ i1 ", "i2~1", "i3 | t2", " i4 ~~ i4"
  ))
  expect_identical(
    parsed$canonical,
    c("f =~ i1", "i2 ~ 1", "i3 | t2", "i4 ~~ i4")
  )
  expect_identical(
    parsed$family,
    c("loadings", "intercepts", "thresholds", "residuals")
  )
  expect_identical(parsed$input_order, 1:4)

  listed <- PsychoMatic:::.ipi_parse_free(list(
    loadings = "f =~ i1",
    thresholds = c("i2 | t1", "i2 | t2")
  ))
  expect_identical(listed$family, c("loadings", "thresholds", "thresholds"))

  expect_error(PsychoMatic:::.ipi_parse_free(c("f =~ i1", " f=~i1 ")), "duplicate")
  expect_error(PsychoMatic:::.ipi_parse_free("i1 ~ 0"), "unsupported|grammar")
  expect_error(PsychoMatic:::.ipi_parse_free("i1 ~~ i2"), "identical|residual")
  expect_error(PsychoMatic:::.ipi_parse_free("i1 | t0"), "unsupported|grammar")
  expect_error(PsychoMatic:::.ipi_parse_free(list(metric = "f =~ i1")), "family|names")
  expect_error(
    PsychoMatic:::.ipi_parse_free(list(intercepts = "f =~ i1")),
    "family"
  )
  expect_error(PsychoMatic:::.ipi_parse_free(c("", NA_character_)), "non-empty|NA")
})

test_that("rationale is required and maps exactly to canonical releases", {
  parsed <- PsychoMatic:::.ipi_parse_free(c("f =~ i1", "i2 ~ 1"))
  global <- PsychoMatic:::.ipi_normalize_rationale("Theory and prior evidence", parsed)
  expect_identical(global$canonical, parsed$canonical)
  expect_identical(global$rationale, rep("Theory and prior evidence", 2L))

  named <- "Item wording"
  names(named) <- " i2~1 "
  named <- c(named, stats::setNames("Prior loading evidence", "f=~i1"))
  mapped <- PsychoMatic:::.ipi_normalize_rationale(named, parsed)
  expect_identical(mapped$rationale, c("Prior loading evidence", "Item wording"))

  for (bad in list(NULL, character(), "", NA_character_)) {
    expect_error(PsychoMatic:::.ipi_normalize_rationale(bad, parsed), "rationale")
  }
  expect_error(
    PsychoMatic:::.ipi_normalize_rationale(c(a = "x", b = "y"), parsed),
    "names|identity|match|Unsupported"
  )
})

test_that("stage resolution respects continuous, ordinal, and binary identification", {
  continuous <- PsychoMatic:::.ipi_resolve_stage("metric", "continuous", "continuous")
  expect_identical(continuous$stage_key, "metric")
  expect_identical(continuous$stage_effective, "metric")
  expect_identical(continuous$prerequisite_key, "configural")

  ordinal <- PsychoMatic:::.ipi_resolve_stage("loadings", "ordinal", "polytomous")
  expect_identical(ordinal$stage_key, "metric")
  expect_identical(ordinal$stage_effective, "loadings")
  alias <- PsychoMatic:::.ipi_resolve_stage("metric", "ordinal", "polytomous")
  expect_identical(alias$stage_requested, "metric")
  expect_identical(alias$stage_effective, "loadings")

  expect_error(
    PsychoMatic:::.ipi_resolve_stage("loadings", "ordinal", "binary"),
    "binary|joint"
  )
  expect_error(
    PsychoMatic:::.ipi_resolve_stage("thresholds", "continuous", "continuous"),
    "stage"
  )
})

test_that("anchors honor factor-specific loadings and item-wide releases", {
  membership <- data.frame(
    factor = c("f1", "f1", "f1", "f2", "f2", "f2"),
    item = c("i1", "i2", "cross", "i3", "i4", "cross"),
    stringsAsFactors = FALSE
  )
  loading <- PsychoMatic:::.ipi_parse_free("f1 =~ cross")
  anchors <- PsychoMatic:::.ipi_anchor_table(
    membership, loading, c("loadings", "intercepts")
  )
  expect_false(anchors$retained[anchors$factor == "f1" & anchors$item == "cross"])
  expect_true(anchors$retained[anchors$factor == "f2" & anchors$item == "cross"])

  intercept <- PsychoMatic:::.ipi_parse_free("cross ~ 1")
  anchors2 <- PsychoMatic:::.ipi_anchor_table(
    membership, intercept, c("loadings", "intercepts")
  )
  expect_true(all(!anchors2$retained[anchors2$item == "cross"]))
})

test_that("partial nesting gate requires exact lower-stage release propagation", {
  parsed <- PsychoMatic:::.ipi_parse_free(c("f =~ i3", "i2 ~ 1"))
  ok <- PsychoMatic:::.ipi_partial_nesting_gate(
    prerequisite_releases = "f =~ i3",
    partial_releases = c("f =~ i3", "i2 ~ 1"),
    requested = parsed,
    prerequisite_constraints = "loadings"
  )
  expect_true(ok$comparable)
  bad <- PsychoMatic:::.ipi_partial_nesting_gate(
    prerequisite_releases = character(),
    partial_releases = c("f =~ i3", "i2 ~ 1"),
    requested = parsed,
    prerequisite_constraints = "loadings"
  )
  expect_false(bad$comparable)
  expect_identical(bad$reason, "mismatch_group_partial_nesting")
})

test_that("fixed reference rows are accepted but ignored releases fail closed", {
  make_record <- function(labels, free) {
    list(parameter_table = data.frame(
      lhs = c("i2", "i2"), op = c("~1", "~1"), rhs = c("", ""),
      group = 1:2, free = free, label = labels,
      plabel = c(".p1.", ".p2."), stringsAsFactors = FALSE
    ))
  }
  parsed <- PsychoMatic:::.ipi_parse_free("i2 ~ 1")
  full <- make_record(c("nu.2", "nu.2"), c(0L, 0L))
  expect_no_error(PsychoMatic:::.ipi_validate_complete_families(
    parsed, full, c("g1", "g2"), "intercepts"
  ))

  released <- make_record(c("nu.2.g1", "nu.2.g2"), c(0L, 1L))
  verified <- PsychoMatic:::.ipi_verify_releases(
    full, released, parsed, c("g1", "g2")
  )
  expect_true(verified$effective_release)

  blank_independent <- make_record(c("", ""), c(1L, 2L))
  verified_blank <- PsychoMatic:::.ipi_verify_releases(
    full, blank_independent, parsed, c("g1", "g2")
  )
  expect_true(verified_blank$effective_release)

  ignored <- PsychoMatic:::.ipi_verify_releases(
    full, full, parsed, c("g1", "g2")
  )
  expect_false(ignored$effective_release)
  expect_identical(ignored$release_reason, "group_partial_was_ignored")

  fixed_blank <- make_record(c("", ""), c(0L, 0L))
  ineffective <- PsychoMatic:::.ipi_verify_releases(
    full, fixed_blank, parsed, c("g1", "g2")
  )
  expect_false(ineffective$effective_release)

  missing_labels <- make_record(c(NA_character_, "nu.2.g2"), c(0L, 1L))
  expect_no_error(verified_missing <- PsychoMatic:::.ipi_verify_releases(
    full, missing_labels, parsed, c("g1", "g2")
  ))
  expect_false(verified_missing$effective_release)

  missing_group <- make_record("nu.2", 0L)
  missing_group$parameter_table$group <- 1L
  expect_error(
    PsychoMatic:::.ipi_validate_complete_families(
      parsed, missing_group, c("g1", "g2"), "intercepts"
    ),
    "every analysis group"
  )
})

test_that("same-stage and inferential metadata mismatches fail independently", {
  measures <- function() list(
    cfi = c(robust = NA, scaled = NA, base = .95),
    rmsea = c(robust = NA, scaled = NA, base = .05),
    srmr = c(robust = NA, scaled = NA, base = .05),
    df = c(robust = NA, scaled = NA, base = 20)
  )
  metadata <- list(
    case_indices = 1:40, variables = paste0("i", 1:3),
    estimator_effective = "ML", estimator_orig = "ML",
    estimand = "means_covariances", ordered = character(), missing = "listwise",
    identification = "std.lv", parameterization = "delta", std_lv = TRUE,
    group_coding = c("g1", "g2"), measurement = "continuous"
  )
  base <- list(
    admissible = TRUE, measures = measures(), fit_metadata = metadata,
    effective_constraints = c("loadings", "intercepts"),
    verified_releases = "i2 ~ 1"
  )
  fields <- names(metadata)
  mutations <- list(
    case_indices = 2:41, variables = rev(metadata$variables),
    estimator_effective = "DWLS", estimator_orig = "MLR",
    estimand = "covariances", ordered = "i1", missing = "fiml",
    identification = "marker", parameterization = "theta", std_lv = FALSE,
    group_coding = rev(metadata$group_coding), measurement = "ordinal"
  )
  for (field in fields) {
    other <- base
    other$fit_metadata[[field]] <- mutations[[field]]
    got <- PsychoMatic:::.ipi_compare_same_stage(base, other, "i2 ~ 1")
    expect_false(got$comparable, info = field)
    expect_identical(got$reason, paste0("mismatch_", field), info = field)
  }
})

test_that("same-stage diagnostic comparator uses the requested direction", {
  measures <- function(cfi, rmsea, srmr, df) list(
    cfi = c(robust = NA, scaled = NA, base = cfi),
    rmsea = c(robust = NA, scaled = NA, base = rmsea),
    srmr = c(robust = NA, scaled = NA, base = srmr),
    df = c(robust = NA, scaled = NA, base = df)
  )
  metadata <- list(
    case_indices = 1:50, variables = paste0("i", 1:6),
    estimator_effective = "ML", estimator_orig = "ML",
    estimand = "means_covariances", ordered = character(), missing = "listwise",
    identification = "std.lv", parameterization = "delta", std_lv = TRUE,
    group_coding = c("A", "B"), measurement = "continuous"
  )
  partial <- list(
    admissible = TRUE, measures = measures(.96, .04, .05, 22),
    effective_constraints = c("loadings", "intercepts"),
    fit_metadata = metadata, verified_releases = "i2 ~ 1",
    group_partial = "i2 ~ 1"
  )
  full <- list(
    admissible = TRUE, measures = measures(.94, .06, .07, 23),
    effective_constraints = c("loadings", "intercepts"),
    fit_metadata = metadata, verified_releases = character(),
    group_partial = character()
  )
  got <- PsychoMatic:::.ipi_compare_same_stage(partial, full, "i2 ~ 1")
  expect_true(got$comparable)
  expect_equal(got$delta_df, 1)
  expect_equal(got$improvement_cfi, .02)
  expect_equal(got$improvement_rmsea, .02)
  expect_equal(got$improvement_srmr, .02)

  extra_release <- partial
  extra_release$group_partial <- c("i2 ~ 1", "i3 ~ 1")
  full$group_partial <- character()
  extra_release$verified_releases <- c("i2 ~ 1", "i3 ~ 1")
  extra <- PsychoMatic:::.ipi_compare_same_stage(
    extra_release, full, "i2 ~ 1"
  )
  expect_false(extra$comparable)
  expect_identical(extra$reason, "mismatch_same_stage_releases")

  unverified_extra <- partial
  unverified_extra$group_partial <- c("i2 ~ 1", "i3 ~ 1")
  extra_argument <- PsychoMatic:::.ipi_compare_same_stage(
    unverified_extra, full, "i2 ~ 1"
  )
  expect_false(extra_argument$comparable)
  expect_identical(extra_argument$reason, "mismatch_same_stage_group_partial")

  no_common_partial <- partial
  no_common_full <- full
  no_common_partial$measures$cfi <- c(robust = .96, scaled = NA, base = NA)
  no_common_full$measures$cfi <- c(robust = NA, scaled = NA, base = .94)
  no_common <- PsychoMatic:::.ipi_compare_same_stage(
    no_common_partial, no_common_full, "i2 ~ 1"
  )
  expect_false(no_common$comparable)
  expect_identical(no_common$reason, "no_common_variant_cfi")

  nonpositive_df <- full
  nonpositive_df$measures <- measures(.94, .06, .07, 22)
  invalid_df <- PsychoMatic:::.ipi_compare_same_stage(
    partial, nonpositive_df, "i2 ~ 1"
  )
  expect_false(invalid_df$comparable)
  expect_identical(invalid_df$reason, "full_minus_partial_df_not_positive")
  expect_equal(invalid_df$delta_df, 0)
})

test_that("configural fit indices cannot override model inadmissibility", {
  measures <- list(
    cfi = c(robust = NA, scaled = NA, base = .95),
    rmsea = c(robust = NA, scaled = NA, base = .05),
    srmr = c(robust = NA, scaled = NA, base = .05)
  )
  inadmissible <- list(
    admissible = FALSE,
    admissibility_reason = "negative_observed_variance",
    measures = measures
  )
  gate <- PsychoMatic:::.ipi_configural_gate_record(
    inadmissible,
    list(configural_cfi = .90, configural_rmsea = .08, configural_srmr = .08)
  )
  expect_false(gate$passed)
  expect_identical(gate$decision, "not_evaluable")
  expect_identical(gate$status, "not_evaluable")
  expect_identical(gate$reason, "negative_observed_variance")

  admissible <- inadmissible
  admissible$admissible <- TRUE
  admissible$admissibility_reason <- "admissible"
  supported <- PsychoMatic:::.ipi_configural_gate_record(
    admissible,
    list(configural_cfi = .90, configural_rmsea = .08, configural_srmr = .08)
  )
  expect_true(supported$passed)
  expect_identical(supported$decision, "supported")
})

test_that("nonconvergence and nonfinite configural output remain not evaluable", {
  measures <- list(
    cfi = c(robust = NA, scaled = NA, base = .95),
    rmsea = c(robust = NA, scaled = NA, base = .05),
    srmr = c(robust = NA, scaled = NA, base = .05)
  )
  reasons <- c("converged", "finite_parameters;nonfinite_cfi")
  for (reason in reasons) {
    gate <- PsychoMatic:::.ipi_configural_gate_record(
      list(
        admissible = FALSE,
        admissibility_reason = reason,
        measures = measures
      ),
      list(configural_cfi = .90, configural_rmsea = .08,
           configural_srmr = .08)
    )
    expect_false(gate$passed, info = reason)
    expect_identical(gate$decision, "not_evaluable", info = reason)
    expect_identical(gate$reason, reason, info = reason)
  }

  nonfinite <- list(
    admissible = TRUE,
    admissibility_reason = "admissible",
    measures = list(
      cfi = c(robust = NA, scaled = NA, base = NA),
      rmsea = measures$rmsea,
      srmr = measures$srmr
    )
  )
  nonfinite_gate <- PsychoMatic:::.ipi_configural_gate_record(
    nonfinite,
    list(configural_cfi = .90, configural_rmsea = .08,
         configural_srmr = .08)
  )
  expect_false(nonfinite_gate$passed)
  expect_identical(nonfinite_gate$decision, "not_evaluable")
  expect_identical(
    nonfinite_gate$reason,
    "configural_indices_missing_or_nonfinite"
  )

  nonfinite_change <- PsychoMatic:::.ipi_decide_change(
    list(
      comparable = TRUE,
      delta_cfi = NA_real_, delta_rmsea = 0, delta_srmr = 0
    ),
    list(cfi = .01, rmsea = .015, srmr_loading = .030,
         srmr_other = .015),
    "scalar"
  )
  expect_identical(nonfinite_change$decision, "not_evaluable")
  expect_identical(
    nonfinite_change$reason,
    "finite_fit_change_indices_are_required"
  )

  expect_true(PsychoMatic:::.ipi_configural_fit_failed(list(
    error = NULL, fit = structure(list(), class = "mock_fit"),
    admissible = FALSE
  )))
  expect_false(PsychoMatic:::.ipi_configural_fit_failed(list(
    error = NULL, fit = structure(list(), class = "mock_fit"),
    admissible = TRUE
  )))
})

test_that("identified ordinal workflow requires semTools before fitting", {
  fit_calls <- 0L
  testthat::local_mocked_bindings(
    .ipi_namespace_available = function(package) {
      !identical(package, "semTools")
    },
    .ipi_fit_record = function(...) {
      fit_calls <<- fit_calls + 1L
      stop("fit attempted", call. = FALSE)
    },
    .package = "PsychoMatic"
  )

  expect_error(
    invariance_partial(
      data.frame(group = c("A", "B"), i1 = c(0, 1)),
      "group", "f =~ i1", free = "f =~ i1", rationale = "Theory.", workflow = "ordinal"
    ),
    "semTools.*required|required.*semTools"
  )
  expect_identical(fit_calls, 0L)
})

test_that("structurally invalid releases fail before any model is fitted", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(n_per_group = 60L, seed = 7315L)
  fit_calls <- 0L
  testthat::local_mocked_bindings(
    .ipi_fit_record = function(...) {
      fit_calls <<- fit_calls + 1L
      stop("A fit must not be attempted for a structurally invalid release.",
           call. = FALSE)
    },
    .package = "PsychoMatic"
  )

  expect_error(
    invariance_partial(
      d, "group", partial_model_6,
      free = "f =~ absent",
      rationale = "A preregistered loading hypothesis.",
      stage = "metric", measurement = "continuous", estimator = "ML",
      ordered = FALSE, language = "eng"
    ),
    "absent.*measurement model|measurement model.*absent"
  )
  expect_identical(fit_calls, 0L)
})

test_that("configural engine failure returns one structured fail-closed result", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(n_per_group = 60L, seed = 7316L)
  fit_calls <- 0L
  engine_error <- simpleError("simulated configural engine failure")
  testthat::local_mocked_bindings(
    .ipi_fit_record = function(prepared, stage_key, stage_detail,
                               constraints, group_partial = character()) {
      fit_calls <<- fit_calls + 1L
      built <- list(
        syntax = prepared$blueprint$model,
        fit = NULL,
        engine = "mock::configural_failure",
        engine_args = list(group.partial = group_partial),
        engine_requests = list(
          estimator = prepared$estimator_info$estimator,
          missing = prepared$missing_request,
          parameterization = prepared$blueprint$parameterization_engine_request,
          ID.fac = prepared$blueprint$ID.fac,
          ID.cat = prepared$blueprint$ID.cat
        ),
        requested_constraints = constraints,
        effective_constraints = constraints,
        group_partial = group_partial,
        parameter_table = data.frame(),
        warnings = list(),
        messages = list(),
        error = engine_error
      )
      PsychoMatic:::.ifa_extract_fit_record(
        built,
        stage = stage_key,
        stage_detail = stage_detail,
        blueprint = prepared$blueprint,
        indicators = prepared$indicators
      )
    },
    .package = "PsychoMatic"
  )

  expect_no_error(out <- invariance_partial(
    d, "group", partial_model_6,
    free = "i3 ~ 1", rationale = "A preregistered intercept hypothesis.",
    stage = "scalar", measurement = "continuous", estimator = "ML",
    ordered = FALSE, language = "eng"
  ))
  expect_s3_class(out, "psychomatic_invariance_partial")
  expect_s3_class(out, "psychomatic_result")
  expect_identical(fit_calls, 1L)
  expect_identical(out$conclusion_status, "not_evaluable")
  expect_true(any(grepl(
    "simulated configural engine failure", out$errors$message,
    fixed = TRUE
  )))
  expect_identical(out$errors$condition[[1L]], engine_error)
  expect_identical(
    out$diagnostics$configural_failure$message,
    "simulated configural engine failure"
  )
  expect_identical(
    names(out$models),
    c(
      "configural", "prerequisite_full", "prerequisite_matched",
      "full", "partial", "prerequisite"
    )
  )
  expect_equal(nrow(out$anchors), 0L)
})

test_that("partial-invariance conclusion states are resolved deterministically", {
  resolve <- PsychoMatic:::.ipi_resolve_conclusion_status
  configural_ok <- list(passed = TRUE, decision = "supported")
  full_not_supported <- list(decision = "not_supported")
  partial_comparable <- list(comparable = TRUE)

  call_resolver <- function(
    releases = TRUE, nesting = TRUE, configural = configural_ok,
    full = full_not_supported, comparison = partial_comparable,
    partial = list(decision = "supported")
  ) {
    resolve(releases, nesting, configural, full, comparison, partial)
  }

  expect_identical(call_resolver(releases = FALSE), "not_evaluable")
  expect_identical(call_resolver(nesting = FALSE), "not_evaluable")
  expect_identical(
    call_resolver(configural = list(passed = FALSE, decision = "not_supported")),
    "configural_not_supported"
  )
  expect_identical(
    call_resolver(configural = list(passed = FALSE, decision = "not_evaluable")),
    "not_evaluable"
  )
  expect_identical(
    call_resolver(full = list(decision = "supported")),
    "full_invariance_supported_partial_not_needed"
  )
  expect_identical(
    call_resolver(comparison = list(comparable = FALSE)),
    "not_evaluable"
  )
  expect_identical(call_resolver(), "partial_invariance_supported")
  expect_identical(
    call_resolver(partial = list(decision = "not_supported")),
    "partial_invariance_not_supported"
  )
  expect_identical(
    call_resolver(partial = list(decision = "not_evaluable")),
    "not_evaluable"
  )
})

test_that("continuous partial invariance forwards releases and freezes rows", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(seed = 7302L)
  original <- d
  expect_silent(
    visible <- withVisible(invariance_partial(
      d, "group", partial_model_6,
      free = "i3 ~ 1",
      rationale = "The item wording differs between groups.",
      stage = "scalar", measurement = "continuous", estimator = "ML",
      ordered = FALSE, report = FALSE, language = "eng", verbose = FALSE
    ))
  )
  out <- visible$value
  expect_false(visible$visible)
  expect_s3_class(out, "psychomatic_invariance_partial")
  expect_s3_class(out, "psychomatic_result")
  expect_identical(class(out)[1:2], c(
    "psychomatic_invariance_partial", "psychomatic_result"
  ))
  expect_identical(out$stage_requested, "scalar")
  expect_identical(out$stage_effective, "scalar")
  expect_identical(out$prerequisite_stage, "metric")
  expect_identical(
    out$models$partial$fit_record$engine_args$group.partial,
    "i3 ~ 1"
  )
  expect_identical(
    out$models$prerequisite_matched$fit_record$engine_args$group.partial,
    character()
  )
  expect_true(out$freed_constraints$effective_release)
  expect_identical(
    out$models$partial$case_idx,
    out$models$prerequisite_matched$case_idx
  )
  expect_identical(d, original)
  expect_null(out$report_path)
  expect_identical(out$effective_settings$ID.fac_engine_request, "marker")
  expect_identical(out$effective_settings$ID.fac_effective, "marker")
  expect_identical(out$effective_settings$ID.cat_engine_request, NA_character_)
  expect_identical(out$effective_settings$ID.cat_effective, NA_character_)
  for (role in names(out$models)) {
    record <- out$models[[role]]$fit_record
    metadata <- record$fit_metadata
    expect_identical(record$case_indices, metadata$case_indices, info = role)
    expect_identical(record$variables, metadata$variables, info = role)
    expect_identical(record$estimator, metadata$estimator_effective, info = role)
    expect_identical(record$missing, metadata$missing, info = role)
    expect_identical(record$parameterization, metadata$parameterization, info = role)
    expect_identical(record$identification, metadata$identification, info = role)
  }
  expect_true(out$conclusion_status %in% c(
    "full_invariance_supported_partial_not_needed",
    "partial_invariance_supported", "partial_invariance_not_supported",
    "configural_not_supported", "not_evaluable"
  ))
})

test_that("matrix and external-group inputs retain original rows and match lavaan", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(n_per_group = 90L, seed = 7314L)
  external_group <- d$group
  x <- as.matrix(d[paste0("i", 1:6)])
  original <- x
  external_group[2L] <- NA
  x[3L, 1L] <- NA_real_
  original[3L, 1L] <- NA_real_

  out <- invariance_partial(
    x, external_group, partial_model_6,
    free = "i3 ~ 1",
    rationale = "A preregistered intercept hypothesis.",
    stage = "scalar", estimator = "ML", ordered = FALSE,
    language = "eng"
  )
  expected_rows <- setdiff(seq_len(nrow(x)), 2L)
  expect_identical(out$analysis_blueprint$group_source, "direct_vector")
  expect_identical(out$models$partial$case_idx, expected_rows)
  expect_identical(out$data_summary$dropped_missing_group, 1L)
  expect_identical(out$data_summary$dropped_missing_items, 0L)
  expect_identical(x, original)

  reference <- lavaan::cfa(
    partial_model_6,
    data = out$analysis_data,
    group = out$analysis_blueprint$group_name,
    estimator = "ML",
    std.lv = FALSE,
    meanstructure = TRUE,
    missing = "fiml",
    group.equal = c("loadings", "intercepts"),
    group.partial = "i3 ~ 1"
  )
  expected_fit <- lavaan::fitMeasures(
    reference, c("df", "cfi", "rmsea", "srmr")
  )
  actual_fit <- lavaan::fitMeasures(
    out$models$partial$fit, c("df", "cfi", "rmsea", "srmr")
  )
  expect_equal(actual_fit, expected_fit, tolerance = 1e-10)
})

test_that("lower-stage releases propagate to the matched prerequisite", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(seed = 7303L)
  out <- invariance_partial(
    d, group, partial_model_6,
    free = "f =~ i3", rationale = "Prior loading evidence.",
    stage = "scalar", estimator = "ML", ordered = FALSE,
    language = "eng"
  )
  expect_identical(
    out$models$prerequisite_matched$fit_record$engine_args$group.partial,
    "f =~ i3"
  )
  expect_identical(
    out$models$partial$fit_record$engine_args$group.partial,
    "f =~ i3"
  )
  expect_true(out$comparison$group_partial_nesting$comparable)
})

test_that("explicit ordinal loadings retain the identified metric key", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(n_per_group = 120L, ordinal = TRUE, seed = 7304L)
  out <- invariance_partial(
    d, "group", partial_model_6,
    free = "f =~ i3", rationale = "Prior ordinal loading evidence.",
    stage = "loadings", measurement = "ordinal", estimator = "WLSMV",
    ordered = TRUE, workflow = "ordinal", language = "eng"
  )
  expect_identical(out$stage_effective, "loadings")
  expect_identical(out$stage, "metric")
  expect_identical(
    out$models$partial$fit_record$engine_args$group.partial,
    "f =~ i3"
  )
  expect_identical(out$effective_settings$parameterization, "theta")

  reference <- semTools::measEq.syntax(
    configural.model = partial_model_6,
    data = out$analysis_data,
    group = out$analysis_blueprint$group_name,
    ordered = out$analysis_blueprint$ordered,
    estimator = "WLSMV",
    parameterization = "theta",
    ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016",
    group.equal = c("thresholds", "loadings"),
    group.partial = "f =~ i3",
    return.fit = TRUE
  )
  expected_fit <- lavaan::fitMeasures(
    reference, c("df", "cfi.scaled", "rmsea.scaled", "srmr")
  )
  actual_fit <- lavaan::fitMeasures(
    out$models$partial$fit,
    c("df", "cfi.scaled", "rmsea.scaled", "srmr")
  )
  expect_equal(actual_fit, expected_fit, tolerance = 1e-10)
})

test_that("ordinal strict release propagates thresholds and frees residuals", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- partial_fixture(n_per_group = 220L, ordinal = TRUE, seed = 7317L)
  releases <- c("i3 | t2", "i4 ~~ i4")
  out <- invariance_partial(
    d, "group", partial_model_6,
    free = releases,
    rationale = "The two preregistered item parameters may differ by group.",
    stage = "strict", measurement = "ordinal", estimator = "WLSMV",
    ordered = TRUE, workflow = "ordinal", language = "eng"
  )

  expect_identical(
    out$models$prerequisite_matched$fit_record$engine_args$group.partial,
    "i3 | t2"
  )
  expect_identical(
    out$models$partial$fit_record$engine_args$group.partial,
    releases
  )
  expect_true(all(out$freed_constraints$effective_release))
  expect_true(out$comparison$group_partial_nesting$comparable)
  expect_true(out$comparison$prerequisite_vs_partial$comparable)
})

test_that("mixed indicators preserve the exact ordered subset", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- simulate_invariance_fixture(
    n_per_group = 220L,
    ordinal = FALSE,
    partial = TRUE,
    seed = 7313L
  )
  items <- paste0("i", 1:6)
  for (item in items[1:3]) {
    d[[item]] <- ordered(
      cut(d[[item]], c(-Inf, -0.7, 0, 0.7, Inf), labels = FALSE),
      levels = 1:4
    )
  }
  d <- d[c("group", items)]
  requested_order <- rev(items[1:3])
  out <- invariance_partial(
    d, "group", partial_model_6,
    free = "f =~ i3",
    rationale = "A preregistered mixed-indicator loading hypothesis.",
    stage = "loadings", measurement = "mixed", estimator = "WLSMV",
    ordered = requested_order, workflow = "ordinal", language = "eng"
  )

  expect_identical(out$measurement, "mixed")
  expect_identical(out$stage_effective, "loadings")
  expect_identical(
    out$models$partial$fit_record$engine_args$ordered,
    items[1:3]
  )
  expect_identical(out$requested_settings$ordered, requested_order)
  expect_identical(out$effective_settings$ordered_engine_request, requested_order)
  expect_identical(
    out$effective_settings$ordered,
    out$models$configural$fit_record$fit_metadata$ordered
  )
  expect_identical(
    out$effective_settings$ordered_effective,
    out$models$configural$fit_record$fit_metadata$ordered
  )
  expect_true(out$freed_constraints$effective_release)
  expect_true(out$comparison$prerequisite_vs_partial$comparable)
  expect_equal(nrow(out$errors), 0L)
})

test_that("input assertions, anchors, higher order, and binary stages fail early", {
  skip_if_not_installed("semTools")
  d <- partial_fixture(seed = 7305L)
  expect_error(
    invariance_partial(
      d, "group", partial_model_6, "i3 ~ 1", "Reason",
      measurement = "ordinal", estimator = "ML", ordered = FALSE
    ),
    "measurement|ordered"
  )
  expect_error(
    invariance_partial(
      d, "group", partial_model_6,
      c("i1 ~ 1", "i2 ~ 1", "i3 ~ 1", "i4 ~ 1", "i5 ~ 1"),
      "Too few anchors", estimator = "ML", ordered = FALSE,
      min_anchor_items = 2L
    ),
    "anchor|retained"
  )
  for (bad in list(1L, 1.5, NA_real_, Inf, c(2, 3))) {
    expect_error(
      invariance_partial(
        d, "group", partial_model_6, "i3 ~ 1", "Reason",
        estimator = "ML", ordered = FALSE, min_anchor_items = bad
      ),
      "min_anchor_items"
    )
  }
  higher <- paste(partial_model_6, "g =~ f", sep = "\n")
  expect_error(
    invariance_partial(
      d, "group", higher, "i3 ~ 1", "Reason",
      estimator = "ML", ordered = FALSE
    ),
    "[Hh]igher-order"
  )

  binary <- d
  binary[paste0("i", 1:6)] <- lapply(
    binary[paste0("i", 1:6)],
    function(x) ordered(as.integer(x > stats::median(x)), levels = 0:1)
  )
  expect_error(
    invariance_partial(
      binary, "group", partial_model_6, "f =~ i3", "Reason",
      stage = "loadings", ordered = TRUE, workflow = "ordinal"
    ),
    "binary|joint"
  )
})

test_that("binary partial invariance uses the jointly identified scalar stage", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  d <- simulate_invariance_fixture(
    n_per_group = 600L,
    ordinal = FALSE,
    partial = TRUE,
    seed = 7312L
  )
  items <- paste0("i", 1:6)
  d[items] <- lapply(d[items], function(pk_score) {
    ordered(as.integer(pk_score > 0), levels = 0:1)
  })
  d <- d[c("group", items)]
  out <- invariance_partial(
    d, "group", partial_model_6,
    free = "i3 | t1",
    rationale = "A preregistered binary-threshold hypothesis.",
    stage = "scalar", measurement = "ordinal", estimator = "WLSMV",
    ordered = TRUE, workflow = "ordinal", language = "eng"
  )

  expect_identical(out$measurement, "ordinal")
  expect_identical(out$effective_settings$category_profile, "binary")
  expect_identical(out$prerequisite_stage, "configural")
  expect_true(out$freed_constraints$effective_release)
  expect_true(out$comparison$prerequisite_vs_partial$comparable)
  expect_gt(out$comparison$prerequisite_vs_partial$delta_df, 0)
  expect_equal(nrow(out$errors), 0L)
})

test_that("language precedence, methods, and reports preserve raw results", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("semTools")
  skip_if_not_installed("openxlsx")
  d <- partial_fixture(seed = 7306L)
  spec <- psychomatic_spec(
    d,
    items = paste0("i", 1:6), group = "group",
    level = "continuous", language = "eng"
  )
  path <- tempfile("partial-invariance-", fileext = ".xlsx")
  out <- invariance_partial(
    spec, 1L, partial_model_6,
    free = "i3 ~ 1", rationale = "Theory.",
    estimator = "ML", ordered = FALSE, report = path, language = NULL
  )
  expect_identical(out$effective_settings$language, "eng")
  expect_identical(out$analysis_blueprint$group_source, "column_index")
  expect_true(file.exists(path))
  expect_identical(out$report_path, normalizePath(path, winslash = "/"))
  required_tables <- c(
    "fit", "freed_constraints", "rationale", "anchors", "anchor_counts",
    "comparisons", "criteria", "warnings", "errors", "conclusion",
    "settings", "references"
  )
  expect_true(all(required_tables %in% names(out$tables)))
  expect_identical(out$tables$references, out$references)
  expect_identical(out$tables$conclusion$status, out$conclusion_status)
  expect_identical(nrow(out$tables$warnings), nrow(out$warnings))
  expect_identical(nrow(out$tables$errors), nrow(out$errors))
  integrated_tables <- PsychoMatic:::.pm_object_tables(out)
  expect_true(all(names(out$tables) %in% names(integrated_tables)))
  for (table_name in names(out$tables)) {
    expect_identical(
      integrated_tables[[table_name]], out$tables[[table_name]],
      info = table_name
    )
  }
  expect_identical(length(openxlsx::getSheetNames(path)), 2L)
  report_cells <- unlist(openxlsx::read.xlsx(path, sheet = 2L, colNames = FALSE))
  expect_true(all(c("Released parameters", "Comparisons", "References") %in% report_cells))
  integrated <- PsychoMatic:::.pm_object_tables(out)
  expect_true(all(required_tables %in% names(integrated)))
  raw <- out$comparison
  expect_silent(capture.output(printed <- print(out)))
  expect_identical(printed, out)
  spanish <- out
  spanish$effective_settings$language <- "esp"
  spanish$conclusion_text <- "Conclusión metodológicamente cauta."
  printed_spanish <- capture.output(print(spanish))
  expect_true(any(grepl(
    "invariancia parcial", printed_spanish,
    ignore.case = TRUE, useBytes = TRUE
  )))
  expect_true(any(grepl("Estado", printed_spanish, useBytes = TRUE)))
  expect_false(any(grepl(
    "Measurement:|Stage:|Status:", printed_spanish,
    useBytes = TRUE
  )))
  summarized <- summary(out)
  expect_true(all(c("fit", "freed_constraints", "anchors", "comparison") %in% names(summarized)))
  expect_identical(out$comparison, raw)
  plot_file <- tempfile("partial-invariance-plot-", fileext = ".pdf")
  grDevices::pdf(plot_file)
  plot_data <- tryCatch(
    expect_silent(plot(out)),
    finally = grDevices::dev.off()
  )
  expect_true(is.data.frame(plot_data))
  expect_error(
    invariance_partial(
      spec, "group", partial_model_6, "i3 ~ 1", "Theory.",
      estimator = "ML", ordered = FALSE, report = path
    ),
    "exists"
  )
})
