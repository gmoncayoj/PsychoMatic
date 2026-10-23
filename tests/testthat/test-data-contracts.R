test_that("language normalization defaults to English and preserves aliases", {
  old <- options(PsychoMatic.language = NULL)
  on.exit(options(old), add = TRUE)
  expect_equal(PsychoMatic:::.pm_normalize_language(NULL), "eng")
  expect_equal(PsychoMatic:::.pm_normalize_language("es"), "esp")
  expect_equal(PsychoMatic:::.pm_normalize_language("english"), "eng")
  expect_error(PsychoMatic:::.pm_normalize_language("fr"), "esp.*eng")
})

test_that("language normalization applies explicit inherited option precedence", {
  old <- options(PsychoMatic.language = "eng")
  on.exit(options(old), add = TRUE)

  expect_identical(PsychoMatic:::.pm_normalize_language(NULL), "eng")
  expect_identical(
    PsychoMatic:::.pm_normalize_language(NULL, inherited = "es"),
    "esp"
  )
  expect_identical(
    PsychoMatic:::.pm_normalize_language("en", inherited = "es"),
    "eng"
  )
})

test_that("result envelopes keep mandatory fields and raw numbers", {
  out <- PsychoMatic:::.pm_result(
    class = "psychomatic_probe",
    call = quote(probe(x)),
    requested_settings = list(method = "auto"),
    effective_settings = list(method = "pearson"),
    data_summary = list(n = 10L, p = 3L),
    results = list(value = 1 / 3),
    warnings = data.frame(code = character(), message = character()),
    errors = data.frame(code = character(), message = character()),
    diagnostics = list(),
    decision_log = data.frame(code = "method", requested = "auto", effective = "pearson", reason = "continuous", severity = "info"),
    conclusion = NA_character_
  )
  expect_s3_class(out, "psychomatic_probe")
  expect_named(out, c("call", "requested_settings", "effective_settings", "data_summary", "results", "warnings", "errors", "decision_log", "diagnostics", "conclusion"))
  expect_identical(out$results$value, 1 / 3)
})

test_that("decisions and non-evaluable states use stable structured fields", {
  decision <- PsychoMatic:::.pm_decision(
    "fallback", "polychoric", "pearson", "sparse cells", "warning"
  )

  expect_s3_class(decision, "data.frame")
  expect_equal(nrow(decision), 1L)
  expect_named(
    decision,
    c("code", "requested", "effective", "reason", "severity")
  )
  expect_identical(decision$severity, "warning")
  expect_identical(
    PsychoMatic:::.pm_not_evaluable("model did not converge"),
    list(status = "not_evaluable", reason = "model did not converge")
  )
})

test_that("decisions reject non-scalar fields instead of recycling rows", {
  valid <- list(
    code = "method",
    requested = "auto",
    effective = "pearson",
    reason = "continuous",
    severity = "info"
  )

  for (field in names(valid)) {
    invalid <- valid
    invalid[[field]] <- c("first", "second")
    expect_error(
      do.call(PsychoMatic:::.pm_decision, invalid),
      paste0("`", field, "`.*single")
    )
  }
  valid["requested"] <- list(NULL)
  null_error <- tryCatch(
    do.call(PsychoMatic:::.pm_decision, valid),
    error = identity
  )
  expect_s3_class(null_error, "error")
  expect_match(conditionMessage(null_error), "`requested`.*single")
})

test_that("scalar-number validation returns valid raw values and rejects invalid ones", {
  expect_identical(
    PsychoMatic:::.pm_validate_scalar_number(1 / 3, "alpha", 0, 1),
    1 / 3
  )
  expect_identical(
    PsychoMatic:::.pm_validate_scalar_number(3L, "iterations", 1, 10, integer = TRUE),
    3L
  )
  expect_error(
    PsychoMatic:::.pm_validate_scalar_number(c(1, 2), "alpha", 0, 1),
    "alpha.*single.*finite"
  )
  expect_error(
    PsychoMatic:::.pm_validate_scalar_number(1.5, "iterations", 1, 10, integer = TRUE),
    "iterations.*integer"
  )
  expect_error(
    PsychoMatic:::.pm_validate_scalar_number(2, "alpha", 0, 1),
    "alpha.*between"
  )
})

test_that("condition capture keeps original classes, messages, and calls", {
  warning_condition <- structure(
    list(message = "fallback used", call = quote(polychoric(x))),
    class = c("psychomatic_warning", "warning", "condition")
  )
  captured <- PsychoMatic:::.pm_capture_conditions({
    warning(warning_condition)
    message("diagnostic note")
    1 / 3
  })

  expect_named(captured, c("value", "warnings", "messages", "error"))
  expect_identical(captured$value, 1 / 3)
  expect_length(captured$warnings, 1L)
  expect_s3_class(captured$warnings[[1L]], "psychomatic_warning")
  expect_identical(conditionMessage(captured$warnings[[1L]]), "fallback used")
  expect_identical(conditionCall(captured$warnings[[1L]]), quote(polychoric(x)))
  expect_length(captured$messages, 1L)
  expect_s3_class(captured$messages[[1L]], "message")
  expect_match(conditionMessage(captured$messages[[1L]]), "diagnostic note")
  expect_null(captured$error)

  error_condition <- structure(
    list(message = "engine failed", call = quote(engine(x))),
    class = c("psychomatic_engine_error", "error", "condition")
  )
  failed <- PsychoMatic:::.pm_capture_conditions(stop(error_condition))
  expect_null(failed$value)
  expect_s3_class(failed$error, "psychomatic_engine_error")
  expect_identical(conditionMessage(failed$error), "engine failed")
  expect_identical(conditionCall(failed$error), quote(engine(x)))
})

test_that("result envelopes structure conditions and reject duplicate fields", {
  warning_condition <- structure(
    list(message = "fallback used", call = quote(polychoric(x))),
    class = c("psychomatic_warning", "warning", "condition")
  )
  error_condition <- structure(
    list(message = "engine failed", call = quote(engine(x))),
    class = c("psychomatic_engine_error", "error", "condition")
  )

  out <- PsychoMatic:::.pm_result(
    class = "psychomatic_probe",
    call = quote(probe(x)),
    requested_settings = list(method = "polychoric"),
    effective_settings = list(method = "pearson"),
    data_summary = list(n = 10L, p = 3L),
    results = list(value = 1 / 3),
    warnings = list(warning_condition),
    errors = list(error_condition),
    diagnostics = list(),
    decision_log = PsychoMatic:::.pm_decision(
      "fallback", "polychoric", "pearson", "sparse cells", "warning"
    ),
    conclusion = PsychoMatic:::.pm_not_evaluable("engine failed")
  )

  expect_s3_class(out, "psychomatic_result")
  expect_s3_class(out$warnings, "data.frame")
  expect_s3_class(out$errors, "data.frame")
  expect_identical(out$warnings$message, "fallback used")
  expect_identical(out$warnings$class[[1L]], class(warning_condition))
  expect_identical(out$warnings$call[[1L]], quote(polychoric(x)))
  expect_s3_class(out$warnings$condition[[1L]], "psychomatic_warning")
  expect_identical(out$errors$message, "engine failed")
  expect_identical(out$errors$class[[1L]], class(error_condition))
  expect_identical(out$errors$call[[1L]], quote(engine(x)))
  expect_s3_class(out$errors$condition[[1L]], "psychomatic_engine_error")

  duplicate_results <- structure(
    list(1, 2),
    names = c("estimate", "estimate")
  )
  expect_error(
    PsychoMatic:::.pm_result(
      class = "psychomatic_probe",
      call = quote(probe(x)),
      requested_settings = list(),
      effective_settings = list(),
      data_summary = list(),
      results = duplicate_results,
      warnings = NULL,
      errors = NULL,
      diagnostics = list(),
      decision_log = data.frame(),
      conclusion = NA_character_
    ),
    "results.*duplicate"
  )
})

test_that("simulation fixtures are deterministic and restore caller RNG state", {
  set.seed(991L)
  original_seed <- .Random.seed

  likert <- simulate_likert_fixture(40L, 111L)
  expect_identical(.Random.seed, original_seed)
  expect_identical(likert, simulate_likert_fixture(40L, 111L))
  expect_equal(dim(likert), c(40L, 12L))
  expect_true(all(unlist(likert, use.names = FALSE) %in% 1:5))

  invariance <- simulate_invariance_fixture(30L, seed = 112L)
  expect_identical(.Random.seed, original_seed)
  expect_equal(dim(invariance), c(60L, 10L))
  expect_equal(levels(invariance$group), c("g1", "g2"))
  expect_true(all(vapply(invariance[paste0("i", 1:9)], is.numeric, logical(1))))

  ordinal <- simulate_invariance_fixture(30L, ordinal = TRUE, seed = 112L)
  expect_identical(.Random.seed, original_seed)
  expect_true(all(vapply(ordinal[paste0("i", 1:9)], is.ordered, logical(1))))
  expect_true(all(vapply(
    ordinal[paste0("i", 1:9)],
    function(x) identical(levels(x), as.character(1:5)),
    logical(1)
  )))

  partial <- simulate_invariance_fixture(30L, partial = TRUE, seed = 112L)
  expect_identical(.Random.seed, original_seed)
  expect_false(identical(invariance$i3, partial$i3))

  bifactor <- simulate_bifactor_fixture(40L, 113L)
  expect_identical(.Random.seed, original_seed)
  expect_identical(bifactor, simulate_bifactor_fixture(40L, 113L))
  expect_equal(dim(bifactor), c(40L, 12L))
  expect_true(all(vapply(bifactor, is.numeric, logical(1))))
})

test_that("simulation fixtures restore an initially absent RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(simulate_likert_fixture(20L, 114L))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("psychomatic_spec excludes declared metadata and classifies items deterministically", {
  d <- data.frame(
    respondent_id = 101:108,
    q1 = c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L),
    q2 = ordered(
      c("low", "mid", "high", "mid", "low", "high", "mid", "low"),
      levels = c("low", "mid", "high")
    ),
    group = rep(c("A", "B"), 4)
  )
  expect_error(psychomatic_spec(d, group = "group"), "possible identifier")
  spec <- psychomatic_spec(
    d,
    items = c("q1", "q2"),
    metadata = "respondent_id",
    group = "group"
  )
  expect_s3_class(spec, "psychomatic_spec")
  expect_identical(spec$items, c("q1", "q2"))
  expect_identical(unname(spec$measurement_level), c("ordinal", "ordinal"))
  expect_identical(spec$group, "group")
})

test_that("automatic item detection checks identifier-like original columns", {
  d <- data.frame(
    respondent_id = paste0("participant-", 1:4),
    q1 = c(1L, 2L, 1L, 2L)
  )

  expect_error(psychomatic_spec(d), "respondent_id")
  declared <- psychomatic_spec(d, metadata = "respondent_id")
  expect_identical(declared$items, "q1")
})

test_that("ambiguous nominal encodings are never alphabetically ordered", {
  d <- data.frame(q1 = c("low", "high", "mid"), q2 = 1:3)
  expect_error(psychomatic_spec(d, items = "q1"), "ambiguous")
  expect_error(
    psychomatic_spec(
      data.frame(q1 = factor(c("low", "high", "mid"))),
      items = "q1"
    ),
    "ambiguous"
  )
  expect_error(
    psychomatic_spec(
      data.frame(q1 = c(TRUE, FALSE, TRUE)),
      items = "q1"
    ),
    "ambiguous"
  )
  spec <- psychomatic_spec(
    d,
    items = "q1",
    ordinal_levels = list(q1 = c("low", "mid", "high"))
  )
  expect_identical(levels(spec$data$q1), c("low", "mid", "high"))
})

test_that("automatic measurement uses the seven-category boundary and explicit overrides", {
  d <- data.frame(
    q7 = rep(1:7, each = 2),
    q8 = rep(1:8, length.out = 14),
    q_fraction = seq(0.1, 1.4, by = 0.1)
  )

  spec <- psychomatic_spec(d, items = c("q8", "q7", "q_fraction"))
  expect_identical(
    unname(spec$measurement_level),
    c("continuous", "ordinal", "continuous")
  )
  overridden <- psychomatic_spec(
    d,
    items = c("q7", "q8"),
    level = c(q7 = "continuous", q8 = "ordinal")
  )
  expect_identical(
    unname(overridden$measurement_level),
    c("continuous", "ordinal")
  )

  mapped <- psychomatic_spec(
    data.frame(q1 = c("second", "first", "third")),
    items = "q1",
    level = "continuous",
    ordinal_levels = list(q1 = c("first", "second", "third"))
  )
  expect_identical(unname(mapped$measurement_level), "ordinal")
  expect_identical(levels(mapped$data$q1), c("first", "second", "third"))
  expect_error(
    psychomatic_spec(
      data.frame(q1 = c("low", "mid", "high")),
      items = "q1",
      ordinal_levels = list(q1 = c("low", "mid"))
    ),
    "cover.*observed"
  )
})

test_that("check_items returns coded issues for invalid psychometric inputs", {
  d <- data.frame(
    a = c(1, 1, 1, 1),
    b = c(1, Inf, NA, 2),
    c = c(NA, NA, NA, NA)
  )
  checked <- check_items(
    psychomatic_spec(d, items = c("a", "b", "c"), level = "continuous")
  )
  expect_s3_class(checked, "psychomatic_check")
  expect_true(all(
    c("zero_variance", "non_finite", "all_missing") %in%
      checked$issues$code
  ))
  expect_true(all(c("severity", "item", "message") %in% names(checked$issues)))
})

test_that("ordinal mappings preserve non-finite numeric issues", {
  spec <- psychomatic_spec(
    data.frame(q1 = c(1, Inf, 1)),
    items = "q1",
    ordinal_levels = list(q1 = c(1, Inf))
  )

  expect_true("non_finite" %in% check_items(spec)$issues$code)
})

test_that("ordinal mappings preserve numeric range violations", {
  spec <- psychomatic_spec(
    data.frame(q1 = c(1, 5, 1)),
    items = "q1",
    ordinal_levels = list(q1 = 1:5),
    ranges = list(q1 = c(1, 4))
  )

  expect_true("range_violation" %in% check_items(spec)$issues$code)
})

test_that("check_items covers categorical, naming, range, and overlap issues", {
  categorical <- ordered(
    c(rep("common", 19), "rare"),
    levels = c("common", "rare", "unused")
  )
  d <- data.frame(
    categorical = categorical,
    other = c(1, rep(NA_real_, 19)),
    check.names = FALSE
  )
  names(d)[2L] <- "bad name"
  spec <- psychomatic_spec(
    d,
    items = 1:2,
    level = c(categorical = "ordinal", "bad name" = "continuous"),
    ranges = list("bad name" = c(2, 5))
  )
  issues <- check_items(spec)$issues
  expect_true(all(c(
    "near_zero_variance", "sparse_category", "empty_category",
    "nonsyntactic_name", "range_violation", "insufficient_overlap"
  ) %in% issues$code))

  duplicated <- data.frame(q1 = 1:4, q2 = c(1, 1, 2, 2), check.names = FALSE)
  names(duplicated) <- c("same", "same")
  duplicate_spec <- psychomatic_spec(
    duplicated,
    items = 1:2,
    level = "continuous"
  )
  expect_true("duplicate_name" %in% check_items(duplicate_spec)$issues$code)
})

test_that("check_items diagnoses correlation-matrix contracts", {
  nonsymmetric <- matrix(c(1, 0.2, 0.4, 1), nrow = 2)
  bad_scale <- matrix(c(2, 0.2, 0.2, 1), nrow = 2)
  indefinite <- matrix(c(1, 2, 2, 1), nrow = 2)

  expect_true(
    "non_symmetric_matrix" %in% check_items(nonsymmetric)$issues$code
  )
  expect_true(
    "not_correlation_matrix" %in% check_items(bad_scale)$issues$code
  )
  expect_true(
    "not_positive_definite" %in% check_items(indefinite)$issues$code
  )
})

test_that("specification and item checks provide stable S3 summaries and plots", {
  spec <- psychomatic_spec(
    data.frame(q1 = c(1, 1, 1, 1), q2 = c(1, 2, 3, 4)),
    items = c("q1", "q2"),
    level = "continuous",
    language = "eng"
  )
  checked <- check_items(spec)

  expect_invisible(print(spec))
  expect_invisible(print(checked))
  spec_summary <- summary(spec)
  check_summary <- summary(checked)
  expect_named(spec_summary, c("data_summary", "items", "decision_log"))
  expect_named(check_summary, c("severity_counts", "issues"))
  expect_identical(checked$effective_settings$language, "eng")
  expect_equal(sum(check_summary$severity_counts$n), nrow(checked$issues))

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  plotted <- plot(checked)
  expect_identical(plotted, check_summary$severity_counts)
})

test_that("index-based reverse keys retain identity with duplicated names", {
  d <- data.frame(
    first = c(1, 2, 3),
    second = c(3, 2, 1),
    check.names = FALSE
  )
  names(d) <- c("same", "same")

  spec <- psychomatic_spec(
    d,
    items = 1:2,
    level = "continuous",
    reverse = 1
  )

  expect_identical(summary(spec)$items$reverse, c(TRUE, FALSE))
})

test_that("existing numeric item preparation consumes a specification selection", {
  d <- data.frame(
    metadata_id = 101:106,
    q1 = c(1, 2, 3, 4, 5, 1),
    q2 = ordered(c(1, 2, 3, 2, 1, 3))
  )
  spec <- psychomatic_spec(
    d,
    items = c("q2", "q1"),
    metadata = "metadata_id"
  )
  screened <- screen_items(spec)
  expect_identical(screened$item, c("q2", "q1"))
})
