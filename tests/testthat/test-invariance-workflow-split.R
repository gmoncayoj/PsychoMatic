test_that("factorial invariance fits four native stages for every estimator", {
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 7120L)
  model <- "F =~ i1 + i2 + i3 + i4 + i5 + i6"
  constraints <- list(configural = character(), metric = "loadings",
    scalar = c("loadings", "intercepts"),
    strict = c("loadings", "intercepts", "residuals"))
  for (est in c("ML", "MLR", "ULS", "WLSMV", "ULSMV")) {
    out <- factorial_invariance_auto(d, "group", model, estimator = est,
      ordered = TRUE, verbose = FALSE)
    expect_identical(out$fit_table$level, names(constraints), info = est)
    categorical <- est %in% c("WLSMV", "ULSMV")
    expect_identical(out$measurement, if (categorical) "ordinal" else "continuous", info = est)
    expect_identical(out$estimator, est)
    expect_identical(out$effective_settings$ordered, if (categorical) paste0("i", 1:6) else character())
    expect_identical(out$effective_settings$source_measurement, "ordinal")
    if (!categorical) {
      expect_true(all(vapply(out$analysis_data[paste0("i", 1:6)], is.numeric, logical(1))))
      expect_equal(out$analysis_data$i1, as.numeric(as.character(d$i1)))
    }
    for (stage in names(constraints)) {
      expect_identical(out$fit_records[[stage]]$effective_constraints, constraints[[stage]])
      expect_identical(any(lavaan::parameterTable(out$fitted_models[[stage]])$op == "|"), categorical)
      native <- .pm_capture_conditions(lavaan::cfa(model, data = out$analysis_data, group = "group",
        estimator = est, ordered = if (categorical) TRUE else FALSE,
        missing = out$analysis_blueprint$missing_engine_request,
        group.equal = constraints[[stage]]))$value
      expect_equal(lavaan::fitMeasures(out$fitted_models[[stage]]),
        lavaan::fitMeasures(native), tolerance = 1e-8)
    }
    if (est == "WLSMV") {
      expect_identical(out$effective_settings$missing, "listwise")
      native <- .pm_capture_conditions(lavaan::cfa(model, data = out$analysis_data, group = "group",
        estimator = est, missing = "listwise", ordered = TRUE,
        group.equal = c("loadings", "intercepts", "residuals")))$value
      expect_equal(lavaan::fitMeasures(out$fitted_models$strict),
        lavaan::fitMeasures(native), tolerance = 1e-8)
    }
    sections <- .pm_analysis_sections(out, language = "eng")
    expect_false("Equal thresholds" %in% unlist(sections$sections))
    expect_match(paste(unlist(sections$sections), collapse = " "), "numeric scores")
  }
})

test_that("categorical encodings have explicit score conversion and missing policy", {
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 7121L)
  model <- "F =~ i1 + i2 + i3 + i4 + i5 + i6"
  levels(d$i1) <- c("very low", "low", "medium", "high", "very high")
  out <- factorial_invariance_auto(d, "group", model, estimator = "MLR", verbose = FALSE)
  expect_identical(out$effective_settings$score_coding$i1, "ordered_level_ranks")
  expect_identical(out$analysis_data$i1, as.numeric(d$i1))
  expect_identical(out$effective_settings$score_coding$i2, "numeric_category_labels")
  expect_error(factorial_invariance_auto(d, "group", model,
    estimator = "WLSMV", missing = "fiml", verbose = FALSE), "FIML|fiml|listwise")
})

test_that("traditional reader reports contain only four stages in both languages", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 7122L)
  out <- factorial_invariance_auto(d, "group", "F =~ i1+i2+i3+i4+i5+i6",
    estimator = "WLSMV", verbose = FALSE)
  directory <- tempfile("factorial-four-stage-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  for (language in c("eng", "esp")) {
    sections <- .pm_analysis_sections(out, language)
    table <- sections$sections[[if (language == "eng") "Fit by stage" else "Ajuste por etapa"]]
    expect_equal(nrow(table), 4L)
    paths <- psychomatic_report(out, format = c("word", "excel"),
      language = language, file_name = file.path(directory, language))
    expect_true(all(file.exists(paths)))
    text <- paste(officer::docx_summary(officer::read_docx(paths[["word"]]))$text, collapse = "\n")
    expect_false(grepl("Equal thresholds|Igualdad de umbrales", text))
    expect_match(text, if (language == "eng") "native lavaan" else "nativas de lavaan", fixed = TRUE)
    expect_match(text, if (language == "eng") "Fit-change heuristics|Configural evidence" else
      "heur.stica|evidencia configural", perl = TRUE)
  }
})
