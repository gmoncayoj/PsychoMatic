test_that("ordinal invariance fits exactly three native identified stages", {
  skip_if_not_installed("semTools")
  d <- simulate_invariance_fixture(200L, ordinal = TRUE, seed = 7110L)
  model <- "F =~ i1 + i2 + i3 + i4 + i5 + i6"
  out <- ordinal_invariance_auto(d, group, model, verbose = FALSE)
  expect_s3_class(out, "ordinal_invariance_auto")
  expect_identical(out$fit_table$level,
    c("configural", "thresholds", "thresholds_loadings"))
  expect_identical(out$estimator, "WLSMV")
  expect_identical(out$measurement, "ordinal")
  expect_identical(out$fit_records$thresholds_loadings$requested_constraints,
    c("thresholds", "loadings"))
  expect_false(any(out$constraint_families$family %in% c("intercepts", "residuals")))
  native <- semTools::measEq.syntax(configural.model = model, data = out$analysis_data,
    group = out$analysis_blueprint$group_name, ordered = out$indicators,
    estimator = "WLSMV", parameterization = "theta", ID.fac = "std.lv",
    ID.cat = "Wu.Estabrook.2016", group.equal = c("thresholds", "loadings"),
    return.fit = TRUE)
  expect_equal(lavaan::fitMeasures(out$fitted_models$thresholds_loadings),
    lavaan::fitMeasures(native), tolerance = 1e-7)
  expect_identical(summary(out)$fit, out$fit_table)
  text <- paste(capture.output(print(out)), collapse = "\n")
  expect_match(text, "Equal thresholds and loadings", fixed = TRUE)
  expect_match(text, "Svetina", fixed = TRUE)
  expect_match(out$conclusion, "thresholds and loadings|configural|threshold", perl = TRUE)
})

test_that("ordinal API rejects nonordinal estimators and continuous indicators", {
  d <- simulate_invariance_fixture(100L, ordinal = TRUE)
  model <- "F =~ i1 + i2 + i3 + i4 + i5 + i6"
  for (est in c("ML", "MLR", "ULS", "DWLS")) {
    expect_error(ordinal_invariance_auto(d, "group", model, estimator = est,
      verbose = FALSE), "WLSMV.*ULSMV")
  }
  expect_error(ordinal_invariance_auto(d, "group", model, ordered = FALSE,
    verbose = FALSE), "ordered")
  cont <- simulate_invariance_fixture(100L)
  expect_error(ordinal_invariance_auto(cont, "group", model, verbose = FALSE),
    "ordinal")
  d$i1 <- as.numeric(d$i1) + seq_len(nrow(d)) / nrow(d)
  expect_error(ordinal_invariance_auto(d, "group", model, verbose = FALSE),
    "ordinal")
})

test_that("ordinal changes use the loading SRMR cutoff and signed differences", {
  cutoffs <- list(cfi = .010, rmsea = .015, srmr_loading = .030, srmr_other = .010)
  expect_identical(.ifa_change_decision(-.011, .001, .020, cutoffs,
    "thresholds_loadings")$decision, "supported")
  expect_identical(.ifa_change_decision(.011, -.020, -.040, cutoffs,
    "thresholds_loadings")$decision, "supported")
})

test_that("unidentified separate ordinal stages are never declared invariant", {
  skip_if_not_installed("semTools")
  d <- simulate_invariance_fixture(300L, ordinal = TRUE, seed = 7112L)
  model <- "F =~ i1 + i2 + i3 + i4 + i5 + i6"
  for (categories in c(2L, 3L)) {
    small <- d
    for (item in paste0("i", 1:9)) {
      small[[item]] <- ordered(if (categories == 2L) as.integer(d[[item]] > 3L) else
        pmin(3L, ceiling(as.integer(d[[item]]) / 2)))
    }
    out <- ordinal_invariance_auto(small, "group", model, verbose = FALSE)
    expect_identical(out$fit_table$level,
      c("configural", "thresholds", "thresholds_loadings"))
    expect_identical(out$fit_table["thresholds", "decision"], "not_evaluable")
    expect_identical(out$fit_table["thresholds", "status"], "skipped_not_identified")
    expect_equal(out$fit_table["thresholds", "delta_df"], 0)
  }
})

test_that("ordinal reports preserve three stages, languages and APA references", {
  skip_if_not_installed("semTools")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  d <- simulate_invariance_fixture(200L, ordinal = TRUE, seed = 7113L)
  group_vector <- d$group
  fit <- ordinal_invariance_auto(d, group_vector,
    "F =~ i1 + i2 + i3 + i4 + i5 + i6", estimator = "ULSMV", verbose = FALSE)
  before <- fit$fit_table
  expect_identical(fit$estimator, "ULSMV")
  directory <- tempfile("ordinal-report-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  for (language in c("eng", "esp")) {
    paths <- psychomatic_report(fit, format = c("word", "excel"),
      language = language, file_name = file.path(directory, language))
    expect_true(all(file.exists(paths)))
    word <- officer::docx_summary(officer::read_docx(paths[["word"]]))
    word_text <- paste(word$text, collapse = "\n")
    expect_match(word_text, if (language == "eng") "Equal thresholds and loadings" else
      "Igualdad de umbrales y cargas", fixed = TRUE)
    expect_match(word_text, "10.1080/10705511.2019.1602776", fixed = TRUE)
    sheets <- openxlsx::getSheetNames(paths[["excel"]])
    excel_text <- paste(unlist(lapply(sheets, function(sheet)
      openxlsx::read.xlsx(paths[["excel"]], sheet = sheet))), collapse = "\n")
    expect_match(excel_text, if (language == "eng") "Equal thresholds and loadings" else
      "Igualdad de umbrales y cargas", fixed = TRUE)
    expect_false(grepl("loadings_intercepts|thresholds_loadings_intercepts", excel_text))
  }
  expect_identical(fit$fit_table, before)
})
