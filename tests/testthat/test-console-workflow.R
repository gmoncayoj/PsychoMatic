test_that("EFA prints the substantive report even when assigned", {
  data <- lavaan::HolzingerSwineford1939[paste0("x", 1:9)]
  output <- capture.output(result <- efa_auto(
    data, refine = "suggest", parallel_iter = 10, seed = 15,
    language = "eng"
  ))
  text <- paste(output, collapse = "\n")
  expect_match(text, "KMO")
  expect_match(text, "Bartlett")
  expect_match(text, "LOADINGS")
  expect_match(text, "REFERENCES")
  expect_false(any(grepl('^\\[1\\] "', output)))
  expect_output(print(result, digits = 4), "[0-9]\\.[0-9]{4}")
})

test_that("CFA defaults to actual lavaan scaled fit measures", {
  data <- lavaan::HolzingerSwineford1939
  model <- "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6\nspeed =~ x7 + x8 + x9"
  output <- capture.output(result <- cfa_auto(data, model, estimator = "MLR"))
  expect_true(length(output) > 5)
  measures <- lavaan::fitMeasures(result$fit)
  for (metric in c("cfi", "tli", "rmsea")) {
    row <- result$fit_indices[result$fit_indices$Metric == metric, ]
    expect_equal(row$Value, unname(measures[paste0(metric, ".scaled")]))
    expect_identical(row$Variant, "scaled")
  }
})

test_that("reliability explains score scope and computation on assignment", {
  data <- lavaan::HolzingerSwineford1939[paste0("x", 1:3)]
  output <- capture.output(result <- reliability_auto(data, level = "continuous", language = "eng"))
  expect_match(paste(output, collapse = "\n"), "covariance")
  expect_match(paste(output, collapse = "\n"), "All selected items")
  expect_false(any(grepl("<missing>|<NA>", output)))
  expect_output(print(result, digits = 4), "[0-9]\\.[0-9]{4}")
})

test_that("invariance prints compact substantive tables on assignment", {
  output <- capture.output(result <- factorial_invariance_auto(
    lavaan::HolzingerSwineford1939, "school",
    "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6\nspeed =~ x7 + x8 + x9",
    estimator = "ML", language = "eng"
  ))
  text <- paste(output, collapse = "\n")
  expect_match(text, "CFI")
  expect_match(text, "REFERENCES")
  expect_match(text, "10.1007/s11336-016-9506-0", fixed = TRUE)
  expect_false(any(grepl("requested_constraints|effective_constraints", output)))
})

test_that("DIF prints wrapped references rather than a raw list", {
  data <- data.frame(i1 = rep(c(0, 1, 1, 0), 30), i2 = rep(c(1, 0, 1), 40))
  output <- capture.output(result <- dif_auto(data, rep(c("A", "B"), each = 60)))
  expect_match(paste(output, collapse = "\n"), "REFERENCES")
  expect_false(any(grepl('^\\$|^attr\\(', output)))
})

test_that("direct and combined reports contain reader-facing sections", {
  skip_if_not_installed("openxlsx")
  data <- lavaan::HolzingerSwineford1939[paste0("x", 1:3)]
  path <- tempfile(fileext = ".xlsx")
  output <- capture.output(rel <- reliability_auto(
    data, level = "continuous", report = TRUE, file_name = path,
    language = "eng", digits = 4
  ))
  expect_true(file.exists(path))
  expect_true("Contents" %in% openxlsx::getSheetNames(path))
  combined <- tempfile(fileext = ".xlsx")
  psychomatic_report(rel, rel, format = "excel", file_name = combined)
  sheets <- openxlsx::getSheetNames(combined)
  expect_false(any(grepl("audit_|result_00", sheets)))
  expect_true(length(sheets) < 16)
  expect_error(reliability_auto(data, report = TRUE, file_name = path), "exist")
})

test_that("descriptives and correlations use console digits without rounding data", {
  data <- data.frame(a = c(1.11, 2.22, 4.87, 3.19), b = c(3.17, 5.23, 6.93, 7.43))
  output <- capture.output(desc <- desc_auto(data, digits = 4, language = "eng"))
  expect_match(paste(output, collapse = "\n"), "2.8475", fixed = TRUE)
  expect_equal(desc$Mean[1], mean(data$a))
  output <- capture.output(cor <- cormat(data, digits = 4))
  expect_match(paste(output, collapse = "\n"), sprintf("%.4f", stats::cor(data)[1, 2]), fixed = TRUE)
  expect_equal(as.numeric(cor), as.numeric(stats::cor(data)))
  expect_length(capture.output(cormat(data, verbose = FALSE)), 0)
})

test_that("failed ordinal invariance explains empty categories visibly", {
  skip_if_not_installed("semTools")
  set.seed(915)
  data <- as.data.frame(matrix(sample(1:5, 2400, TRUE), ncol = 6))
  names(data) <- paste0("i", 1:6)
  data$group <- rep(c("A", "B"), each = 200)
  data$i1[data$group == "A" & data$i1 == 5] <- 4
  output <- capture.output(result <- ordinal_invariance_auto(
    data, "group", "f =~ i1 + i2 + i3 + i4 + i5 + i6", estimator = "WLSMV"
  ))
  expect_match(paste(output, collapse = "\n"), "empty|vacia")
  expect_false(any(result$fit_table$decision == "supported"))
})

test_that("partial invariance prints and exports directly", {
  .pm_skip_extended()
  skip_if_not_installed("semTools")
  skip_if_not_installed("openxlsx")
  path <- tempfile(fileext = ".xlsx")
  output <- capture.output(result <- invariance_partial(
    lavaan::HolzingerSwineford1939, "school",
    "visual =~ x1 + x2 + x3\ntextual =~ x4 + x5 + x6\nspeed =~ x7 + x8 + x9",
    free = "x2 ~ 1", rationale = "Teaching example of a prespecified intercept release.",
    estimator = "ML", measurement = "continuous", report = TRUE,
    file_name = path, language = "eng"
  ))
  expect_true(file.exists(path))
  expect_match(paste(output, collapse = "\n"), "PARTIAL INVARIANCE")
})

test_that("Spanish EFA reports have unique headings and export to Word", {
  .pm_skip_extended()
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  result <- efa_auto(lavaan::HolzingerSwineford1939[paste0("x", 1:9)],
                    refine = "suggest", parallel_iter = 10, seed = 916, verbose = FALSE)
  result$sugerencias <- c(result$sugerencias, "Heywood sentinel for methodological review")
  report <- PsychoMatic:::.pm_analysis_sections(result, "esp")
  expect_true(any(grepl("Heywood sentinel", unlist(report$sections))))
  expect_true(all(vapply(report$sections, function(x) anyDuplicated(names(x)) == 0L, logical(1))))
  path <- tempfile(fileext = ".docx")
  expect_no_error(psychomatic_report(result, format = "word", language = "esp", file_name = path))
  expect_true(file.exists(path))
  text <- officer::docx_summary(officer::read_docx(path))$text
  expect_true(any(grepl("N items", text, fixed = TRUE)))
})

test_that("reader reports retain comparison recommendations and alignment item summaries", {
  comparison <- structure(list(language = "eng", fit = data.frame(model = "M1", cfi = .95),
    recommendation = list(model = "M1", basis = "fit_indices", message = "Model M1 is the descriptive candidate.")),
    class = "psychomatic_model_comparison")
  report <- PsychoMatic:::.pm_analysis_sections(comparison)
  expect_match(report$sections$Recommendation$Text, "Model M1")
  alignment <- structure(list(language = "eng", item_results = list(
    loadings = list(summary = data.frame(item = "i1", noninvariant = TRUE)),
    intercepts = list(summary = data.frame(item = "i1", noninvariant = FALSE)))),
    class = "inv_align_auto_result")
  report <- PsychoMatic:::.pm_analysis_sections(alignment)
  expect_identical(report$sections$Loadings$noninvariant, TRUE)
  expect_identical(report$sections$Intercepts$noninvariant, FALSE)
})

test_that("reader workbooks do not reference nonexistent drawing parts", {
  skip_if_not_installed("openxlsx")
  path <- tempfile(fileext = ".xlsx")
  psychomatic_report(data.frame(Item = "i1", value = pi), file_name = path)
  directory <- tempfile("report-parts-")
  dir.create(directory)
  utils::unzip(path, exdir = directory)
  relations <- list.files(file.path(directory, "xl", "worksheets", "_rels"), full.names = TRUE)
  text <- unlist(lapply(relations, readLines, warn = FALSE))
  expect_false(any(grepl("/drawing\"|/vmlDrawing\"", text)))
})
