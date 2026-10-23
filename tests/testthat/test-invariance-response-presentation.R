test_that("traditional invariance shows observed response options without a coding table", {
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 7130L)
  fit <- factorial_invariance_auto(d, "group", "F =~ i1+i2+i3+i4+i5+i6",
    estimator = "WLSMV", verbose = FALSE)
  for (language in c("eng", "esp")) {
    sections <- .pm_analysis_sections(fit, language)$sections
    context <- sections[[if (language == "eng") "Data and method" else "Datos y metodo"]]
    label <- if (language == "eng") "Response options (observed)" else "Opciones de respuesta observadas"
    expect_true(label %in% context[[1L]])
    expect_identical(context[[2L]][context[[1L]] == label], "5")
    expect_false(any(context[[1L]] %in% c("Measurement", "Nivel de medicion")))
    expect_false(any(names(sections) %in% c("Score coding", "Codificacion de puntuaciones")))
    display <- fit
    display$effective_settings$language <- language
    text <- paste(capture.output(print(display)), collapse = "\n")
    expect_match(text, label, fixed = TRUE)
    expect_false(grepl("SCORE CODING|CODIFICACION DE PUNTUACIONES", text))
  }
  expect_identical(fit$measurement, "ordinal")
  expect_identical(fit$effective_settings$score_coding, list())
})

test_that("observed response counts handle different item categories and old results", {
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 7131L)
  d$i1 <- ordered(pmin(as.numeric(as.character(d$i1)), 4), levels = 1:5)
  fit <- factorial_invariance_auto(d, "group", "F =~ i1+i2+i3+i4+i5+i6",
    estimator = "WLSMV", verbose = FALSE)
  expect_equal(fit$item_info$category_counts[["i1"]], 4)
  for (language in c("eng", "esp")) {
    sections <- .pm_analysis_sections(fit, language)$sections
    context <- sections[[if (language == "eng") "Data and method" else "Datos y metodo"]]
    expect_true(any(grepl("4, 5", context[[2L]], fixed = TRUE)))
    expect_true(any(grepl(if (language == "eng") "varies by item" else "segun el item",
      context[[2L]], fixed = TRUE)))
  }
  old <- fit
  old$item_info$category_counts <- NULL
  context <- .pm_analysis_sections(old, "eng")$sections[["Data and method"]]
  expect_true("Not available" %in% context[[2L]])
})

test_that("invariance reports share the concise response summary", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  d <- simulate_invariance_fixture(150L, ordinal = TRUE, seed = 7132L)
  fit <- factorial_invariance_auto(d, "group", "F =~ i1+i2+i3+i4+i5+i6",
    estimator = "WLSMV", verbose = FALSE)
  directory <- tempfile("response-summary-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  for (language in c("eng", "esp")) {
    paths <- psychomatic_report(fit, format = c("word", "excel"),
      language = language, file_name = file.path(directory, language))
    word <- paste(officer::docx_summary(officer::read_docx(paths[["word"]]))$text, collapse = " ")
    excel <- paste(unlist(openxlsx::read.xlsx(paths[["excel"]], sheet = 2, colNames = FALSE)), collapse = " ")
    label <- if (language == "eng") "Response options (observed)" else "Opciones de respuesta observadas"
    for (text in c(word, excel)) {
      expect_match(text, label, fixed = TRUE)
      expect_false(grepl("Score coding|Codificacion de puntuaciones", text))
    }
  }
})
