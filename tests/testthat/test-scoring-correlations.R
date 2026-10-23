test_that("reverse scoring deduplicates keys and refuses suffix collisions", {
  d <- data.frame(q1 = 1:5, q1_r = 5:1)

  once <- reverse_items(
    d["q1"],
    items = c("q1", "q1"),
    min_score = 1,
    max_score = 5,
    overwrite = TRUE
  )

  expect_identical(once$q1, 5:1)
  expect_identical(attr(once, "reversed_items"), "q1")
  expect_error(
    reverse_items(d, items = "q1", min_score = 1, max_score = 5),
    "already exists"
  )
})

test_that("reverse scoring rejects responses outside the declared range", {
  d <- data.frame(q1 = c(0, 1, 5))

  expect_error(
    reverse_items(
      d,
      items = "q1",
      min_score = 1,
      max_score = 5,
      overwrite = TRUE
    ),
    "outside.*range"
  )
})

test_that("reverse scoring rejects non-finite responses but preserves NA", {
  for (bad in list(Inf, -Inf, NaN)) {
    expect_error(
      reverse_items(
        data.frame(q1 = c(1, bad, 5)),
        items = "q1",
        min_score = 1,
        max_score = 5,
        overwrite = TRUE
      ),
      "finite"
    )
  }

  out <- reverse_items(
    data.frame(q1 = c(1, NA, 5)),
    items = "q1",
    min_score = 1,
    max_score = 5,
    overwrite = TRUE
  )
  expect_equal(out$q1, c(5, NA, 1))
})

test_that("min_valid prorates eligible incomplete rows for every scoring policy", {
  d <- data.frame(
    q1 = c(1, 1, 1),
    q2 = c(2, NA, NA),
    q3 = c(3, 3, NA)
  )

  for (impute in c("none", "person_mean")) {
    mean_score <- score_scale(
      d,
      items = names(d),
      method = "mean",
      min_valid = 2 / 3,
      impute = impute
    )
    sum_score <- score_scale(
      d,
      items = names(d),
      method = "sum",
      min_valid = 2 / 3,
      impute = impute
    )

    expect_equal(mean_score$scores$total, c(2, 2, NA))
    expect_equal(sum_score$scores$total, c(6, 6, NA))
  }
})

test_that("scale scoring rejects non-finite values before scoring", {
  for (bad in list(Inf, -Inf, NaN)) {
    expect_error(
      score_scale(
        data.frame(q1 = c(1, bad), q2 = c(2, 3)),
        items = c("q1", "q2"),
        method = "mean",
        min_valid = 1 / 2,
        impute = "none"
      ),
      "finite"
    )
  }
})

test_that("duplicate scale labels preserve every result under unique keys", {
  keys <- structure(
    list(c("q1", "q2"), c("q2", "q3")),
    names = c("scale", "scale")
  )

  out <- score_scale(
    data.frame(q1 = 1:4, q2 = 2:5, q3 = 3:6),
    keys = keys
  )

  expect_equal(ncol(out$scores), 2L)
  expect_equal(nrow(out$reliability), 2L)
  expect_equal(length(unique(names(out$scores))), 2L)
  expect_equal(out$scale_map$label, c("scale", "scale"))
  expect_identical(out$scale_map$key, names(out$scores))
})

test_that("correlations expose pairwise diagnostics and reject ambiguous characters", {
  d <- data.frame(
    a = c(1, 2, NA, NA),
    b = c(1, 2, 3, NA),
    c = c(NA, NA, 1, 2)
  )

  r <- cormat(d, type = "pearson")

  expect_equal(attr(r, "pairwise_n")["a", "b"], 2L)
  expect_true("insufficient_overlap" %in% attr(r, "diagnostics")$code)
  expect_identical(attr(r, "requested_type"), "pearson")
  expect_identical(attr(r, "effective_type"), "pearson")
  expect_error(
    cormat(data.frame(a = c("low", "high"), b = 1:2), type = "poly"),
    "ordered"
  )
})

test_that("explicit ordinal levels preserve declared category order", {
  prepared <- PsychoMatic:::.cormat_prepare_data(
    data.frame(
      a = c("high", "low", "medium", NA),
      b = c(3, 1, 2, 1),
      stringsAsFactors = FALSE
    ),
    type = "poly",
    ordinal_levels = list(a = c("low", "medium", "high"))
  )

  expect_identical(prepared$a, c(3L, 1L, 2L, NA_integer_))
  expect_identical(prepared$b, c(3, 1, 2, 1))
})

test_that("continuous descriptives have fixed width and raw internal values", {
  set.seed(340)
  d <- as.data.frame(
    replicate(20, seq_len(5000) + stats::rnorm(5000), simplify = FALSE)
  )

  out <- desc_auto(d, digits = 1, language = "eng", report = FALSE)

  expect_lte(ncol(out), 12L)
  expect_equal(nrow(out), 20L)
  expect_true(is.data.frame(attr(out, "response_percentages")))
  expect_equal(nrow(attr(out, "response_percentages")), 0L)

  exact <- desc_auto(data.frame(q = c(0, 0, 1)), digits = 1, language = "eng")
  expect_equal(exact$Mean, 1 / 3, tolerance = 1e-14)
})

test_that("response percentages are long and limited to ordinal items", {
  d <- data.frame(
    q = c(1, 1, 2, NA, 2, 3, 3, 3),
    z = seq(0.1, 0.8, by = 0.1)
  )

  out <- desc_auto(d, language = "eng", report = FALSE)
  percentages <- attr(out, "response_percentages")

  expect_identical(
    names(percentages),
    c("item", "category", "n", "percent")
  )
  expect_identical(percentages$item, rep("q", 3L))
  expect_identical(percentages$category, c("1", "2", "3"))
  expect_identical(percentages$n, c(2L, 2L, 3L))
  expect_equal(percentages$percent, c(200 / 7, 200 / 7, 300 / 7))
})

test_that("Excel sheet names are unique without case sensitivity", {
  existing <- c("Results", "results_2")

  expect_identical(
    PsychoMatic:::.pm_excel_sheet_name("results", existing),
    "results_3"
  )
  long <- PsychoMatic:::.pm_excel_sheet_name(strrep("a", 40), character())
  expect_lte(nchar(long), 31L)
  expect_identical(
    PsychoMatic:::.pm_excel_sheet_name("a/b:c", character()),
    "a_b_c"
  )
})

test_that("Excel writers retain requested-to-effective sheet mappings", {
  skip_if_not_installed("openxlsx")
  target <- tempfile(fileext = ".xlsx")
  on.exit(unlink(target), add = TRUE)
  tables <- structure(
    list(data.frame(x = 1), data.frame(x = 2)),
    names = c("Results", "results")
  )

  written <- PsychoMatic:::.pm_write_excel_report(tables, target)
  mapping <- attr(written, "sheet_map")

  expect_identical(
    names(mapping),
    c("object_index", "display_label", "requested_sheet", "effective_sheet")
  )
  expect_identical(mapping$object_index, 1:2)
  expect_identical(mapping$display_label, c("Results", "results"))
  expect_identical(mapping$requested_sheet, c("Results", "results"))
  expect_identical(mapping$effective_sheet, c("Results", "results_2"))
})

test_that("safe paths add only the requested extension and never overwrite by default", {
  target <- tempfile()
  xlsx <- PsychoMatic:::.pm_safe_output_path(target, extension = ".xlsx")

  expect_true(endsWith(xlsx, ".xlsx"))
  expect_identical(
    xlsx,
    normalizePath(paste0(target, ".xlsx"), winslash = "/", mustWork = FALSE)
  )
  writeLines("existing", xlsx)
  on.exit(unlink(xlsx), add = TRUE)
  expect_error(
    PsychoMatic:::.pm_safe_output_path(target, extension = ".xlsx"),
    "exists"
  )
  expect_identical(
    PsychoMatic:::.pm_safe_output_path(target, extension = ".xlsx", overwrite = TRUE),
    xlsx
  )
})

test_that("cormat refuses to overwrite an existing explicit path", {
  skip_if_not_installed("openxlsx")
  target <- tempfile(fileext = ".xlsx")
  writeLines("sentinel", target)
  on.exit(unlink(target), add = TRUE)
  d <- data.frame(
    a = c(1, 2, 3, 4, 5, 4, 3, 2),
    b = c(2, 2, 3, 4, 4, 5, 3, 1)
  )

  expect_error(
    cormat(d, report = TRUE, file_name = target),
    "exists"
  )
  expect_identical(readLines(target), "sentinel")
})

test_that("desc_auto refuses to overwrite an existing explicit path", {
  skip_if_not_installed("writexl")
  target <- tempfile(fileext = ".xlsx")
  writeLines("sentinel", target)
  on.exit(unlink(target), add = TRUE)

  expect_error(
    desc_auto(
      data.frame(a = 1:5, b = 2:6),
      report = TRUE,
      file_name = target
    ),
    "exists"
  )
  expect_identical(readLines(target), "sentinel")
})

test_that("psychomatic_report refuses to overwrite an existing explicit path", {
  skip_if_not_installed("openxlsx")
  target <- tempfile(fileext = ".xlsx")
  writeLines("sentinel", target)
  on.exit(unlink(target), add = TRUE)

  expect_error(
    psychomatic_report(
      descriptives = data.frame(item = "q1", mean = 3),
      format = "excel",
      file_name = target
    ),
    "exists"
  )
  expect_identical(readLines(target), "sentinel")
})
