test_that("bundled example datasets have stable CRAN-friendly shapes", {
  data("psychomatic_continuous", package = "PsychoMatic")
  data("psychomatic_ordinal", package = "PsychoMatic")
  data("psychomatic_multigroup", package = "PsychoMatic")

  expect_s3_class(psychomatic_continuous, "data.frame")
  expect_equal(dim(psychomatic_continuous), c(120L, 6L))
  expect_true(all(vapply(psychomatic_continuous, is.numeric, logical(1))))

  expect_s3_class(psychomatic_ordinal, "data.frame")
  expect_equal(dim(psychomatic_ordinal), c(120L, 6L))
  expect_true(all(unlist(psychomatic_ordinal, use.names = FALSE) %in% 1:5))

  expect_s3_class(psychomatic_multigroup, "data.frame")
  expect_equal(nrow(psychomatic_multigroup), 140L)
  expect_true("group" %in% names(psychomatic_multigroup))
  expect_equal(length(unique(psychomatic_multigroup$group)), 2L)
  expect_true(all(vapply(psychomatic_multigroup[paste0("mg", 1:6)], is.numeric, logical(1))))
})

test_that("screen_items flags basic item-quality problems", {
  items <- data.frame(
    item1 = c(1, 2, 3, 4, NA, 5),
    item2 = c(5, 4, 3, 2, 1, 1),
    item3 = c(1, 1, 1, 1, 1, 1)
  )

  screened <- screen_items(items)

  expect_s3_class(screened, "psychomatic_item_screen")
  expect_true(all(c(
    "item", "n", "missing_pct", "sd", "skewness", "kurtosis",
    "zero_variance", "floor_pct", "ceiling_pct", "item_rest_correlation",
    "probable_reverse"
  ) %in% names(screened)))
  expect_true(screened$zero_variance[screened$item == "item3"])
  expect_true(screened$probable_reverse[screened$item == "item2"])
  removed_alias <- paste0("check", "_items")
  expect_false(removed_alias %in% getNamespaceExports("PsychoMatic"))
})

test_that("reverse_items recodes selected and detected reverse-keyed items", {
  items <- data.frame(
    item1 = 1:5,
    item2 = 5:1,
    item3 = c(1, 2, 3, 4, 5)
  )

  manual <- reverse_items(items, items = "item2", min_score = 1, max_score = 5, overwrite = TRUE)
  expect_equal(manual$item2, 1:5)
  expect_equal(attr(manual, "reversed_items"), "item2")

  detected <- reverse_items(items, detect = TRUE, min_score = 1, max_score = 5, overwrite = TRUE)
  expect_true("item2" %in% attr(detected, "reversed_items"))
})

test_that("score_scale computes scale scores with missing-data rules", {
  items <- data.frame(
    item1 = c(1, 2, NA, 4),
    item2 = c(2, 3, 4, 5),
    item3 = c(3, 4, 5, 5)
  )

  scored <- score_scale(
    items,
    keys = list(total = c("item1", "item2", "item3")),
    method = "mean",
    min_valid = 2 / 3,
    impute = "person_mean"
  )

  expect_s3_class(scored, "psychomatic_scores")
  expect_true("total" %in% names(scored$scores))
  expect_equal(nrow(scored$scores), nrow(items))
  expect_true(is.data.frame(scored$reliability))
  expect_true(all(c("scale", "n_items", "alpha") %in% names(scored$reliability)))
})

test_that("psychomatic_report writes an Excel report to a temporary path", {
  skip_if_not_installed("openxlsx")

  items <- data.frame(
    item1 = c(1, 2, 3, 4, 5),
    item2 = c(2, 3, 4, 5, 5),
    item3 = c(1, 2, 2, 4, 5)
  )
  screened <- screen_items(items)
  scored <- score_scale(items, items = names(items))
  path <- tempfile("psychomatic-report-", fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)

  out <- psychomatic_report(
    screening = screened,
    scores = scored,
    format = "excel",
    file_name = path,
    language = "eng"
  )

  expect_equal(normalizePath(out, winslash = "/", mustWork = FALSE),
               normalizePath(path, winslash = "/", mustWork = FALSE))
  expect_true(file.exists(path))
})

test_that("dif_auto returns logistic DIF screening results", {
  set.seed(11)
  group <- rep(c("A", "B"), each = 40)
  theta <- rnorm(80)
  items <- data.frame(
    item1 = as.integer(theta + rnorm(80) > 0),
    item2 = as.integer(theta + rnorm(80) + ifelse(group == "B", 0.4, 0) > 0),
    item3 = as.integer(theta + rnorm(80) > 0)
  )

  dif <- dif_auto(items, group = group, alpha = .05)

  expect_s3_class(dif, "psychomatic_dif")
  expect_true(all(c("item", "p_uniform", "p_nonuniform", "flagged") %in% names(dif$results)))
  expect_equal(nrow(dif$results), 3L)
  expect_true(is.data.frame(dif$interpretation))
  expect_true(all(c(
    "item", "dif_type", "decision", "severity", "recommendation"
  ) %in% names(dif$interpretation)))
  expect_equal(nrow(dif$interpretation), 3L)
  expect_true(is.data.frame(dif$summary))
  expect_true(all(c(
    "n_items", "flagged_items", "uniform_flags", "nonuniform_flags",
    "flagged_pct", "recommendation"
  ) %in% names(dif$summary)))
  expect_true(is.data.frame(dif$references))
  expect_equal(nrow(dif$references), 5L)
  expect_true(all(c("method", "reference") %in% names(dif$references)))
  expect_true(any(grepl("logistic regression", dif$references$method)))
  expect_true(all(grepl("https://doi.org/", dif$references$reference)))
})

test_that("compare_models_auto compares candidate CFA models", {
  skip_if_not_installed("lavaan")

  models <- list(
    one_factor = "general =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9",
    three_factor = paste(
      "visual =~ x1 + x2 + x3",
      "textual =~ x4 + x5 + x6",
      "speed =~ x7 + x8 + x9",
      sep = "\n"
    )
  )
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]

  comparison <- compare_models_auto(hs, models = models, estimator = "ML", language = "eng")

  expect_s3_class(comparison, "psychomatic_model_comparison")
  expect_equal(nrow(comparison$fit), 2L)
  expect_true(all(c("model", "converged", "cfi", "rmsea", "srmr", "aic", "bic") %in% names(comparison$fit)))
  expect_true(comparison$recommendation$model %in% names(models))
})

test_that("compare_models_auto supports WLSMV with scaled fit-index recommendations", {
  skip_if_not_installed("lavaan")

  data("psychomatic_ordinal", package = "PsychoMatic")
  items <- psychomatic_ordinal[1:100, ]
  models <- list(
    one_factor = "f =~ ord1 + ord2 + ord3 + ord4 + ord5 + ord6",
    two_factor = paste(
      "f1 =~ ord1 + ord2 + ord3",
      "f2 =~ ord4 + ord5 + ord6",
      sep = "\n"
    )
  )

  comparison <- compare_models_auto(items, models = models, estimator = "WLSMV", language = "eng")

  expect_s3_class(comparison, "psychomatic_model_comparison")
  expect_equal(comparison$estimator, "WLSMV")
  expect_true(all(c(
    "cfi_scaled", "tli_scaled", "rmsea_scaled",
    "comparison_cfi", "comparison_tli", "comparison_rmsea", "fit_score"
  ) %in% names(comparison$fit)))
  expect_true(all(is.na(comparison$fit$aic)))
  expect_true(all(is.na(comparison$fit$bic)))
  expect_equal(comparison$recommendation$basis, "fit_indices")
  expect_true(comparison$recommendation$model %in% names(models))
})

test_that("irt_auto fits a small graded response model when mirt is available", {
  skip_if_not_installed("mirt")

  data("psychomatic_ordinal", package = "PsychoMatic")
  result <- irt_auto(psychomatic_ordinal[1:80, 1:4], model = "grm")

  expect_s3_class(result, "psychomatic_irt")
  expect_equal(result$itemtype, "graded")
  expect_true(is.data.frame(result$item_parameters))
  expect_true(is.data.frame(result$interpretation))
  expect_true(all(c(
    "item", "discrimination", "discrimination_level", "difficulty_center",
    "threshold_min", "threshold_max", "threshold_ordered", "recommendation"
  ) %in% names(result$interpretation)))
  expect_true(is.data.frame(result$model_summary))
  expect_true(all(c(
    "model", "itemtype", "n_items", "items_to_review",
    "mean_discrimination", "recommendation"
  ) %in% names(result$model_summary)))
  expect_true(is.data.frame(result$references))
  expect_equal(nrow(result$references), 5L)
  expect_true(all(c("method", "reference") %in% names(result$references)))
  expect_true(any(grepl("graded response", result$references$method)))
  expect_true(all(grepl("https://doi.org/", result$references$reference)))
})
