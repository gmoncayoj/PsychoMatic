hs_model <- paste(
  "visual =~ x1 + x2 + x3",
  "textual =~ x4 + x5 + x6",
  "speed =~ x7 + x8 + x9",
  sep = "\n"
)

quiet_run <- function(expr) {
  result <- NULL

  invisible(capture.output(
    result <- suppressWarnings(
      suppressMessages(
        eval.parent(substitute(expr))
      )
    )
  ))

  result
}

test_that("desc_auto devuelve descriptivos para variables numericas", {
  data <- data.frame(
    item1 = c(1, 2, 3, 4, 5),
    item2 = c(2, 2, 3, 4, 5),
    item3 = c(5, 4, 3, 2, 1)
  )

  res <- quiet_run(desc_auto(data, digits = 2, language = "eng", report = FALSE))

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 3L)
  expect_true(all(c("Item", "Mean", "SD", "Skewness", "Kurtosis") %in% names(res)))
})

test_that("cormat calcula matrices y exporta Excel", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("psych")

  data <- data.frame(
    item1 = c(1, 2, 3, 4, 5, 4, 3, 2),
    item2 = c(2, 2, 3, 4, 4, 5, 3, 1),
    item3 = c(5, 4, 3, 2, 1, 2, 3, 4)
  )

  pearson <- cormat(data, type = "pearson")
  expect_true(is.matrix(pearson))
  expect_equal(dim(pearson), c(3L, 3L))
  expect_equal(unname(pearson), unname(stats::cor(data, use = "pairwise.complete.obs")))

  poly <- quiet_run(cormat(data, type = "poly"))
  expect_true(is.matrix(poly))
  expect_equal(dim(poly), c(3L, 3L))

  binary_data <- data.frame(
    item1 = c(0, 1, 1, 0, 1, 0, 1, 0),
    item2 = c(0, 1, 0, 0, 1, 1, 1, 0),
    item3 = c(1, 1, 0, 0, 1, 0, 1, 0)
  )

  tetra <- quiet_run(cormat(binary_data, type = "tetra"))
  expect_true(is.matrix(tetra))
  expect_equal(dim(tetra), c(3L, 3L))

  tmp <- tempfile("psychomatic-cormat-")
  on.exit(unlink(paste0(tmp, ".xlsx")), add = TRUE)

  quiet_run(cormat(data, type = "pearson", report = TRUE, file_name = tmp))

  expect_true(file.exists(paste0(tmp, ".xlsx")))
})

test_that("efa_auto genera un resultado estructurado y export_efa crea Excel", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("openxlsx")

  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3", "x4", "x5", "x6")]
  res <- quiet_run(efa_auto(hs, rotation = "oblique", verbose = FALSE, language = "eng"))

  expect_s3_class(res, "efa_auto")
  expect_gte(res$n_factores, 1)
  expect_true(!is.null(res$cargas))
  expect_true("rotation" %in% names(res))

  tmp <- tempfile("psychomatic-efa-")
  on.exit(unlink(paste0(tmp, ".xlsx")), add = TRUE)

  quiet_run(export_efa(res, format = "excel", file_name = tmp))

  expect_true(file.exists(paste0(tmp, ".xlsx")))
})

test_that("efa_auto normaliza valores de rotacion legados", {
  expect_equal(PsychoMatic:::.efa_normalize_rotation("oblicua"), "oblique")
  expect_equal(PsychoMatic:::.efa_normalize_rotation("ortogonal"), "orthogonal")
})

test_that("cfa_auto ajusta el ejemplo de Holzinger-Swineford", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("openxlsx")

  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  res <- quiet_run(cfa_auto(hs, model = hs_model, language = "eng"))

  expect_s3_class(res, "cfa_auto")
  expect_true(isTRUE(res$converged))
  expect_equal(sort(res$factors), sort(c("speed", "textual", "visual")))
  expect_true(!is.null(res$fit))
  expect_true(all(c("Index", "Value", "Evaluation") %in% names(res$fit_indices)))
  expect_true(all(c("Loading", "p_value", "CI_lower", "CI_upper", "Alert") %in% names(res$factor_loadings)))

  printed <- capture.output(print(res))
  expect_true(any(grepl("95% CI", printed, fixed = TRUE)))

  tmp <- tempfile("psychomatic-cfa-")
  on.exit(unlink(paste0(tmp, ".xlsx")), add = TRUE)

  quiet_run(export_cfa(res, format = "excel", file_name = tmp))

  expect_true(file.exists(paste0(tmp, ".xlsx")))
})

test_that("inv_align_auto en modo matrices devuelve un resumen valido", {
  res <- quiet_run(inv_align_auto(
    lambda = matrix(c(0.70, 0.80, 0.75, 0.72, 0.78, 0.74), nrow = 2, byrow = TRUE),
    nu = matrix(c(0.00, 0.00, 0.00, 0.10, 0.05, 0.08), nrow = 2, byrow = TRUE),
    language = "eng"
  ))

  expect_s3_class(res, "inv_align_auto_result")
  expect_true(is.data.frame(summary(res)))
  expect_equal(nrow(summary(res)), 2L)
  expect_true(is.character(res$report))
})

test_that("factorial_invariance_auto devuelve la tabla esperada de niveles", {
  skip_if_not_installed("lavaan")

  hs <- lavaan::HolzingerSwineford1939[, c("school", paste0("x", 1:9))]
  res <- quiet_run(factorial_invariance_auto(
    hs,
    group = "school",
    model = hs_model,
    language = "eng",
    report = FALSE
  ))

  expect_s3_class(res, "factorial_invariance_auto")
  expect_true(is.data.frame(res$fit_table))
  expect_equal(rownames(res$fit_table), c("configural", "metric", "scalar", "strict"))
  expect_true(is.character(res$conclusion) && nzchar(res$conclusion))
})
