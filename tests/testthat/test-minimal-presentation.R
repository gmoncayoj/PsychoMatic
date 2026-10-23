test_that("correlations keep pairwise N internally without printing its matrix", {
  d <- data.frame(a = c(1, 2, 3, 4), b = c(2, NA, 5, 8))
  txt <- capture.output(r <- cormat(d, type = "pearson"))
  expect_false(any(grepl("N POR PAR", txt, fixed = TRUE)))
  expect_equal(attr(r, "pairwise_n")[1, 2], 3)
})

test_that("EFA variance uses percentages with explicit denominators", {
  x <- efa_auto(lavaan::HolzingerSwineford1939[paste0("x", 1:9)],
    parallel_iter = 10, seed = 15, refine = "suggest", verbose = FALSE)
  sections <- PsychoMatic:::.pm_analysis_sections(x, "esp")$sections
  table <- sections[["Varianza explicada"]]
  expect_true("Varianza total (%)" %in% names(table))
  expect_equal(table[["Varianza total (%)"]],
    unname(100 * x$varianza_explicada["Proportion Var", ]))
  expect_true(any(grepl("percentil", unlist(sections), fixed = TRUE)))
})

test_that("ordinal CFA reports observed-code Mardia and hides redundant hierarchy", {
  set.seed(234)
  f <- rnorm(350)
  d <- as.data.frame(replicate(4, as.integer(cut(f + rnorm(350),
    breaks = c(-Inf, -1, -.3, .3, 1, Inf)))))
  names(d) <- paste0("i", 1:4)
  x <- cfa_auto(d, "f =~ i1 + i2 + i3 + i4", verbose = FALSE)
  s <- PsychoMatic:::.pm_analysis_sections(x, "esp")$sections
  expect_identical(x$mardia$status, "ok")
  expect_true(any(grepl("ordinal", s[["Nota sobre Mardia"]]$Text)))
  loads <- s[["Cargas factoriales estandarizadas (IC 95%)"]]
  expect_true(is.data.frame(loads))
  expect_false(any(c("Etapa", "Nivel factorial") %in% names(loads)))
  expect_false("Etapa" %in% names(s$Confiabilidad))
  x$conclusion <- list(reason = "The CFA solution is inadmissible: post_check_failed, negative_variance.")
  x$diagnostics$admissibility$admissible <- FALSE
  x$diagnostics$admissibility$reason_codes <- c("post_check_failed", "negative_variance")
  text <- PsychoMatic:::.pm_analysis_sections(x, "esp")$sections$Conclusion$Text
  expect_false(any(grepl("The CFA|negative_variance", text)))
  expect_true(any(grepl("varianza negativa", text)))
})

test_that("invariance delta labels identify the difference without altering values", {
  d <- data.frame(delta_cfi = -.01, delta_rmsea = .005, delta_srmr = .004)
  es <- PsychoMatic:::.pm_label_table(d, "esp")
  en <- PsychoMatic:::.pm_label_table(d, "eng")
  expect_identical(names(es), paste0("\u0394", c("CFI", "RMSEA", "SRMR")))
  expect_identical(names(en), names(es))
  expect_equal(unname(unlist(es)), c(-.01, .005, .004))
})
