test_that("reader adequacy uses inclusive boundaries without hiding failures", {
  raw <- data.frame(metric = c("omega", "omega", "ave", "ave", "omega"),
    value = c(.70, .699, .50, .499, .95),
    status = c(rep("ok", 4), "not_evaluable"),
    reason = c(rep(NA_character_, 4), "Invalid solution"))
  out <- .pm_coefficient_evaluation(raw, "value", "metric", "eng")
  expect_equal(out$status, c("Adequate", "Inadequate", "Adequate", "Inadequate", "not_evaluable"))
  expect_equal(out$reason, c(">= 0.70", "< 0.70", ">= 0.50", "< 0.50", "Invalid solution"))
  expect_equal(raw$status[1], "ok")
  es <- .pm_coefficient_evaluation(raw, "value", "metric", "esp")
  expect_equal(es$status[1:2], c("Adecuado", "No adecuado"))
})

test_that("CFA note identifiers preserve suggestion order in both languages", {
  suggestions <- lapply(1:2, function(i) list(type="residual_covariance",
    syntax=paste0("R1 ~~ R", i + 1), mi=10 + i, epc=.1,
    same_factor=TRUE, note=paste("Explanation", i)))
  x <- structure(list(respecification_suggestions=suggestions), class="cfa_auto")
  en <- .pm_analysis_sections(x, "eng")$sections
  es <- .pm_analysis_sections(x, "esp")$sections
  expect_equal(en[["Suggestions for review"]]$No, 1:2)
  expect_equal(es[["Sugerencias para revisar"]][["Nro."]], 1:2)
  expect_equal(en$Notes$Text, c("1. Explanation 1", "2. Explanation 2"))
  expect_false("Note" %in% names(en[["Suggestions for review"]]))
  exported <- .pm_analysis_sections(x, "eng", destination = "report")$sections
  expect_equal(exported[["Suggestions for review"]]$Note, c("Explanation 1", "Explanation 2"))
  expect_false("Notes" %in% names(exported))
  expect_false("No" %in% names(exported[["Suggestions for review"]]))
})

test_that("Spanish reader translates metadata and missing-case warnings", {
  tab <- data.frame(Information=c("Missing", "Criteria", "Keying"),
    Value=c("listwise", "categorical_descriptive_extrapolation", "none"))
  es <- .pm_label_table(tab, "esp")
  expect_equal(es$Valor, c("Casos completos", "Orientaciones descriptivas para items ordinales", "Ninguna"))
  x <- structure(list(warnings=data.frame(message=
    "161 cases were removed by the item missing-data policy.")), class="cfa_auto")
  s <- .pm_analysis_sections(x, "esp")$sections
  expect_match(s$Advertencias$Text, "Se excluyeron 161 casos", fixed=TRUE)
  expect_equal(x$warnings$message, "161 cases were removed by the item missing-data policy.")
})
