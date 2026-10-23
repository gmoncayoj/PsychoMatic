test_that("compact bilingual Word and Excel exports run in the distribution suite", {
  skip_if_not_installed("officer")
  skip_if_not_installed("openxlsx")
  folder <- tempfile("compact-reports-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  for (language in c("esp", "eng")) {
    paths <- psychomatic_report(data.frame(Probe = "precision", value = pi),
      format = c("word", "excel"), language = language,
      file_name = file.path(folder, language))
    expect_true(all(file.exists(paths)))
    text <- officer::docx_summary(officer::read_docx(paths[["word"]]))$text
    expect_true(any(grepl("precision", text, fixed = TRUE)))
    sheets <- openxlsx::getSheetNames(paths[["excel"]])
    # Title rows make read.xlsx infer text columns. Select the actual data row
    # by its stable label before parsing its numeric cell, without fixed rows.
    tab <- openxlsx::read.xlsx(paths[["excel"]], sheet = sheets[[2L]])
    cell <- tab[[2L]][which(tab[[1L]] == "precision")]
    expect_equal(as.numeric(cell), pi, tolerance = 1e-14)
  }
})

test_that("native network bootstrap and permutation stages run with a small budget", {
  skip_if_not_installed("EGAnet")
  skip_if_not_installed("bootnet")
  skip_if_not_installed("NetworkComparisonTest")
  d <- .pm_fixture_with_seed(493, {
    n <- 300L
    z <- matrix(rnorm(n * 3L), n, 3L)
    x <- as.data.frame(sapply(1:12, function(j) .75 * z[, ceiling(j / 4)] + rnorm(n, sd = .66)))
    names(x) <- paste0("q", 1:12)
    x$sex <- rep(c("F", "M"), length.out = n)
    x
  })
  x <- network_auto(d, items = 1:12, group = "sex", uva = NULL,
    wording = FALSE, hierarchy = FALSE, n_boot = 5L, n_perm = 3L,
    consensus_iter = 5L, n_cores = 1L, seed = 493L, verbose = FALSE)
  active <- x$results$stages$status != "skipped"
  expect_true(all(x$results$stages$status[active] == "completed"),
    info = paste(x$results$stages$reason, collapse = "; "))
  expect_equal(x$results$bootstrap$structural$valid, 5)
  expect_equal(x$results$invariance$p_resolution, .25)
  loadings <- x$results$invariance$loadings
  expect_equal(nrow(loadings), 12L)
  expect_true("valid_permutations" %in% names(loadings))
  expect_length(loadings$valid_permutations, 12L)
  expect_true(all(is.finite(loadings$valid_permutations)))
  expect_true(all(loadings$valid_permutations == 3L))
  global <- x$results$invariance$global
  expect_equal(nrow(global), 2L)
  expect_true(all(c("p_raw", "p_BH") %in% names(global)))
  expect_length(global$p_raw, 2L)
  expect_length(global$p_BH, 2L)
  expect_true(all(is.finite(global$p_raw)))
  expect_true(all(is.finite(global$p_BH)))
  expect_equal(global$p_BH, stats::p.adjust(global$p_raw, "BH"))
  diagnostics <- network_auto(d, items = 1:12, uva = NULL,
    bootstrap = FALSE, invariance = FALSE, consensus_iter = 5L,
    n_cores = 1L, seed = 493L, verbose = FALSE)
  stages <- diagnostics$results$stages
  selected <- stages$stage %in% c("riEGA", "hierarchy")
  expect_equal(sum(selected), 2L)
  expect_true(all(stages$status[selected] == "completed"),
    info = paste(stages$reason[selected], collapse = "; "))
})
