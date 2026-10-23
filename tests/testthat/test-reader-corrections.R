test_that("correlation title appears once in direct console output", {
  txt <- capture.output(cormat(data.frame(a=1:6, b=c(2,1,3,6,4,5)), type="pearson"))
  expect_equal(sum(trimws(txt) == "PEARSON CORRELATIONS"), 1L)
})

test_that("screening supports bilingual console and safe Excel Word reports", {
  d <- data.frame(a=c(1,1,2,3,NA), b=c(3,2,2,1,1))
  en <- capture.output(x <- screen_items(d))
  es <- capture.output(y <- screen_items(d, language="esp"))
  expect_true(any(grepl("ITEM SCREENING", en, fixed=TRUE)))
  expect_true(any(grepl("REVISION DE ITEMS", es, fixed=TRUE)))
  expect_identical(names(x), names(y))
  expect_equal(x$floor_pct, c(.5,.4))
  expect_equal(y$floor_pct, x$floor_pct)
  expect_output(screen_items(d, verbose=FALSE), NA)
  expect_error(screen_items(d, language="invalid"))
  path <- tempfile("screen-")
  suppressMessages(z <- screen_items(d, language="esp", report=TRUE,
    format=c("excel","word"), file_name=path, verbose=FALSE))
  expect_true(all(file.exists(attr(z, "report_paths"))))
  book <- lapply(openxlsx::getSheetNames(paste0(path,".xlsx")), function(sheet) {
    t <- openxlsx::read.xlsx(paste0(path,".xlsx"), sheet=sheet)
    c(names(t), unlist(t))
  })
  expect_true(any(grepl("Suelo", unlist(book), fixed=TRUE)))
  doc <- officer::docx_summary(officer::read_docx(paste0(path,".docx")))
  expect_true(any(grepl("Suelo", doc$text, fixed=TRUE)))
  expect_error(screen_items(d, report=TRUE, file_name=path, verbose=FALSE), "exist")
})

test_that("invariance successful conclusion is concise and failure remains explicit", {
  tab <- data.frame(level=c("configural","strict"), decision="supported")
  expect_identical(PsychoMatic:::.ifa_conclusion_text(tab,"eng"),
    "Fit-change heuristics were sequentially compatible through the 'strict' level.")
  tab$decision <- "not_evaluable"
  expect_match(PsychoMatic:::.ifa_conclusion_text(tab,"eng"), "must not be interpreted")
})

test_that("alignment reader keeps group identities and separates flagged groups", {
  items <- data.frame(item=c("i1","i2"), joint_parameter=c(.7,.8),
    distinct_parameters=c(1,2), noninvariant_groups=c(0,1),
    percent_noninvariant=c(0,50), flagged_groups=c("","Peru"), invariant=c(TRUE,FALSE))
  x <- structure(list(summary=data.frame(parameter="loadings", R2=.95,
    percent_noninvariant=25, cutoff=25, invariant=TRUE),
    item_results=list(loadings=list(summary=items), intercepts=list(summary=items)),
    group_parameters=data.frame(alpha0=c(0,.3),psi0=c(1,1.2),row.names=c("Chile","Peru")),
    dimensions=list(group_names=c("Chile","Peru")), references=NULL), class="inv_align_auto_result")
  s <- PsychoMatic:::.pm_analysis_sections(x,"eng")$sections
  expect_identical(s[["Group parameters"]]$Group, c("Chile","Peru"))
  expect_equal(s[["Group parameters"]][["Latent SD"]], c(1,1.2))
  expect_false(any(grepl("_",names(s$Loadings),fixed=TRUE)))
  expect_true("Flagged groups: loadings" %in% names(s))
  expect_true(any(grepl("Peru",unlist(s[["Flagged groups: loadings"]]),fixed=TRUE)))
  expect_equal(x$item_results$loadings$summary$joint_parameter,c(.7,.8))
})
