# Fail early with the native load error if a runner cannot support the reports.
for (package in c("officer", "flextable", "openxlsx")) {
  tryCatch(loadNamespace(package), error = function(error) {
    stop("CI cannot load ", package, ": ", conditionMessage(error), call. = FALSE)
  })
}

word <- tempfile(fileext = ".docx")
excel <- tempfile(fileext = ".xlsx")
svg <- tempfile(fileext = ".svg")
document <- officer::read_docx()
document <- flextable::body_add_flextable(document,
  flextable::flextable(data.frame(value = pi)))
print(document, target = word)
openxlsx::write.xlsx(data.frame(value = pi), excel)
stopifnot(nrow(officer::docx_summary(officer::read_docx(word))) > 0L,
  isTRUE(all.equal(openxlsx::read.xlsx(excel)$value, pi)))
grDevices::svg(svg, width = 3, height = 3)
tryCatch(graphics::plot(1:3, 1:3, main = "CI graphics"),
  finally = grDevices::dev.off())
stopifnot(all(file.exists(c(word, excel, svg))),
  all(file.info(c(word, excel, svg))$size > 0))
unlink(c(word, excel, svg))
cat("PASS: Word, Excel and SVG exports are supported by this runner\n")
