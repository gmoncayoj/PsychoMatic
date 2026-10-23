test_that('manuscript interpretation is the final section in either language', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(), verbose = FALSE,
    rotation_args = list(rstarts = 3), language = 'esp')
  esp <- PsychoMatic:::.pm_analysis_sections(out, 'esp')
  eng <- PsychoMatic:::.pm_analysis_sections(out, 'eng')
  expect_identical(tail(names(esp$sections), 2), c('Referencias', 'Interpretaci\u00f3n'))
  expect_identical(tail(names(eng$sections), 2), c('References', 'Interpretation'))
  expect_match(tail(esp$sections, 1)[[1]]$Text, '400')
  expect_match(tail(eng$sections, 1)[[1]]$Text, '400')
  text <- capture.output(print(out))
  expect_match(paste(text, collapse = ' '), 'INTERPRETACI\u00d3N')
  expect_false(grepl('Este borrador describe', paste(text, collapse = ' ')))
  audit <- PsychoMatic:::audit_report(out)
  expect_true(is.data.frame(audit$settings))
})

test_that('reader reports record actual rotation controls', {
  out <- esem_auto(esem_fixture(), factors = esem_factors(), rotation = 'geomin',
    rotation_args = list(rstarts = 3, geomin_epsilon = .01), verbose = FALSE)
  sections <- PsychoMatic:::.pm_analysis_sections(out, 'eng')
  meta <- sections$sections[['Data and decisions']]
  expect_true(all(c('Random starts', 'Seed', 'Geomin epsilon') %in% meta$Information))
  expect_identical(meta$Value[meta$Information == 'Geomin epsilon'], '0.01')
})

test_that('ESEM report paths are protected and draft survives Word and Excel export', {
  .pm_skip_extended()
  skip_if_not_installed('officer'); skip_if_not_installed('openxlsx')
  out <- esem_auto(esem_fixture(), factors = esem_factors(), verbose = FALSE,
    rotation_args = list(rstarts = 3))
  folder <- tempfile('esem-report-'); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  before <- list.files(folder, all.files = TRUE)
  old <- setwd(folder); on.exit(setwd(old), add = TRUE)
  invisible(esem_auto(esem_fixture(), factors = esem_factors(), verbose = FALSE,
    rotation_args = list(rstarts = 3)))
  expect_identical(list.files(folder, all.files = TRUE), before)
  for (lang in c('eng', 'esp')) {
    paths <- psychomatic_report(out, format = c('word', 'excel'),
      file_name = file.path(folder, lang), language = lang)
    expect_true(all(file.exists(paths)))
    sheets <- openxlsx::getSheetNames(paths[['excel']])
    tab <- openxlsx::read.xlsx(paths[['excel']], sheet = tail(sheets, 1))
    cells <- as.character(unlist(tab, use.names = FALSE))
    expect_true(all((if (lang == 'eng') c('Preparation','Item map','Factor map','Versions') else
      c('Preparaci\u00f3n','Mapa de \u00edtems','Mapa de factores','Versiones')) %in% cells))
    interpretation <- if (lang == 'eng') 'Interpretation' else 'Interpretaci\u00f3n'
    reference <- if (lang == 'eng') 'References' else 'Referencias'
    expect_true(match(reference, tab[[1]]) < match(interpretation, tab[[1]]))
    present <- rowSums(!is.na(tab) & as.matrix(tab) != '') > 0
    last_row <- as.character(unlist(tail(tab[present, , drop = FALSE], 1)))
    expect_match(paste(last_row, collapse = ' '), '400')
    expect_match(paste(unlist(tab), collapse = ' '), '400')
    expect_false(grepl('Este borrador describe|This draft describes', paste(unlist(tab), collapse = ' ')))
    wd <- officer::read_docx(paths[['word']])
    expect_match(paste(officer::docx_summary(wd)$text, collapse = ' '), '400')
    body <- officer::docx_body_xml(wd)
    paragraphs <- xml2::xml_text(xml2::xml_find_all(body, './/w:p',
      ns = xml2::xml_ns(body)))
    nonempty <- paragraphs[nzchar(trimws(paragraphs))]
    expect_match(tail(nonempty, 1), '400')
    expect_false(grepl('Este borrador describe|This draft describes', tail(nonempty, 1)))
    expect_error(psychomatic_report(out, file_name = file.path(folder, lang)), 'exist|overwrite')
  }
})
