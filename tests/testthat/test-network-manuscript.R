manuscript_network_fixture <- function() {
  x <- list(language='eng',data_summary=list(n=314L,p_retained=3L,p_original=3L),
    effective_settings=list(corr='spearman',gamma=.5,consensus_iter=100L,seed=2026,
      seed_source='explicit',n_boot=100L,n_perm=499L,uva=NULL),
    conclusion=list(status='partial'),warnings=data.frame(stage='riEGA',message='empty network'),
    errors=data.frame(stage=character(),message=character()),references=.pm_network_references(),
    diagnostics=list(limitations=.pm_network_limitations()))
  x$results <- list(ega=list(n.dim=2L,wc=c(R1=1,R2=1,R3=2)),
    retained_items=c('R1','R2','R3'),loadings=list(primary=data.frame(item=c('R1','R2','R3'),
      dimension=c(1L,1L,2L),loading=c(.6,.4,.5),magnitude=c('large','moderate','large'))),
    stages=data.frame(stage=c('UVA','EGA','riEGA','bootstrap','invariance'),
      status=c('skipped','completed','not_evaluable','completed','completed'),
      elapsed_seconds=0,reason=c('UVA was not requested (uva=NULL); redundancy was not assessed.','','Empty residual partition.','','')),
    bootstrap=list(counts=data.frame(stage='bootEGA',requested=100L,valid=100L),
      structural=list(frequency=data.frame(dimensions=2L,frequency=.95),
        dimensions=data.frame(dimension=1:2,structural_consistency=c(.92,.6),average_item_stability=c(.97,.65)),
        items=data.frame(item=c('R1','R2','R3'),stability=c(.99,.95,.65))),
      cs=data.frame(metric=c('strength','betweenness'),CS=c(.596,.051),interpretation=c('preferable','unstable'))),
    invariance=list(groups=data.frame(group=c('F','M'),n=c(171L,143L)),p_resolution=1/500,
      configural_supported=TRUE,
      global=data.frame(comparison='F vs M',test=c('network_structure','global_strength'),
        statistic=c(.256,.413),p_raw=c(.259,.007),p_BH=c(.259,.014)),
      loadings=data.frame(comparison='F vs M',item=c('R1','R2','R3'),dimension=c(1L,1L,2L),
        difference=c(.1,.03,-.1),p_raw=c(.04,.6,.7),p_BH=c(.12,.8,.9),
        evaluable=TRUE,valid_permutations=499L)))
  class(x) <- c('psychomatic_network','psychomatic_result');x
}

test_that('reader exposes signed cross-loadings and every item stability in both languages', {
  x <- manuscript_network_fixture()
  x$results$loadings$matrix <- matrix(c(.6,.4,-.07,-.1,.12,.5),3,
    dimnames=list(c('R1','R2','R3'),c('1','2')))
  before <- x
  for(lang in c('eng','esp')) {
    s <- .pm_network_sections(x,lang)$sections
    name <- if(lang=='eng') 'Revised cross-loadings' else 'Cargas cruzadas revisadas'
    expect_true(name %in% names(s))
    if(!name %in% names(s)) next
    t <- s[[name]]
    expect_equal(t$D1,c(NA_real_,NA_real_,-.07))
    expect_equal(t$D2,c(-.1,.12,NA_real_))
    stability <- s[[if(lang=='eng') 'Item stability' else 'Estabilidad de \u00edtems']]
    expect_equal(stability[[1]],c('R1','R2','R3'))
    expect_equal(stability[[ncol(stability)]],c(.99,.95,.65))
    expect_lt(match(name,names(s)),match(if(lang=='eng')'Loading interpretation' else 'Interpretaci\u00f3n de cargas',names(s)))
  }
  expect_identical(x,before)
})

test_that('bootstrap narrative distinguishes dimensionality count from exact community recovery', {
  x <- manuscript_network_fixture()
  sections <- .pm_network_sections(x)$sections
  expect_match(sections[['Bootstrap interpretation']]$Text,'number of communities')
  expect_equal(sections[['Structural consistency']]$structural_consistency,c(.92,.6))
})

test_that('default network reader is compact and retains scientific qualifications', {
  x <- manuscript_network_fixture();before <- x
  for(lang in c('eng','esp')) {
    s <- .pm_network_sections(x,lang)$sections
    expect_false(any(c('Stage explanations','Warnings','Stages and conditions',
      'Advertencias','Etapas y condiciones','Edge differences','Original configural structure') %in% names(s)))
    expect_true(any(c('Detected communities','Comunidades detectadas') %in% names(s)))
    text <- paste(unlist(s),collapse=' ')
    expect_match(text,'R1, R2')
    expect_match(text,'R3')
    expect_match(text,if(lang=='eng')'global strength.*difference'else'diferencia.*fuerza global')
    expect_match(text,if(lang=='eng')'No item.*BH'else'Ning.n .tem.*BH')
    expect_match(text,if(lang=='eng')'does not establish equivalence'else'no demuestra equivalencia')
    expect_match(text,if(lang=='eng')'riEGA.*not evaluable'else'riEGA.*no evaluable')
  }
  expect_identical(x,before)
  audit <- .pm_network_sections(x,detail='audit')$sections
  expect_true(all(c('Stage explanations','Warnings','Stages and conditions') %in% names(audit)))
  expect_false(any(grepl('WARNINGS|STAGE EXPLANATIONS',capture.output(print(x)))))
  expect_true(any(grepl('WARNINGS',capture.output(print(x,detail='audit')))))
})

test_that('invariance interpretation uses BH and never declares missing tests invariant', {
  x <- manuscript_network_fixture()
  x$results$invariance$loadings$p_BH <- c(.03,.05,NA_real_)
  x$results$invariance$loadings$evaluable[3] <- FALSE
  s <- .pm_network_sections(x)$sections
  tab <- s[['Items with loading differences (BH)']]
  expect_identical(tab$item,'R1')
  expect_match(paste(unlist(s),collapse=' '),'not evaluable|incomplete')
  x$results$invariance$configural_supported <- FALSE
  expect_match(paste(unlist(.pm_network_sections(x)$sections),collapse=' '),'exploratory.*configural|configural.*not supported')
  x$results$invariance <- NULL
  x$results$stages$status[5] <- 'interrupted'
  expect_match(paste(unlist(.pm_network_sections(x)$sections),collapse=' '),'invariance.*interrupted')
})

test_that('APA formatting preserves references and italicizes only proper spans', {
  article <- 'Smith, A. (2025). A title. Journal of Tests, 57(4), Article 114. https://doi.org/10.1000/test'
  book <- 'Brown, T. A. (2015). Confirmatory factor analysis for applied research (2nd ed.). Guilford Press.'
  for(ref in c(article,book,'Unstructured reference')) {
    runs <- .pm_reference_runs(ref)
    expect_identical(paste(runs$text,collapse=''),ref)
    expect_false(grepl('\033',.pm_reference_console(ref,ansi=FALSE)))
  }
  a <- .pm_reference_runs(article)
  expect_identical(a$text[a$italic],c('Journal of Tests','57'))
  b <- .pm_reference_runs(book)
  expect_identical(b$text[b$italic],'Confirmatory factor analysis for applied research')
  expect_match(.pm_reference_console(article,ansi=TRUE),'\033\\[3mJournal of Tests\033\\[23m')
})

test_that('shared APA formatter covers the bibliographies of all analysis families', {
  cfa <- .build_references(list(estimator='MLR',mardia=list(p_skew=NA_real_),
    reliability=list(type='second_order'),respecifications=NULL))
  refs <- unique(c(.efa_references(),cfa,.pm_network_references(),alignment_apa_references(),
    .pm_irt_references()$reference,.pm_dif_references()$reference,
    .pm_reference_text(.ifa_reference_table()),.cormat_references('poly')))
  for(ref in refs) {
    runs <- .pm_reference_runs(ref)
    expect_identical(paste(runs$text,collapse=''),ref)
    if(grepl(', [0-9]+(\\(|,)',ref)) expect_true(any(runs$italic),info=ref)
  }
  software <- .efa_references()[5]
  expect_false(any(.pm_reference_runs(software)$italic))
})

test_that('Word and Excel preserve partial italics and full plain reference text', {
  skip_if_not_installed('officer');skip_if_not_installed('openxlsx')
  refs <- c('Smith, A. (2025). A title. Journal of Tests, 57(4), Article 114. https://doi.org/10.1000/test',
    'Brown, T. A. (2015). Confirmatory factor analysis for applied research (2nd ed.). Guilford Press.',
    'A & B (2025). A <title>. Quality & Quantity, 44, 1-10. https://doi.org/10.1000/escape')
  report <- list(title='APA reference test',sections=list(References=data.frame(Text=refs)),figures=list())
  dir <- tempfile('apa-report-');dir.create(dir)
  on.exit(unlink(dir,recursive=TRUE),add=TRUE)
  word <- file.path(dir,'test.docx');excel <- file.path(dir,'test.xlsx')
  .pm_write_analysis_word(list(report),word)
  .pm_write_analysis_excel(list(report),excel)
  expect_true(all(refs %in% officer::docx_summary(officer::read_docx(word))$text))
  expect_true(all(refs %in% unlist(openxlsx::read.xlsx(excel,sheet=2,colNames=FALSE))))
  utils::unzip(word,files='word/document.xml',exdir=file.path(dir,'word'))
  utils::unzip(excel,files='xl/sharedStrings.xml',exdir=file.path(dir,'excel'))
  w <- paste(readLines(file.path(dir,'word','word/document.xml'),warn=FALSE),collapse='')
  e <- paste(readLines(file.path(dir,'excel','xl/sharedStrings.xml'),warn=FALSE),collapse='')
  expect_match(w,'<w:i[^>]*>.*Journal of Tests')
  expect_match(e,'<i/>.*Journal of Tests')
  expect_match(e,'<i/>.*Confirmatory factor analysis for applied research')
  expect_match(e,'Quality &amp; Quantity')
  expect_false(grepl('\033',paste(w,e)))
  audit_word <- file.path(dir,'audit.docx');audit_excel <- file.path(dir,'audit.xlsx')
  tables <- list(references=data.frame(reference=refs))
  .pm_write_word_report(tables,audit_word)
  .pm_write_excel_report(tables,audit_excel)
  expect_true(all(refs %in% officer::docx_summary(officer::read_docx(audit_word))$text))
  expect_true(all(refs %in% openxlsx::read.xlsx(audit_excel)$reference))
  utils::unzip(audit_word,files='word/document.xml',exdir=file.path(dir,'audit-word'))
  utils::unzip(audit_excel,files='xl/sharedStrings.xml',exdir=file.path(dir,'audit-excel'))
  aw <- paste(readLines(file.path(dir,'audit-word','word/document.xml'),warn=FALSE),collapse='')
  ae <- paste(readLines(file.path(dir,'audit-excel','xl/sharedStrings.xml'),warn=FALSE),collapse='')
  expect_match(aw,'<w:i[^>]*>.*Journal of Tests')
  expect_match(ae,'<i/>.*Journal of Tests')
})

test_that('console references can wrap italic spans without formatting log files', {
  ref <- .pm_network_references()[1]
  formatted <- .pm_reference_console(ref,ansi=TRUE,width=60L)
  plain <- .pm_reference_console(ref,ansi=FALSE,width=60L)
  expect_identical(gsub('\033\\[[0-9]+m','',formatted),plain)
  expect_true(all(nchar(strsplit(plain,'\n',fixed=TRUE)[[1]]) <= 60L))
  old <- options(PsychoMatic.ansi=FALSE);on.exit(options(old),add=TRUE)
  expect_false(any(grepl('\033',capture.output(print(manuscript_network_fixture())))))
  x <- structure(list(report=ref,references=ref),class='inv_align_auto_result')
  options(PsychoMatic.ansi=TRUE)
  expect_true(any(grepl('\033',capture.output(print(x)))))
  expect_identical(x$report,ref)
})
