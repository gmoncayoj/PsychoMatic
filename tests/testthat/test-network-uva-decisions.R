.pm_uva_decision_fixture <- function() {
  m <- matrix(c(0,.4,.08,.06, .4,0,.18,.1, .08,.18,0,.03,
    .06,.1,.03,0),4,dimnames=list(letters[1:4],letters[1:4]))
  list(wto=list(matrix=m,pairwise=data.frame(node_i='a',node_j='b',wto=.4)),
    redundant=list(a='b'),keep_remove=list(keep='a',remove='b'))
}

test_that('UVA modes are explicit and never prompt', {
  expect_null(formals(network_auto)$uva)
  local_mocked_bindings(readline=function(...)stop('Unexpected interactive prompt'),.package='base')
  u <- .pm_uva_decision_fixture()
  for (mode in c('diagnostic','removed')) {
    output <- capture.output(s <- .pm_network_uva_select(u,letters[1:4],
      mode,.25))
    expect_length(output,0)
    expect_identical(s$removed,if(mode=='removed')'b'else character())
    expect_identical(s$effective_mode,mode)
    expect_equal(s$audit$removed_score,.18)
    expect_equal(s$audit$retained_score,.08)
    expect_match(s$audit$reason,'lowest maximum wTO')
  }
})

test_that('no UVA candidates preserves the requested mode and every item', {
  u <- .pm_uva_decision_fixture()
  u$wto$pairwise$wto <- .2
  u$redundant <- list();u$keep_remove <- list(keep=character(),remove=character())
  s <- .pm_network_uva_select(u,letters[1:4],'removed',.25)
  expect_identical(s$decision,'no_candidates')
  expect_identical(s$effective_mode,'removed')
  expect_identical(s$retained,letters[1:4])
})

test_that('UVA audits larger sets and ties using the native representative', {
  u <- .pm_uva_decision_fixture()
  u$redundant <- list(a=c('b','c'));u$keep_remove <- list(keep='b',remove=c('a','c'))
  s <- .pm_network_uva_select(u,letters[1:4],'removed',.25)
  expect_identical(s$audit$representative,c('b','b'))
  expect_equal(s$audit$retained_score,rep((.4+.18)/3,2))
  expect_true(all(grepl('highest mean within-set wTO',s$audit$reason)))
  u <- .pm_uva_decision_fixture();u$wto$matrix['b',c('c','d')] <- c(.08,.06)
  u$wto$matrix[c('c','d'),'b'] <- c(.08,.06)
  s <- .pm_network_uva_select(u,letters[1:4],'removed',.25)
  expect_true(s$audit$tied)
  expect_match(s$audit$reason,'tie.*native')
})

test_that('network weights and item-dimension loadings occupy different edges', {
  skip_if_not_installed('qgraph')
  graph <- matrix(c(0,.42,.42,0),2,dimnames=list(c('b','a'),c('b','a')))
  loads <- list(matrix=matrix(c(.31,.2,.1,-.5),2,
    dimnames=list(c('a','b'),c('1','2'))),
    primary=data.frame(item=c('a','b'),dimension=c('1','2'),loading=c(.31,-.5)))
  x <- list(language='eng',results=list(ega=list(network=graph,wc=c(b=2,a=1)),
    loadings=loads))
  class(x) <- 'psychomatic_network'
  path <- tempfile(fileext='.pdf');grDevices::pdf(path, encoding='WinAnsi.enc')
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  n <- plot(x,type='network')
  expect_equal(n$graphAttributes$Nodes$labels,c('b','a'))
  expect_equal(n$Edgelist$weight,.42)
  expect_equal(as.numeric(n$graphAttributes$Edges$labels),.42)
  l <- plot(x,type='loadings')
  expect_equal(sort(l$Edgelist$weight),c(-.5,.31))
  expect_equal(sort(as.numeric(l$graphAttributes$Edges$labels)),c(-.5,.31))
  expect_true('loadings_matrix' %in% .pm_network_plot_types(x))
  expect_match(.pm_network_figure_caption('network',FALSE),'partial correlations')
  expect_match(.pm_network_figure_caption('loadings',FALSE),'Item.*dimension')
})

test_that('NULL truly skips UVA and does not invent redundancy results', {
  skip_if_not_installed('EGAnet')
  local_mocked_bindings(UVA=function(...)stop('UVA must not run'),.package='EGAnet')
  d <- simulate_likert_fixture(n=180,seed=391)
  x <- network_auto(d,wording=FALSE,hierarchy=FALSE,bootstrap=FALSE,
    invariance=FALSE,consensus_iter=5L,seed=391,verbose=FALSE)
  expect_null(x$results$uva)
  expect_null(x$requested_settings$uva)
  expect_null(x$effective_settings$uva)
  expect_identical(x$results$retained_items,names(d))
  expect_identical(x$results$removed_items,character())
  expect_identical(x$conclusion$status,'completed')
  expect_identical(x$results$stages$status[x$results$stages$stage=='UVA'],'skipped')
  expect_false('uva' %in% network_plots(x))
  for(lang in c('eng','esp')) {
    sections <- .pm_network_sections(x,lang,'report')
    expect_false(any(c('Redundant pairs','Pares redundantes','UVA selection','Selecci\u00f3n UVA') %in% names(sections$sections)))
    text <- sections$sections[[if(lang=='eng')'Method'else'M\u00e9todo']]$Text
    expect_match(text,if(lang=='eng')'not performed'else'no se ejecut\u00f3')
    expect_false('uva' %in% vapply(sections$figures,`[[`,character(1),'type'))
  }
  for(mode in list('ask','reduce','invalid',NA_character_,c('diagnostic','removed'),TRUE))
    expect_error(network_auto(d,uva=mode,bootstrap=FALSE,invariance=FALSE),'uva')
})

test_that('explicit native UVA decisions reach final bilingual results and reports', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=220,seed=137)
  d$duplicate <- d[[1]]
  local_mocked_bindings(readline=function(...)stop('Unexpected prompt'),.package='base')
  run <- function(mode)network_auto(d,uva=mode,wording=FALSE,hierarchy=FALSE,bootstrap=FALSE,
    invariance=FALSE,consensus_iter=5L,seed=137,verbose=FALSE)
  capture.output(kept <- run('diagnostic'))
  expect_identical(kept$requested_settings$uva,'diagnostic')
  expect_identical(kept$effective_settings$uva,'diagnostic')
  expect_identical(kept$results$removed_items,character())
  expect_identical(kept$results$retained_items,names(d))
  expect_identical(kept$results$uva_selection$decision,'explicit_retain')
  capture.output(reduced <- run('removed'))
  expect_identical(reduced$effective_settings$uva,'removed')
  expect_setequal(reduced$results$removed_items,reduced$results$uva$keep_remove$remove)
  expect_identical(rownames(reduced$results$ega$network),reduced$results$retained_items)
  expect_true(all(reduced$results$uva_selection$audit$removed))
  for(es in c(FALSE,TRUE)) {
    diagnostic_sections <- .pm_network_sections(kept,if(es)'esp'else'eng')$sections
    proposed <- diagnostic_sections[[if(es)'Propuesta UVA y motivos de selecci\u00f3n'else'UVA proposal and selection reasons']]
    expect_true(all(grepl(if(es)'Se propone'else'Propose',proposed[[if(es)'Motivo'else'reason']])))
    sections <- .pm_network_sections(reduced,if(es)'esp'else'eng')$sections
    audit <- sections[[if(es)'Propuesta UVA y motivos de selecci\u00f3n'else'UVA proposal and selection reasons']]
    expect_true(all(grepl(if(es)'Se elimin\u00f3'else'Removed',audit[[if(es)'Motivo'else'reason']])))
    expect_setequal(sections[[if(es)'\u00cdtems eliminados por UVA'else'Items removed by UVA']][[1]],
      reduced$results$removed_items)
  }
})
