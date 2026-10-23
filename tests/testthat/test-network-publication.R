test_that('publication selection omits unavailable diagnostics without deleting them', {
  x <- list(results=list(ega=list(wc=c(a=1,b=2)),loadings=list(matrix=diag(2)),
    wording=list(status='not_evaluable'),hierarchy=NULL,bootstrap=list(),
    invariance=list(loadings=data.frame(evaluable=TRUE,p_BH=.5,valid_permutations=99L))),
    effective_settings=list(n_perm=99L))
  types <- .pm_network_publication_types(x)
  expect_true(all(c('ega_stability','loadings','dimension_stability','group_networks','loading_differences') %in% types))
  expect_false(any(c('wording','tefi','centrality','edge_differences') %in% types))
  expect_true(all(c('network','item_stability','wording','tefi','edge_differences') %in% .pm_network_plot_types(x)))
})

test_that('combined EGA stability figure is exportable and preserves signed estimates', {
  skip_if_not_installed('qgraph')
  graph <- matrix(c(0,-.42,-.42,0),2,dimnames=list(c('a','b'),c('a','b')))
  x <- structure(list(language='eng',results=list(ega=list(network=graph,wc=c(a=1,b=2)),
    bootstrap=list(structural=list(items=data.frame(item=c('b','a'),stability=c(.8,1)),
      dimensions=data.frame(dimension=1:2,structural_consistency=c(1,.8)))),
    loadings=list(primary=data.frame(item=c('a','b'),dimension=1:2,loading=c(.4,.5))))),
    class='psychomatic_network')
  before <- x
  expect_true('ega_stability' %in% network_plots(x))
  if(!'ega_stability' %in% network_plots(x)) return(invisible(NULL))
  directory <- tempfile();dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  paths <- network_plots(x,directory,types=c('ega_stability','item_stability'),formats='pdf')
  expect_length(paths,2L)
  expect_true(all(file.info(paths)$size>1000))
  expect_identical(x,before)
})

test_that('network item flags use complete BH tests for each comparison', {
  t <- data.frame(comparison=rep(c('F vs M','F vs X'),each=4),item=rep(c('a','b','c','d'),2),
    difference=0,p_BH=c(.049,.05,NA,.001,.7,.01,.2,.4),evaluable=TRUE,
    valid_permutations=c(499L,499L,499L,8L,499L,499L,499L,499L))
  x <- list(results=list(invariance=list(loadings=t)),effective_settings=list(n_perm=499L))
  expect_equal(.pm_network_invariance_status(x,'F vs M',c('a','b','c','d','e')),
    c('difference','no_difference','unavailable','unavailable','unavailable'))
  expect_equal(.pm_network_invariance_status(x,'F vs X',c('a','b')),
    c('no_difference','difference'))
})

test_that('replication points preserve near-one, zero and unavailable item values', {
  path <- tempfile(fileext='.pdf');grDevices::pdf(path)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  result <- .pm_network_stability_points(c(.998,NA_real_,0),c('b','a','c'),c(2,1,1),'Stability')
  expect_equal(result$item,c('a','c','b'))
  expect_equal(result$stability,c(NA_real_,0,.998))
  unavailable <- .pm_network_stability_points(c(NA_real_,NA_real_),c('a','b'),c(1,1),'Stability')
  expect_true(all(is.na(unavailable$stability)))
})

test_that('loading comparison plot distinguishes BH flags, incomplete tests and zero', {
  t <- data.frame(comparison='F vs M',item=c('R1','R2','R3'),dimension=c('1','1','2'),
    difference=c(-.1,.2,.3),p_BH=c(.02,.07,NA_real_),evaluable=c(TRUE,TRUE,FALSE),
    valid_permutations=c(499L,499L,0L))
  points <- .pm_network_loading_plot_data(t,499L)
  expect_identical(points$significant,c(TRUE,FALSE,FALSE))
  expect_identical(points$available,c(TRUE,TRUE,FALSE))
  expect_equal(points$difference,t$difference)
  t$valid_permutations[1] <- 10L
  expect_false(.pm_network_loading_plot_data(t,499L)$significant[1])
})

test_that('community colors are stable across reordered item labels', {
  a <- .pm_network_community_colors(c(R1=1,R2=2,R3=3))
  b <- .pm_network_community_colors(c(R3=3,R1=1,R2=2))
  expect_identical(a,b)
  expect_length(unique(unname(a)),3L)
})

test_that('network plots expose a native community legend and preserve edge weights', {
  skip_if_not_installed('qgraph')
  graph <- matrix(c(0,-.42,-.42,0),2,dimnames=list(c('a','b'),c('a','b')))
  x <- structure(list(language='eng',results=list(ega=list(network=graph,wc=c(a=1,b=2)))),class='psychomatic_network')
  path <- tempfile(fileext='.pdf');grDevices::pdf(path)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  result <- plot(x)
  expect_true(result$plotOptions$legend)
  expect_identical(result$plotOptions$legend.mode,'groups')
  expect_equal(result$Edgelist$weight,-.42)
  expect_identical(x$results$ega$network,graph)
})
