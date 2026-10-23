test_that('polychoric uses native ordinal estimation and rejects continuous items', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=150,seed=893)
  x <- network_auto(d,corr='polychoric',wording=FALSE,hierarchy=FALSE,
    bootstrap=FALSE,invariance=FALSE,consensus_iter=5,seed=893,verbose=FALSE)
  native <- EGAnet::EGA(as.matrix(d),corr='auto',ordinal.categories=7,
    na.data='listwise',model='glasso',algorithm='louvain',gamma=.5,
    consensus.iter=5,seed=893,plot.EGA=FALSE)
  expect_equal(x$results$ega$network,native$network)
  expect_equal(x$results$ega$correlation,native$correlation)
  expect_identical(x$effective_settings$corr,'polychoric')
  d[[1]] <- seq_len(nrow(d))/nrow(d)
  expect_error(network_auto(d,items=seq_along(d),corr='polychoric',bootstrap=FALSE),'ordinal')
})

test_that('fixed-partition loading permutations preserve the native test', {
  skip_if_not_installed('EGAnet')
  d <- .pm_fixture_with_seed(814,{
    z <- matrix(rnorm(300*3),300,3)
    as.data.frame(sapply(1:12,function(j)as.numeric(cut(.8*z[,ceiling(j/4)]+rnorm(300,sd=.5),
      breaks=c(-Inf,-1,-.3,.3,1,Inf)))))
  })
  groups <- rep(c('F','M'),each=150)
  for(corr in c('spearman','pearson','polychoric')) {
    settings <- list(corr=corr,gamma=.5,missing='listwise',consensus_iter=5,
      n_boot=3,n_perm=3,n_cores=1,seed=814,ordinal_categories=7,verbose=FALSE)
    ega <- .pm_network_ega(d,settings)
    a <- .pm_network_metric_invariance(d,groups,ega$wc,settings)
    b <- EGAnet::invariance(as.matrix(d),groups,structure=ega$wc,
      corr=if(corr=='polychoric')'auto' else corr,ordinal.categories=7,
      na.data='listwise',model='glasso',algorithm='louvain',gamma=.5,
      consensus.iter=5,iter=3,ncores=1,seed=814,verbose=FALSE,
      loading.method='revised')
    expect_equal(a$groups$loadingsDifference,b$groups$loadingsDifference)
    expect_equal(a$permutation$groups,b$permutation$groups)
    expect_equal(a$permutation$loadingsDifference,b$permutation$loadingsDifference)
  }
})

test_that('hierarchical bootstrap respects consensus and correlation arguments', {
  skip_if_not_installed('EGAnet')
  local_mocked_bindings(bootEGA=function(...)list(...),.package='EGAnet')
  settings <- list(corr='polychoric',ordinal_categories=7,missing='listwise',gamma=.5,
    consensus_iter=12,seed=814,n_boot=3,n_cores=2)
  args <- .pm_network_boot_ega(data.frame(a=1:5,b=5:1,c=1:5),settings,'hierEGA')
  expect_null(args$algorithm)
  expect_identical(args$consensus.iter,12)
  expect_identical(args$lower.algorithm,'louvain')
  expect_identical(args$corr,'auto')
  expect_identical(args$ordinal.categories,7)
})

test_that('interrupts preserve earlier stages and stop subsequent work', {
  skip_if_not_installed('EGAnet')
  local_mocked_bindings(.pm_network_wording=function(...)stop('Interrupted\nSource: community/louvain.c:669'),
    .pm_network_hierarchy=function(...)stop('Must not run'),.package='PsychoMatic')
  d <- simulate_likert_fixture(n=100,seed=83)
  x <- network_auto(d,bootstrap=FALSE,invariance=FALSE,consensus_iter=5,seed=83,verbose=FALSE)
  expect_equal(x$results$stages$status[x$results$stages$stage=='riEGA'],'interrupted')
  expect_equal(x$conclusion$status,'partial')
  expect_null(x$results$hierarchy)
  expect_true(all(is.finite(x$results$stages$elapsed_seconds)))
  expect_false(any(grepl('Must not run',x$errors$message)))
})

test_that('R interrupt conditions are also recorded in partial results', {
  skip_if_not_installed('EGAnet')
  local_mocked_bindings(.pm_network_wording=function(...)stop(structure(
    list(message='User interrupt',call=NULL),class=c('interrupt','condition'))),.package='PsychoMatic')
  d <- simulate_likert_fixture(n=100,seed=83)
  x <- network_auto(d,bootstrap=FALSE,invariance=FALSE,consensus_iter=5,seed=83,verbose=FALSE)
  expect_equal(x$results$stages$status[x$results$stages$stage=='riEGA'],'interrupted')
  expect_null(x$results$hierarchy)
  expect_match(x$errors$message,'User interrupt')
})

test_that('ordinal limits and riEGA threshold are explicit rather than silent method changes', {
  expect_error(.pm_network_corr_args('polychoric',data.frame(a=1:12)),'11 ordinal categories')
  settings <- list(corr='polychoric')
  w <- .pm_network_wording(data.frame(a=rep(1:8,10)),NULL,NULL,settings)
  expect_identical(w$status,'not_evaluable')
  expect_match(w$reason,'seven-category')
})

test_that('bootnet estimator survives workers with another PsychoMatic version', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=150,seed=893)
  settings <- list(corr='polychoric',ordinal_categories=7,missing='listwise',gamma=.5)
  fun <- .pm_network_estimator_function(settings)
  expect_identical(parent.env(environment(fun)),baseenv())
  expected <- .pm_network_estimator(d,corr='polychoric')
  cluster <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cluster),add=TRUE)
  result <- parallel::clusterCall(cluster,function(fun,d,libs) {
    .libPaths(libs)
    fun(d)
  },fun,d,.libPaths())
  expect_equal(result[[1]],expected)
  expect_equal(result[[2]],expected)
})

test_that('progress and hierarchy separate network evidence from supplementary CFA', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=150,seed=82)
  text <- capture.output(x <- network_auto(d,wording=FALSE,bootstrap=FALSE,invariance=FALSE,
    consensus_iter=5,seed=82,verbose=TRUE))
  expect_true(any(grepl('Starting EGA',text,fixed=TRUE)))
  expect_true(any(grepl('Finished EGA',text,fixed=TRUE)))
  h <- x$results$hierarchy
  expect_true(all(c('order','dimensions','TEFI') %in% names(h$network_structure)))
  sections <- .pm_network_sections(x)$sections
  expect_true('Network hierarchy (hierEGA)' %in% names(sections))
  expect_true('Hierarchy interpretation' %in% names(sections))
  audit <- .pm_network_sections(x,detail='audit')$sections
  expect_true('Hierarchy method' %in% names(audit))
  if(nrow(h$candidates))expect_true('Supplementary CFA candidates' %in% names(audit))
})
