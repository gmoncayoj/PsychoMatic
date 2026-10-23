test_that('network method choices reproduce native EGA', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=160,seed=609)
  for(model in c('glasso','TMFG')) for(algorithm in c('louvain','leiden','walktrap')) {
    x <- network_auto(d,model=model,algorithm=algorithm,wording=FALSE,
      hierarchy=FALSE,bootstrap=FALSE,invariance=FALSE,consensus_iter=5,
      seed=609,verbose=FALSE)
    native <- EGAnet::EGA(as.matrix(d),corr='spearman',na.data='listwise',
      model=model,algorithm=algorithm,gamma=.5,consensus.iter=5,seed=609,
      objective_function='modularity',
      uni.method=if(algorithm=='louvain')'louvain' else 'LE',plot.EGA=FALSE)
    expect_equal(x$results$ega$network,native$network)
    expect_equal(x$results$ega$wc,native$wc)
    expect_identical(x$effective_settings$model,model)
    expect_identical(x$effective_settings$algorithm,algorithm)
    estimator <- .pm_network_estimator_function(x$effective_settings)
    expect_equal(estimator(d),native$network,ignore_attr=TRUE)
    method <- .pm_network_sections(x)$sections$Method$Text
    expect_true(any(grepl(model,method,fixed=TRUE)))
    if(model=='TMFG')expect_false(any(grepl('gamma =',method,fixed=TRUE)))
  }
})

test_that('chosen bootstrap and hierarchy algorithms are propagated', {
  skip_if_not_installed('EGAnet')
  local_mocked_bindings(bootEGA=function(...)list(...),.package='EGAnet')
  settings <- list(corr='spearman',missing='listwise',gamma=.5,model='TMFG',
    algorithm='leiden',bootstrap_type='parametric',consensus_iter=5,
    seed=609,n_boot=5,n_cores=1)
  for(variant in c('EGA','riEGA','hierEGA')) {
    args <- .pm_network_boot_ega(data.frame(a=1:5,b=5:1,c=1:5),settings,variant)
    expect_identical(args$model,'TMFG')
    expect_identical(args$type,'parametric')
    expect_null(args$gamma)
    expect_null(args$consensus.iter)
    if(variant=='hierEGA') {
      expect_null(args$algorithm)
      expect_identical(args$lower.algorithm,'leiden')
      expect_identical(args$higher.algorithm,'leiden')
    } else expect_identical(args$algorithm,'leiden')
  }
})

test_that('parametric bootstrap is opt-in and preserves native replications', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=160,seed=609)
  settings <- list(corr='spearman',missing='listwise',gamma=.5,model='glasso',
    algorithm='louvain',bootstrap_type='parametric',consensus_iter=5,
    seed=609,n_boot=5,n_cores=1)
  result <- .pm_network_seed(609,.pm_network_boot_ega(d,settings))
  native <- .pm_network_seed(609,EGAnet::bootEGA(as.matrix(d),corr='spearman',na.data='listwise',
    model='glasso',algorithm='louvain',gamma=.5,consensus.iter=5,seed=609,
    iter=5,ncores=1,type='parametric',plot.itemStability=FALSE,verbose=FALSE))
  expect_identical(result$type,'parametric')
  expect_equal(result$bootGraphs,native$bootGraphs)
  expect_equal(result$boot.wc,native$boot.wc)
  expect_identical(eval(formals(network_auto)$bootstrap_type),c('resampling','parametric'))
  expect_error(network_auto(d,model='invalid'),'arg')
  expect_error(network_auto(d,algorithm='invalid'),'arg')
  expect_error(network_auto(d,bootstrap_type='auto'),'arg')
})

test_that('TMFG loading permutations reproduce the fixed-partition native test', {
  skip_if_not_installed('EGAnet')
  d <- simulate_likert_fixture(n=160,seed=609)
  settings <- list(corr='spearman',missing='listwise',gamma=.5,model='TMFG',
    algorithm='leiden',bootstrap_type='parametric',consensus_iter=5,
    seed=609,n_boot=5,n_perm=3,n_cores=1)
  ega <- .pm_network_ega(d,settings)
  groups <- rep(c('F','M'),each=80)
  result <- .pm_network_seed(609,.pm_network_metric_invariance(d,groups,ega$wc,settings))
  native <- .pm_network_seed(609,EGAnet::invariance(as.matrix(d),groups,
    structure=ega$wc,corr='spearman',na.data='listwise',model='TMFG',
    algorithm='walktrap',uni.method='LE',iter=3,ncores=1,seed=609,
    loading.method='revised',verbose=FALSE))
  expect_equal(result$groups$loadingsDifference,native$groups$loadingsDifference)
  expect_equal(result$permutation$loadingsDifference,native$permutation$loadingsDifference)
  expect_equal(result$permutation$groups,native$permutation$groups)
})

test_that('disabled BGGM fails explicitly instead of falling back to GLASSO', {
  skip_if_not_installed('EGAnet')
  if(!any(grepl('not available at the moment',deparse(EGAnet::network.estimation),fixed=TRUE)))
    skip('Installed EGAnet does not disable BGGM')
  d <- simulate_likert_fixture(n=100,seed=609)
  expect_error(network_auto(d,model='BGGM',bootstrap=FALSE,invariance=FALSE,
    wording=FALSE,hierarchy=FALSE,verbose=FALSE),'BGGM.*not available')
})

test_that('native variant limitations retain explicit failures without changing methods', {
  skip_if_not_installed('EGAnet')
  settings <- list(corr='spearman',missing='listwise',model='TMFG',algorithm='leiden',
    bootstrap_type='parametric',n_boot=3,n_cores=1,seed=609)
  d <- data.frame(a=1:5,b=5:1,c=1:5,d=5:1)
  local_mocked_bindings(bootEGA=function(...)stop("invalid 'use' argument"),.package='EGAnet')
  expect_error(.pm_network_boot_ega(d,settings,'riEGA'),'invalid native correlation missing-data')
  local_mocked_bindings(hierEGA=function(...)stop("invalid 'nrow' value (< 0)"),.package='EGAnet')
  expect_error(.pm_network_hierarchy(d,NULL,NULL,settings),'at least four nodes')
})
