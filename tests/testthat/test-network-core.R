test_that("network input protects raw wording and validates keys", {
  d <- data.frame(a=1:20, b=20:1, c=rep(1:5, 4), sex=rep(c("F", "M"),10))
  p <- .pm_network_input(d, c("a","b","c"), "sex", "b", list(b=c(1,20)), "listwise")
  expect_equal(p$raw$b, d$b)
  expect_equal(p$keyed$b, 21-d$b)
  expect_error(.pm_network_input(d, c("a","b","c"), "sex", "b", NULL, "listwise"), "ranges")
  expect_error(network_auto(d, items=c("a","b","c"), n_boot=-1), "n_boot")
})

test_that("network core reproduces native Spearman glasso and revised loadings", {
  skip_if_not_installed("EGAnet")
  d <- simulate_likert_fixture(n=180, seed=391)
  x <- network_auto(d, uva="diagnostic", wording=FALSE, hierarchy=FALSE,
    bootstrap=FALSE, invariance=FALSE, consensus_iter=10, seed=391, verbose=FALSE)
  e <- EGAnet::EGA(as.matrix(d), corr="spearman", model="glasso", algorithm="louvain",
    consensus.iter=10, gamma=.5, na.data="listwise", seed=391, plot.EGA=FALSE)
  expect_s3_class(x, "psychomatic_network")
  expect_s3_class(x, "psychomatic_result")
  expect_equal(x$results$ega$network, e$network)
  expect_equal(x$results$loadings$matrix,
    EGAnet::net.loads(e, loading.method="revised", ordered="variable")$std)
  expect_equal(x$results$centrality$strength, unname(rowSums(abs(e$network))))
  expect_equal(x$results$tefi$value[1], as.numeric(e$TEFI))
  expect_identical(x$effective_settings$language, "eng")
})

test_that("BH is applied once within each declared family", {
  p <- data.frame(family=c("global","global","edges","edges"), p_raw=c(.01,.04,.01,.2))
  out <- .pm_network_bh(p)
  expect_equal(out$p_BH, c(.02,.04,.02,.2))
  expect_equal(.pm_network_loading_magnitude(c(.19,.20,-.35,.50)),
    c("below_small","small","moderate","large"))
})

test_that("NCT edge differences retain signs and precision from the graphs", {
  a<-matrix(c(0,.123456789,.123456789,0),2,dimnames=list(c("x","y"),c("x","y")))
  b<-matrix(c(0,.234567891,.234567891,0),2,dimnames=dimnames(a))
  native<-list(nw1=a,nw2=b,einv.pvals=data.frame(Var1=1,Var2=2,p=.4321,stat=.111))
  t<-.pm_network_nct_edges(native,"F vs M",c("x","y"))
  expect_equal(t$difference,a[1,2]-b[1,2])
  expect_equal(t$item1,"x");expect_equal(t$item2,"y")
  expect_equal(t$p_raw,.4321)
})

test_that("centrality matches native unnormalized shortest-path conventions", {
  skip_if_not_installed("qgraph")
  a<-matrix(c(0,.5,0,.5,0,.5,0,.5,0),3,dimnames=list(letters[1:3],letters[1:3]))
  native<-qgraph::centrality_auto(a)$node.centrality
  actual<-.pm_network_centrality(a)
  expect_equal(actual$closeness,unname(native$Closeness))
  expect_equal(actual$betweenness,unname(native$Betweenness))
})

test_that("Spanish presentation preserves arbitrary item and group identifiers", {
  t<-data.frame(item=c("small","strength"),group=c("edges","centrality"),
    dimension=c("large","small"),magnitude=c("small","large"))
  translated<-.pm_network_translate_table(t,TRUE)
  expect_identical(translated[["Ítem"]],t$item)
  expect_identical(translated[["Grupo"]],t$group)
  expect_identical(translated[["Dimensión"]],t$dimension)
  expect_identical(translated[["Magnitud"]],c("pequeña","grande"))
})

test_that("mixed hierarchical candidate encoding keeps every ordinal indicator", {
  d<-data.frame(a=seq(.01,.8,length.out=80),b=rep(1:4,20),c=rep(1:5,16))
  s<-psychomatic_spec(d,items=names(d))
  e<-.pm_network_cfa_encoding(names(d),s)
  expect_equal(e$ordered,c("x2","x3"))
  expect_equal(e$estimator,"WLSMV")
})

test_that("undefined observed loading differences cannot become significant", {
  perm<-matrix(rep(c(.1,.2,.3),4),4,3,byrow=TRUE)
  actual<-c(NA_real_,NaN,Inf,.25)
  p<-.pm_network_perm_p(actual,perm,3)
  expect_true(all(is.na(p[1:3])))
  expect_equal(p[4],.5)
  perm[4,1]<-NA_real_
  expect_true(is.na(.pm_network_perm_p(actual,perm,3)[4]))
})

test_that("network permutation p values never display as zero", {
  t<-data.frame(p_raw=1/5001,p_BH=1/5001)
  displayed<-.pm_display_table(t,3)
  expect_equal(displayed$p_raw,"<0.001")
  expect_equal(displayed$p_BH,"<0.001")
})

test_that("RI diagnostics require a converged admissible native fit", {
  skip_if_not_installed("EGAnet")
  d<-simulate_likert_fixture(n=180,seed=991)
  fit<-suppressWarnings(lavaan::cfa(paste("RI =~",paste(paste0("1*",names(d)),collapse=" + ")),
    data=d,control=list(iter.max=1)))
  expect_false(isTRUE(lavaan::lavInspect(fit,"converged")))
  expect_error(.pm_network_validate_ri(list(RI=list(fit=fit))),"converge")
  expect_error(.pm_network_validate_ri(list(RI=list(fit=structure("failed",class="try-error")))),"converge")
})
