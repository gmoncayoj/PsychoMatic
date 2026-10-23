.pm_network_fixture <- function(n=300L, seed=493L) {
  .pm_fixture_with_seed(seed, {
    z <- matrix(rnorm(n*3L),n,3L)
    x <- as.data.frame(sapply(1:12,function(j).75*z[,ceiling(j/4)]+rnorm(n,sd=.66)))
    names(x)<-paste0("q",1:12);x$sex<-rep(c("F","M"),length.out=n);x
  })
}

test_that("requested network stages run against native packages", {
  .pm_skip_extended()
  skip_if_not_installed("EGAnet");skip_if_not_installed("bootnet");skip_if_not_installed("NetworkComparisonTest")
  d<-.pm_network_fixture()
  x<-network_auto(d,items=1:12,group="sex",uva="diagnostic",n_boot=20L,n_perm=9L,consensus_iter=5L,seed=493L,verbose=FALSE)
  expect_equal(x$results$stages$status,rep("completed",nrow(x$results$stages)),info=paste(x$results$stages$reason,collapse="; "))
  if(any(x$results$stages$status=="failed"))return(invisible(NULL))
  expect_true(all(x$warnings$stage %in% x$results$stages$stage))
  expect_equal(nrow(x$errors),length(x$errors$stage))
  expect_equal(x$results$bootstrap$estimated$graph,unclass(x$results$ega$network),ignore_attr=TRUE)
  expect_equal(x$results$bootstrap$structural$valid,20)
  i<-x$results$invariance
  expect_equal(i$p_resolution,.1)
  expect_equal(i$original_item_set,paste0("q",1:12))
  expect_equal(i$global$p_BH,stats::p.adjust(i$global$p_raw,"BH"))
  expect_equal(i$loadings$p_BH,stats::p.adjust(i$loadings$p_raw,"BH"))
  native<-i$native_metric
  actual<-native$groups$loadingsDifference[[1]]
  perm<-native$permutation$loadingsDifference[[1]]
  expect_equal(i$loadings$p_raw,unname((1+rowSums(abs(perm)>=abs(actual)))/10))
  expect_equal(i$loadings$difference,unname(actual))
  expect_true(all(i$loadings$valid_permutations==9))
  expect_equal(x$results$wording$original$correlation,
    EGAnet::EGA(as.matrix(d[1:12]),corr="spearman",na.data="listwise",algorithm="louvain",
      consensus.iter=5,seed=493,plot.EGA=FALSE)$correlation)
  skip_if_not_installed("callr");skip_if_not_installed("pkgload")
  directory<-tempfile("cached-network-");dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  path<-file.path(directory,"network.rds");saveRDS(x,path)
  dev<-pkgload::is_dev_package("PsychoMatic")
  root<-if(dev)pkgload::pkg_path()else NULL
  callr::r(function(path,root,libs) {
    .libPaths(libs)
    if(!is.null(root))pkgload::load_all(root,quiet=TRUE)else loadNamespace("PsychoMatic")
    cached<-readRDS(path)
    PsychoMatic::network_plots(cached,dirname(path),types=c("edge_accuracy","centrality_stability"),formats="png")
  },args=list(path=path,root=root,libs=.libPaths()))
  expect_true(all(file.exists(file.path(directory,c("edge_accuracy.png","centrality_stability.png")))))
})

test_that("UVA, disconnected centrality and specifications have explicit behavior", {
  skip_if_not_installed("EGAnet")
  d<-.pm_network_fixture();d$duplicate<-d$q1
  x<-network_auto(d,items=c(1:12,14),group="sex",uva="removed",wording=FALSE,hierarchy=FALSE,bootstrap=FALSE,invariance=FALSE,consensus_iter=5,verbose=FALSE)
  expect_true(length(x$results$removed_items)>=1)
  expect_false(all(c("q1","duplicate") %in% x$results$retained_items))
  a<-matrix(0,3,3,dimnames=list(letters[1:3],letters[1:3]));a[1,2]<-a[2,1]<-.3
  c<-.pm_network_centrality(a)
  expect_true(all(is.na(c$closeness)))
  expect_false(any(c$connected))
  s<-psychomatic_spec(d,items=1:12,group="sex")
  z<-network_auto(s,wording=FALSE,hierarchy=FALSE,bootstrap=FALSE,invariance=FALSE,uva="diagnostic",consensus_iter=5,verbose=FALSE)
  expect_equal(z$data_summary$group,"sex")
  expect_error(network_auto(s,items=1:3,bootstrap=FALSE),"psychomatic_spec")
})

test_that("session-derived seeds control native bootstrap and permutation RNG", {
  .pm_skip_extended()
  skip_if_not_installed("EGAnet");skip_if_not_installed("bootnet")
  skip_if_not_installed("NetworkComparisonTest")
  d <- .pm_network_fixture()
  run <- function() network_auto(d,items=1:12,group="sex",uva="diagnostic",
    wording=FALSE,hierarchy=FALSE,n_boot=5L,n_perm=3L,consensus_iter=5L,verbose=FALSE)
  set.seed(79);a <- run();after_a <- .Random.seed
  set.seed(79);b <- run();after_b <- .Random.seed
  expect_true(all(a$results$stages$status=="completed"),
    info=paste(a$results$stages$reason,collapse="; "))
  expect_true(all(b$results$stages$status=="completed"),
    info=paste(b$results$stages$reason,collapse="; "))
  expect_identical(a$effective_settings$seed,b$effective_settings$seed)
  expect_identical(after_a,after_b)
  expect_equal(a$results$bootstrap$structural$native$boot.wc,
    b$results$bootstrap$structural$native$boot.wc)
  expect_equal(a$results$invariance$native_metric$permutation$loadingsDifference,
    b$results$invariance$native_metric$permutation$loadingsDifference)
  expect_equal(a$results$invariance$global,b$results$invariance$global)
})

test_that("network exports preserve values, figures and caller RNG", {
  .pm_skip_extended()
  skip_if_not_installed("EGAnet");skip_if_not_installed("officer");skip_if_not_installed("openxlsx")
  d<-.pm_network_fixture()
  set.seed(18);before<-.Random.seed
  x<-network_auto(d,items=1:12,group="sex",uva="diagnostic",wording=FALSE,hierarchy=FALSE,bootstrap=FALSE,invariance=FALSE,consensus_iter=5,seed=18L,verbose=FALSE)
  expect_identical(.Random.seed,before)
  directory<-tempfile("net-report-");dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  paths<-network_plots(x,directory,types="network",formats=c("pdf","svg","png"))
  expect_true(all(file.exists(paths)))
  expect_error(network_plots(x,directory,types="network"),"exist")
  eng<-.pm_network_sections(x,"eng","report");esp<-.pm_network_sections(x,"esp","report")
  expect_equal(eng$sections[["Revised primary loadings"]]$loading,
    esp$sections[["Cargas primarias revisadas"]][["Carga"]])
  paths<-psychomatic_report(x,format=c("word","excel"),file_name=file.path(directory,"network"),language="esp")
  word<-file.path(directory,"network.docx");excel<-file.path(directory,"network.xlsx")
  expect_true(file.exists(word));expect_true(file.exists(excel))
  expect_true(any(grepl("word/media/",utils::unzip(word,list=TRUE)$Name)))
  expect_true(any(grepl("xl/media/",utils::unzip(excel,list=TRUE)$Name)))
  utils::unzip(excel,exdir=file.path(directory,"xlsx-parts"))
  rels<-readLines(file.path(directory,"xlsx-parts","xl","worksheets","_rels","sheet2.xml.rels"),warn=FALSE)
  expect_false(any(grepl('/vmlDrawing"',rels)),info="No phantom VML relationship in a report with images")
  expect_true(any(grepl('/drawing"',rels)),info="Actual figure relationships must survive cleanup")
  expect_false(any(grepl("raw_data|row_index",officer::docx_summary(officer::read_docx(word))$text)))
})
