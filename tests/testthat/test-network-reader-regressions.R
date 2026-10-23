test_that("network seeds are optional and respect the caller's RNG", {
  expect_null(formals(network_auto)$seed)
  set.seed(79)
  before <- .Random.seed
  generated <- .pm_network_seed(NULL, stats::runif(3))
  expect_false(identical(before, .Random.seed))
  set.seed(79)
  expect_equal(generated, stats::runif(3))
  before <- .Random.seed
  expect_equal(.pm_network_seed(2401L, stats::runif(3)),
    .pm_network_seed(2401, stats::runif(3)))
  expect_identical(.Random.seed, before)
})

test_that("reader output restricts redundant pairs without deleting diagnostic data", {
  skip_if_not_installed("EGAnet")
  d <- simulate_likert_fixture(n=180,seed=392)
  x <- network_auto(d, uva="diagnostic", wording=FALSE, hierarchy=FALSE,
    bootstrap=FALSE, invariance=FALSE, consensus_iter=5, seed=NULL, verbose=FALSE)
  expect_type(x$effective_settings$seed,"double")
  expect_identical(x$effective_settings$seed_source,"session_generated")
  expect_null(x$requested_settings$seed)
  pairs <- data.frame(node_i=c("a","a","b","c"),node_j=c("b","c","c","d"),
    wto=c(.249,.25,.251,NA_real_))
  x$results$uva$wto$pairwise <- pairs
  x$results$stages$reason <- c("", "", "", "Check the empirical graph.")
  for (language in c("eng","esp")) {
    sections <- .pm_network_sections(x, language)$sections
    tab <- sections[[if(language=="eng")"Redundant pairs"else"Pares redundantes"]]
    expect_equal(tab$wto,c(.25,.251))
    expect_false(any(c("Versions","Versiones") %in% names(sections)))
    method <- sections[[if(language=="eng")"Method"else"M\u00e9todo"]]$Text
    expect_false(any(grepl("seed 2401|semilla 2401",method)))
    expect_false(any(c("Stages and conditions","Etapas y condiciones") %in% names(sections)))
    stage_table <- .pm_network_sections(x, language, detail="audit")$sections[[if(language=="eng")"Stages and conditions"else"Etapas y condiciones"]]
    expect_false(any(c("reason","Motivo") %in% names(stage_table)))
    expect_equal(nrow(stage_table),nrow(x$results$stages))
  }
  expect_identical(x$results$uva$wto$pairwise,pairs)
  expect_true(nrow(x$diagnostics$versions)>0)
  printed <- capture.output(print(x))
  expect_false(any(grepl("^[1-3]\\. reason:",printed)))
  x$results$uva$wto$pairwise <- pairs[1,,drop=FALSE]
  expect_match(.pm_network_sections(x,"eng")$sections[["Redundant pairs"]]$Text,
    "No pairs.*0.25")
  x$results$tefi$value <- NA_real_
  path <- tempfile(fileext=".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  expect_warning(plot(x,type="tefi"),NA)
})

test_that("an empty RI graph is reported as not evaluable without TEFI warnings", {
  wc <- setNames(c(1,1,2),letters[1:3])
  graph <- matrix(0,3,3,dimnames=list(letters[1:3],letters[1:3]))
  original <- list(wc=wc,network=graph,n.dim=2,TEFI=-1,correlation=diag(3))
  adjusted <- list(EGA=list(wc=setNames(rep(NA_real_,3),letters[1:3]),
    network=graph,n.dim=0),TEFI=NaN)
  local_mocked_bindings(.pm_network_ega=function(data,settings,variant="EGA")
    if(variant=="riEGA")adjusted else original,
    .pm_network_validate_ri=function(adjusted)list(converged=TRUE,admissible=TRUE))
  d <- data.frame(a=1:20,b=20:1,c=rep(1:5,4))
  settings <- list(corr="spearman",missing="listwise",reverse=character())
  expect_warning(w <- .pm_network_wording(d,d,original,settings),NA)
  expect_identical(w$status,"not_evaluable")
  expect_match(w$reason,"no estimable communities")
  expect_true(all(is.na(w$tefi$value[2:3])))
  expect_false(any(w$tefi$comparable[2:3]))
  adjusted$EGA$wc <- setNames(c(1,NA_real_,1),letters[1:3])
  adjusted$EGA$n.dim <- 1
  adjusted$EGA$network[1,3] <- adjusted$EGA$network[3,1] <- .4
  expect_warning(partial <- .pm_network_wording(d,d,original,settings),NA)
  expect_identical(partial$status,"not_evaluable")
  expect_true(all(is.na(partial$tefi$value[2:3])))
})

test_that("unavailable RI diagnostics do not block main bootstrap or invariance", {
  skip_if_not_installed("EGAnet");skip_if_not_installed("bootnet")
  skip_if_not_installed("NetworkComparisonTest")
  d <- simulate_likert_fixture(n=180,seed=397)
  d$sex <- rep(c("F","M"),90)
  local_mocked_bindings(
    .pm_network_wording=function(raw,keyed,ega,settings)
      list(status="not_evaluable",reason="Residual partition unavailable.",tefi=data.frame()),
    .pm_network_boot_ega=function(...)stop("RI bootstrap must not be called"),
    .pm_network_bootstrap=function(...)list(main_bootstrap_completed=TRUE),
    .pm_network_invariance=function(...)list(invariance_completed=TRUE))
  x <- network_auto(d,items=1:12,group="sex",uva="diagnostic",hierarchy=FALSE,
    n_boot=2L,n_perm=2L,consensus_iter=5L,verbose=FALSE)
  stages <- setNames(x$results$stages$status,x$results$stages$stage)
  expect_equal(unname(stages[c("riEGA","riEGA stability")]),rep("not_evaluable",2))
  expect_equal(unname(stages[c("bootstrap","invariance")]),rep("completed",2))
  expect_true(x$results$bootstrap$main_bootstrap_completed)
  expect_true(x$results$invariance$invariance_completed)
  expect_equal(nrow(x$errors),0)
  expect_identical(x$conclusion$status,"partial")
})
