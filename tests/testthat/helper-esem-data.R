esem_fixture <- function(ordinal = FALSE, missing = FALSE) {
  d <- PsychoMatic:::.efa_with_seed(101, lavaan::simulateData('
    a =~ .8*x1 + .7*x2 + .75*x3 + .12*x4
    b =~ .8*x4 + .7*x5 + .75*x6 + .12*x7
    c =~ .8*x7 + .7*x8 + .75*x9 + .12*x1
    a ~~ .25*b + .2*c
    b ~~ .3*c', sample.nobs = 400))
  if (ordinal) for (nm in names(d)) d[[nm]] <- cut(d[[nm]],
    c(-Inf, stats::quantile(d[[nm]], c(.2, .4, .6, .8)), Inf), labels = FALSE)
  if (missing) { d[1:4, 1] <- NA; d[5:8, 5] <- NA }
  d
}
esem_factors <- function() list(A = paste0('x', 1:3), B = paste0('x', 4:6), C = paste0('x', 7:9))
esem_prepare <- function(d = esem_fixture(), factors = esem_factors(),
                         n_factors = NULL, items = NULL, rotation = 'target',
                         ordered = NULL, estimator = NULL, missing = 'listwise',
                         rotation_args = list()) {
  PsychoMatic:::.pm_esem_prepare(d, factors, n_factors, items, rotation,
                               ordered, estimator, missing, rotation_args, 'eng')
}
