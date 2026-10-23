test_that("EFA omega shares standalone reliability estimation for each dimension", {
  skip_if_not_installed("semTools")
  set.seed(260926)
  f <- rnorm(400)
  continuous <- as.data.frame(sapply(c(.5, .65, .75, .85), function(a)
    a * f + sqrt(1 - a^2) * rnorm(400)))
  names(continuous) <- paste0("i", 1:4)
  ordinal <- as.data.frame(lapply(continuous, function(x)
    as.numeric(cut(x, c(-Inf, -.8, -.2, .5, 1.2, Inf)))))
  distribution <- data.frame(Item = names(continuous), Dimension = "F1")
  incomplete <- ordinal
  incomplete[1:5, 1] <- NA
  for (d in list(continuous, ordinal, incomplete)) {
    standalone <- reliability_auto(d, verbose = FALSE)
    efa <- PsychoMatic:::.efa_reliability(d, distribution)$F1
    expected <- standalone$coefficients$estimate[standalone$coefficients$coefficient == "omega"]
    expect_true(is.finite(expected))
    expect_equal(efa$omega, expected, tolerance = 1e-10)
    expect_identical(efa$omega_status, "ok")
    expect_equal(efa$omega_diagnostics$n, sum(complete.cases(d)))
  }
})
