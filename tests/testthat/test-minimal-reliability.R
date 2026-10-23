test_that("automatic omega matches congeneric SEM and item deletion refits", {
  skip_if_not_installed("semTools")
  set.seed(817)
  z <- rnorm(400)
  x <- as.data.frame(sapply(c(.6, .7, .8, .9), function(l) l*z+rnorm(400)))
  names(x) <- paste0("i", 1:4)
  out <- reliability_auto(x, level = "continuous", verbose = FALSE)
  expect_true("omega" %in% out$coefficients$coefficient)
  fit <- lavaan::cfa("f =~ i1+i2+i3+i4", x, std.lv = TRUE)
  expect_equal(out$coefficients$estimate[out$coefficients$coefficient == "omega"],
               as.numeric(unlist(semTools::compRelSEM(fit))), tolerance = 1e-6)
  expect_equal(nrow(out$item_stats), 4L)
  expect_equal(out$item_stats$item_rest_r[1], cor(x$i1, rowSums(x[-1])))
  reduced <- lavaan::cfa("f =~ i2+i3+i4", x, std.lv = TRUE)
  expect_equal(out$item_stats$omega_if_deleted[1], as.numeric(unlist(semTools::compRelSEM(reduced))), tolerance = 1e-6)
})

test_that("omega closes on unidentified deletion and missing dependency", {
  x <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")]
  out <- reliability_auto(x, level = "continuous", verbose = FALSE)
  expect_true(all(is.na(out$item_stats$omega_if_deleted)))
  expect_true(all(out$item_stats$omega_if_deleted_reason == "insufficient_items"))
  local_mocked_bindings(.pm_semtools_available = function() FALSE, .package = "PsychoMatic")
  out <- reliability_auto(x, level = "continuous", verbose = FALSE)
  omega <- subset(out$coefficients, coefficient == "omega")
  expect_equal(omega$reason, "semtools_not_available")
  expect_true(is.na(omega$estimate))
})

test_that("CFA factor alpha uses its own observed indicators", {
  fit <- lavaan::cfa("a =~ x1+x2+x3\nb =~ x4+x5+x6", lavaan::HolzingerSwineford1939)
  rows <- PsychoMatic:::.pm_standard_reliability_values(fit)$rows
  observed <- subset(rows, coefficient == "alpha_observed" & factor == "a")
  s <- cov(lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")])
  expect_equal(observed$value, 3/2*(1-sum(diag(s))/sum(s)))
})

test_that("ordinal omega targets observed categories and reports missingness", {
  skip_if_not_installed("semTools")
  set.seed(721)
  z <- rnorm(500)
  x <- as.data.frame(replicate(4, as.integer(cut(z + rnorm(500), c(-Inf, -.8, .1, .9, Inf)))))
  names(x) <- paste0("i", 1:4)
  x[1, 1] <- NA
  out <- reliability_auto(x, level = "ordinal", verbose = FALSE)
  fit <- lavaan::cfa("f =~ i1+i2+i3+i4", x, std.lv = TRUE, ordered = names(x))
  omega <- subset(out$coefficients, coefficient == "omega")
  expect_equal(omega$estimate, as.numeric(unlist(semTools::compRelSEM(fit,
    obs.var = TRUE, tau.eq = FALSE, ord.scale = TRUE))), tolerance = 1e-6)
  expect_equal(omega$estimand, "observed_score_congeneric")
  expect_equal(out$item_stats$n, rep(499L, 4))
  expect_equal(out$item_stats$omega_if_deleted_n, c(500L, 499L, 499L, 499L))
  table <- PsychoMatic:::.pm_analysis_sections(out, "esp")$sections[["Analisis de items: confiabilidad al eliminar"]]
  expect_equal(table[["N casos"]], rep(499L, 4))
  expect_false("N para omega" %in% names(table))
  rows <- PsychoMatic:::.pm_standard_reliability_values(fit)$rows
  expect_true(all(c("alpha_observed", "alpha_ordinal") %in% rows$coefficient))
})

test_that("automatic omega fails closed for constant items", {
  constant <- data.frame(a = rep(1, 30), b = 1:30, c = 30:1)
  out <- suppressWarnings(reliability_auto(constant, level = "continuous", verbose = FALSE))
  expect_equal(subset(out$coefficients, coefficient == "omega")$reason, "constant_item")
})

test_that("automatic omega fails closed for a native improper model", {
  .pm_skip_extended()
  set.seed(33)
  a <- rnorm(300); b <- rnorm(300)
  x <- data.frame(a = a, b = b, c = a + b + rnorm(300, sd = .1))
  out <- reliability_auto(x, level = "continuous", verbose = FALSE)
  expect_equal(subset(out$coefficients, coefficient == "omega")$status, "not_evaluable")
  expect_true(length(out$diagnostics$omega$warnings) > 0L)
})
