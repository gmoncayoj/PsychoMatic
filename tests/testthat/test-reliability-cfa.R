make_bifactor_nine_fixture <- function(n = 500L, seed = 2505L, groups = FALSE) {
  .pm_fixture_with_seed(seed, {
    group <- if (isTRUE(groups)) factor(rep(c("A", "B"), length.out = n)) else NULL
    g <- stats::rnorm(n)
    s <- matrix(stats::rnorm(n * 3L), nrow = n)
    y <- vapply(seq_len(9L), function(j) {
      group_shift <- if (isTRUE(groups) && j == 1L) ifelse(group == "B", 0.25, 0) else 0
      0.65 * g + 0.40 * s[, ceiling(j / 3)] + group_shift +
        stats::rnorm(n, sd = sqrt(1 - 0.65^2 - 0.40^2))
    }, numeric(n))
    out <- as.data.frame(y, stringsAsFactors = FALSE)
    names(out) <- paste0("i", seq_len(9L))
    if (isTRUE(groups)) out <- data.frame(group = group, out, check.names = FALSE)
    out
  })
}

bifactor_nine_model <- function(correlated = FALSE) {
  syntax <- c(
    paste("g =~", paste0("i", 1:9, collapse = " + ")),
    "s1 =~ i1 + i2 + i3",
    "s2 =~ i4 + i5 + i6",
    "s3 =~ i7 + i8 + i9"
  )
  if (isTRUE(correlated)) syntax <- c(syntax, "g ~~ 0.25*s1")
  paste(syntax, collapse = "\n")
}

test_that("Spearman-Brown keeps its sign and fails at a zero denominator", {
  expect_equal(PsychoMatic:::.pm_spearman_brown(-0.50), -2)
  expect_true(is.na(PsychoMatic:::.pm_spearman_brown(-1)))
})

test_that("EFA two-item reliability no longer makes negative correlations positive", {
  d <- data.frame(i1 = 1:10, i2 = c(10, 8, 9, 7, 6, 5, 4, 3, 1, 2))
  distribution <- data.frame(
    Item = c("i1", "i2"), Dimension = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  r <- stats::cor(d$i1, d$i2)
  out <- PsychoMatic:::.efa_reliability(d, distribution)
  expect_equal(out$F1$spearman_brown, 2 * r / (1 + r), tolerance = 1e-15)
  expect_lt(out$F1$spearman_brown, 0)
})

test_that("omega uses the full residual covariance matrix", {
  lambda <- matrix(c(0.8, 0.7, 0.6), ncol = 1)
  psi <- matrix(1, nrow = 1, ncol = 1)
  theta <- matrix(c(0.36, 0.08, 0, 0.08, 0.51, 0, 0, 0, 0.64),
                  3, 3, byrow = TRUE)
  common <- lambda %*% psi %*% t(lambda)
  expected <- sum(common) / sum(common + theta)
  expect_equal(
    PsychoMatic:::.pm_omega_from_matrices(lambda, psi, theta),
    expected,
    tolerance = 1e-12
  )
})

test_that("fit selection uses robust, then scaled, then base by index", {
  fit <- c(cfi = 0.90, cfi.scaled = 0.91, cfi.robust = 0.92)
  picked <- PsychoMatic:::.pm_pick_fit_variant(fit, "cfi")
  expect_equal(picked$value, 0.92)
  expect_identical(picked$variant, "robust")

  fit["cfi.robust"] <- NA_real_
  expect_identical(PsychoMatic:::.pm_pick_fit_variant(fit, "cfi")$variant, "scaled")
  fit["cfi.scaled"] <- Inf
  expect_identical(PsychoMatic:::.pm_pick_fit_variant(fit, "cfi")$variant, "base")
})

test_that("PUC requires complete, disjoint, orthogonal bifactor structure", {
  membership <- list(
    g = paste0("i", 1:9),
    s1 = paste0("i", 1:3),
    s2 = paste0("i", 4:6),
    s3 = paste0("i", 7:9)
  )
  ok <- PsychoMatic:::.pm_validate_puc_structure(membership, general = "g")
  expect_true(ok$supported)
  expect_equal(ok$puc, 27 / 36, tolerance = 1e-15)

  overlap <- membership
  overlap$s2 <- c("i3", "i4", "i5")
  bad_overlap <- PsychoMatic:::.pm_validate_puc_structure(overlap, general = "g")
  expect_false(bad_overlap$supported)
  expect_identical(bad_overlap$code, "unsupported_overlap")
  expect_true(is.na(bad_overlap$puc))

  missing_general <- membership
  missing_general$g <- paste0("i", 1:8)
  bad_general <- PsychoMatic:::.pm_validate_puc_structure(missing_general, general = "g")
  expect_false(bad_general$supported)
  expect_identical(bad_general$code, "missing_general_loadings")
  expect_true(is.na(bad_general$puc))

  nonorthogonal <- PsychoMatic:::.pm_validate_puc_structure(
    membership, general = "g", orthogonal = FALSE
  )
  expect_false(nonorthogonal$supported)
  expect_identical(nonorthogonal$code, "nonorthogonal_factors")
  expect_true(is.na(nonorthogonal$puc))
})

test_that("MVN Mardia uses 6.3 arguments and preserves raw precision", {
  received <- NULL
  fake <- list(multivariate_normality = data.frame(
    Test = c("Mardia Kurtosis", "Mardia Skewness"),
    Statistic = c(1.234567890123, 17.987654321098),
    p.value = c(0.219876543210, 0.031234567890),
    Method = "asymptotic",
    check.names = FALSE
  ))
  testthat::local_mocked_bindings(
    .pm_mvn_available = function() TRUE,
    .pm_mvn_call = function(data, mvn_test, use_population, alpha, tidy, descriptives) {
      received <<- list(
        data = data,
        mvn_test = mvn_test,
        use_population = use_population,
        alpha = alpha,
        tidy = tidy,
        descriptives = descriptives
      )
      fake
    },
    .package = "PsychoMatic"
  )
  d <- data.frame(x = 1:8, y = c(2, 4, 3, 8, 7, 9, 11, 10))
  out <- PsychoMatic:::.pm_mvn_mardia(d, alpha = 0.037)
  expect_identical(received$mvn_test, "mardia")
  expect_identical(received$use_population, TRUE)
  expect_identical(received$alpha, 0.037)
  expect_identical(received$tidy, FALSE)
  expect_identical(received$descriptives, FALSE)
  expect_identical(out$status, "ok")
  expect_identical(out$skew_stat, 17.987654321098)
  expect_identical(out$kurt_stat, 1.234567890123)
  expect_identical(out$p_skew, 0.031234567890)
  expect_identical(out$p_kurt, 0.219876543210)
  expect_identical(out$normal, FALSE)
})

test_that("MVN normality decision applies alpha inclusively", {
  fake <- list(multivariate_normality = data.frame(
    Test = c("Mardia Skewness", "Mardia Kurtosis"),
    Statistic = c(1, 2),
    p.value = c(0.05, 0.05),
    check.names = FALSE
  ))
  testthat::local_mocked_bindings(
    .pm_mvn_available = function() TRUE,
    .pm_mvn_call = function(...) fake,
    .package = "PsychoMatic"
  )
  out <- PsychoMatic:::.pm_mvn_mardia(data.frame(x = 1:8, y = 2:9), 0.05)
  expect_identical(out$normal, TRUE)
})

test_that("MVN absence and errors are explicit and preserve the condition", {
  testthat::local_mocked_bindings(
    .pm_mvn_available = function() FALSE,
    .package = "PsychoMatic"
  )
  absent <- PsychoMatic:::.pm_mvn_mardia(data.frame(x = 1:5, y = 2:6), 0.05)
  expect_identical(absent$status, "not_available")

  sentinel <- structure(
    list(message = "sentinel MVN failure", call = quote(MVN::mvn(probe))),
    class = c("psychomatic_mvn_sentinel", "error", "condition")
  )
  testthat::local_mocked_bindings(
    .pm_mvn_available = function() TRUE,
    .pm_mvn_call = function(...) stop(sentinel),
    .package = "PsychoMatic"
  )
  failed <- PsychoMatic:::.pm_mvn_mardia(data.frame(x = 1:5, y = 2:6), 0.05)
  expect_identical(failed$status, "not_evaluable")
  expect_s3_class(failed$condition, "psychomatic_mvn_sentinel")
  expect_identical(conditionMessage(failed$condition), "sentinel MVN failure")
  expect_identical(conditionCall(failed$condition), quote(MVN::mvn(probe)))
})

test_that("MVN Mardia agrees with the installed MVN 6.3 raw table", {
  set.seed(2506)
  d <- as.data.frame(matrix(stats::rnorm(360), ncol = 3))
  names(d) <- c("x", "y", "z")
  direct <- MVN::mvn(
    d,
    mvn_test = "mardia",
    use_population = TRUE,
    alpha = 0.05,
    tidy = FALSE,
    descriptives = FALSE
  )$multivariate_normality
  out <- PsychoMatic:::.pm_mvn_mardia(d, 0.05)
  expect_identical(out$status, "ok")
  expect_equal(
    c(out$skew_stat, out$kurt_stat),
    direct$Statistic[match(c("Mardia Skewness", "Mardia Kurtosis"), direct$Test)],
    tolerance = 1e-14
  )
  expect_equal(
    c(out$p_skew, out$p_kurt),
    direct$p.value[match(c("Mardia Skewness", "Mardia Kurtosis"), direct$Test)],
    tolerance = 1e-14
  )
})

test_that("bounded Mardia agrees with MVN including its small-sample correction", {
  x <- matrix(
    c(1, 2, 2, 1, 3, 4, 4, 3, 5, 7, 7, 5, 8, 6, 9, 10),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("x", "y"))
  )
  direct <- MVN::mvn(
    as.data.frame(x),
    mvn_test = "mardia",
    use_population = TRUE,
    alpha = 0.05,
    tidy = FALSE,
    descriptives = FALSE
  )$multivariate_normality
  bounded <- PsychoMatic:::.pm_mardia_blockwise(x, block_size = 2L)
  expect_equal(
    unname(bounded[c("skew_statistic", "kurtosis_z")]),
    direct$Statistic[match(c("Mardia Skewness", "Mardia Kurtosis"), direct$Test)],
    tolerance = 1e-8
  )
  expect_equal(
    unname(bounded[c("skew_p", "kurtosis_p")]),
    direct$p.value[match(c("Mardia Skewness", "Mardia Kurtosis"), direct$Test)],
    tolerance = 1e-8
  )
})

test_that("reliability_auto identifies estimand, matrix, keying, and status", {
  d <- simulate_likert_fixture(220L, 801L)
  out <- reliability_auto(d)
  expect_s3_class(out, "psychomatic_reliability")
  required <- c("coefficient", "estimate", "estimand", "matrix", "keying", "status", "reason")
  expect_true(all(required %in% names(out$results$coefficients)))
  ok <- out$results$coefficients$status == "ok"
  expect_true(all(is.finite(out$results$coefficients$estimate[ok])))
  expect_identical(out$effective_settings$level, "ordinal")
  expect_identical(out$diagnostics$correlation_matrix |> rownames(), names(d))
  expect_identical(out$diagnostics$correlation_matrix |> colnames(), names(d))
})

test_that("observed alpha uses raw data and matches psych", {
  set.seed(2507)
  d <- as.data.frame(matrix(stats::rnorm(800), ncol = 4))
  names(d) <- paste0("x", 1:4)
  d$x2 <- 0.5 * d$x1 + stats::rnorm(nrow(d), sd = 0.7)
  d$x3 <- 0.4 * d$x1 + stats::rnorm(nrow(d), sd = 0.8)
  d$x4 <- 0.3 * d$x1 + stats::rnorm(nrow(d), sd = 0.9)
  expected_fit <- NULL
  invisible(utils::capture.output({
    expected_fit <- suppressWarnings(psych::alpha(d, check.keys = FALSE, warnings = FALSE))
  }))
  expected <- expected_fit$total$raw_alpha
  out <- reliability_auto(d, level = "continuous")
  row <- out$results$coefficients[out$results$coefficients$coefficient == "alpha_observed", ]
  expect_identical(row$matrix, "raw_covariance")
  expect_equal(row$estimate, unname(expected), tolerance = 1e-12)
})

test_that("ordinal alpha is calculated from an explicitly named ordinal matrix", {
  d <- simulate_likert_fixture(300L, 2508L)[1:6]
  out <- reliability_auto(d, level = "ordinal")
  R <- out$diagnostics$correlation_matrix
  p <- ncol(R)
  expected <- p / (p - 1) * (1 - sum(diag(R)) / sum(R))
  row <- out$results$coefficients[out$results$coefficients$coefficient == "alpha_ordinal", ]
  expect_identical(row$matrix, "polychoric")
  expect_identical(rownames(R), names(d))
  expect_identical(colnames(R), names(d))
  expect_equal(row$estimate, expected, tolerance = 1e-12)
})

test_that("reliability fail-closes edge cases without clamping", {
  one <- reliability_auto(data.frame(x = 1:10), level = "continuous")
  expect_true(all(one$results$coefficients$status == "not_evaluable"))
  expect_true(any(one$results$coefficients$reason == "insufficient_items"))

  zero_total <- reliability_auto(
    data.frame(x = rep(c(-1, 1), 10), y = rep(c(1, -1), 10)),
    level = "continuous"
  )
  expect_true(any(zero_total$results$coefficients$reason == "zero_total_score_variance"))

  nonfinite <- reliability_auto(
    data.frame(x = c(1, 2, Inf, 4), y = c(1, 2, 3, 4)),
    level = "continuous"
  )
  expect_true(any(nonfinite$results$coefficients$reason == "non_finite_input"))
  expect_true(all(nonfinite$results$coefficients$status == "not_evaluable"))

  invalid <- PsychoMatic:::.pm_reliability_row(
    "omega_total", 1.04, "model_based", "model_implied", "none",
    domain = c(0, 1)
  )
  expect_identical(invalid$status, "not_evaluable")
  expect_identical(invalid$reason, "outside_mathematical_domain")
  expect_identical(invalid$diagnostic_value[[1]], 1.04)
  expect_true(is.na(invalid$estimate))
})

test_that("two-item reliability returns signed Spearman-Brown", {
  d <- data.frame(x = 1:8, y = 8:1)
  out <- reliability_auto(d, level = "continuous")
  row <- out$results$coefficients[
    out$results$coefficients$coefficient == "spearman_brown",
  ]
  expect_identical(row$status, "not_evaluable")
  expect_identical(row$reason, "zero_denominator")

  d2 <- data.frame(x = 1:10, y = c(10, 8, 9, 7, 6, 5, 4, 3, 1, 2))
  r <- stats::cor(d2$x, d2$y)
  expected <- 2 * r / (1 + r)
  out2 <- reliability_auto(d2, level = "continuous")
  row2 <- out2$results$coefficients[
    out2$results$coefficients$coefficient == "spearman_brown",
  ]
  expect_equal(row2$estimate, expected, tolerance = 1e-12)
  expect_lt(row2$estimate, 0)
})

test_that("CFA model reliability agrees with current semTools APIs", {
  skip_if_not_installed("semTools")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    sep = "\n"
  )
  fit <- lavaan::cfa(model, data = hs, std.lv = TRUE)
  out <- reliability_auto(hs, level = "continuous", model = fit)
  rows <- out$results$coefficients
  cr <- semTools::compRelSEM(fit)
  ave <- semTools::AVE(fit)
  for (factor in names(cr)) {
    observed <- rows$estimate[rows$coefficient == "composite_reliability" & rows$factor == factor]
    expect_equal(observed, as.numeric(cr[[factor]]), tolerance = 1e-12)
  }
  for (factor in names(ave)) {
    observed <- rows$estimate[rows$coefficient == "ave" & rows$factor == factor]
    expect_equal(observed, unname(ave[[factor]]), tolerance = 1e-12)
  }
})

test_that("CFA runs when MVN and semTools namespaces are unavailable", {
  testthat::local_mocked_bindings(
    .pm_mvn_available = function() FALSE,
    .pm_semtools_available = function() FALSE,
    .package = "PsychoMatic"
  )
  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")]
  out <- suppressWarnings(cfa_auto(
    hs,
    model = "visual =~ x1 + x2 + x3",
    estimator = "ML",
    language = "eng"
  ))
  expect_s3_class(out, "cfa_auto")
  expect_true(out$converged)
  expect_identical(out$reliability$status, "not_available")
})

test_that("CFA promotes captured dependency conditions to structured fields", {
  rel_warning <- structure(
    list(message = "sentinel reliability warning", call = quote(compRelSEM(probe))),
    class = c("psychomatic_rel_warning", "warning", "condition")
  )
  testthat::local_mocked_bindings(
    .compute_standard_reliability = function(...) list(
      type = "standard",
      status = "ok",
      reason = NA_character_,
      values = matrix(0.8, nrow = 1L, dimnames = list("omega", "visual")),
      warnings = list(rel_warning),
      messages = list()
    ),
    .package = "PsychoMatic"
  )
  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")]
  out <- suppressWarnings(cfa_auto(
    hs, model = "visual =~ x1 + x2 + x3", estimator = "ML", language = "eng"
  ))
  row <- out$warnings[out$warnings$message == "sentinel reliability warning", ]
  expect_equal(nrow(row), 1L)
  expect_true(inherits(row$condition[[1L]], "psychomatic_rel_warning"))

  mvn_error <- structure(
    list(message = "sentinel diagnostic failure", call = quote(MVN::mvn(probe))),
    class = c("psychomatic_mvn_error", "error", "condition")
  )
  testthat::local_mocked_bindings(
    .pm_mvn_mardia = function(...) list(
      status = "not_evaluable", reason = "mvn_error", result = NULL,
      condition = mvn_error, normal = NA, warnings = list(), messages = list(),
      skew_stat = NA_real_, kurt_stat = NA_real_,
      p_skew = NA_real_, p_kurt = NA_real_
    ),
    .package = "PsychoMatic"
  )
  out_error <- suppressWarnings(cfa_auto(
    lavaan::HolzingerSwineford1939[, paste0("x", 1:9)],
    model = paste(
      "visual =~ x1 + x2 + x3",
      "textual =~ x4 + x5 + x6",
      "speed =~ x7 + x8 + x9",
      sep = "\n"
    ),
    language = "eng"
  ))
  row_error <- out_error$errors[out_error$errors$message == "sentinel diagnostic failure", ]
  expect_equal(nrow(row_error), 1L)
  expect_true(inherits(row_error$condition[[1L]], "psychomatic_mvn_error"))
})

test_that("CFA preserves an unexpected reliability condition", {
  sentinel <- structure(
    list(message = "sentinel reliability failure", call = quote(compRelSEM(probe))),
    class = c("psychomatic_reliability_error", "error", "condition")
  )
  testthat::local_mocked_bindings(
    .compute_standard_reliability = function(...) stop(sentinel),
    .package = "PsychoMatic"
  )
  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")]
  out <- suppressWarnings(cfa_auto(
    hs, model = "visual =~ x1 + x2 + x3", estimator = "ML", language = "eng"
  ))
  row <- out$errors[out$errors$message == "sentinel reliability failure", ]
  expect_equal(nrow(row), 1L)
  expect_true(inherits(row$condition[[1L]], "psychomatic_reliability_error"))
  expect_identical(conditionCall(row$condition[[1L]]), quote(compRelSEM(probe)))
})

test_that("reliability_auto exposes captured semTools warnings structurally", {
  sentinel <- structure(
    list(message = "sentinel semTools warning", call = quote(compRelSEM(probe))),
    class = c("psychomatic_semtools_warning", "warning", "condition")
  )
  testthat::local_mocked_bindings(
    .pm_semtools_coefficients = function(...) list(
      rows = PsychoMatic:::.pm_reliability_row(
        "composite_reliability", 0.8, "model_based", "model_implied", "none",
        factor = "f", domain = c(0, 1)
      ),
      diagnostics = list(status = "ok", warnings = list(sentinel), messages = list())
    ),
    .package = "PsychoMatic"
  )
  out <- reliability_auto(
    data.frame(x = 1:10, y = 1:10),
    level = "continuous",
    model = structure(list(), class = "lavaan")
  )
  row <- out$warnings[out$warnings$message == "sentinel semTools warning", ]
  expect_equal(nrow(row), 1L)
  expect_true(inherits(row$condition[[1L]], "psychomatic_semtools_warning"))
})

test_that("the bifactor orthogonality control is a trailing compatible argument", {
  expected_prefix <- c(
    "data", "model", "ordered", "estimator", "std.lv", "mi_threshold",
    "n_mi", "alpha_norm", "language"
  )
  actual <- names(formals(cfa_auto))
  expect_identical(actual[seq_along(expected_prefix)], expected_prefix)
  expect_identical(actual[10L], "bifactor_orthogonal")
  expect_identical(formals(cfa_auto)$bifactor_orthogonal, TRUE)

  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")]
  positional <- suppressWarnings(cfa_auto(
    hs, "visual =~ x1 + x2 + x3", NULL, "ML", FALSE, 10, 10, 0.05, "eng"
  ))
  expect_s3_class(positional, "cfa_auto")
  expect_identical(positional$estimator, "ML")
  expect_identical(positional$language, "eng")
})

test_that("cfa_auto records each effective robust fit variant", {
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    sep = "\n"
  )
  out <- suppressWarnings(cfa_auto(hs, model, estimator = "MLR", language = "eng", fit_variant = "robust"))
  expect_true(all(c("Base", "Scaled", "Robust", "Variant") %in% names(out$fit_indices)))
  for (metric in c("cfi", "tli", "rmsea")) {
    row <- out$fit_indices[out$fit_indices$Metric == metric, ]
    direct <- PsychoMatic:::.pm_pick_fit_variant(out$all_fit_measures, metric)
    expect_identical(row$Variant, direct$variant)
    expect_equal(row$Value, direct$value, tolerance = 1e-14)
  }
})

test_that("bifactor omega total and hierarchical follow direct matrices", {
  skip_if_not_installed("semTools")
  d <- make_bifactor_nine_fixture()
  fit <- lavaan::cfa(
    bifactor_nine_model(), data = d, std.lv = TRUE, orthogonal = TRUE
  )
  rel <- PsychoMatic:::.compute_bifactor_reliability(
    fit,
    general_factor = "g",
    specific_factors = c("s1", "s2", "s3"),
    observed_items = paste0("i", 1:9)
  )
  matrices <- lavaan::lavInspect(fit, "std")
  lambda <- matrices$lambda
  psi <- matrices$psi
  theta <- matrices$theta
  common <- lambda %*% psi %*% t(lambda)
  denominator <- sum(common + theta)
  general <- lambda[, "g", drop = FALSE] %*%
    psi["g", "g", drop = FALSE] %*%
    t(lambda[, "g", drop = FALSE])
  expect_equal(rel$omega_total, sum(common) / denominator, tolerance = 1e-12)
  expect_equal(rel$omega_h, sum(general) / denominator, tolerance = 1e-12)
  expect_equal(rel$puc, 27 / 36, tolerance = 1e-15)
  expect_equal(
    rel$omega_total,
    as.numeric(semTools::compRelSEM(fit, obs.var = FALSE)[["g"]]),
    tolerance = 0.01
  )
})

test_that("PUC is unavailable when fitted bifactor factors are not orthogonal", {
  d <- make_bifactor_nine_fixture()
  fit <- suppressWarnings(lavaan::cfa(
    bifactor_nine_model(correlated = TRUE), data = d, std.lv = TRUE
  ))
  rel <- PsychoMatic:::.compute_bifactor_reliability(
    fit, "g", c("s1", "s2", "s3"), paste0("i", 1:9)
  )
  expect_true(is.na(rel$puc))
  expect_identical(rel$puc_diagnostic$code, "nonorthogonal_factors")
})

test_that("bifactor reliability keeps each multigroup matrix set separate", {
  d <- make_bifactor_nine_fixture(n = 600L, groups = TRUE)
  fit <- lavaan::cfa(
    bifactor_nine_model(), data = d, group = "group", std.lv = TRUE,
    orthogonal = TRUE
  )
  rel <- PsychoMatic:::.compute_bifactor_reliability(
    fit, "g", c("s1", "s2", "s3"), paste0("i", 1:9)
  )
  matrices <- lavaan::lavInspect(fit, "std")
  expect_length(rel$groups, 2L)
  expect_identical(names(rel$groups), names(matrices))
  for (group in names(matrices)) {
    expected <- PsychoMatic:::.pm_omega_from_matrices(
      matrices[[group]]$lambda,
      matrices[[group]]$psi,
      matrices[[group]]$theta
    )
    expect_equal(rel$groups[[group]]$omega_total, expected, tolerance = 1e-12)
  }
  expect_false(isTRUE(all.equal(
    rel$groups[[1]]$lambda,
    rel$groups[[2]]$lambda,
    tolerance = 1e-12
  )))
})

test_that("modification suggestions preserve exact lavaan syntax", {
  mi <- data.frame(
    lhs = c("factor", "i1", "i1", "i1"),
    op = c("=~", "~~", "~", "|"),
    rhs = c("item", "i2", "1", "t1"),
    mi = c(20, 19, 18, 17),
    epc = c(0.31, 0.22, 0.13, 0.04),
    stringsAsFactors = FALSE
  )
  out <- PsychoMatic:::.build_respec_suggestions(
    mi, list(factor = c("i1", "i2")), language = "eng"
  )
  expect_identical(
    vapply(out, `[[`, character(1), "syntax"),
    c("factor =~ item", "i1 ~~ i2", "i1 ~ 1", "i1 | t1")
  )
  expect_identical(vapply(out, `[[`, character(1), "lhs"), mi$lhs)
  expect_identical(vapply(out, `[[`, character(1), "op"), mi$op)
  expect_identical(vapply(out, `[[`, character(1), "rhs"), mi$rhs)
  required <- c("type", "syntax", "mi", "epc", "same_factor", "theory_warning")
  expect_true(all(vapply(out, function(x) all(required %in% names(x)), logical(1))))
})

test_that("second-order CFA tidy output includes L1/L2 loadings and reliability", {
  skip_if_not_installed("semTools")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    "g =~ visual + textual + speed",
    sep = "\n"
  )
  out <- suppressWarnings(cfa_auto(hs, model, estimator = "ML", language = "eng"))
  tidy <- summary(out)
  expect_true(all(c("level", "Factor", "Item", "Loading") %in% names(tidy$loadings)))
  expect_true(all(c("L1", "L2") %in% unique(tidy$loadings$level)))
  expect_true(all(c("level", "factor", "metric", "value", "status", "reason") %in%
                    names(tidy$reliability)))
  expect_true(all(c("L1", "L2") %in% unique(tidy$reliability$level)))
  expect_true(any(tidy$reliability$metric == "composite_reliability" &
                    tidy$reliability$level == "L1"))
  expect_true(any(tidy$reliability$metric == "composite_reliability" &
                    tidy$reliability$level == "L2"))

  report_tables <- PsychoMatic:::.pm_object_tables(out)
  expect_true(all(c("L1", "L2") %in% unique(report_tables$loadings$level)))
  expect_true(all(c("L1", "L2") %in% unique(report_tables$reliability$level)))
  expect_output(print(out), "Second-order|Segundo Orden")

  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  PsychoMatic:::.cfa_export_excel(out, path)
  expect_true("Reliability" %in% openxlsx::getSheetNames(path))
  exported <- openxlsx::read.xlsx(path, sheet = "Reliability")
  expect_true(all(c("L1", "L2") %in% unique(exported$level)))
})

test_that("public reliability S3 methods expose reusable objects", {
  expect_silent(
    out <- reliability_auto(simulate_likert_fixture(180L, 2510L)[1:5], verbose = FALSE)
  )
  expect_invisible(print(out))
  s <- summary(out)
  expect_true(is.list(s))
  expect_true(is.data.frame(s$coefficients))
  plot_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_path)
  tryCatch(
    expect_invisible(plot(out)),
    finally = {
      grDevices::dev.off()
      unlink(plot_path)
    }
  )
})

test_that("reliability_auto and its S3 methods are exported and registered", {
  expect_true("reliability_auto" %in% getNamespaceExports("PsychoMatic"))
  methods <- getNamespaceInfo(asNamespace("PsychoMatic"), "S3methods")
  registered <- paste(methods[, 1L], methods[, 2L], sep = ".")
  expect_true(all(c(
    "print.psychomatic_reliability",
    "summary.psychomatic_reliability",
    "plot.psychomatic_reliability"
  ) %in% registered))
})

test_that("CFA has a single non-divergent printer", {
  expect_false(exists(
    ".print_cfa_auto_legacy",
    envir = asNamespace("PsychoMatic"),
    inherits = FALSE
  ))
})

test_that("two-item ordinal reliability fail-closes when polychoric estimation fails", {
  sentinel <- structure(
    list(
      message = "sentinel polychoric failure",
      call = quote(psych::polychoric(probe))
    ),
    class = c("psychomatic_polychoric_sentinel", "error", "condition")
  )
  testthat::local_mocked_bindings(
    .pm_capture_conditions = function(...) {
      list(value = NULL, error = sentinel, warnings = list(), messages = list())
    },
    .package = "PsychoMatic"
  )

  expect_no_error(
    out <- reliability_auto(
      data.frame(i1 = rep(1:5, 4), i2 = rep(5:1, 4)),
      level = "ordinal"
    )
  )
  sb <- out$results$coefficients[
    out$results$coefficients$coefficient == "spearman_brown",
  ]
  expect_identical(sb$status, "not_evaluable")
  expect_identical(sb$reason, "polychoric_not_evaluable")
  expect_true(is.na(sb$estimate))
  expect_equal(nrow(out$errors), 1L)
  expect_s3_class(out$errors$condition[[1L]], "psychomatic_polychoric_sentinel")
  expect_identical(
    conditionCall(out$errors$condition[[1L]]),
    quote(psych::polychoric(probe))
  )
})

test_that("bifactor detection requires exactly one general-factor candidate", {
  items <- paste0("i", 1:6)
  ambiguous <- PsychoMatic:::.pm_detect_bifactor_structure(
    list(
      g1 = items,
      g2 = items,
      s1 = items[1:3],
      s2 = items[4:6]
    ),
    observed_items = items,
    is_second_order = FALSE
  )
  expect_false(ambiguous$is_bifactor)
  expect_identical(ambiguous$code, "ambiguous_general_structure")
  expect_identical(ambiguous$general_candidates, c("g1", "g2"))
  expect_null(ambiguous$general_factor)
  expect_false(ambiguous$specialized_indices_supported)

  absent <- PsychoMatic:::.pm_detect_bifactor_structure(
    list(s1 = items[1:3], s2 = items[4:6], s3 = items[c(1, 4, 6)]),
    observed_items = items,
    is_second_order = FALSE
  )
  expect_false(absent$is_bifactor)
  expect_identical(absent$code, "ambiguous_general_structure")
  expect_length(absent$general_candidates, 0L)
  expect_false(absent$specialized_indices_supported)
})

test_that("second-order reliability targets an explicit unit-weighted item composite", {
  lambda <- matrix(
    c(.80, 0, .70, 0, 0, .75, 0, .65),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("i", 1:4), c("f1", "f2"))
  )
  beta <- matrix(c(.90, .60), ncol = 1, dimnames = list(c("f1", "f2"), "g"))
  phi <- matrix(1.20, 1, 1, dimnames = list("g", "g"))
  psi_u <- matrix(c(.20, .04, .04, .30), 2, 2, byrow = TRUE,
                  dimnames = list(c("f1", "f2"), c("f1", "f2")))
  theta <- matrix(
    c(.40, .05, 0, 0,
      .05, .45, 0, 0,
      0, 0, .50, -.03,
      0, 0, -.03, .55),
    4, 4, byrow = TRUE,
    dimnames = list(paste0("i", 1:4), paste0("i", 1:4))
  )
  weights <- c(1, 1, 1, 1)
  higher <- lambda %*% beta %*% phi %*% t(beta) %*% t(lambda)
  lower <- lambda %*% psi_u %*% t(lambda)
  total <- higher + lower + theta
  expected_omega <- as.numeric(
    crossprod(weights, higher %*% weights) /
      crossprod(weights, total %*% weights)
  )
  expected_ave <- sum(diag(higher)) / sum(diag(total))

  out <- PsychoMatic:::.pm_second_order_composite(
    lambda, beta, phi, psi_u, theta, weights
  )
  expect_identical(out$status, "ok")
  expect_identical(out$estimand, "unit_weighted_observed_composite")
  expect_equal(out$weights, weights)
  expect_equal(out$omega, expected_omega, tolerance = 1e-14)
  expect_equal(out$ave, expected_ave, tolerance = 1e-14)
  expect_equal(out$components$higher_order, higher, tolerance = 1e-14)
  expect_equal(out$components$first_order_disturbance, lower, tolerance = 1e-14)
  expect_equal(out$components$residual, theta, tolerance = 1e-14)
})

test_that("second-order CFA rows use the downstream observed-item composite", {
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    "g =~ visual + textual + speed",
    sep = "\n"
  )
  fit <- lavaan::cfa(model, data = hs, std.lv = TRUE)
  m <- lavaan::lavInspect(fit, "std")
  first <- c("visual", "textual", "speed")
  lambda <- m$lambda[, first, drop = FALSE]
  beta <- m$beta[first, "g", drop = FALSE]
  phi <- m$psi["g", "g", drop = FALSE]
  psi_u <- m$psi[first, first, drop = FALSE]
  higher <- lambda %*% beta %*% phi %*% t(beta) %*% t(lambda)
  lower <- lambda %*% psi_u %*% t(lambda)
  weights <- rep(1, nrow(lambda))
  expected <- as.numeric(
    crossprod(weights, higher %*% weights) /
      crossprod(weights, (higher + lower + m$theta) %*% weights)
  )

  rows <- PsychoMatic:::.pm_second_order_matrix_rows(fit, "g")
  observed <- rows$value[
    rows$factor == "g" & rows$metric == "composite_reliability"
  ]
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_identical(
    rows$estimand[rows$factor == "g" & rows$metric == "composite_reliability"],
    "unit_weighted_observed_composite"
  )
})

test_that("semTools outputs normalize vector, data-frame, and list forms to long", {
  vector_rows <- PsychoMatic:::.pm_semtools_long(
    c(f1 = .70, f2 = .80),
    coefficient = "composite_reliability",
    level = "L1"
  )
  expect_identical(
    names(vector_rows),
    c("group", "level", "factor", "coefficient", "value")
  )
  expect_identical(vector_rows$group, c("group1", "group1"))
  expect_identical(vector_rows$factor, c("f1", "f2"))
  expect_equal(vector_rows$value, c(.70, .80))

  frame_rows <- PsychoMatic:::.pm_semtools_long(
    data.frame(group = c("A", "B"), f1 = c(.40, .50), f2 = c(.60, .70)),
    coefficient = "ave",
    group_labels = c("A", "B"),
    level = "L1"
  )
  expect_equal(nrow(frame_rows), 4L)
  expect_identical(unique(frame_rows$group), c("A", "B"))
  expect_identical(unique(frame_rows$factor), c("f1", "f2"))

  list_rows <- PsychoMatic:::.pm_semtools_long(
    list(f1 = c(A = .71, B = .72), f2 = c(A = .81, B = .82)),
    coefficient = "composite_reliability",
    group_labels = c("A", "B"),
    level = "L1"
  )
  expect_equal(nrow(list_rows), 4L)
  expect_equal(
    list_rows$value[list_rows$factor == "f2" & list_rows$group == "B"],
    .82
  )
})

test_that("real semTools multigroup values remain distinct in both CFA adapters", {
  skip_if_not_installed("semTools")
  hs <- lavaan::HolzingerSwineford1939[, c(paste0("x", 1:9), "school")]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    sep = "\n"
  )
  fit <- lavaan::cfa(model, data = hs, group = "school", std.lv = TRUE)
  group_labels <- lavaan::lavInspect(fit, "group.label")
  direct_cr <- semTools::compRelSEM(fit)
  direct_ave <- semTools::AVE(fit)

  standard <- PsychoMatic:::.pm_standard_reliability_values(fit)
  public <- PsychoMatic:::.pm_semtools_coefficients(fit, "none")
  expect_identical(standard$status, "ok")
  expect_equal(nrow(standard$rows), 18L)
  expect_equal(nrow(public$rows), 12L)
  expect_identical(sort(unique(standard$rows$group)), sort(group_labels))
  expect_identical(sort(unique(public$rows$group)), sort(group_labels))

  for (factor in names(direct_cr)) {
    for (group in group_labels) {
      expected_cr <- unname(direct_cr[[factor]][[group]])
      expect_equal(
        standard$rows$value[
          standard$rows$coefficient == "composite_reliability" &
            standard$rows$factor == factor & standard$rows$group == group
        ],
        expected_cr,
        tolerance = 1e-12
      )
      expect_equal(
        public$rows$estimate[
          public$rows$coefficient == "composite_reliability" &
            public$rows$factor == factor & public$rows$group == group
        ],
        expected_cr,
        tolerance = 1e-12
      )
    }
  }
  for (factor in setdiff(names(direct_ave), "group")) {
    for (group in group_labels) {
      expected_ave <- direct_ave[direct_ave$group == group, factor]
      expect_equal(
        standard$rows$value[
          standard$rows$coefficient == "ave" & standard$rows$factor == factor &
            standard$rows$group == group
        ],
        expected_ave,
        tolerance = 1e-12
      )
    }
  }
})

test_that("a converged Heywood CFA is inadmissible and interpretations fail closed", {
  target <- matrix(
    c(1, .8, .8, .8, 1, .5, .8, .5, 1),
    3, 3, byrow = TRUE
  )
  n <- 500L
  raw <- outer(
    seq_len(n), seq_len(3),
    function(i, j) sin(i * (j + .37)) + cos(i * (j + .11))
  )
  orthonormal <- qr.Q(qr(scale(raw, center = TRUE, scale = FALSE)))[, 1:3, drop = FALSE]
  data <- as.data.frame((orthonormal * sqrt(n - 1)) %*% chol(target))
  names(data) <- paste0("x", 1:3)

  out <- suppressWarnings(cfa_auto(
    data,
    model = "f =~ x1 + x2 + x3",
    estimator = "ML",
    language = "eng"
  ))
  expect_true(out$converged)
  expect_false(out$admissible)
  expect_identical(out$admissibility$status, "not_evaluable")
  expect_true(any(out$admissibility$reason_codes %in%
                    c("post_check_failed", "negative_variance", "heywood_loading")))
  expect_true(is.finite(out$fit_indices$Value[out$fit_indices$Metric == "cfi"]))
  evaluable_metrics <- out$fit_indices$Metric %in%
    c("chisq_ratio", "cfi", "tli", "rmsea", "srmr")
  expect_true(all(out$fit_indices$Evaluation[evaluable_metrics] == "Not evaluable"))
  expect_identical(out$reliability$status, "not_evaluable")
  expect_identical(out$reliability$reason, "inadmissible_cfa_solution")
  expect_identical(out$conclusion$status, "not_evaluable")
})

test_that("CFA printer uses the effective Mardia alpha", {
  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3")]
  out <- suppressWarnings(cfa_auto(
    hs,
    model = "visual =~ x1 + x2 + x3",
    estimator = "ML",
    alpha_norm = .10,
    language = "eng"
  ))
  out$mardia <- list(
    skew_stat = 1,
    kurt_stat = 2,
    p_skew = .075,
    p_kurt = .075,
    status = "ok",
    reason = NA_character_
  )
  printed <- capture.output(print(out))
  decisions <- grep("Skewness|Kurtosis", printed, value = TRUE)
  expect_length(decisions, 2L)
  expect_true(all(grepl("\\[Rejected\\]", decisions)))
  expect_identical(out$effective_settings$alpha_norm, .10)
})

test_that("second-order multigroup aliases and CFA printing preserve groups", {
  skip_if_not_installed("semTools")
  hs <- lavaan::HolzingerSwineford1939[, c(paste0("x", 1:9), "school")]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    "g =~ visual + textual + speed",
    sep = "\n"
  )
  fit <- lavaan::cfa(model, data = hs, group = "school", std.lv = TRUE)
  groups <- lavaan::lavInspect(fit, "group.label")
  rel <- PsychoMatic:::.compute_second_order_reliability(fit, "g")
  l2 <- rel$rows[rel$rows$level == "L2", , drop = FALSE]

  expect_setequal(unique(l2$group), groups)
  expect_named(rel$second_order, groups, ignore.order = TRUE)
  for (group in groups) {
    expect_named(rel$second_order[[group]], "g")
    alias <- rel$second_order[[group]][["g"]]
    canonical <- l2[l2$group == group & l2$factor == "g", , drop = FALSE]
    expect_identical(alias$metric, canonical$metric)
    expect_equal(alias$value, canonical$value, tolerance = 1e-12)
  }

  single_fit <- lavaan::cfa(model, data = hs, std.lv = TRUE)
  single_rel <- PsychoMatic:::.compute_second_order_reliability(single_fit, "g")
  expect_named(single_rel$second_order, "g")
  expect_s3_class(single_rel$second_order[["g"]], "data.frame")

  cfa_result <- suppressWarnings(cfa_auto(
    hs[paste0("x", 1:9)], model = model, estimator = "ML", language = "eng"
  ))
  cfa_result$reliability <- rel
  printed <- capture.output(print(cfa_result))
  for (group in groups) {
    expect_true(any(grepl(group, printed, fixed = TRUE)))
  }
})

test_that("reliability print and plot label every multigroup coefficient", {
  skip_if_not_installed("semTools")
  hs <- lavaan::HolzingerSwineford1939[, c(paste0("x", 1:9), "school")]
  model <- paste(
    "visual =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed =~ x7 + x8 + x9",
    sep = "\n"
  )
  fit <- lavaan::cfa(model, data = hs, group = "school", std.lv = TRUE)
  groups <- lavaan::lavInspect(fit, "group.label")
  out <- reliability_auto(
    hs, items = paste0("x", 1:9), level = "continuous", model = fit
  )

  printed <- capture.output(print(out))
  for (group in groups) {
    expect_true(any(grepl(group, printed, fixed = TRUE)))
  }

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  plotted <- plot(out)
  grDevices::dev.off()
  unlink(plot_file)

  expect_true("plot_label" %in% names(plotted))
  grouped <- !is.na(plotted$group) & nzchar(plotted$group)
  expect_true(all(vapply(
    which(grouped),
    function(i) grepl(plotted$group[[i]], plotted$plot_label[[i]], fixed = TRUE),
    logical(1)
  )))
  expect_equal(anyDuplicated(plotted$plot_label[grouped]), 0L)
})

test_that("second-order aliases distinguish labelled, missing, and empty groups", {
  l2_rows <- data.frame(
    group = c("A", "A", NA_character_, NA_character_, "", ""),
    level = "L2",
    factor = "g",
    metric = rep(c("composite_reliability", "ave"), 3L),
    value = c(.111, .112, .221, .222, .331, .332),
    status = "ok",
    reason = NA_character_,
    estimand = "unit_weighted_observed_composite",
    stringsAsFactors = FALSE
  )
  current_rows <- l2_rows
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(status = "not_evaluable"),
    .pm_second_order_matrix_rows = function(...) current_rows,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(NULL, "g")
  expect_identical(out$rows$group, l2_rows$group)
  expect_named(out$second_order, c("A", "<missing>", "<empty>"))
  expect_equal(out$second_order$A$g$value, c(.111, .112))
  expect_equal(out$second_order[["<missing>"]]$g$value, c(.221, .222))
  expect_equal(out$second_order[["<empty>"]]$g$value, c(.331, .332))
  expect_equal(
    sum(vapply(out$second_order, function(group) nrow(group$g), integer(1))),
    nrow(l2_rows)
  )

  current_rows <- l2_rows[FALSE, , drop = FALSE]
  expect_no_error(empty <- PsychoMatic:::.compute_second_order_reliability(NULL, "g"))
  expect_identical(empty$second_order, list())
  expect_identical(empty$status, "not_evaluable")

  current_rows <- l2_rows[
    !is.na(l2_rows$group) & l2_rows$group == "A", , drop = FALSE
  ]
  single <- PsychoMatic:::.compute_second_order_reliability(NULL, "g")
  expect_named(single$second_order, "g")
  expect_equal(single$second_order$g$value, c(.111, .112))
})

test_that("CFA second-order printing attributes mixed group identities exactly once", {
  rows <- data.frame(
    group = c("A", "A", NA_character_, NA_character_, "", ""),
    level = "L2",
    factor = "g",
    metric = rep(c("composite_reliability", "ave"), 3L),
    value = c(.111, .112, .221, .222, .331, .332),
    status = "ok",
    reason = NA_character_,
    estimand = "unit_weighted_observed_composite",
    stringsAsFactors = FALSE
  )
  out <- structure(list(
    language = "eng",
    effective_settings = list(alpha_norm = .05),
    model_type = "second_order",
    factors = "g",
    n_items = 1L,
    n = 3L,
    pct_missing = 0,
    n_complete = 3L,
    n_categories = 5L,
    estimator = "ML",
    estimator_reason = "synthetic printer fixture",
    converged = TRUE,
    mardia = list(p_skew = NA_real_),
    fi_variant_note = "synthetic printer fixture",
    fit_indices = data.frame(
      Metric = character(), Indice = character(), Valor = numeric(),
      Evaluacion = character(), stringsAsFactors = FALSE
    ),
    loadings = data.frame(
      Factor = character(), Item = character(), Lambda = numeric(), SE = numeric(),
      z = numeric(), p = numeric(), IC_inf = numeric(), IC_sup = numeric(),
      Alerta = character(), stringsAsFactors = FALSE
    ),
    so_loadings = NULL,
    r_squared = NULL,
    factor_correlations = NULL,
    reliability = list(type = "second_order", status = "ok", rows = rows),
    respecifications = list()
  ), class = "cfa_auto")
  canonical_groups <- out$reliability$rows$group

  printed <- capture.output(print(out))
  for (group in c("A", "<missing>", "<empty>")) {
    expect_true(any(grepl(group, printed, fixed = TRUE)))
  }
  for (value in rows$value) {
    expect_equal(
      sum(grepl(sprintf("%.3f", value), printed, fixed = TRUE)),
      1L
    )
  }
  expect_identical(out$reliability$rows$group, canonical_groups)
})

test_that("reliability presentation disambiguates mixed group identities", {
  rows <- data.frame(
    group = c("A", NA_character_, ""),
    factor = "f",
    coefficient = "alpha",
    estimate = c(.711, .722, .733),
    status = "ok",
    reason = NA_character_,
    stringsAsFactors = FALSE
  )
  out <- structure(list(
    effective_settings = list(language = "eng"),
    results = list(coefficients = rows),
    diagnostics = list(),
    conclusion = list()
  ), class = "psychomatic_reliability")
  canonical_groups <- out$results$coefficients$group

  printed <- capture.output(print(out))
  expect_true(any(grepl("<missing>", printed, fixed = TRUE)))
  expect_true(any(grepl("<empty>", printed, fixed = TRUE)))
  expect_identical(out$results$coefficients$group, canonical_groups)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  plotted <- plot(out)
  grDevices::dev.off()
  unlink(plot_file)
  expect_identical(
    plotted$plot_label,
    c("A | f: alpha", "<missing> | f: alpha", "<empty> | f: alpha")
  )
  expect_equal(anyDuplicated(plotted$plot_label), 0L)
  expect_identical(plotted$group, canonical_groups)

  empty <- out
  empty$results$coefficients <- rows[FALSE, , drop = FALSE]
  expect_no_error(capture.output(print(empty)))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_no_error(plotted_empty <- plot(empty))
  grDevices::dev.off()
  unlink(plot_file)
  expect_equal(nrow(plotted_empty), 0L)

  single <- out
  single$results$coefficients <- rows[1L, , drop = FALSE]
  expect_no_error(capture.output(print(single)))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_no_error(plotted_single <- plot(single))
  grDevices::dev.off()
  unlink(plot_file)
  expect_identical(plotted_single$plot_label, "A | f: alpha")
})

test_that("second-order reliability retains groups with missing L2 results", {
  l1_rows <- data.frame(
    group = c("A", "B"),
    level = "L1",
    factor = "f1",
    coefficient = "composite_reliability",
    value = c(.71, .72),
    status = "ok",
    reason = NA_character_,
    stringsAsFactors = FALSE
  )
  l2_rows <- data.frame(
    group = c("A", "A"),
    level = "L2",
    factor = "g",
    metric = c("composite_reliability", "ave"),
    value = c(.81, .61),
    status = "ok",
    reason = NA_character_,
    estimand = c(
      "unit_weighted_observed_composite",
      "mean_observed_variance_due_to_higher_order_factor"
    ),
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "ok", rows = l1_rows, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) l2_rows,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(NULL, "g")
  expect_identical(out$status, "partially_evaluable")
  expect_identical(out$reason, "second_order_partially_evaluable")
  expect_named(out$second_order, c("A", "B"))
  expect_named(out$second_order$B, "g")
  missing_l2 <- out$second_order$B$g
  expect_identical(
    missing_l2$metric,
    c("composite_reliability", "ave")
  )
  expect_true(all(is.na(missing_l2$value)))
  expect_true(all(missing_l2$status == "not_evaluable"))
  expect_true(all(missing_l2$reason == "second_order_group_not_evaluable"))
  canonical_b <- out$rows[
    out$rows$level == "L2" & out$rows$group == "B", , drop = FALSE
  ]
  expect_identical(canonical_b$metric, missing_l2$metric)
  expect_identical(canonical_b$value, missing_l2$value)
})

test_that("reserved-looking literal group values have stable escaped labels", {
  raw_groups <- c(NA_character_, "", "<missing>", "<empty>", "A")
  expect_identical(
    PsychoMatic:::.pm_group_label(raw_groups),
    c("<missing>", "<empty>", "\\<missing>", "\\<empty>", "A")
  )

  l2_rows <- data.frame(
    group = rep(raw_groups[1:4], each = 2L),
    level = "L2",
    factor = "g",
    metric = rep(c("composite_reliability", "ave"), 4L),
    value = seq(.11, .18, by = .01),
    status = "ok",
    reason = NA_character_,
    estimand = "unit_weighted_observed_composite",
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(status = "not_evaluable"),
    .pm_second_order_matrix_rows = function(...) l2_rows,
    .package = "PsychoMatic"
  )
  out <- PsychoMatic:::.compute_second_order_reliability(NULL, "g")
  expect_named(
    out$second_order,
    c("<missing>", "<empty>", "\\<missing>", "\\<empty>")
  )
  expect_equal(out$second_order[["<missing>"]]$g$value, c(.11, .12))
  expect_equal(out$second_order[["\\<missing>"]]$g$value, c(.15, .16))
  expect_equal(out$second_order[["<empty>"]]$g$value, c(.13, .14))
  expect_equal(out$second_order[["\\<empty>"]]$g$value, c(.17, .18))
  expect_identical(out$rows$group, l2_rows$group)
})

test_that("reliability presentation escapes reserved-looking literal groups", {
  rows <- data.frame(
    group = c(NA_character_, "", "<missing>", "<empty>"),
    factor = "f",
    coefficient = "alpha",
    estimate = c(.71, .72, .73, .74),
    status = "ok",
    reason = NA_character_,
    stringsAsFactors = FALSE
  )
  out <- structure(list(
    effective_settings = list(language = "eng"),
    results = list(coefficients = rows),
    diagnostics = list(),
    conclusion = list()
  ), class = "psychomatic_reliability")
  canonical_groups <- out$results$coefficients$group

  printed <- capture.output(print(out))
  printed_rows <- trimws(printed)
  printed_labels <- c(
    "<missing>", "<empty>", "\\\\<missing>", "\\\\<empty>"
  )
  for (label in printed_labels) {
    expect_equal(sum(startsWith(printed_rows, paste0(label, " "))), 1L)
  }
  expect_identical(out$results$coefficients$group, canonical_groups)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  plotted <- plot(out)
  grDevices::dev.off()
  unlink(plot_file)
  expect_identical(
    plotted$plot_label,
    c(
      "<missing> | f: alpha", "<empty> | f: alpha",
      "\\<missing> | f: alpha", "\\<empty> | f: alpha"
    )
  )
  expect_equal(anyDuplicated(plotted$plot_label), 0L)
  expect_identical(plotted$group, canonical_groups)
})

test_that("second-order reliability completes the group-factor-metric grid", {
  l1_rows <- data.frame(
    group = c("A", "B"), level = "L1", factor = "f1",
    coefficient = "composite_reliability", value = c(.71, .72),
    status = "ok", reason = NA_character_, stringsAsFactors = FALSE
  )
  l2_rows <- data.frame(
    group = rep(c("A", "B"), each = 2L), level = "L2", factor = "g1",
    metric = rep(c("composite_reliability", "ave"), 2L),
    value = c(.81, .61, .82, .62), status = "ok", reason = NA_character_,
    estimand = rep(c(
      "unit_weighted_observed_composite",
      "mean_observed_variance_due_to_higher_order_factor"
    ), 2L),
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "ok", rows = l1_rows, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) l2_rows,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(
    NULL, c("g1", "g2", "g1")
  )
  l2 <- out$rows[out$rows$level == "L2", , drop = FALSE]
  expect_identical(
    paste(l2$group, l2$factor, l2$metric, sep = "/"),
    c(
      "A/g1/composite_reliability", "A/g1/ave",
      "A/g2/composite_reliability", "A/g2/ave",
      "B/g1/composite_reliability", "B/g1/ave",
      "B/g2/composite_reliability", "B/g2/ave"
    )
  )
  expect_identical(out$status, "partially_evaluable")
  expect_identical(out$reason, "second_order_partially_evaluable")
  expect_named(out$second_order, c("A", "B"))
  expect_named(out$second_order$A, c("g1", "g2"))
  expect_named(out$second_order$B, c("g1", "g2"))
  missing_factor <- l2[l2$factor == "g2", , drop = FALSE]
  expect_true(all(is.na(missing_factor$value)))
  expect_true(all(missing_factor$status == "not_evaluable"))
  expect_true(all(missing_factor$reason == "second_order_factor_not_evaluable"))
})

test_that("second-order reliability distinguishes a group-specific missing factor", {
  l1_rows <- data.frame(
    group = c("A", "B"), level = "L1", factor = "f1",
    coefficient = "composite_reliability", value = c(.71, .72),
    status = "ok", reason = NA_character_, stringsAsFactors = FALSE
  )
  l2_rows <- data.frame(
    group = c(rep("A", 4L), rep("B", 2L)), level = "L2",
    factor = c(rep(c("g1", "g2"), each = 2L), rep("g1", 2L)),
    metric = rep(c("composite_reliability", "ave"), 3L),
    value = c(.81, .61, .73, .53, .82, .62),
    status = "ok", reason = NA_character_,
    estimand = rep(c(
      "unit_weighted_observed_composite",
      "mean_observed_variance_due_to_higher_order_factor"
    ), 3L),
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "ok", rows = l1_rows, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) l2_rows,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(NULL, c("g1", "g2"))
  missing_factor <- out$rows[
    out$rows$level == "L2" & out$rows$group == "B" &
      out$rows$factor == "g2", , drop = FALSE
  ]
  expect_identical(missing_factor$metric, c("composite_reliability", "ave"))
  expect_true(all(is.na(missing_factor$value)))
  expect_true(all(missing_factor$status == "not_evaluable"))
  expect_true(all(missing_factor$reason == "second_order_factor_not_evaluable"))
  alias <- out$second_order$B$g2
  canonical <- missing_factor[c(
    "metric", "value", "status", "reason"
  )]
  rownames(alias) <- NULL
  rownames(canonical) <- NULL
  expect_equal(alias, canonical)
})

test_that("second-order reliability identifies a missing metric cell", {
  l1_rows <- data.frame(
    group = c("A", "B"), level = "L1", factor = "f1",
    coefficient = "composite_reliability", value = c(.71, .72),
    status = "ok", reason = NA_character_, stringsAsFactors = FALSE
  )
  l2_rows <- data.frame(
    group = c("A", "B", "B"), level = "L2", factor = "g1",
    metric = c("composite_reliability", "composite_reliability", "ave"),
    value = c(.81, .82, .62), status = "ok", reason = NA_character_,
    estimand = c(
      "unit_weighted_observed_composite",
      "unit_weighted_observed_composite",
      "mean_observed_variance_due_to_higher_order_factor"
    ),
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "ok", rows = l1_rows, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) l2_rows,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(NULL, "g1")
  missing_metric <- out$rows[
    out$rows$level == "L2" & out$rows$group == "A" &
      out$rows$factor == "g1" & out$rows$metric == "ave", , drop = FALSE
  ]
  expect_equal(nrow(missing_metric), 1L)
  expect_true(is.na(missing_metric$value))
  expect_identical(missing_metric$status, "not_evaluable")
  expect_identical(missing_metric$reason, "second_order_metric_not_evaluable")
  expect_identical(out$status, "partially_evaluable")
})

test_that("empty L2 output still exposes the complete multigroup grid", {
  l1_rows <- data.frame(
    group = c("A", "B"), level = "L1", factor = "f1",
    coefficient = "composite_reliability", value = c(.71, .72),
    status = "ok", reason = NA_character_, stringsAsFactors = FALSE
  )
  empty_l2 <- data.frame(
    group = character(), level = character(), factor = character(),
    metric = character(), value = numeric(), status = character(),
    reason = character(), estimand = character(), stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "ok", rows = l1_rows, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) empty_l2,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(NULL, c("g1", "g2"))
  l2 <- out$rows[out$rows$level == "L2", , drop = FALSE]
  expect_equal(nrow(l2), 8L)
  expect_identical(unique(l2$group), c("A", "B"))
  expect_identical(unique(l2$factor), c("g1", "g2"))
  expect_identical(unique(l2$metric), c("composite_reliability", "ave"))
  expect_true(all(is.na(l2$value)))
  expect_true(all(l2$status == "not_evaluable"))
  expect_true(all(l2$reason == "second_order_group_not_evaluable"))
  expect_identical(out$status, "not_evaluable")
  expect_identical(out$reason, "second_order_not_evaluable")
  expect_named(out$second_order, c("A", "B"))
  expect_named(out$second_order$A, c("g1", "g2"))
  expect_named(out$second_order$B, c("g1", "g2"))
})

test_that("lavaan labels alone establish the multigroup L2 grid", {
  hs <- lavaan::HolzingerSwineford1939[, c("x1", "x2", "x3", "school")]
  fit <- lavaan::cfa(
    "visual =~ x1 + x2 + x3", data = hs, group = "school", std.lv = TRUE
  )
  groups <- as.character(lavaan::lavInspect(fit, "group.label"))
  empty_l2 <- data.frame(
    group = character(), level = character(), factor = character(),
    metric = character(), value = numeric(), status = character(),
    reason = character(), estimand = character(), stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "not_evaluable", rows = NULL, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) empty_l2,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(fit, c("g1", "g2"))
  l2 <- out$rows[out$rows$level == "L2", , drop = FALSE]
  expect_identical(unique(l2$group), groups)
  expect_equal(nrow(l2), length(groups) * 2L * 2L)
  expect_named(out$second_order, groups)
  expect_named(out$second_order[[groups[[1L]]]], c("g1", "g2"))
  expect_named(out$second_order[[groups[[2L]]]], c("g1", "g2"))
  expect_identical(out$status, "not_evaluable")
})

test_that("an empty one-group L2 grid keeps the flat legacy alias", {
  l1_rows <- data.frame(
    group = "group1", level = "L1", factor = "f1",
    coefficient = "composite_reliability", value = .71,
    status = "ok", reason = NA_character_, stringsAsFactors = FALSE
  )
  empty_l2 <- data.frame(
    group = character(), level = character(), factor = character(),
    metric = character(), value = numeric(), status = character(),
    reason = character(), estimand = character(), stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .pm_standard_reliability_values = function(...) list(
      status = "ok", rows = l1_rows, values = NULL,
      compRelSEM = NULL, AVE = NULL, warnings = list(), messages = list(),
      condition = NULL
    ),
    .pm_second_order_matrix_rows = function(...) empty_l2,
    .package = "PsychoMatic"
  )

  out <- PsychoMatic:::.compute_second_order_reliability(NULL, c("g1", "g2"))
  expect_named(out$second_order, c("g1", "g2"))
  expect_s3_class(out$second_order$g1, "data.frame")
  expect_s3_class(out$second_order$g2, "data.frame")
  expect_identical(out$second_order$g1$metric, c("composite_reliability", "ave"))
  expect_true(all(out$second_order$g2$status == "not_evaluable"))
  expect_true(all(out$second_order$g2$reason == "second_order_group_not_evaluable"))
})
