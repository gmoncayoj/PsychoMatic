test_that("model-tool public signatures and legacy classes remain stable", {
  expect_identical(
    head(names(formals(irt_auto)), 7L),
    c("data", "model", "factors", "itemtype", "technical", "language", "digits")
  )
  expect_identical(
    head(names(formals(compare_models_auto)), 5L),
    c("data", "models", "ordered", "estimator", "language")
  )
  expect_identical(
    names(formals(dif_auto))[1:6],
    c("data", "group", "items", "alpha", "language", "digits")
  )
  expect_identical(names(formals(dif_auto))[7], "polytomous")
  expect_identical(eval(formals(dif_auto)$polytomous), c("not_evaluable", "median"))
})

test_that("IRT automatic selection never silently assumes Rasch", {
  expect_identical(
    PsychoMatic:::.pm_auto_irt_model(c(i1 = 2L, i2 = 2L)),
    list(model = "2pl", itemtype = c(i1 = "2PL", i2 = "2PL"))
  )
  expect_identical(
    PsychoMatic:::.pm_auto_irt_model(c(i1 = 4L, i2 = 5L)),
    list(model = "grm", itemtype = c(i1 = "graded", i2 = "graded"))
  )
  expect_identical(
    PsychoMatic:::.pm_auto_irt_model(c(i1 = 2L, i2 = 4L)),
    list(model = "mixed", itemtype = c(i1 = "2PL", i2 = "graded"))
  )
})

test_that("IRT retains all signed slopes and a separate magnitude", {
  row <- c(a1 = -0.3, a2 = 0.4, a3 = 0, a4 = 0.5, d = 1)
  out <- PsychoMatic:::.pm_discrimination_components(row)
  expect_identical(out$dimensions, c(a1 = -0.3, a2 = 0.4, a3 = 0, a4 = 0.5))
  expect_equal(out$magnitude, sqrt(0.3^2 + 0.4^2 + 0.5^2), tolerance = 1e-15)

  interpreted <- data.frame(
    item = "i1", a1 = -0.3, a2 = 0.4, a3 = 0, a4 = 0.5,
    discrimination_magnitude = out$magnitude,
    stringsAsFactors = FALSE
  )
  expect_equal(
    PsychoMatic:::.pm_item_discrimination(interpreted),
    out$magnitude,
    tolerance = 1e-15
  )
})

test_that("IRT validates dimensions, items, and explicit item types before fitting", {
  binary <- data.frame(i1 = rep(0:1, 20), i2 = rep(1:0, 20))
  for (bad in list(NA_real_, Inf, 0, 1.5, c(1, 2))) {
    expect_error(irt_auto(binary, factors = bad), "factors")
  }
  expect_error(irt_auto(binary, technical = 1), "technical")
  expect_error(irt_auto(binary, digits = -1), "digits")
  expect_error(irt_auto(binary, model = "grm"), "incompatible")
  expect_error(
    irt_auto(binary, itemtype = c(i1 = "2PL", wrong = "2PL")),
    "names"
  )
  expect_error(irt_auto(binary, itemtype = "2PL"), "one value per item|exactly")

  duplicate <- binary
  names(duplicate) <- c("i", "i")
  expect_error(irt_auto(duplicate), "duplicate")
  expect_error(irt_auto(data.frame(i1 = rep(1, 40), i2 = rep(0:1, 20))), "constant")
  expect_error(irt_auto(data.frame(i1 = rep(NA_real_, 40), i2 = rep(0:1, 20))), "missing")
  expect_error(irt_auto(data.frame(i1 = c(0:1, rep(NaN, 38)), i2 = rep(0:1, 20))), "finite")
  expect_error(irt_auto(data.frame(i1 = rep(c(0, 1.5), 20), i2 = rep(0:1, 20))), "integer")
  expect_error(irt_auto(data.frame(i1 = factor(rep(c("b", "a"), 20)), i2 = rep(0:1, 20))), "unordered")
})

test_that("IRT fit failures preserve explicit NULL legacy slots and the original condition", {
  skip_if_not_installed("mirt")
  too_small <- data.frame(
    i1 = rep(0:1, 20),
    i2 = rep(1:0, 20)
  )

  got <- irt_auto(too_small)

  expect_true(all(c("fit", "coefficients", "raw_coefficients") %in% names(got)))
  expect_null(got[["fit"]])
  expect_null(got$fit)
  expect_null(got[["coefficients"]])
  expect_null(got[["raw_coefficients"]])
  expect_false(identical(got$fit, got$fit_statistics))
  expect_true(nrow(got$errors) >= 1L)
  fit_error <- which(vapply(got$errors$condition, inherits, logical(1), what = "error"))[1L]
  expect_true(is.finite(fit_error))
  expect_identical(
    conditionMessage(got$errors$condition[[fit_error]]),
    got$reason
  )
})

test_that("IRT real fits match public raw mirt coefficients and preserve raw precision", {
  skip_if_not_installed("mirt")
  set.seed(909)
  theta <- stats::rnorm(240)
  binary <- data.frame(
    i1 = stats::rbinom(240, 1, stats::plogis(0.8 * theta - 0.2)),
    i2 = stats::rbinom(240, 1, stats::plogis(1.1 * theta + 0.1)),
    i3 = stats::rbinom(240, 1, stats::plogis(0.7 * theta - 0.4)),
    i4 = stats::rbinom(240, 1, stats::plogis(1.3 * theta + 0.3))
  )

  set.seed(99)
  got3 <- irt_auto(binary, digits = 3, technical = list(NCYCLES = 300))
  set.seed(99)
  got7 <- irt_auto(binary, digits = 7, technical = list(NCYCLES = 300))

  expect_s3_class(got3, "psychomatic_irt")
  expect_s3_class(got3, "psychomatic_result")
  expect_identical(got3$model, "2pl")
  expect_true(got3$converged)
  expect_true(got3$admissible)
  expect_identical(got3$raw_coefficients, mirt::coef(got3$fit, IRTpars = FALSE, simplify = TRUE))
  expect_equal(got3$raw_item_parameters$a1,
               unname(mirt::coef(got3$fit, IRTpars = FALSE, simplify = TRUE)$items[, "a1"]),
               tolerance = 1e-12)
  expect_identical(got3$raw_item_parameters, got7$raw_item_parameters)
  expect_identical(got3$fit_statistics, got7$fit_statistics)
  expect_false(identical(got3$interpretation, got7$interpretation))
  if (identical(got3$fit_statistics$status, "available")) {
    expect_equal(got3$fit_statistics$M2, mirt::M2(got3$fit), tolerance = 1e-12)
  }
})

test_that("IRT solution assessment fails closed for nonconvergence and nonfinite parameters", {
  finite <- data.frame(item = "i1", a1 = 1, d = 0)
  nonfinite <- data.frame(item = "i1", a1 = Inf, d = 0)
  expect_false(PsychoMatic:::.pm_assess_irt_solution(FALSE, finite)$admissible)
  expect_match(
    PsychoMatic:::.pm_assess_irt_solution(FALSE, finite)$reason,
    "converge",
    ignore.case = TRUE
  )
  expect_false(PsychoMatic:::.pm_assess_irt_solution(TRUE, nonfinite)$admissible)
  expect_match(
    PsychoMatic:::.pm_assess_irt_solution(TRUE, nonfinite)$reason,
    "nonfinite",
    ignore.case = TRUE
  )

  mixed <- data.frame(
    item = c("binary", "graded"),
    a1 = c(0.8, 1.1),
    d = c(-0.2, NA_real_),
    g = c(0, NA_real_),
    u = c(1, NA_real_),
    d1 = c(NA_real_, 0.9),
    d2 = c(NA_real_, 0.1),
    d3 = c(NA_real_, -0.8),
    stringsAsFactors = FALSE
  )
  applicable <- PsychoMatic:::.pm_assess_irt_solution(
    TRUE,
    mixed,
    itemtype = c(binary = "2PL", graded = "graded"),
    n_categories = c(binary = 2L, graded = 4L)
  )
  expect_true(applicable$admissible)

  mixed$d2[2] <- Inf
  expect_false(PsychoMatic:::.pm_assess_irt_solution(
    TRUE,
    mixed,
    itemtype = c(binary = "2PL", graded = "graded"),
    n_categories = c(binary = 2L, graded = 4L)
  )$admissible)
})

test_that("IRT mixed fits treat structural NA parameters as inapplicable", {
  skip_if_not_installed("mirt")
  set.seed(904)
  n <- 600
  theta <- stats::rnorm(n)
  binary <- vapply(seq_len(4), function(j) {
    stats::rbinom(n, 1, stats::plogis((0.8 + 0.1 * j) * theta - 0.3 + 0.15 * j))
  }, numeric(n))
  graded <- vapply(seq_len(4), function(j) {
    as.integer(cut(
      (0.8 + 0.1 * j) * theta + stats::rlogis(n),
      breaks = c(-Inf, -0.8, 0, 0.8, Inf),
      labels = FALSE
    )) - 1L
  }, numeric(n))
  dat <- as.data.frame(cbind(binary, graded))
  names(dat) <- paste0("i", seq_len(8))

  got <- irt_auto(dat, technical = list(NCYCLES = 500))
  expect_true(got$converged)
  expect_true(got$admissible)
  expect_identical(got$conclusion$status, "evaluable")
  expect_identical(got$model, "mixed")
  expect_identical(unname(got$effective_itemtype), c(rep("2PL", 4), rep("graded", 4)))
  expect_equal(nrow(got$model_summary), 1L)
  expect_identical(got$model_summary$itemtype, "2PL, graded")
  expect_true(all(is.finite(got$raw_item_parameters$discrimination_magnitude)))
})

test_that("IRT summary decisions use unrounded discrimination evidence", {
  item_parameters <- data.frame(
    item = "i1",
    a1 = 0.6496,
    d = 0,
    discrimination_magnitude = 0.6496,
    stringsAsFactors = FALSE
  )
  interpretation_2 <- PsychoMatic:::.pm_irt_interpretation(
    item_parameters, model = "2pl", itemtype = "2PL", language = "eng", digits = 2
  )
  interpretation_4 <- PsychoMatic:::.pm_irt_interpretation(
    item_parameters, model = "2pl", itemtype = "2PL", language = "eng", digits = 4
  )
  summary_2 <- PsychoMatic:::.pm_irt_model_summary(
    interpretation_2, "2pl", "2PL", 100, "eng", 2,
    raw_discrimination = item_parameters$discrimination_magnitude
  )
  summary_4 <- PsychoMatic:::.pm_irt_model_summary(
    interpretation_4, "2pl", "2PL", 100, "eng", 4,
    raw_discrimination = item_parameters$discrimination_magnitude
  )
  expect_identical(interpretation_2$recommendation, interpretation_4$recommendation)
  expect_identical(summary_2$recommendation, summary_4$recommendation)
  expect_false(identical(summary_2$mean_discrimination, summary_4$mean_discrimination))
})

test_that("DIF excludes the focal item and freezes identical rows for every nested fit", {
  set.seed(912)
  n <- 360
  group <- rep(c("A", "B"), each = n / 2)
  theta <- stats::rnorm(n)
  dat <- data.frame(
    focal = stats::rbinom(n, 1, stats::plogis(theta + 0.45 * (group == "B"))),
    r1 = stats::rbinom(n, 1, stats::plogis(0.9 * theta)),
    r2 = stats::rbinom(n, 1, stats::plogis(1.1 * theta - 0.2)),
    r3 = stats::rbinom(n, 1, stats::plogis(0.8 * theta + 0.2))
  )
  dat$r2[c(2, 17, 201)] <- NA
  original <- dat

  got <- dif_auto(dat, group = group, items = names(dat), alpha = 0.05)
  row <- got$results[got$results$item == "focal", , drop = FALSE]
  score <- rowMeans(dat[c("r1", "r2", "r3")], na.rm = TRUE)
  score[!is.finite(score)] <- NA_real_
  keep <- stats::complete.cases(dat$focal, score, group)
  direct <- data.frame(y = dat$focal[keep], score = score[keep], group = factor(group[keep]))
  base <- stats::glm(y ~ score, direct, family = stats::binomial())
  uniform <- stats::glm(y ~ score + group, direct, family = stats::binomial())
  nonuniform <- stats::glm(y ~ score * group, direct, family = stats::binomial())
  lr_u <- stats::anova(base, uniform, test = "Chisq")
  lr_n <- stats::anova(uniform, nonuniform, test = "Chisq")

  expect_identical(row$status, "evaluable")
  expect_identical(row$original_indices[[1]], which(keep))
  expect_identical(row$base_indices[[1]], row$original_indices[[1]])
  expect_identical(row$uniform_indices[[1]], row$original_indices[[1]])
  expect_identical(row$nonuniform_indices[[1]], row$original_indices[[1]])
  expect_equal(row$p_uniform, lr_u$`Pr(>Chi)`[2], tolerance = 1e-10)
  expect_equal(row$p_nonuniform, lr_n$`Pr(>Chi)`[2], tolerance = 1e-10)
  expect_equal(row$lr_uniform, lr_u$Deviance[2], tolerance = 1e-10)
  expect_equal(row$lr_nonuniform, lr_n$Deviance[2], tolerance = 1e-10)
  expect_equal(row$df_uniform, lr_u$Df[2], tolerance = 1e-12)
  expect_equal(row$df_nonuniform, lr_n$Df[2], tolerance = 1e-12)
  expect_equal(row$coefficients[[1]], stats::coef(nonuniform), tolerance = 1e-10)
  expect_equal(row$odds_ratios[[1]], exp(stats::coef(nonuniform)), tolerance = 1e-10)
  direct_r2 <- c(
    base = PsychoMatic:::.pm_nagelkerke_r2(base),
    uniform = PsychoMatic:::.pm_nagelkerke_r2(uniform),
    nonuniform = PsychoMatic:::.pm_nagelkerke_r2(nonuniform)
  )
  expect_equal(row$nagelkerke_base, unname(direct_r2["base"]), tolerance = 1e-12)
  expect_equal(row$nagelkerke_uniform, unname(direct_r2["uniform"]), tolerance = 1e-12)
  expect_equal(row$nagelkerke_nonuniform, unname(direct_r2["nonuniform"]), tolerance = 1e-12)
  expect_equal(
    row$delta_r2_uniform,
    unname(direct_r2["uniform"] - direct_r2["base"]),
    tolerance = 1e-12
  )
  expect_equal(
    row$delta_r2_nonuniform,
    unname(direct_r2["nonuniform"] - direct_r2["uniform"]),
    tolerance = 1e-12
  )
  expect_equal(nrow(got$warnings), 0L)
  expect_identical(dat, original)
})

test_that("DIF distinguishes non-evaluable items from an evaluable no-DIF result", {
  group <- rep(c("A", "B"), each = 20)
  dat <- data.frame(
    constant = rep(0, 40),
    varying = rep(c(0, 1), 20),
    poly = rep(0:3, 10)
  )
  got <- dif_auto(dat, group = group, alpha = 0.05)

  expect_identical(got$results$status[got$results$item == "constant"], "not_evaluable")
  expect_identical(got$results$status[got$results$item == "poly"], "not_evaluable")
  expect_match(got$results$reason[got$results$item == "poly"], "ordinal|polytomous|IRT", ignore.case = TRUE)
  expect_gte(got$summary$not_evaluable_items, 2L)
  expect_false(grepl("no DIF", got$summary$recommendation, ignore.case = TRUE))

  median <- dif_auto(dat[c("poly", "varying")], group = group, polytomous = "median")
  expect_identical(median$results$method[median$results$item == "poly"], "median_dichotomized")
  expect_match(median$results$note[median$results$item == "poly"], "warning|advertencia", ignore.case = TRUE)
  expect_true(any(median$warnings$code == "dif_median_dichotomization"))
  expect_true(any(grepl("median|polytomous", median$warnings$message, ignore.case = TRUE)))
})

test_that("DIF matching scores exclude numeric metadata outside the requested item set", {
  set.seed(914)
  n <- 240
  group <- rep(c("A", "B"), each = n / 2)
  z <- stats::rnorm(n)
  items <- data.frame(
    i1 = stats::rbinom(n, 1, stats::plogis(z)),
    i2 = stats::rbinom(n, 1, stats::plogis(0.8 * z)),
    i3 = stats::rbinom(n, 1, stats::plogis(1.2 * z))
  )
  with_metadata <- items
  with_metadata$respondent_id <- seq_len(n) * 1000000
  with_metadata$age <- rep(18:77, length.out = n)

  reference <- dif_auto(items, group, items = names(items))
  got <- dif_auto(with_metadata, group, items = names(items))
  expect_identical(got$results, reference$results)
  expect_identical(
    lapply(got$diagnostics$items, `[[`, "score"),
    lapply(reference$diagnostics$items, `[[`, "score")
  )
})

test_that("DIF fails closed for unusable matching scores and operationally small groups", {
  group <- rep(c("A", "B"), each = 20)
  one <- dif_auto(data.frame(i1 = rep(c(0, 1), 20)), group = group)
  expect_identical(one$results$status, "not_evaluable")
  expect_match(one$results$reason, "matching|rest", ignore.case = TRUE)

  small_group <- c(rep("A", 35), rep("B", 5))
  dat <- data.frame(i1 = rep(c(0, 1), 20), i2 = rep(c(1, 0), 20))
  small <- dif_auto(dat, group = small_group)
  expect_true(all(small$results$status == "not_evaluable"))
  expect_true(all(grepl("group", small$results$reason, ignore.case = TRUE)))

  expect_error(dif_auto(dat, group = rep("A", 40)), "two groups")
  expect_error(dif_auto(dat, group = group, items = c("i1", "i1")), "unique")
  expect_error(dif_auto(dat, group = group, alpha = 1), "alpha")
  expect_error(dif_auto(dat, group = group, digits = 1.5), "digits")
})

test_that("DIF fitting diagnostics reject separation instead of declaring no DIF", {
  group <- rep(c("A", "B"), each = 30)
  dat <- data.frame(
    focal = rep(c(0, 1), each = 30),
    r1 = rep(c(0, 1), each = 30),
    r2 = rep(c(0, 1), each = 30)
  )
  got <- dif_auto(dat, group = group, items = names(dat))
  focal <- got$results[got$results$item == "focal", , drop = FALSE]
  expect_identical(focal$status, "not_evaluable")
  expect_match(focal$reason, "separation|degenerate|finite", ignore.case = TRUE)
  expect_true(is.na(focal$flagged))
})

test_that("DIF raw evidence does not depend on presentation digits", {
  set.seed(919)
  n <- 300
  g <- rep(c("A", "B"), each = n / 2)
  z <- stats::rnorm(n)
  dat <- data.frame(
    i1 = stats::rbinom(n, 1, stats::plogis(z)),
    i2 = stats::rbinom(n, 1, stats::plogis(0.7 * z)),
    i3 = stats::rbinom(n, 1, stats::plogis(1.2 * z))
  )
  a <- dif_auto(dat, g, digits = 2)
  b <- dif_auto(dat, g, digits = 7)
  expect_identical(a$results, b$results)
  expect_false(identical(a$interpretation, b$interpretation))
})

test_that("model comparison preserves duplicate labels under stable positional IDs", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  models <- list(
    duplicate = "f =~ x1 + x2 + x3 + x4 + x5 + x6",
    duplicate = "f1 =~ x1 + x2 + x3\nf2 =~ x4 + x5 + x6"
  )
  got <- compare_models_auto(hs, models, estimator = "ML", language = "eng")
  expect_identical(got$fit$model_id, c("model_001", "model_002"))
  expect_identical(got$fit$model, c("duplicate", "duplicate"))
  expect_identical(names(got$fitted_models), c("model_001", "model_002"))
  expect_s3_class(got, "psychomatic_model_comparison")
  expect_s3_class(got, "psychomatic_result")
})

test_that("model comparison preserves NULL fit slots by positional model ID", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  models <- list(
    invalid_fit = "f =~ absent1 + absent2 + absent3",
    valid_fit = "f =~ x1 + x2 + x3 + x4 + x5 + x6"
  )
  got <- compare_models_auto(hs, models, estimator = "ML", language = "eng")
  expect_identical(names(got$fitted_models), c("model_001", "model_002"))
  expect_length(got$fitted_models, 2L)
  expect_null(got$fitted_models[["model_001"]])
  expect_s4_class(got$fitted_models[["model_002"]], "lavaan")
})

test_that("model comparison rejects candidates using different variables or actual rows", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  vars <- compare_models_auto(
    hs,
    list(a = "f =~ x1 + x2 + x3", b = "f =~ x4 + x5 + x6"),
    estimator = "ML",
    language = "eng"
  )
  expect_identical(vars$recommendation$basis, "none")
  expect_true(is.na(vars$recommendation$model))
  expect_true(all(!vars$fit$comparable))
  expect_true(all(c(
    "information_criteria_status", "information_criteria_reason"
  ) %in% names(vars$fit)))
  expect_true(all(is.na(vars$fit$aic)))
  expect_true(all(is.na(vars$fit$bic)))
  expect_true(all(vars$fit$information_criteria_status == "not_comparable"))
  expect_true(all(nzchar(vars$fit$information_criteria_reason)))

  hs$x1[1:12] <- NA
  rows <- compare_models_auto(
    hs,
    list(a = "f =~ x1 + x2 + x3", b = "f =~ x2 + x3 + x4"),
    estimator = "ML",
    language = "eng"
  )
  expect_identical(rows$recommendation$basis, "none")
  expect_true(is.na(rows$recommendation$model_id))
  expect_false(identical(rows$diagnostics$models[[1]]$case_indices,
                         rows$diagnostics$models[[2]]$case_indices))
  expect_true(all(c(
    "information_criteria_status", "information_criteria_reason"
  ) %in% names(rows$fit)))
  expect_true(all(is.na(rows$fit$aic)))
  expect_true(all(is.na(rows$fit$bic)))
  expect_true(all(rows$fit$information_criteria_status == "not_comparable"))
  expect_true(all(nzchar(rows$fit$information_criteria_reason)))
})

test_that("model comparison retains a valid compatible subset and ignores variable order", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  got <- compare_models_auto(
    hs,
    list(
      first = "f =~ x1 + x2 + x3",
      reordered = "f =~ x3 + x2 + x1",
      other = "f =~ x4 + x5 + x6"
    ),
    estimator = "ML",
    language = "eng"
  )
  expect_identical(got$fit$comparable, c(TRUE, TRUE, FALSE))
  expect_identical(
    got$fit$information_criteria_status,
    c("comparable", "comparable", "not_comparable")
  )
  expect_true(is.na(got$fit$aic[3L]))
  expect_true(is.na(got$fit$bic[3L]))
  expect_identical(
    sort(got$diagnostics$models[[1]]$observed_variables),
    sort(got$diagnostics$models[[2]]$observed_variables)
  )
  expect_identical(got$recommendation$basis, "none")
  expect_match(got$recommendation$message, "tie", ignore.case = TRUE)
})

test_that("model comparison distinguishes requested and effective ordered settings", {
  skip_if_not_installed("lavaan")
  set.seed(921)
  dat <- as.data.frame(replicate(6, sample(0:3, 240, replace = TRUE)))
  names(dat) <- paste0("i", 1:6)
  models <- list(
    one = "f =~ i1 + i2 + i3 + i4 + i5 + i6",
    two = "f1 =~ i1 + i2 + i3\nf2 =~ i4 + i5 + i6"
  )
  expect_no_warning(
    got <- compare_models_auto(dat, models, estimator = "WLSMV", language = "eng")
  )
  expect_null(got$requested_settings$ordered)
  expect_identical(got$effective_settings$ordered, names(dat))
  diagnostic_original <- unique(vapply(
    got$diagnostics$models,
    `[[`,
    character(1),
    "estimator_original"
  ))
  diagnostic_computational <- unique(vapply(
    got$diagnostics$models,
    `[[`,
    character(1),
    "estimator_effective"
  ))
  expect_identical(got$effective_settings$estimator_original, unname(diagnostic_original))
  expect_identical(got$effective_settings$estimator_computational, unname(diagnostic_computational))
  expect_identical(got$effective_settings$estimator, got$effective_settings$estimator_original)
  expect_true(all(is.na(got$fit$aic)))
  expect_true(all(is.na(got$fit$bic)))
  expect_true(nrow(got$warnings) > 0L)
  expect_true(any(grepl("variance", got$warnings$message, ignore.case = TRUE)))
})

test_that("ULS remains a continuous estimator unless ordered variables are explicit", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  models <- list(
    one = "f =~ x1 + x2 + x3 + x4 + x5 + x6",
    two = "f1 =~ x1 + x2 + x3\nf2 =~ x4 + x5 + x6"
  )
  got <- compare_models_auto(hs, models, estimator = "ULS", language = "eng")
  expect_true("ordered" %in% names(got))
  expect_null(got[["ordered", exact = TRUE]])
  expect_null(got$requested_settings$ordered)
  expect_identical(got$effective_settings$ordered, character())
  expect_true(all(vapply(
    got$diagnostics$models,
    function(x) identical(x$ordered_variables, character()),
    logical(1)
  )))
  expect_true(all(vapply(
    got$diagnostics$models,
    function(x) identical(x$estimand, "covariances"),
    logical(1)
  )))
  expect_true(all(!vapply(got$diagnostics$models, `[[`, logical(1), "meanstructure")))
  expect_identical(got$effective_settings$estimator_original, "ULS")
  expect_identical(got$effective_settings$estimator_computational, "ULS")
})

test_that("model comparison does not rank distinct effective estimands", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  covariance_syntax <- "f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
  got <- compare_models_auto(
    hs,
    list(
      covariance_only = covariance_syntax,
      with_means = paste(covariance_syntax, "x1 ~ 1", sep = "\n")
    ),
    estimator = "ML",
    language = "eng"
  )

  expect_identical(
    vapply(got$diagnostics$models, `[[`, logical(1), "meanstructure"),
    c(model_001 = FALSE, model_002 = TRUE)
  )
  expect_identical(
    vapply(got$diagnostics$models, `[[`, character(1), "estimand"),
    c(model_001 = "covariances", model_002 = "means_covariances")
  )
  expect_true(all(!got$fit$comparable))
  expect_identical(got$recommendation$basis, "none")
  expect_true(is.na(got$recommendation$model_id))
  expect_true(all(is.na(got$fit$aic)))
  expect_true(all(is.na(got$fit$bic)))
  expect_true(all(got$fit$information_criteria_status == "not_comparable"))
  expect_true(all(nzchar(got$fit$information_criteria_reason)))
})

test_that("common CFA index variants keep labels paired to exact values", {
  fit <- data.frame(
    model_id = c("model_001", "model_002"),
    cfi = c(0.90, 0.91), cfi_scaled = c(0.92, 0.93), cfi_robust = c(0.94, 0.95),
    tli = c(0.89, 0.90), tli_scaled = c(0.91, 0.92), tli_robust = c(NA, 0.94),
    rmsea = c(0.08, 0.07), rmsea_scaled = c(0.07, 0.06), rmsea_robust = c(0.06, 0.05),
    srmr = c(0.07, 0.06),
    stringsAsFactors = FALSE
  )
  out <- PsychoMatic:::.pm_select_common_fit_variants(fit)
  expect_identical(out$variants[c("cfi", "tli", "rmsea", "srmr")],
                   c(cfi = "scaled", tli = "scaled", rmsea = "scaled", srmr = "base"))
  expect_identical(out$values$cfi, c(0.92, 0.93))
  expect_identical(out$values$tli, c(0.91, 0.92))
  expect_identical(out$values$rmsea, c(0.07, 0.06))
  expect_identical(out$values$srmr, c(0.07, 0.06))
})

test_that("real CFA evidence matches public lavaan metadata and fit measures", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
  models <- list(
    one = "g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9",
    three = paste("v =~ x1 + x2 + x3", "t =~ x4 + x5 + x6", "s =~ x7 + x8 + x9", sep = "\n")
  )
  got <- compare_models_auto(hs, models, estimator = "ML", language = "eng")
  direct <- lavaan::cfa(models[[1]], hs, estimator = "ML")
  expect_identical(got$diagnostics$models[[1]]$observed_variables, lavaan::lavNames(direct, "ov"))
  expect_identical(got$diagnostics$models[[1]]$case_indices,
                   sort(unique(as.integer(unlist(lavaan::lavInspect(direct, "case.idx"))))))
  expect_equal(got$fit$cfi[1], unname(lavaan::fitMeasures(direct, "cfi")), tolerance = 1e-10)
  expect_equal(got$fit$rmsea[1], unname(lavaan::fitMeasures(direct, "rmsea")), tolerance = 1e-10)
})

test_that("model recommendation fails closed for one model, all NA, errors, and exact ties", {
  skip_if_not_installed("lavaan")
  hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:6)]
  one <- compare_models_auto(
    hs,
    list(only = "f =~ x1 + x2 + x3 + x4 + x5 + x6"),
    estimator = "ML",
    language = "eng"
  )
  expect_identical(one$recommendation$basis, "none")
  expect_true(is.na(one$recommendation$model))

  tied <- compare_models_auto(
    hs,
    list(a = "f =~ x1 + x2 + x3 + x4 + x5 + x6",
         b = "f =~ x1 + x2 + x3 + x4 + x5 + x6"),
    estimator = "ML",
    language = "eng"
  )
  expect_identical(tied$recommendation$basis, "none")
  expect_match(tied$recommendation$message, "tie|empate", ignore.case = TRUE)

  bad <- compare_models_auto(
    hs,
    list(a = "not valid lavaan syntax", b = "also ??? invalid"),
    estimator = "ML",
    language = "eng"
  )
  expect_true(all(!bad$fit$admissible))
  expect_true(all(nzchar(bad$fit$reason)))
  expect_identical(bad$recommendation$basis, "none")

  synthetic <- data.frame(
    model_id = c("model_001", "model_002"), model = c("a", "b"),
    converged = TRUE, admissible = TRUE, comparable = TRUE,
    comparison_cfi = NA_real_, comparison_tli = NA_real_,
    comparison_rmsea = NA_real_, comparison_srmr = NA_real_,
    fit_score = NA_real_, stringsAsFactors = FALSE
  )
  none <- PsychoMatic:::.pm_model_recommendation(synthetic, "eng")
  expect_identical(none$basis, "none")
  expect_true(is.na(none$model_id))
})
