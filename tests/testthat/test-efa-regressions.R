test_that("efa_auto refits after the final removal", {
  set.seed(701)
  f <- rnorm(250)
  d <- data.frame(
    q1 = f + rnorm(250, 0, .2),
    q2 = f + rnorm(250, 0, .2),
    q3 = f + rnorm(250, 0, .2),
    weak = rnorm(250)
  )

  out <- suppressMessages(suppressWarnings(efa_auto(
    d,
    min_loading = .70,
    max_iter = 1,
    parallel_iter = 10,
    verbose = FALSE,
    seed = 11
  )))

  expect_equal(rownames(out$cargas), colnames(out$final_data))
  expect_false(any(out$items_eliminados$Item %in% rownames(out$cargas)))
})

test_that("parallel count can retain zero factors", {
  expect_identical(
    PsychoMatic:::.pm_parallel_count(c(.40, .30), c(.50, .35)),
    0L
  )
  expect_identical(
    PsychoMatic:::.pm_parallel_count(c(.60, .30), c(.50, .35)),
    1L
  )
})

test_that("blockwise Mardia equals the literal small-sample fixture", {
  x <- matrix(
    c(1, 2, 2, 1, 3, 4, 4, 3, 5, 7, 7, 5, 8, 6),
    ncol = 2,
    byrow = TRUE
  )
  expected <- c(
    b1p = 1.186540390827097,
    skew_statistic = 2.636756424060216,
    skew_df = 4,
    skew_p = 0.620325885451584,
    b2p = 5.349113459690803,
    kurtosis_z = -0.876698317438317,
    kurtosis_p = 0.380650520910857
  )

  expect_equal(
    PsychoMatic:::.pm_mardia_blockwise(x, block_size = 2L),
    expected,
    tolerance = 1e-8
  )
})

test_that("Mardia caps requested blocks below n without changing the fixture", {
  x <- matrix(
    c(1, 2, 2, 1, 3, 4, 4, 3, 5, 7, 7, 5, 8, 6),
    ncol = 2,
    byrow = TRUE
  )
  expected <- c(
    b1p = 1.186540390827097,
    skew_statistic = 2.636756424060216,
    skew_df = 4,
    skew_p = 0.620325885451584,
    b2p = 5.349113459690803,
    kurtosis_z = -0.876698317438317,
    kurtosis_p = 0.380650520910857
  )
  block_ranges <- get0(
    ".pm_mardia_block_ranges",
    envir = asNamespace("PsychoMatic"),
    inherits = FALSE
  )

  expect_true(is.function(block_ranges))
  if (!is.function(block_ranges)) return(invisible(NULL))

  for (requested in c(7L, 100L)) {
    ranges <- block_ranges(nrow(x), requested)
    rows_per_block <- vapply(ranges, length, integer(1))

    expect_lte(max(rows_per_block), 6L)
    expect_equal(
      PsychoMatic:::.pm_mardia_blockwise(x, block_size = requested),
      expected,
      tolerance = 1e-8
    )
  }
})

test_that("singular Mardia covariance is a traceable condition", {
  singular <- cbind(q1 = 1:8, q2 = 1:8)
  condition <- tryCatch(
    PsychoMatic:::.pm_mardia_blockwise(singular),
    error = identity
  )

  expect_s3_class(condition, "psychomatic_mardia_singular")
  expect_identical(condition$code, "mardia_singular_covariance")
  expect_match(conditionMessage(condition), "singular")
})

test_that("suggest mode never deletes an item", {
  d <- simulate_likert_fixture(180L, 703L)
  out <- suppressMessages(suppressWarnings(efa_auto(
    d,
    refine = "suggest",
    max_iter = 2,
    parallel_iter = 10,
    verbose = FALSE,
    seed = 12
  )))

  expect_identical(colnames(out$final_data), colnames(d))
  expect_identical(out$final_data, out$datos_finales)
  expect_true(is.data.frame(out$screening_suggestions))
  expect_true(nrow(out$screening_suggestions) >= 0L)
})

test_that("polychoric fallback preserves the original condition and provenance", {
  d <- data.frame(
    q1 = rep(1:5, 8),
    q2 = rep(c(1, 2, 3, 5, 4), 8),
    q3 = rep(c(5, 4, 3, 2, 1), 8)
  )
  original <- structure(
    list(
      message = "forced sparse-cell failure",
      call = quote(psych::polychoric(d))
    ),
    class = c("forced_polychoric_error", "error", "condition")
  )
  fail_polychoric <- function(...) stop(original)

  out <- PsychoMatic:::.efa_estimate_correlations(
    d,
    method = "polychoric",
    repair = "none",
    polychoric_fun = fail_polychoric
  )

  expect_identical(out$requested_method, "polychoric")
  expect_identical(out$effective_method, "pearson")
  expect_identical(out$warnings[[1L]], original)
  expect_identical(conditionCall(out$warnings[[1L]]), quote(psych::polychoric(d)))
  expect_true(any(out$decision_log$code == "correlation_fallback"))
})

test_that("polychoric fallback is restricted to numeric ordered codes", {
  continuous <- data.frame(
    q1 = seq(0, 1, length.out = 30),
    q2 = seq(1, 3, length.out = 30)^2,
    q3 = sin(seq(0, 3, length.out = 30))
  )
  fail_polychoric <- function(...) stop("forced failure")

  expect_error(
    PsychoMatic:::.efa_estimate_correlations(
      continuous,
      method = "polychoric",
      repair = "none",
      polychoric_fun = fail_polychoric
    ),
    "numeric ordered codes"
  )
})

test_that("nearPD repair is opt-in and records the pre-repair eigenvalue", {
  d <- data.frame(
    q1 = 1:20,
    q2 = 1:20,
    q3 = c(1:10, 20:11)
  )

  raw <- PsychoMatic:::.efa_estimate_correlations(
    d,
    method = "pearson",
    repair = "none"
  )
  repaired <- PsychoMatic:::.efa_estimate_correlations(
    d,
    method = "pearson",
    repair = "nearPD"
  )

  expect_false(raw$repaired)
  expect_true(repaired$repaired)
  expect_true(is.finite(repaired$pre_repair_min_eigenvalue))
  expect_lte(repaired$pre_repair_min_eigenvalue, sqrt(.Machine$double.eps))
  expect_true(any(repaired$decision_log$code == "correlation_repair"))
  expect_gt(min(eigen(repaired$matrix, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("EFA exposes exact requested and effective parallel settings", {
  expect_identical(formals(efa_auto)$parallel_iter, 500L)

  set.seed(704)
  f <- rnorm(120)
  d <- data.frame(
    q1 = f + rnorm(120, sd = .3),
    q2 = f + rnorm(120, sd = .3),
    q3 = f + rnorm(120, sd = .3),
    q4 = f + rnorm(120, sd = .3)
  )

  for (mode in c("remove", "suggest")) {
    out <- suppressMessages(suppressWarnings(efa_auto(
      d,
      refine = mode,
      max_iter = 1,
      parallel_iter = 9L,
      parallel_quantile = .90,
      verbose = FALSE,
      seed = 705
    )))

    expect_identical(out$requested_settings$refine, mode)
    expect_identical(out$effective_settings$refine, mode)
    expect_identical(out$effective_settings$parallel_analysis$n_iter, 9L)
    expect_identical(out$effective_settings$parallel_analysis$quantile, .90)
    expect_identical(out$effective_settings$parallel_analysis$seed, 705L)
  }
})

test_that("parallel engine uses factor-eigenvalue quantiles and restores RNG", {
  set.seed(799)
  d <- matrix(rnorm(360), ncol = 4)
  R <- cor(d)
  before <- .Random.seed
  out <- suppressMessages(suppressWarnings(PsychoMatic:::.efa_parallel_analysis(
    R,
    n_obs = nrow(d),
    n_iter = 7,
    quantile = .90,
    seed = 19
  )))

  expect_identical(.Random.seed, before)
  expect_identical(out$target, "factor_eigenvalues")
  expect_equal(
    out$simulated,
    apply(out$simulated_eigenvalues, 2, quantile, probs = .90, names = FALSE)
  )
  expect_identical(
    out$n_factors,
    PsychoMatic:::.pm_parallel_count(out$observed, out$simulated)
  )
})

test_that("mixed numeric EFA records Pearson classification and a fixed sample", {
  set.seed(711)
  f <- rnorm(100)
  d <- data.frame(
    ordinal1 = as.integer(cut(f + rnorm(100, sd = .4), 5, labels = FALSE)),
    ordinal2 = as.integer(cut(f + rnorm(100, sd = .4), 5, labels = FALSE)),
    continuous1 = f + rnorm(100, sd = .4),
    continuous2 = f + rnorm(100, sd = .4)
  )
  d$continuous1[c(2, 9)] <- NA
  d$ordinal1[5] <- NA

  out <- suppressMessages(suppressWarnings(efa_auto(
    d,
    refine = "suggest",
    parallel_iter = 5,
    verbose = FALSE,
    seed = 51
  )))

  expect_identical(out$data_summary$measurement, "mixed_numeric")
  expect_identical(out$correlation$requested_method, "pearson")
  expect_identical(out$effective_settings$missing, "complete_cases")
  expect_true(any(out$decision_log$code == "missing_data"))
  expect_true(any(out$decision_log$code == "measurement_classification"))
  expect_true(all(out$correlation$pairwise_n == nrow(out$final_data)))
})

test_that("EFA resolves measurement settings for every current item set", {
  ordinal <- simulate_likert_fixture(120L, 719L)[1:3]
  mixed <- ordinal
  mixed$continuous <- seq_len(nrow(mixed)) / nrow(mixed)
  resolver <- get0(
    ".efa_resolve_current_settings",
    envir = asNamespace("PsychoMatic"),
    inherits = FALSE
  )

  expect_true(is.function(resolver))
  if (!is.function(resolver)) return(invisible(NULL))

  initial <- resolver(mixed, normality_ok = FALSE, language = "eng")
  final <- resolver(ordinal, normality_ok = FALSE, language = "eng")

  expect_identical(initial$measurement, "mixed_numeric")
  expect_identical(initial$requested_method, "pearson")
  expect_identical(final$measurement, "ordinal")
  expect_identical(final$requested_method, "polychoric")
  expect_identical(initial$estimator$fm, "minres")
  expect_identical(final$estimator$fm, "minres")
})

test_that("EFA adapts from mixed to ordinal after screening its only continuous item", {
  ordinal <- simulate_likert_fixture(180L, 720L)[1:3]
  mixed <- ordinal
  mixed$continuous <- seq_len(nrow(mixed)) / nrow(mixed)
  screen_calls <- 0L
  deterministic_screen <- function(...) {
    screen_calls <<- screen_calls + 1L
    if (screen_calls == 1L) {
      return(data.frame(
        Item = "continuous",
        Razon = "Deterministic transition fixture",
        Prioridad = Inf,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      Item = character(), Razon = character(), Prioridad = numeric(),
      stringsAsFactors = FALSE
    )
  }
  testthat::local_mocked_bindings(
    .efa_evaluar_items = deterministic_screen,
    .package = "PsychoMatic"
  )

  out <- suppressMessages(suppressWarnings(efa_auto(
    mixed,
    refine = "remove",
    max_iter = 1L,
    parallel_iter = 8L,
    verbose = FALSE,
    language = "eng",
    seed = 721L
  )))

  expect_identical(out$data_summary$measurement_initial, "mixed_numeric")
  expect_identical(out$data_summary$measurement_final, "ordinal")
  expect_identical(out$correlation$requested_method, "polychoric")
  expect_identical(out$requested_settings$correlation_policy, "adaptive_by_current_items")
  expect_identical(out$effective_settings$measurement, "ordinal")
  expect_identical(out$estimador$fm, "minres")
  expect_identical(out$parallel_analysis$extraction, "minres")
  expect_equal(out$efa$r, out$correlation$matrix, tolerance = 1e-12)
  expect_true(any(out$decision_log$code == "measurement_transition"))
  expect_identical(nrow(out$final_data), nrow(mixed))
})

test_that("EFA preserves diagnostic conditions but flattens them for presentation", {
  d <- data.frame(
    q1 = 1:40,
    q2 = 1:40,
    q3 = c(1:20, 40:21),
    q4 = rep(c(1, 2, 4, 3), 10)
  )
  out <- suppressMessages(suppressWarnings(efa_auto(
    d,
    refine = "suggest",
    repair = "nearPD",
    parallel_iter = 3,
    verbose = FALSE,
    seed = 54
  )))

  singular_row <- which(out$warnings$code == "mardia_singular_covariance")
  expect_length(singular_row, 1L)
  expect_s3_class(out$warnings$condition[[singular_row]], "psychomatic_mardia_singular")
  display <- summary(out)$warnings
  expect_true(all(vapply(display, is.atomic, logical(1))))
  expect_false("condition" %in% names(display))
  expect_true(all(c("code", "message", "class", "call") %in% names(display)))
})

test_that("screening language is descriptive and legacy alias metadata is explicit", {
  set.seed(706)
  f <- rnorm(120)
  d <- data.frame(
    q1 = f + rnorm(120, sd = .25),
    q2 = f + rnorm(120, sd = .25),
    q3 = f + rnorm(120, sd = .25),
    q4 = f + rnorm(120, sd = .25)
  )
  out <- suppressMessages(suppressWarnings(efa_auto(
    d,
    refine = "suggest",
    parallel_iter = 8,
    verbose = FALSE,
    language = "eng",
    seed = 707
  )))

  output <- capture.output(print(out))
  stored_text <- unlist(list(out$conclusion, out$decisions$rationale, out$decisions$decision))
  expect_false(any(grepl("optimal", c(output, stored_text), ignore.case = TRUE)))
  expect_identical(out$screening_passed, nrow(out$screening_suggestions) == 0L)
  expect_true(out$compatibility$aliases$solucion_optima$deprecated)
  expect_warning(
    legacy <- PsychoMatic:::.efa_legacy(out, "solucion_optima"),
    "deprecated"
  )
  expect_identical(legacy, out$solucion_optima)
})

test_that("EFA export fails closed on an existing path", {
  skip_if_not_installed("openxlsx")
  set.seed(712)
  f <- rnorm(60)
  d <- data.frame(
    q1 = f + rnorm(60, sd = .4), q2 = f + rnorm(60, sd = .4),
    q3 = f + rnorm(60, sd = .4), q4 = f + rnorm(60, sd = .4)
  )
  out <- suppressMessages(suppressWarnings(efa_auto(
    d, refine = "suggest", parallel_iter = 3, verbose = FALSE, seed = 53
  )))
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)

  suppressMessages(export_efa(out, format = "excel", file_name = path))
  expect_error(
    suppressMessages(export_efa(out, format = "excel", file_name = path)),
    "already exists"
  )
})

test_that("dimensionality sensitivity returns a transparent deterministic grid", {
  set.seed(708)
  f <- rnorm(90)
  d <- data.frame(
    q1 = f + rnorm(90, sd = .5),
    q2 = f + rnorm(90, sd = .5),
    q3 = f + rnorm(90, sd = .5),
    q4 = f + rnorm(90, sd = .5)
  )

  one <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = "auto",
    quantiles = c(.90, .95),
    n_iter = 5,
    seed = c(31, 32),
    language = "eng"
  )))
  two <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = "auto",
    quantiles = c(.90, .95),
    n_iter = 5,
    seed = c(31, 32),
    language = "eng"
  )))

  expect_s3_class(one, "psychomatic_dimensionality")
  expect_equal(one$configurations, two$configurations)
  expect_setequal(unique(one$configurations$requested_method), c("pearson", "spearman"))
  expect_equal(nrow(one$configurations), 8L)
  expect_true(all(c(
    "requested_method", "effective_method", "quantile", "seed", "n_iter",
    "n_factors", "status", "reason"
  ) %in% names(one$configurations)))
  expect_true(all(c("minimum", "maximum") %in% names(one$solution_range)))
  expect_true(is.data.frame(summary(one)$configurations))
  expect_true(is.data.frame(summary(one)$agreement))
  expect_false(any(grepl(
    "vote|winner|recommended solution",
    unlist(list(one$conclusion, one$agreement)),
    ignore.case = TRUE
  )))
})

test_that("dimensionality sensitivity uses one fixed complete-case sample", {
  set.seed(713)
  d <- as.data.frame(matrix(rnorm(320), ncol = 4))
  names(d) <- paste0("q", 1:4)
  d[1, 1] <- NA_real_
  d[2, 2] <- NA_real_
  d[3, 3] <- NA_real_
  complete <- d[stats::complete.cases(d), , drop = FALSE]

  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = "pearson",
    quantiles = .90,
    n_iter = 4,
    seed = 33,
    language = "eng"
  )))
  reference <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    complete,
    correlations = "pearson",
    quantiles = .90,
    n_iter = 4,
    seed = 33,
    language = "eng"
  )))

  expect_identical(out$data_summary$n_original, 80L)
  expect_identical(out$data_summary$n_complete, 77L)
  expect_identical(out$data_summary$n, 77L)
  expect_identical(out$effective_settings$missing, "complete_cases")
  expect_true(any(out$decision_log$code == "missing_data"))
  expect_equal(out$configurations, reference$configurations)
})

test_that("dimensionality sensitivity inherits specification language", {
  d <- simulate_likert_fixture(60L, 714L)[1:4]
  specification <- psychomatic_spec(
    d,
    items = names(d),
    language = "eng"
  )

  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    specification,
    correlations = "spearman",
    n_iter = 3,
    seed = 34
  )))

  expect_identical(out$language, "eng")
  expect_identical(out$effective_settings$language, "eng")
})

test_that("dimensionality sensitivity honors specification measurement levels", {
  d <- simulate_likert_fixture(60L, 717L)[1:4]
  specification <- psychomatic_spec(
    d,
    items = names(d),
    level = "continuous",
    language = "eng"
  )

  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    specification,
    correlations = "auto",
    n_iter = 3,
    seed = 37
  )))

  expect_identical(out$data_summary$measurement, "continuous")
  expect_setequal(
    unique(out$configurations$requested_method),
    c("pearson", "spearman")
  )
})

test_that("ordinal sensitivity evaluates ordered-factor Spearman cells", {
  ordinal <- simulate_likert_fixture(70L, 715L)[1:4]
  ordered_data <- as.data.frame(lapply(
    ordinal,
    function(x) ordered(x, levels = 1:5)
  ))

  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    ordered_data,
    correlations = "auto",
    n_iter = 3,
    seed = 35,
    language = "eng"
  )))
  spearman <- out$configurations[
    out$configurations$requested_method == "spearman",
    ,
    drop = FALSE
  ]

  expect_identical(nrow(spearman), 1L)
  expect_identical(spearman$status, "ok")
  expect_identical(spearman$effective_method, "spearman")
})

test_that("dimensionality sensitivity returns its effective matrices", {
  set.seed(716)
  d <- as.data.frame(matrix(rnorm(240), ncol = 4))
  names(d) <- paste0("q", 1:4)

  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = c("pearson", "spearman"),
    n_iter = 3,
    seed = 36,
    language = "eng"
  )))

  expect_setequal(names(out$effective_matrices), c("pearson", "spearman"))
  expect_true(all(vapply(out$effective_matrices, is.matrix, logical(1))))
  expect_equal(
    out$effective_matrices$pearson,
    stats::cor(d, method = "pearson"),
    tolerance = 1e-12
  )
  expect_identical(out$results$effective_matrices, out$effective_matrices)
})

test_that("dimensionality sensitivity reports methods effective after fallback", {
  original_estimator <- getFromNamespace(
    ".efa_estimate_correlations",
    "PsychoMatic"
  )
  fallback_estimator <- function(data, method, repair) {
    result <- original_estimator(data, method = "pearson", repair = repair)
    result$requested_method <- method
    result$effective_method <- "pearson"
    result
  }
  testthat::local_mocked_bindings(
    .efa_estimate_correlations = fallback_estimator,
    .package = "PsychoMatic"
  )
  d <- simulate_likert_fixture(60L, 718L)[1:4]

  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = "polychoric",
    n_iter = 3,
    seed = 38,
    language = "eng"
  )))

  expect_identical(out$configurations$effective_method, "pearson")
  expect_identical(
    unname(out$effective_settings$correlations[["polychoric"]]),
    "pearson"
  )
})

test_that("automatic dimensionality methods distinguish ordinal and binary items", {
  ordinal <- simulate_likert_fixture(80L, 709L)[1:5]
  ordinal_out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    ordinal,
    n_iter = 3,
    seed = 41,
    language = "eng"
  )))
  expect_setequal(
    unique(ordinal_out$configurations$requested_method),
    c("polychoric", "spearman")
  )

  binary <- as.data.frame(lapply(ordinal[1:4], function(x) as.integer(x >= 3)))
  binary_out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    binary,
    n_iter = 3,
    seed = 42,
    language = "eng"
  )))
  expect_setequal(
    unique(binary_out$configurations$requested_method),
    c("tetrachoric", "pearson")
  )
  expect_true(any(grepl("sensitivity contrast", binary_out$configurations$method_role)))
})

test_that("dimensionality methods print, summarize, and plot reusable data", {
  set.seed(710)
  d <- as.data.frame(matrix(rnorm(200), ncol = 4))
  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = "pearson",
    quantiles = .95,
    n_iter = 3,
    seed = 43,
    language = "eng"
  )))

  printed <- capture.output(print(out))
  expect_true(any(grepl("requested", printed, ignore.case = TRUE)))
  expect_true(any(grepl("range", printed, ignore.case = TRUE)))
  expect_true(any(grepl("quantiles.*0.95", printed, ignore.case = TRUE)))
  expect_true(any(grepl("seeds.*43", printed, ignore.case = TRUE)))
  expect_true(any(grepl("n_iter.*3", printed, ignore.case = TRUE)))

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  plotted <- plot(out)
  grDevices::dev.off()
  unlink(path)
  expect_equal(plotted, out$configurations[out$configurations$status == "ok", ])
})

test_that("invalid sensitivity configurations fail early", {
  d <- as.data.frame(matrix(rnorm(120), ncol = 4))
  expect_error(dimensionality_sensitivity(d, quantiles = 1), "quantiles")
  expect_error(dimensionality_sensitivity(d, n_iter = 0), "n_iter")
  expect_error(dimensionality_sensitivity(d, seed = NA), "seed")
  expect_error(dimensionality_sensitivity(d, report = NA), "report")
})

test_that("an entirely skipped sensitivity grid has an NA range", {
  d <- data.frame(q1 = 1:30, q2 = 1:30, q3 = 30:1)
  out <- suppressMessages(suppressWarnings(dimensionality_sensitivity(
    d,
    correlations = "pearson",
    n_iter = 3,
    seed = 52,
    language = "eng"
  )))

  expect_true(all(out$configurations$status == "skipped"))
  expect_true(all(is.na(out$solution_range)))
  expect_match(out$conclusion, "No configuration was evaluable")
})
