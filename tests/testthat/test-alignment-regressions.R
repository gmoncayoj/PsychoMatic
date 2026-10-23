alignment_fixture <- function() {
  lambda <- matrix(
    c(0.70, 0.80, 0.75, 0.72, 0.78, 0.74),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("B", "A"), c("i2", "i1", "i3"))
  )
  nu <- matrix(
    c(0.00, 0.00, 0.00, 0.10, 0.05, 0.08),
    nrow = 2L,
    byrow = TRUE,
    dimnames = dimnames(lambda)
  )
  list(lambda = lambda, nu = nu)
}

alignment_real_outputs <- function() {
  x <- alignment_fixture()
  fit <- suppressWarnings(sirt::invariance.alignment(
    lambda = x$lambda,
    nu = x$nu,
    wgt = NULL,
    align.scale = c(lambda = 0.40, nu = 0.20),
    align.pow = c(lambda = 0.25, nu = 0.25),
    eps = 1e-3,
    center = FALSE,
    optimizer = "optim",
    fixed = NULL,
    meth = 1,
    vcov = NULL
  ))
  constraints <- sirt::invariance_alignment_constraints(
    model = fit,
    lambda_parm_tol = 0.40,
    nu_parm_tol = 0.20
  )
  list(fit = fit, constraints = constraints)
}

test_that("inv_align_auto preserves its public and legacy contracts", {
  skip_if_not_installed("sirt")

  expect_identical(
    head(names(formals(inv_align_auto)), 21L),
    c(
      "lambda", "nu", "data", "group", "sampling_weights", "wgt",
      "language", "model", "align.scale", "align.pow", "parm_tol",
      "noninvariance_cutoff", "center", "optimizer", "fixed", "meth",
      "eps", "vcov", "digits", "config_args", "alignment_args"
    )
  )
  defaults <- formals(inv_align_auto)
  expect_identical(defaults$lambda, NULL)
  expect_identical(defaults$nu, NULL)
  expect_identical(deparse(defaults$language), "c(\"eng\", \"esp\")")
  expect_identical(deparse(defaults$model), "c(\"2PM\", \"1PM\")")
  expect_identical(eval(defaults$align.scale), c(lambda = 0.40, nu = 0.20))
  expect_identical(eval(defaults$align.pow), c(lambda = 0.25, nu = 0.25))

  x <- alignment_fixture()
  positional <- unname(list(
    x$lambda, x$nu, NULL, NULL, NULL, NULL, "eng", "2PM",
    c(lambda = 0.40, nu = 0.20), c(lambda = 0.25, nu = 0.25), NULL,
    25, FALSE, "optim", NULL, 1, 1e-3, NULL, 3, list(), list()
  ))
  out <- suppressWarnings(do.call(inv_align_auto, positional))

  expect_identical(class(out)[1:2], c("inv_align_auto_result", "psychomatic_result"))
  expect_true(all(c(
    "call", "input_mode", "settings", "dimensions", "configural",
    "alignment", "constraints", "summary", "item_results",
    "group_parameters", "aligned_item_parameters", "noninvariance",
    "references", "report"
  ) %in% names(out)))
  expect_s3_class(summary(out), "data.frame")
  expect_identical(nrow(summary(out)), 2L)
  expect_null(out$requested_settings$fixed)
  expect_null(out$effective_settings$fixed_engine_request)
  expect_true(out$effective_settings$fixed_effective)
  expect_identical(out$effective_settings$fixed, out$alignment$fixed)
  printed <- NULL
  capture.output(printed <- print(out))
  expect_identical(printed, out)
})

test_that("named matrices are matched by identity and agree with direct sirt", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  nu_scrambled <- x$nu[c("A", "B"), c("i3", "i2", "i1"), drop = FALSE]
  wgt <- c(A = 3, B = 1)

  before_lambda <- serialize(x$lambda, NULL)
  before_nu <- serialize(nu_scrambled, NULL)
  before_wgt <- serialize(wgt, NULL)

  out <- expect_silent(inv_align_auto(
    lambda = x$lambda,
    nu = nu_scrambled,
    wgt = wgt,
    language = "eng", verbose = FALSE
  ))
  expected_wgt <- matrix(
    rep(c(1, 3), ncol(x$lambda)),
    nrow = nrow(x$lambda),
    dimnames = dimnames(x$lambda)
  )
  direct <- suppressWarnings(sirt::invariance.alignment(
    lambda = x$lambda,
    nu = x$nu,
    wgt = expected_wgt,
    align.scale = c(lambda = 0.40, nu = 0.20),
    align.pow = c(lambda = 0.25, nu = 0.25),
    eps = 1e-3,
    center = FALSE,
    optimizer = "optim",
    fixed = NULL,
    meth = 1,
    vcov = NULL
  ))

  expect_identical(out$effective_settings$lambda, x$lambda)
  expect_identical(out$effective_settings$nu, x$nu)
  expect_identical(out$effective_settings$wgt, expected_wgt)
  expect_equal(out$group_parameters, direct$pars, tolerance = 1e-10)
  expect_equal(out$aligned_item_parameters, direct$itempars.aligned, tolerance = 1e-10)
  expect_identical(serialize(x$lambda, NULL), before_lambda)
  expect_identical(serialize(nu_scrambled, NULL), before_nu)
  expect_identical(serialize(wgt, NULL), before_wgt)
})

test_that("matrix identities fail early when incomplete or ambiguous", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  bad <- list()

  z <- x$lambda
  rownames(z)[2] <- rownames(z)[1]
  bad$duplicate <- list(z, x$nu)
  z <- x$lambda
  rownames(z)[2] <- NA_character_
  bad$missing <- list(z, x$nu)
  z <- x$lambda
  colnames(z)[2] <- ""
  bad$empty <- list(z, x$nu)
  z <- x$lambda
  colnames(z) <- NULL
  bad$partial <- list(z, x$nu)
  z <- unname(x$nu)
  dim(z) <- dim(x$nu)
  bad$one_sided <- list(x$lambda, z)
  z <- x$nu
  rownames(z)[2] <- "C"
  bad$mismatch <- list(x$lambda, z)

  for (case in bad) {
    expect_error(
      inv_align_auto(lambda = case[[1L]], nu = case[[2L]]),
      "name|identif",
      ignore.case = TRUE
    )
  }
})

test_that("wholly unnamed matrices receive stable recorded identities", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  lambda <- unname(x$lambda)
  dim(lambda) <- dim(x$lambda)
  nu <- unname(x$nu)
  dim(nu) <- dim(x$nu)

  out <- suppressWarnings(inv_align_auto(lambda = lambda, nu = nu))

  expect_identical(out$dimensions$group_names, c("G1", "G2"))
  expect_identical(out$dimensions$item_names, c("item1", "item2", "item3"))
  expect_identical(out$dimensions$names_generated, TRUE)
  expect_identical(out$effective_settings$identity_matching, "generated")
})

test_that("alignment weights are validated and matched by identity", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  wm <- matrix(
    c(6, 3, 5, 2, 4, 1),
    nrow = 2L,
    dimnames = list(c("A", "B"), c("i3", "i2", "i1"))
  )
  out <- suppressWarnings(inv_align_auto(lambda = x$lambda, nu = x$nu, wgt = wm))
  expect_identical(
    out$effective_settings$wgt,
    wm[c("B", "A"), c("i2", "i1", "i3"), drop = FALSE]
  )
  expect_identical(out$effective_settings$weight_matching, "named")

  bad <- list(
    1,
    c(B = 1, C = 2),
    c(B = 1, A = Inf),
    c(B = 1, A = -1),
    c(B = 0, A = 0),
    matrix(1, 3, 3),
    `rownames<-`(matrix(1, 2, 3), c("B", "A"))
  )
  for (value in bad) {
    expect_error(
      inv_align_auto(lambda = x$lambda, nu = x$nu, wgt = value),
      "wgt|weight|name|dimension|zero|finite|negative",
      ignore.case = TRUE
    )
  }
})

test_that("data-mode sampling weights fail before the configural engine", {
  skip_if_not_installed("sirt")
  dat <- data.frame(i1 = 1:8, i2 = 2:9)
  grp <- rep(c("A", "B"), each = 4)
  bad <- list(
    1:7,
    c(rep(1, 7), Inf),
    c(rep(1, 7), -1),
    rep(0, 8)
  )
  for (value in bad) {
    expect_error(
      inv_align_auto(data = dat, group = grp, sampling_weights = value),
      "sampling_weights|finite|negative|zero",
      ignore.case = TRUE
    )
  }
})

test_that("entirely missing group-item cells fail before the configural engine", {
  skip_if_not_installed("sirt")
  dat <- data.frame(
    i1 = c(rep(NA_real_, 4), 1:4),
    i2 = 1:8
  )
  grp <- rep(c("A", "B"), each = 4)
  called <- FALSE
  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) {
      called <<- TRUE
      stop("configural engine must not run")
    },
    .package = "PsychoMatic"
  )

  expect_error(
    inv_align_auto(data = dat, group = grp),
    "group-item|entirely missing|completely missing",
    ignore.case = TRUE
  )
  expect_false(called)
})

test_that("data-mode records exact sampling weights and benign config options", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()
  dat <- data.frame(
    i2 = c(1, 3, 2, 4, 6, 5, 7, 8),
    i1 = c(2, 1, 4, 3, 5, 7, 6, 9),
    i3 = c(3, 2, 1, 5, 4, 6, 8, 7)
  )
  grp <- rep(c("B", "A"), each = 4)
  sampling <- c(1, 2, 1, 2, 3, 2, 3, 2)
  config_extra <- list(estimator = "ML", verbose = TRUE)

  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) {
      list(lambda = x$lambda, nu = x$nu, vcov = NULL)
    },
    .pm_alignment_engine = function(args) real$fit,
    .pm_alignment_constraints_engine = function(args) real$constraints,
    .package = "PsychoMatic"
  )
  out <- inv_align_auto(
    data = dat,
    group = grp,
    sampling_weights = sampling,
    config_args = config_extra,
    language = "eng"
  )

  expect_identical(out$conclusion_status, "interpretable")
  expect_identical(out$requested_settings$sampling_weights, sampling)
  expect_identical(out$effective_settings$sampling_weights, sampling)
  expect_identical(out$requested_settings$config_args, config_extra)
  expect_identical(out$effective_settings$config_args, config_extra)
  expect_identical(out$engine_calls$configural$arguments$weights, sampling)
  expect_identical(out$engine_calls$configural$arguments$estimator, "ML")
  expect_identical(out$engine_calls$configural$arguments$verbose, TRUE)
  expect_identical(
    names(out$engine_calls$configural$arguments),
    c("dat", "group", "weights", "model", "verbose", "estimator")
  )
})

test_that("requested settings retain raw inputs and matrix-only settings are honest", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  raw_scale <- c(nu = 0.20, lambda = 0.40)
  raw_power <- c(nu = 0.25, lambda = 0.30)
  raw_sampling <- c(4, 2, 1)
  raw_config <- list(estimator = "ML")

  out <- suppressWarnings(inv_align_auto(
    lambda = x$lambda,
    nu = x$nu,
    sampling_weights = raw_sampling,
    language = "en",
    model = "1PM",
    align.scale = raw_scale,
    align.pow = raw_power,
    parm_tol = NULL,
    config_args = raw_config
  ))

  expect_identical(out$requested_settings$lambda, x$lambda)
  expect_identical(out$requested_settings$nu, x$nu)
  expect_identical(out$requested_settings$language, "en")
  expect_identical(out$requested_settings$model, "1PM")
  expect_identical(out$requested_settings$align.scale, raw_scale)
  expect_identical(out$requested_settings$align.pow, raw_power)
  expect_null(out$requested_settings$parm_tol)
  expect_identical(out$requested_settings$sampling_weights, raw_sampling)
  expect_identical(out$requested_settings$config_args, raw_config)

  expect_identical(out$effective_settings$language, "eng")
  expect_identical(out$effective_settings$align.scale, c(lambda = 0.40, nu = 0.20))
  expect_identical(out$effective_settings$align.pow, c(lambda = 0.30, nu = 0.25))
  expect_identical(out$effective_settings$parm_tol, c(lambda = 0.40, nu = 0.20))
  expect_null(out$effective_settings$model)
  expect_null(out$effective_settings$sampling_weights)
  expect_null(out$effective_settings$config_args)
  expect_identical(out$effective_settings$configural_status, "not_applicable")
  expect_identical(out$effective_settings$sampling_weights_status, "not_applicable")
  expect_identical(out$effective_settings$config_args_status, "not_applicable")
})

test_that("data-mode suppresses text printed by the configural engine", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()
  dat <- data.frame(
    i1 = c(2, 1, 4, 3, 5, 7, 6, 9),
    i2 = c(1, 3, 2, 4, 6, 5, 7, 8),
    i3 = c(3, 2, 1, 5, 4, 6, 8, 7)
  )
  grp <- rep(c("B", "A"), each = 4)

  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) {
      cat("Compute CFA for group B\n")
      list(lambda = x$lambda, nu = x$nu, vcov = NULL)
    },
    .pm_alignment_engine = function(args) real$fit,
    .pm_alignment_constraints_engine = function(args) real$constraints,
    .package = "PsychoMatic"
  )

  out <- expect_silent(inv_align_auto(data = dat, group = grp, language = "eng", verbose = FALSE))
  expect_identical(out$conclusion_status, "interpretable")
})

test_that("stdout from every alignment engine is bounded and never escapes", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()
  dat <- data.frame(
    i1 = c(2, 1, 4, 3, 5, 7, 6, 9),
    i2 = c(1, 3, 2, 4, 6, 5, 7, 8),
    i3 = c(3, 2, 1, 5, 4, 6, 8, 7)
  )
  grp <- rep(c("B", "A"), each = 4)

  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) {
      cat("configural trace\n")
      list(lambda = x$lambda, nu = x$nu, vcov = NULL)
    },
    .pm_alignment_engine = function(args) {
      cat(paste0("alignment trace ", seq_len(80L), collapse = "\n"), "\n")
      real$fit
    },
    .pm_alignment_constraints_engine = function(args) {
      cat("constraints trace\n")
      real$constraints
    },
    .package = "PsychoMatic"
  )

  out <- expect_silent(inv_align_auto(data = dat, group = grp, language = "eng", verbose = FALSE))
  expect_identical(out$conclusion_status, "interpretable")
  expect_match(out$diagnostics$engine_stdout$configural, "configural trace")
  expect_match(out$diagnostics$engine_stdout$constraints, "constraints trace")
  expect_lte(length(out$diagnostics$engine_stdout$alignment), 50L)
  expect_lte(sum(nchar(out$diagnostics$engine_stdout$alignment)), 4000L)
  expect_true(out$diagnostics$engine_stdout_truncated$alignment)
})

test_that("data mode accepts NA, rejects NaN and Inf, and matches public sirt", {
  skip_if_not_installed("sirt")
  set.seed(804)
  n <- 240L
  grp <- rep(c("A", "B"), each = n / 2L)
  latent <- stats::rnorm(n)
  dat <- data.frame(
    i1 = 0.80 * latent + stats::rnorm(n, 0, 0.50),
    i2 = 0.75 * latent + stats::rnorm(n, 0, 0.55),
    i3 = 0.70 * latent + stats::rnorm(n, 0, 0.60),
    i4 = 0.65 * latent + stats::rnorm(n, 0, 0.65)
  )
  dat$i1[c(3, 125)] <- NA_real_
  dat$i2[c(8, 140, 201)] <- NA_real_
  dat$i3[17] <- NA_real_

  invisible(utils::capture.output(
    direct_config <- sirt::invariance_alignment_cfa_config(
      dat = dat,
      group = grp,
      weights = NULL,
      model = "2PM",
      verbose = FALSE
    )
  ))
  direct_alignment <- suppressWarnings(sirt::invariance.alignment(
    lambda = direct_config$lambda,
    nu = direct_config$nu,
    wgt = NULL,
    align.scale = c(lambda = 0.40, nu = 0.20),
    align.pow = c(lambda = 0.25, nu = 0.25),
    eps = 1e-3,
    center = FALSE,
    optimizer = "optim",
    fixed = NULL,
    meth = 1,
    vcov = direct_config$vcov
  ))

  out <- expect_silent(inv_align_auto(data = dat, group = grp, language = "eng", verbose = FALSE))
  expect_identical(out$conclusion_status, "interpretable")
  expect_equal(out$configural$lambda, direct_config$lambda, tolerance = 1e-10)
  expect_equal(out$configural$nu, direct_config$nu, tolerance = 1e-10)
  expect_equal(out$group_parameters, direct_alignment$pars, tolerance = 1e-10)
  expect_identical(out$data_summary$missing_values, 6L)
  expect_identical(
    out$data_summary$missing_by_item,
    c(i1 = 2L, i2 = 3L, i3 = 1L, i4 = 0L)
  )
  expect_identical(out$data_summary$missing_rows, 6L)
  expect_identical(out$data_summary$complete_rows, 234L)
  expect_identical(
    out$data_summary$missing_by_row,
    as.integer(rowSums(is.na(dat)))
  )
  expect_identical(out$data_summary$n_input_by_group, c(A = 120L, B = 120L))
  expect_identical(out$data_summary$n_complete_by_group, c(A = 117L, B = 117L))
  expect_null(out$requested_settings$missing_policy)
  expect_null(out$effective_settings$missing_engine_request)
  expect_identical(out$data_summary$missing_policy, "listwise")
  expect_identical(out$data_summary$missing_policy_effective, "listwise")
  expect_identical(out$data_summary$configural_n_by_group, direct_config$N)
  expect_identical(out$data_summary$n_effective_by_group, direct_config$N)
  expect_identical(out$effective_settings$missing_policy, "listwise")

  invisible(utils::capture.output(
    direct_config_fiml <- sirt::invariance_alignment_cfa_config(
      dat = dat,
      group = grp,
      weights = NULL,
      model = "2PM",
      verbose = FALSE,
      missing = "fiml"
    )
  ))
  fiml <- expect_silent(inv_align_auto(
    data = dat,
    group = grp,
    language = "eng",
    config_args = list(missing = "fiml"), verbose = FALSE
  ))
  expect_identical(fiml$requested_settings$missing_policy, "fiml")
  expect_identical(fiml$effective_settings$missing_engine_request, "fiml")
  expect_identical(fiml$effective_settings$missing_policy, "fiml")
  expect_identical(fiml$data_summary$missing_policy_effective, "fiml")
  expect_identical(fiml$engine_calls$configural$arguments$missing, "fiml")
  expect_identical(fiml$data_summary$configural_n_by_group, direct_config_fiml$N)
  expect_identical(fiml$data_summary$n_effective_by_group, direct_config_fiml$N)
  expect_equal(unname(fiml$data_summary$n_effective_by_group), c(120, 120))

  bad_nan <- dat
  bad_nan$i1[1] <- NaN
  bad_inf <- dat
  bad_inf$i2[1] <- Inf
  expect_error(inv_align_auto(data = bad_nan, group = grp), "NaN|infinite|finite")
  expect_error(inv_align_auto(data = bad_inf, group = grp), "NaN|infinite|finite")
})

test_that("data-mode group counts follow the effective configural order", {
  skip_if_not_installed("sirt")
  set.seed(818)
  group <- c(rep("B", 220L), rep("A", 180L))
  latent <- stats::rnorm(length(group))
  data <- data.frame(
    i1 = 0.80 * latent + stats::rnorm(length(group), 0, 0.50),
    i2 = 0.75 * latent + stats::rnorm(length(group), 0, 0.55),
    i3 = 0.70 * latent + stats::rnorm(length(group), 0, 0.60),
    i4 = 0.65 * latent + stats::rnorm(length(group), 0, 0.65)
  )
  data$i1[1:3] <- NA_real_
  data$i2[221:225] <- NA_real_

  out <- suppressWarnings(inv_align_auto(
    data = data,
    group = group,
    language = "eng"
  ))

  expect_identical(out$conclusion_status, "interpretable")
  expect_identical(out$data_summary$group_names, c("A", "B"))
  expect_identical(out$data_summary$n_input_by_group, c(A = 180L, B = 220L))
  expect_identical(out$data_summary$n_complete_by_group, c(A = 175L, B = 217L))
  expect_identical(
    names(out$data_summary$n_effective_by_group),
    out$data_summary$group_names
  )
})

test_that("data-mode rejects configural group identifier drift", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  lambda <- x$lambda
  nu <- x$nu
  rownames(lambda) <- c("A", "C")
  rownames(nu) <- c("A", "C")
  data <- data.frame(
    i1 = c(2, 1, 4, 3, 5, 7, 6, 9),
    i2 = c(1, 3, 2, 4, 6, 5, 7, 8),
    i3 = c(3, 2, 1, 5, 4, 6, 8, 7)
  )
  group <- rep(c("A", "B"), each = 4L)
  alignment_called <- FALSE

  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) {
      list(lambda = lambda, nu = nu, vcov = NULL, N = c(A = 4L, C = 4L))
    },
    .pm_alignment_engine = function(args) {
      alignment_called <<- TRUE
      stop("alignment engine must not run")
    },
    .package = "PsychoMatic"
  )
  out <- inv_align_auto(data = data, group = group, language = "eng")

  expect_identical(out$conclusion_status, "not_evaluable")
  expect_identical(
    out$diagnostics$reason,
    "configural_parameter_matrices_invalid"
  )
  expect_identical(out$effective_settings$configural_status, "invalid_output")
  expect_false(alignment_called)
})

test_that("all protected pass-through arguments are rejected", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  for (key in c("dat", "group", "weights", "model")) {
    expect_error(
      inv_align_auto(
        lambda = x$lambda,
        nu = x$nu,
        config_args = stats::setNames(list(1), key)
      ),
      "protected",
      ignore.case = TRUE
    )
  }
  for (key in c(
    "lambda", "nu", "wgt", "align.scale", "align.pow", "eps",
    "center", "optimizer", "fixed", "meth", "vcov"
  )) {
    expect_error(
      inv_align_auto(
        lambda = x$lambda,
        nu = x$nu,
        alignment_args = stats::setNames(list(1), key)
      ),
      "protected",
      ignore.case = TRUE
    )
  }
})

test_that("pass-through lists require unique complete names and retain benign options", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  expect_error(
    inv_align_auto(lambda = x$lambda, nu = x$nu, alignment_args = list(FALSE)),
    "named",
    ignore.case = TRUE
  )
  duplicate <- structure(list(FALSE, FALSE), names = c("num_deriv", "num_deriv"))
  expect_error(
    inv_align_auto(lambda = x$lambda, nu = x$nu, alignment_args = duplicate),
    "duplicate",
    ignore.case = TRUE
  )

  extra <- list(eps_grid = c(0, -1), num_deriv = FALSE)
  out <- suppressWarnings(inv_align_auto(
    lambda = x$lambda,
    nu = x$nu,
    alignment_args = extra
  ))
  expect_identical(out$engine_calls$alignment$arguments$eps_grid, c(0, -1))
  expect_identical(out$engine_calls$alignment$arguments$num_deriv, FALSE)
  expect_identical(out$requested_settings$alignment_args, extra)
  expect_identical(out$effective_settings$alignment_args, extra)
})

test_that("semantic tuning pairs reorder by name and invalid scalars fail", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  out <- suppressWarnings(inv_align_auto(
    lambda = x$lambda,
    nu = x$nu,
    align.scale = c(nu = 0.2, lambda = 0.4),
    align.pow = c(nu = 0.25, lambda = 0.25),
    parm_tol = c(nu = 0.2, lambda = 0.4),
    optimizer = "nlminb",
    fixed = FALSE,
    meth = 2
  ))
  expect_identical(out$effective_settings$align.scale, c(lambda = 0.4, nu = 0.2))
  expect_identical(out$effective_settings$parm_tol, c(lambda = 0.4, nu = 0.2))
  expect_identical(out$engine_calls$alignment$arguments$optimizer, "nlminb")
  expect_identical(out$engine_calls$alignment$arguments$fixed, FALSE)
  expect_identical(out$effective_settings$fixed_engine_request, FALSE)
  expect_identical(out$effective_settings$fixed_effective, FALSE)
  expect_identical(out$effective_settings$fixed, FALSE)
  expect_identical(out$engine_calls$alignment$arguments$meth, 2)

  bad_pairs <- list(
    c(nu = 0.2, other = 0.4),
    structure(c(0.2, 0.4), names = c("nu", "nu")),
    c(lambda = 0.4, nu = Inf)
  )
  for (value in bad_pairs) {
    expect_error(
      inv_align_auto(lambda = x$lambda, nu = x$nu, align.scale = value),
      "align.scale|lambda|nu|finite|duplicate",
      ignore.case = TRUE
    )
  }
  expect_error(inv_align_auto(lambda = x$lambda, nu = x$nu, eps = 0), "eps")
  expect_error(inv_align_auto(lambda = x$lambda, nu = x$nu, optimizer = "bad"), "optimizer")
  expect_error(inv_align_auto(lambda = x$lambda, nu = x$nu, meth = 99), "meth")
  expect_error(inv_align_auto(lambda = x$lambda, nu = x$nu, digits = 2.5), "digits")
})

test_that("alignment vcov must match the parameter vector and be a covariance matrix", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  called <- 0L
  testthat::local_mocked_bindings(
    .pm_alignment_engine = function(args) {
      called <<- called + 1L
      stop("alignment engine must not run")
    },
    .package = "PsychoMatic"
  )

  asymmetric <- diag(12)
  asymmetric[1, 2] <- 0.1
  indefinite <- diag(12)
  indefinite[1, 1] <- -1
  mismatched_names <- diag(12)
  parameter_names <- unlist(lapply(seq_len(2L), function(group_index) {
    c(
      paste0(colnames(x$lambda), "_lam_Gr", group_index),
      paste0(colnames(x$lambda), "_nu_Gr", group_index)
    )
  }), use.names = FALSE)
  rownames(mismatched_names) <- parameter_names
  colnames(mismatched_names) <- rev(parameter_names)

  invalid <- list(
    wrong_dimension = diag(11),
    asymmetric = asymmetric,
    indefinite = indefinite,
    mismatched_names = mismatched_names
  )
  for (candidate in invalid) {
    out <- inv_align_auto(
      lambda = x$lambda,
      nu = x$nu,
      vcov = candidate,
      language = "eng"
    )
    expect_identical(out$conclusion_status, "not_evaluable")
    expect_identical(out$diagnostics$reason, "alignment_vcov_invalid")
    expect_null(out$engine_calls$alignment)
    expect_true(all(is.na(out$summary$R2)))
  }
  expect_identical(called, 0L)
})

test_that("engine and output failures are non-evaluable without favorable prose", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()

  run_with_fit <- function(fit = real$fit, constraints = real$constraints) {
    testthat::local_mocked_bindings(
      .pm_alignment_engine = function(args) fit,
      .pm_alignment_constraints_engine = function(args) constraints,
      .package = "PsychoMatic"
    )
    inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
  }

  cases <- list()
  z <- real$fit
  z$pars[,] <- NA_real_
  cases$all_na_pars <- z
  z <- real$fit
  z$itempars.aligned[1, 1] <- Inf
  cases$nonfinite_item <- z
  z <- real$fit
  z$lambda.aligned <- z$lambda.aligned[-1, , drop = FALSE]
  cases$wrong_dimension <- z
  z <- real$fit
  z$res_optim$converged <- FALSE
  cases$nonconvergence <- z
  cases$empty_output <- list()

  for (fit in cases) {
    out <- suppressWarnings(run_with_fit(fit = fit))
    expect_identical(out$conclusion_status, "not_evaluable")
    expect_true(all(is.na(out$summary$R2)))
    expect_false(grepl("meets|greater invariance", out$report, ignore.case = TRUE))
    expect_false(grepl("meets|greater invariance", out$conclusion, ignore.case = TRUE))
  }

  sentinel <- structure(
    list(message = "alignment sentinel", call = quote(sirt::invariance.alignment())),
    class = c("alignment_sentinel", "error", "condition")
  )
  testthat::local_mocked_bindings(
    .pm_alignment_engine = function(args) stop(sentinel),
    .package = "PsychoMatic"
  )
  failed <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
  expect_identical(failed$conclusion_status, "not_evaluable")
  expect_equal(nrow(failed$errors), 1L)
  expect_s3_class(failed$errors$condition[[1L]], "alignment_sentinel")

  testthat::local_mocked_bindings(
    .pm_alignment_engine = function(args) real$fit,
    .pm_alignment_constraints_engine = function(args) stop("constraint sentinel"),
    .package = "PsychoMatic"
  )
  failed_constraints <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
  expect_identical(failed_constraints$conclusion_status, "not_evaluable")
  expect_match(failed_constraints$diagnostics$reason, "constraint")

  malformed_constraints <- real$constraints
  malformed_constraints$lambda_list$parm_dif[,] <- NA_real_
  testthat::local_mocked_bindings(
    .pm_alignment_engine = function(args) real$fit,
    .pm_alignment_constraints_engine = function(args) malformed_constraints,
    .package = "PsychoMatic"
  )
  failed_constraint_output <- inv_align_auto(
    lambda = x$lambda,
    nu = x$nu,
    language = "eng"
  )
  expect_identical(failed_constraint_output$conclusion_status, "not_evaluable")
  expect_match(failed_constraint_output$diagnostics$reason, "lambda_list")
})

test_that("non-finite alignment objectives fail closed before constraints", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()

  cases <- list()
  z <- real$fit
  z$fopt <- Inf
  cases$nonfinite_fopt <- z
  z <- real$fit
  z$res_optim$value <- NaN
  cases$nonfinite_optimizer_value <- z
  z <- real$fit
  z$fopt <- NULL
  z$res_optim$value <- NULL
  cases$missing_objective <- z

  for (fit in cases) {
    testthat::local_mocked_bindings(
      .pm_alignment_engine = function(args) fit,
      .pm_alignment_constraints_engine = function(args) {
        stop("constraint engine must not run")
      },
      .package = "PsychoMatic"
    )
    out <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
    expect_identical(out$conclusion_status, "not_evaluable")
    expect_identical(
      out$diagnostics$reason,
      "alignment_objective_not_evaluable"
    )
    expect_null(out$engine_calls$constraints)
    expect_identical(nrow(out$errors), 0L)
    expect_true(all(is.na(out$summary$R2)))
    expect_false(grepl("meets|greater invariance", out$report, ignore.case = TRUE))
  }
})

test_that("negative alignment objectives fail closed", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()

  negative <- real$fit
  negative$fopt <- -1
  negative$res_optim$value <- -1
  cases <- list(
    negative = list(
      fit = negative,
      reason = "alignment_objective_out_of_range"
    )
  )
  for (case in cases) {
    testthat::local_mocked_bindings(
      .pm_alignment_engine = function(args) case$fit,
      .pm_alignment_constraints_engine = function(args) {
        stop("constraint engine must not run")
      },
      .package = "PsychoMatic"
    )
    out <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
    expect_identical(out$conclusion_status, "not_evaluable")
    expect_identical(out$diagnostics$reason, case$reason)
    expect_null(out$engine_calls$constraints)
    expect_true(all(is.na(out$summary$R2)))
  }
})

test_that("distinct finite non-negative alignment objectives are accepted", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()
  distinct <- real$fit
  distinct$fopt <- 1
  distinct$res_optim$value <- 2

  testthat::local_mocked_bindings(
    .pm_alignment_engine = function(args) distinct,
    .pm_alignment_constraints_engine = function(args) real$constraints,
    .package = "PsychoMatic"
  )
  out <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")

  expect_identical(out$conclusion_status, "interpretable")
})

test_that("valid sirt objectives from different stages need not be equal", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()

  direct <- suppressWarnings(sirt::invariance.alignment(
    lambda = x$lambda,
    nu = x$nu,
    align.scale = c(lambda = 0.40, nu = 0.20),
    align.pow = c(lambda = 0.25, nu = 0.25),
    eps = 1e-3,
    center = FALSE,
    optimizer = "optim",
    fixed = TRUE,
    meth = 2
  ))
  expect_identical(direct$res_optim$convergence, 0L)
  expect_false(isTRUE(all.equal(
    direct$fopt,
    direct$res_optim$value,
    tolerance = sqrt(.Machine$double.eps)
  )))

  out <- suppressWarnings(inv_align_auto(
    lambda = x$lambda,
    nu = x$nu,
    language = "eng",
    align.scale = c(lambda = 0.40, nu = 0.20),
    align.pow = c(lambda = 0.25, nu = 0.25),
    eps = 1e-3,
    center = FALSE,
    optimizer = "optim",
    fixed = TRUE,
    meth = 2
  ))

  expect_identical(out$conclusion_status, "interpretable")
  expect_identical(out$diagnostics$reason, "alignment_outputs_valid")
})

test_that("constraint percentages must be bounded and coherent with sirt counts", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()

  run_with_constraints <- function(constraints) {
    testthat::local_mocked_bindings(
      .pm_alignment_engine = function(args) real$fit,
      .pm_alignment_constraints_engine = function(args) constraints,
      .package = "PsychoMatic"
    )
    inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
  }

  for (value in c(-0.1, 101)) {
    malformed <- real$constraints
    malformed$lambda_list$prop_noninvariance <- value
    out <- run_with_constraints(malformed)
    expect_identical(out$conclusion_status, "not_evaluable")
    expect_identical(
      out$diagnostics$reason,
      "lambda_list_prop_noninvariance_out_of_range"
    )
    expect_true(all(is.na(out$summary$percent_noninvariant)))
    expect_true(all(is.na(out$summary$invariant)))
    expect_false(grepl("meets|greater invariance", out$report, ignore.case = TRUE))
  }

  inconsistent <- real$constraints
  inconsistent$lambda_list$prop_noninvariance <- 10
  out <- run_with_constraints(inconsistent)
  expect_identical(out$conclusion_status, "not_evaluable")
  expect_identical(
    out$diagnostics$reason,
    "lambda_list_prop_noninvariance_inconsistent"
  )
  expect_true(all(is.na(out$summary$percent_noninvariant)))
})

test_that("constraint availability indicators must be strictly binary", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()
  for (invalid_value in c(0.5, 2)) {
    malformed <- real$constraints
    malformed$lambda_list$miss_items[,] <- invalid_value
    malformed$lambda_list$N_total <- sum(malformed$lambda_list$miss_items)
    malformed$lambda_list$prop_noninvariance <-
      100 * sum(malformed$lambda_list$N_unique_parm_groups) /
      malformed$lambda_list$N_total

    testthat::local_mocked_bindings(
      .pm_alignment_engine = function(args) real$fit,
      .pm_alignment_constraints_engine = function(args) malformed,
      .package = "PsychoMatic"
    )
    out <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")

    expect_identical(out$conclusion_status, "not_evaluable")
    expect_identical(
      out$diagnostics$reason,
      "lambda_list_miss_items_not_binary"
    )
    expect_true(all(is.na(out$summary$percent_noninvariant)))
  }
})

test_that("configural failures are retained as non-evaluable", {
  skip_if_not_installed("sirt")
  dat <- data.frame(i1 = 1:8, i2 = c(2:8, 10))
  grp <- rep(c("A", "B"), each = 4)
  sentinel <- structure(
    list(message = "config sentinel", call = quote(sirt::invariance_alignment_cfa_config())),
    class = c("config_sentinel", "error", "condition")
  )
  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) stop(sentinel),
    .package = "PsychoMatic"
  )
  out <- inv_align_auto(data = dat, group = grp, language = "eng")
  expect_identical(out$conclusion_status, "not_evaluable")
  expect_match(out$diagnostics$reason, "configural")
  expect_s3_class(out$errors$condition[[1L]], "config_sentinel")
  expect_identical(out$effective_settings$configural_status, "failed")
  expect_null(out$effective_settings$fixed_engine_request)
  expect_null(out$effective_settings$fixed_effective)
  expect_null(out$effective_settings$fixed)
  expect_null(out$engine_calls$alignment)
})

test_that("configural early exits record the effective terminal status", {
  skip_if_not_installed("sirt")
  dat <- data.frame(i1 = 1:8, i2 = c(2:8, 10))
  grp <- rep(c("A", "B"), each = 4)

  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) list(converged = FALSE),
    .package = "PsychoMatic"
  )
  nonconverged <- inv_align_auto(data = dat, group = grp, language = "eng")
  expect_identical(nonconverged$conclusion_status, "not_evaluable")
  expect_identical(
    nonconverged$effective_settings$configural_status,
    "nonconverged"
  )
  expect_null(nonconverged$engine_calls$alignment)

  x <- alignment_fixture()
  testthat::local_mocked_bindings(
    .pm_alignment_config_engine = function(args) {
      list(lambda = x$lambda, nu = x$nu, vcov = -diag(12))
    },
    .package = "PsychoMatic"
  )
  invalid <- inv_align_auto(data = dat, group = grp, language = "eng")
  expect_identical(invalid$conclusion_status, "not_evaluable")
  expect_identical(invalid$effective_settings$configural_status, "invalid_output")
  expect_identical(invalid$diagnostics$reason, "alignment_vcov_invalid")
  expect_null(invalid$engine_calls$alignment)
})

test_that("unreported convergence does not override otherwise valid output", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  real <- alignment_real_outputs()
  real$fit$res_optim$converged <- NULL
  real$fit$res_optim$convergence <- NULL
  testthat::local_mocked_bindings(
    .pm_alignment_engine = function(args) real$fit,
    .pm_alignment_constraints_engine = function(args) real$constraints,
    .package = "PsychoMatic"
  )
  out <- inv_align_auto(lambda = x$lambda, nu = x$nu, language = "eng")
  expect_identical(out$conclusion_status, "interpretable")
  expect_identical(out$diagnostics$convergence_status, "not_reported")
})

test_that("digits changes presentation only and defaults are recreated per call", {
  skip_if_not_installed("sirt")
  x <- alignment_fixture()
  low <- suppressWarnings(inv_align_auto(lambda = x$lambda, nu = x$nu, digits = 2))
  high <- suppressWarnings(inv_align_auto(lambda = x$lambda, nu = x$nu, digits = 6))
  expect_identical(low$summary, high$summary)
  expect_identical(low$group_parameters, high$group_parameters)
  expect_identical(low$aligned_item_parameters, high$aligned_item_parameters)
  expect_false(identical(low$report, high$report))

  custom <- suppressWarnings(inv_align_auto(
    lambda = x$lambda,
    nu = x$nu,
    alignment_args = list(eps_grid = c(0, -1))
  ))
  default_again <- suppressWarnings(inv_align_auto(lambda = x$lambda, nu = x$nu))
  expect_true("eps_grid" %in% names(custom$engine_calls$alignment$arguments))
  expect_false("eps_grid" %in% names(default_again$engine_calls$alignment$arguments))
})
