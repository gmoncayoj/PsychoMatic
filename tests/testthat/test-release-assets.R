.release_root <- function() {
  normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
}

.release_path <- function(...) file.path(.release_root(), ...)

.read_release_file <- function(path) {
  if (!grepl("^(?:[A-Za-z]:|/)", path, perl = TRUE)) path <- .release_path(path)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

test_that("DESCRIPTION declares the local 1.0.0 candidate and dependency floors", {
  description <- read.dcf(.release_path("DESCRIPTION"))
  expect_identical(unname(description[1L, "Version"]), "1.0.0")
  expect_match(description[1L, "Imports"], "lavaan \\(>= 0\\.7-3\\)")
  suggests <- description[1L, "Suggests"]
  expect_match(suggests, "MVN \\(>= 6\\.3\\)")
  expect_match(suggests, "semTools \\(>= 0\\.5-9\\)")
  expect_match(suggests, "sirt \\(>= 4\\.2-133\\)")
  expect_match(suggests, "mirt")
  expect_match(suggests, "openxlsx")
  expect_match(suggests, "officer")
  expect_match(suggests, "flextable")
})

test_that("demo generators are deterministic and preserve caller RNG state", {
  expect_true(exists(".pm_generate_likert_demo", asNamespace("PsychoMatic")))
  expect_true(exists(".pm_generate_invariance_demo", asNamespace("PsychoMatic")))
  likert_generator <- get(".pm_generate_likert_demo", asNamespace("PsychoMatic"))
  invariance_generator <- get(".pm_generate_invariance_demo", asNamespace("PsychoMatic"))

  old_kind <- RNGkind()
  set.seed(991L)
  old_seed <- .Random.seed
  first <- likert_generator()
  expect_identical(RNGkind(), old_kind)
  expect_identical(.Random.seed, old_seed)
  second <- likert_generator()
  expect_identical(serialize(first, NULL), serialize(second, NULL))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  invariance_generator()
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("likert_demo has its documented reproducible structure", {
  expect_true(file.exists(.release_path("tests", "fixtures", "likert_demo.rda")))
  environment <- new.env(parent = emptyenv())
  load(.release_path("tests", "fixtures", "likert_demo.rda"), envir = environment)
  observed <- environment$likert_demo
  expected <- PsychoMatic:::.pm_generate_likert_demo()

  expect_identical(observed, expected)
  expect_identical(dim(observed), c(300L, 12L))
  expect_identical(names(observed), sprintf("item%02d", 1:12))
  values <- unlist(observed, use.names = FALSE)
  expect_true(all(is.na(values) | values %in% 1:5))
  expect_true(anyNA(observed))
  expect_identical(attr(observed, "reverse_keys"), c("item03", "item07", "item11"))
  expect_identical(attr(observed, "real_participant_data"), FALSE)
  expect_identical(attr(observed, "generator_seed"), 2401L)
  expect_type(attr(observed, "factor_structure"), "list")
})

test_that("invariance_demo contains ordered groups/items and the generated threshold signal", {
  expect_true(file.exists(.release_path("tests", "fixtures", "invariance_demo.rda")))
  environment <- new.env(parent = emptyenv())
  load(.release_path("tests", "fixtures", "invariance_demo.rda"), envir = environment)
  observed <- environment$invariance_demo
  expected <- PsychoMatic:::.pm_generate_invariance_demo()

  expect_identical(observed, expected)
  expect_identical(dim(observed), c(600L, 10L))
  expect_identical(names(observed), c("group", paste0("i", 1:9)))
  expect_true(is.ordered(observed$group))
  expect_identical(levels(observed$group), c("A", "B"))
  expect_true(all(vapply(observed[paste0("i", 1:9)], is.ordered, logical(1L))))
  expect_true(all(vapply(
    observed[paste0("i", 1:9)],
    function(item) identical(levels(item), as.character(1:5)),
    logical(1L)
  )))
  expect_identical(attr(observed, "known_noninvariance")$parameter, "i3 | t3")
  expect_identical(attr(observed, "real_participant_data"), FALSE)
  signal <- with(observed, {
    mean(as.integer(i3[group == "B"]) <= 3L, na.rm = TRUE) -
      mean(as.integer(i3[group == "A"]) <= 3L, na.rm = TRUE)
  })
  controls <- vapply(setdiff(paste0("i", 1:9), "i3"), function(item) {
    values <- as.integer(observed[[item]])
    mean(values[observed$group == "B"] <= 3L, na.rm = TRUE) -
      mean(values[observed$group == "A"] <= 3L, na.rm = TRUE)
  }, numeric(1L))
  expect_gt(signal, 0.12)
  expect_gt(signal, max(abs(controls)) + 0.04)
})

test_that("all five datasets load and use appropriate rda compression", {
  expected <- c(
    "psychomatic_continuous", "psychomatic_ordinal", "psychomatic_multigroup",
    "likert_demo", "invariance_demo"
  )
  for (name in expected) {
    path <- .release_path("tests", "fixtures", paste0(name, ".rda"))
    expect_true(file.exists(path), info = name)
    environment <- new.env(parent = emptyenv())
    expect_identical(load(path, envir = environment), name)
  }
  checked <- tools::checkRdaFiles(.release_path("tests", "fixtures"))
  new_rows <- checked[basename(rownames(checked)) %in%
                        c("likert_demo.rda", "invariance_demo.rda"), , drop = FALSE]
  expect_identical(nrow(new_rows), 2L)
  expect_true(all(new_rows$compress == "xz"))
})

test_that("the twelve shipped scripts are concise, assigned, and portable", {
  assigned <- c(
    EXAMPLE_DESC_AUTO.R = "desc_auto",
    EXAMPLE_EFA_AUTO.R = "efa_auto",
    EXAMPLE_ESEM_AUTO.R = "esem_auto",
    EXAMPLE_CFA_AUTO.R = "cfa_auto",
    EXAMPLE_FACTORIAL_INVARIANCE_AUTO.R = "factorial_invariance_auto",
    EXAMPLE_ORDINAL_INVARIANCE_AUTO.R = "ordinal_invariance_auto",
    EXAMPLE_INVARIANCE_ALIGNMENT_AUTO.R = "inv_align_auto",
    EXAMPLE_FULL_PSYCHOMATIC_WORKFLOW.R = "psychomatic_report",
    EXAMPLE_SPEC_CHECK_SCORING.R = "psychomatic_spec",
    EXAMPLE_RELIABILITY_AUTO.R = "reliability_auto",
    EXAMPLE_INVARIANCE_PARTIAL.R = "invariance_partial",
    EXAMPLE_IRT_DIF.R = "irt_auto"
  )
  paths <- .release_path("inst", "examples", names(assigned))
  expect_true(all(file.exists(paths)))
  expect_identical(length(list.files(
    .release_path("inst", "examples"), pattern = "^EXAMPLE_.*\\.R$"
  )), length(assigned))

  for (i in seq_along(paths)) {
    text <- .read_release_file(paths[[i]])
    expect_match(text, paste0("\\b", assigned[[i]], "\\s*\\("), info = paths[[i]])
    expect_false(grepl("\\blibrary\\s*\\(|\\brequire\\s*\\(", text), info = paths[[i]])
    expect_false(grepl(":::", text, fixed = TRUE), info = paths[[i]])
    expect_false(grepl("(?<![[:alnum:]])[A-Za-z]:[/\\\\]", text, perl = TRUE), info = paths[[i]])
    expect_false(grepl("download\\.file|browseURL|url\\s*\\(|setwd\\s*\\(", text), info = paths[[i]])
  }
  irt_text <- .read_release_file("inst/examples/EXAMPLE_IRT_DIF.R")
  expect_match(irt_text, "dif_auto\\s*\\(")
})

test_that("release documentation covers public additions and methodological cautions", {
  files <- c("README.md", "NEWS.md", "R/PsychoMatic-package.R", list.files(
    .release_path("vignettes"), pattern = "\\.Rmd$", full.names = TRUE
  ))
  text <- paste(vapply(files, .read_release_file, character(1L)), collapse = "\n")
  for (function_name in c(
    "psychomatic_spec", "reliability_auto",
    "invariance_partial", "psychomatic_report"
  )) {
    expect_match(text, paste0(function_name, "\\(\\)"), info = function_name)
  }
  expect_match(.read_release_file("README.md"), "1\\.0\\.0")
  expect_false(grepl("Version 0\\.3\\.0 focuses", .read_release_file("README.md")))
  expect_match(text, "heuristic", ignore.case = TRUE)
  expect_match(text, "not[_ -]?evaluable|no evaluable", ignore.case = TRUE)
  expect_match(text, "anchor|ancla", ignore.case = TRUE)
  expect_match(text, "rationale|justificaci", ignore.case = TRUE)
  expect_match(text, "partial invariance|invariancia parcial", ignore.case = TRUE)
  expect_match(text, "overwrite", ignore.case = TRUE)
  expect_match(text, "full precision|raw value|precisi.n completa|valor bruto", ignore.case = TRUE)
})

test_that("shipped R Markdown and data scripts avoid attached packages and internals", {
  paths <- c(
    list.files(.release_path("vignettes"), pattern = "\\.Rmd$", full.names = TRUE),
    .release_path("data-raw", "generate_demo_data.R")
  )
  expect_true(all(file.exists(paths)))
  for (path in paths) {
    text <- .read_release_file(path)
    expect_false(grepl("\\blibrary\\s*\\(|\\brequire\\s*\\(", text), info = path)
    expect_false(grepl(":::", text, fixed = TRUE), info = path)
    expect_false(grepl("(?<![[:alnum:]])[A-Za-z]:[/\\\\]", text, perl = TRUE), info = path)
  }
})

test_that("portability checks distinguish Windows paths from reference URLs", {
  pattern <- "(?<![[:alnum:]])[A-Za-z]:[/\\\\]"
  expect_true(grepl(pattern, 'C:/Users/example/data.csv', perl = TRUE))
  expect_true(grepl(pattern, '"D:\\data\\input.csv"', perl = TRUE))
  expect_false(grepl(pattern, 'https://doi.org/10.1111/bmsp.12336', perl = TRUE))
})

test_that("ESEM runtime strings use portable Unicode escapes", {
  for (file in c('esem_methods.R', 'esem_results.R')) {
    lines <- readLines(.release_path('R', file), encoding = 'UTF-8', warn = FALSE)
    lines <- lines[!grepl('^\\s*#', lines)]
    expect_false(any(grepl('[^\\x01-\\x7F]', lines, perl = TRUE)), info = file)
  }
})

test_that("build exclusions cover all development-only roots", {
  ignore <- readLines(.release_path(".Rbuildignore"), warn = FALSE)
  expected <- c(
    "^tests/testthat/test-release-assets[.]R$",
    "^docs$", "^docs/", "^data-raw$", "^data-raw/",
    "^\\.github$", "^\\.github/",
    "^\\.superpowers$", "^\\.superpowers/",
    "^\\.agents$", "^\\.agents/",
    "^\\.codex$", "^\\.codex/"
  )
  expect_true(all(expected %in% ignore), info = paste(setdiff(expected, ignore), collapse = ", "))
})

test_that("CI matrix covers CRAN-relevant platforms and R versions", {
  skip_if_not_installed("yaml")
  path <- .release_path(".github", "workflows", "R-CMD-check.yaml")
  expect_true(file.exists(path))
  workflow <- yaml::read_yaml(path)
  text <- .read_release_file(path)
  expect_type(workflow, "list")
  for (value in c(
    "ubuntu-latest", "windows-latest", "macos-latest",
    "release", "devel", "oldrel-1"
  )) {
    expect_match(text, value, fixed = TRUE)
  }
  expect_match(text, "--as-cran", fixed = TRUE)
  expect_match(text, "error-on:.*warning", ignore.case = TRUE)
  expect_false(grepl("secret|token|password", text, ignore.case = TRUE))
})
