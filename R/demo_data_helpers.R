.pm_demo_with_seed <- function(seed, generator) {
  old_kind <- RNGkind()
  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (seed_existed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (seed_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  set.seed(seed)
  generator()
}

.pm_demo_latent_factors <- function(n, correlation_loading = 0.45) {
  general <- stats::rnorm(n)
  residual_scale <- sqrt(1 - correlation_loading^2)
  cbind(
    factor1 = correlation_loading * general + residual_scale * stats::rnorm(n),
    factor2 = correlation_loading * general + residual_scale * stats::rnorm(n),
    factor3 = correlation_loading * general + residual_scale * stats::rnorm(n)
  )
}

.pm_demo_discretize <- function(x, thresholds) {
  as.integer(findInterval(x, thresholds, rightmost.closed = FALSE) + 1L)
}

.pm_generate_likert_demo <- function() {
  seed <- 2401L
  .pm_demo_with_seed(seed, function() {
    n <- 300L
    item_names <- sprintf("item%02d", seq_len(12L))
    domains <- list(
      factor1 = item_names[1:4],
      factor2 = item_names[5:8],
      factor3 = item_names[9:12]
    )
    reverse_keys <- c("item03", "item07", "item11")
    factor_index <- rep(seq_len(3L), each = 4L)
    loadings <- c(
      0.78, 0.72, 0.70, 0.76,
      0.74, 0.79, 0.68, 0.73,
      0.77, 0.71, 0.75, 0.69
    )
    latent <- .pm_demo_latent_factors(n)
    continuous <- vapply(seq_along(item_names), function(i) {
      loadings[[i]] * latent[, factor_index[[i]]] +
        sqrt(1 - loadings[[i]]^2) * stats::rnorm(n)
    }, numeric(n))
    colnames(continuous) <- item_names
    continuous[, reverse_keys] <- -continuous[, reverse_keys, drop = FALSE]

    threshold_shifts <- seq(-0.12, 0.12, length.out = length(item_names))
    responses <- vapply(seq_along(item_names), function(i) {
      .pm_demo_discretize(
        continuous[, i],
        c(-1.10, -0.35, 0.35, 1.10) + threshold_shifts[[i]]
      )
    }, integer(n))
    colnames(responses) <- item_names
    missing <- matrix(
      stats::runif(length(responses)) < 0.025,
      nrow = n,
      ncol = length(item_names)
    )
    missing[1L, 1L] <- TRUE
    responses[missing] <- NA_integer_
    out <- as.data.frame(responses, check.names = FALSE)

    attr(out, "generator_seed") <- seed
    attr(out, "simulation_model") <- paste(
      "Three correlated latent factors with four indicators each,",
      "item-specific thresholds, reverse-keyed indicators, and MCAR missingness."
    )
    attr(out, "factor_structure") <- domains
    attr(out, "reverse_keys") <- reverse_keys
    attr(out, "scoring_keys") <- stats::setNames(
      ifelse(item_names %in% reverse_keys, -1L, 1L),
      item_names
    )
    attr(out, "missing_rate") <- mean(is.na(responses))
    attr(out, "real_participant_data") <- FALSE
    out
  })
}

.pm_generate_invariance_demo <- function() {
  seed <- 2402L
  .pm_demo_with_seed(seed, function() {
    n_per_group <- 300L
    n <- 2L * n_per_group
    group <- ordered(
      rep(c("A", "B"), each = n_per_group),
      levels = c("A", "B")
    )
    item_names <- paste0("i", seq_len(9L))
    domains <- list(
      factor1 = item_names[1:3],
      factor2 = item_names[4:6],
      factor3 = item_names[7:9]
    )
    factor_index <- rep(seq_len(3L), each = 3L)
    loadings <- c(0.78, 0.73, 0.76, 0.72, 0.80, 0.69, 0.77, 0.74, 0.71)
    latent <- .pm_demo_latent_factors(n)
    continuous <- vapply(seq_along(item_names), function(i) {
      loadings[[i]] * latent[, factor_index[[i]]] +
        sqrt(1 - loadings[[i]]^2) * stats::rnorm(n)
    }, numeric(n))
    colnames(continuous) <- item_names

    base_thresholds <- c(-1.10, -0.35, 0.25, 1.15)
    threshold_shift <- 0.60
    responses <- matrix(NA_integer_, nrow = n, ncol = length(item_names))
    for (i in seq_along(item_names)) {
      if (identical(item_names[[i]], "i3")) {
        group_a <- group == "A"
        shifted <- base_thresholds
        shifted[[3L]] <- shifted[[3L]] + threshold_shift
        responses[group_a, i] <- .pm_demo_discretize(
          continuous[group_a, i],
          base_thresholds
        )
        responses[!group_a, i] <- .pm_demo_discretize(
          continuous[!group_a, i],
          shifted
        )
      } else {
        responses[, i] <- .pm_demo_discretize(
          continuous[, i],
          base_thresholds
        )
      }
    }
    colnames(responses) <- item_names
    missing <- matrix(
      stats::runif(length(responses)) < 0.01,
      nrow = n,
      ncol = length(item_names)
    )
    responses[missing] <- NA_integer_
    items <- lapply(seq_along(item_names), function(i) {
      ordered(as.character(responses[, i]), levels = as.character(1:5))
    })
    names(items) <- item_names
    out <- data.frame(group = group, items, check.names = FALSE)

    attr(out, "generator_seed") <- seed
    attr(out, "simulation_model") <- paste(
      "Three correlated latent factors with invariant loadings and thresholds,",
      "except a generated group-B shift in the third threshold of i3."
    )
    attr(out, "factor_structure") <- domains
    attr(out, "known_noninvariance") <- list(
      parameter = "i3 | t3",
      group = "B",
      threshold_index = 3L,
      threshold_shift = threshold_shift,
      direction = "higher_threshold"
    )
    attr(out, "group_sizes") <- c(A = n_per_group, B = n_per_group)
    attr(out, "missing_rate") <- mean(is.na(responses))
    attr(out, "real_participant_data") <- FALSE
    out
  })
}
