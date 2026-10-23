.pm_fixture_with_seed <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != trunc(seed)) {
    stop("`seed` must be a non-negative integer.", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  eval.parent(substitute(code))
}

.pm_fixture_size <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 2 || x > .Machine$integer.max || x != trunc(x)) {
    stop(sprintf("`%s` must be an integer of at least 2.", name), call. = FALSE)
  }
  as.integer(x)
}

simulate_likert_fixture <- function(n = 300L, seed = 2401L) {
  n <- .pm_fixture_size(n, "n")

  .pm_fixture_with_seed(seed, {
    general <- stats::rnorm(n)
    specific <- matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
    loadings <- seq(0.62, 0.78, length.out = 12L)
    specific_loading <- 0.25
    continuous <- matrix(NA_real_, nrow = n, ncol = 12L)

    for (item in seq_len(12L)) {
      domain <- ((item - 1L) %/% 4L) + 1L
      residual_sd <- sqrt(1 - loadings[item]^2 - specific_loading^2)
      continuous[, item] <-
        loadings[item] * general +
        specific_loading * specific[, domain] +
        stats::rnorm(n, sd = residual_sd)
    }

    ordinal <- apply(
      continuous,
      2L,
      function(x) as.integer(cut(
        x,
        breaks = c(-Inf, -1.1, -0.35, 0.35, 1.1, Inf),
        labels = FALSE
      ))
    )
    ordinal <- as.data.frame(ordinal, stringsAsFactors = FALSE)
    names(ordinal) <- paste0("i", seq_len(12L))
    ordinal
  })
}

simulate_invariance_fixture <- function(n_per_group = 150L,
                                        ordinal = FALSE,
                                        partial = FALSE,
                                        seed = 2402L) {
  n_per_group <- .pm_fixture_size(n_per_group, "n_per_group")
  if (!is.logical(ordinal) || length(ordinal) != 1L || is.na(ordinal)) {
    stop("`ordinal` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(partial) || length(partial) != 1L || is.na(partial)) {
    stop("`partial` must be TRUE or FALSE.", call. = FALSE)
  }

  .pm_fixture_with_seed(seed, {
    group <- factor(
      rep(c("g1", "g2"), each = n_per_group),
      levels = c("g1", "g2")
    )
    n <- length(group)
    latent <- stats::rnorm(n, mean = ifelse(group == "g2", 0.15, 0))
    loadings <- seq(0.65, 0.81, length.out = 9L)
    continuous <- matrix(NA_real_, nrow = n, ncol = 9L)

    for (item in seq_len(9L)) {
      item_loading <- rep(loadings[item], n)
      item_intercept <- rep(0, n)
      if (isTRUE(partial) && item == 3L) {
        item_loading[group == "g2"] <- 0.35
        item_intercept[group == "g2"] <- 0.45
      }
      continuous[, item] <-
        item_intercept + item_loading * latent +
        stats::rnorm(n, sd = sqrt(1 - loadings[item]^2))
    }

    items <- as.data.frame(continuous, stringsAsFactors = FALSE)
    names(items) <- paste0("i", seq_len(9L))
    if (isTRUE(ordinal)) {
      items[] <- lapply(items, function(x) {
        ordered(
          cut(
            x,
            breaks = c(-Inf, -1.1, -0.35, 0.35, 1.1, Inf),
            labels = FALSE
          ),
          levels = 1:5
        )
      })
    }

    data.frame(group = group, items, check.names = FALSE)
  })
}

simulate_bifactor_fixture <- function(n = 400L, seed = 2403L) {
  n <- .pm_fixture_size(n, "n")

  .pm_fixture_with_seed(seed, {
    general <- stats::rnorm(n)
    specific <- matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
    general_loading <- 0.65
    specific_loading <- 0.40
    residual_sd <- sqrt(1 - general_loading^2 - specific_loading^2)
    items <- matrix(NA_real_, nrow = n, ncol = 12L)

    for (item in seq_len(12L)) {
      domain <- ((item - 1L) %/% 4L) + 1L
      items[, item] <-
        general_loading * general +
        specific_loading * specific[, domain] +
        stats::rnorm(n, sd = residual_sd)
    }

    items <- as.data.frame(items, stringsAsFactors = FALSE)
    names(items) <- paste0("i", seq_len(12L))
    items
  })
}
