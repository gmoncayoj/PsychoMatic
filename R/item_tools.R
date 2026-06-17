#' Screen item-level data quality
#'
#' Computes a compact item-screening table for numeric psychometric items:
#' missingness, variance, distribution shape, floor and ceiling proportions,
#' item-rest correlations, and probable reverse-keyed items.
#'
#' @param data Data frame or matrix containing item responses.
#' @param digits Number of decimal places used in the returned table.
#' @param floor_cutoff Proportion at the minimum score used to flag floor effects.
#' @param ceiling_cutoff Proportion at the maximum score used to flag ceiling effects.
#' @param min_item_rest Item-rest correlation below which an item is flagged as
#'   probably reverse keyed or inconsistent with the remaining items.
#'
#' @return A data frame of class `psychomatic_item_screen`.
#'
#' @examples
#' items <- data.frame(
#'   item1 = c(1, 2, 3, 4, 5),
#'   item2 = c(5, 4, 3, 2, 1),
#'   item3 = c(1, 1, 2, 2, 3)
#' )
#' screen_items(items)
#'
#' @export
screen_items <- function(data,
                         digits = 3,
                         floor_cutoff = 0.80,
                         ceiling_cutoff = 0.80,
                         min_item_rest = -0.10) {
  items <- .pm_prepare_numeric_items(data)
  item_names <- names(items)

  stats_rows <- lapply(item_names, function(item) {
    x <- items[[item]]
    observed <- x[!is.na(x)]
    n_obs <- length(observed)
    item_sd <- stats::sd(observed)
    item_min <- if (n_obs > 0L) min(observed) else NA_real_
    item_max <- if (n_obs > 0L) max(observed) else NA_real_
    item_mean <- if (n_obs > 0L) mean(observed) else NA_real_

    rest_names <- setdiff(item_names, item)
    rest_data <- items[rest_names]
    rest_data <- rest_data[vapply(rest_data, .pm_has_variation, logical(1))]
    item_rest <- NA_real_
    if (.pm_has_variation(x) && length(rest_data) > 0L) {
      rest_score <- rowMeans(rest_data, na.rm = TRUE)
      rest_score[!is.finite(rest_score)] <- NA_real_
      item_rest <- suppressWarnings(stats::cor(x, rest_score, use = "complete.obs"))
    }

    floor_pct <- if (n_obs > 0L) mean(observed == item_min) else NA_real_
    ceiling_pct <- if (n_obs > 0L) mean(observed == item_max) else NA_real_

    data.frame(
      item = item,
      n = n_obs,
      missing = sum(is.na(x)),
      missing_pct = mean(is.na(x)),
      mean = item_mean,
      sd = item_sd,
      min = item_min,
      max = item_max,
      skewness = .pm_skewness(observed),
      kurtosis = .pm_kurtosis(observed),
      n_categories = length(unique(observed)),
      zero_variance = !.pm_has_variation(x),
      floor_pct = floor_pct,
      ceiling_pct = ceiling_pct,
      floor_flag = isTRUE(!is.na(floor_pct) && floor_pct >= floor_cutoff),
      ceiling_flag = isTRUE(!is.na(ceiling_pct) && ceiling_pct >= ceiling_cutoff),
      item_rest_correlation = item_rest,
      probable_reverse = isTRUE(!is.na(item_rest) && item_rest < min_item_rest),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, stats_rows)
  numeric_cols <- vapply(out, is.numeric, logical(1))
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) round(x, digits))
  class(out) <- c("psychomatic_item_screen", "data.frame")
  out
}

#' Reverse-score selected items
#'
#' Reverses item scores using either a user-defined score range or the observed
#' range of the selected items. Items can be selected explicitly or detected
#' from negative item-rest correlations.
#'
#' @param data Data frame or matrix containing item responses.
#' @param items Character vector of item names to reverse. Required unless
#'   `detect = TRUE`.
#' @param min_score,max_score Optional theoretical score range. If omitted, the
#'   observed range among selected items is used.
#' @param suffix Suffix used for reversed columns when `overwrite = FALSE`.
#' @param overwrite If `TRUE`, replace the original columns. If `FALSE`, append
#'   reversed columns.
#' @param detect If `TRUE`, reverse items flagged by `screen_items()` as probable
#'   reverse-keyed items.
#'
#' @return A data frame with attributes `reversed_items` and `reverse_range`.
#'
#' @examples
#' items <- data.frame(item1 = 1:5, item2 = 5:1)
#' reverse_items(items, items = "item2", min_score = 1, max_score = 5)
#'
#' @export
reverse_items <- function(data,
                          items = NULL,
                          min_score = NULL,
                          max_score = NULL,
                          suffix = "_r",
                          overwrite = FALSE,
                          detect = FALSE) {
  out <- as.data.frame(data)
  numeric_items <- .pm_prepare_numeric_items(out)

  selected <- character(0)
  if (!is.null(items)) {
    selected <- as.character(items)
  }
  if (isTRUE(detect)) {
    screened <- screen_items(numeric_items)
    selected <- unique(c(selected, screened$item[screened$probable_reverse]))
  }

  if (length(selected) == 0L) {
    stop("No items were selected for reverse scoring.", call. = FALSE)
  }

  missing_items <- setdiff(selected, names(out))
  if (length(missing_items) > 0L) {
    stop("Items not found in `data`: ", paste(missing_items, collapse = ", "), call. = FALSE)
  }
  non_numeric <- selected[!vapply(out[selected], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop("Items selected for reverse scoring must be numeric: ",
         paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  observed <- unlist(out[selected], use.names = FALSE)
  observed <- observed[!is.na(observed)]
  if (length(observed) == 0L && (is.null(min_score) || is.null(max_score))) {
    stop("`min_score` and `max_score` are required when selected items are all missing.",
         call. = FALSE)
  }
  if (is.null(min_score)) min_score <- min(observed)
  if (is.null(max_score)) max_score <- max(observed)
  if (!is.numeric(min_score) || !is.numeric(max_score) || min_score >= max_score) {
    stop("`min_score` and `max_score` must define an increasing numeric range.",
         call. = FALSE)
  }

  for (item in selected) {
    reversed <- min_score + max_score - out[[item]]
    if (isTRUE(overwrite)) {
      out[[item]] <- reversed
    } else {
      out[[paste0(item, suffix)]] <- reversed
    }
  }

  attr(out, "reversed_items") <- selected
  attr(out, "reverse_range") <- c(min = min_score, max = max_score)
  out
}

#' Score one or more psychometric scales
#'
#' Computes scale scores using sum or mean scoring, optional reverse scoring,
#' optional person-mean imputation, and a minimum valid-response rule.
#'
#' @param data Data frame or matrix containing item responses.
#' @param keys Named list where each element contains the item names for one
#'   scale. If `NULL`, `items` are scored as a single scale named `total`.
#' @param items Optional item names used when `keys = NULL`.
#' @param method Scoring method: `"mean"` or `"sum"`.
#' @param min_valid Minimum number or proportion of valid item responses required
#'   to compute a scale score.
#' @param impute Missing-data handling: `"none"` or `"person_mean"`.
#' @param reverse Optional character vector of items to reverse before scoring.
#' @param min_score,max_score Optional score range used when `reverse` is provided.
#'
#' @return A list of class `psychomatic_scores` with score and reliability tables.
#'
#' @examples
#' items <- data.frame(item1 = c(1, 2, NA), item2 = c(2, 3, 4), item3 = c(3, 4, 5))
#' score_scale(items, items = names(items), min_valid = 2 / 3)
#'
#' @export
score_scale <- function(data,
                        keys = NULL,
                        items = NULL,
                        method = c("mean", "sum"),
                        min_valid = 0.80,
                        impute = c("none", "person_mean"),
                        reverse = NULL,
                        min_score = NULL,
                        max_score = NULL) {
  method <- match.arg(method)
  impute <- match.arg(impute)
  item_data <- as.data.frame(data)

  if (!is.null(reverse)) {
    item_data <- reverse_items(
      item_data,
      items = reverse,
      min_score = min_score,
      max_score = max_score,
      overwrite = TRUE
    )
  }

  if (is.null(keys)) {
    if (is.null(items)) {
      items <- names(.pm_prepare_numeric_items(item_data))
    }
    keys <- list(total = as.character(items))
  }
  if (!is.list(keys) || length(keys) == 0L || is.null(names(keys)) ||
      any(!nzchar(names(keys)))) {
    stop("`keys` must be a named list of item names.", call. = FALSE)
  }

  scores <- data.frame(row = seq_len(nrow(item_data)))
  reliability <- lapply(names(keys), function(scale) {
    scale_items <- as.character(keys[[scale]])
    missing_items <- setdiff(scale_items, names(item_data))
    if (length(missing_items) > 0L) {
      stop("Items not found for scale `", scale, "`: ",
           paste(missing_items, collapse = ", "), call. = FALSE)
    }
    if (!all(vapply(item_data[scale_items], is.numeric, logical(1)))) {
      stop("All scoring items must be numeric.", call. = FALSE)
    }

    values <- as.matrix(item_data[scale_items])
    valid_count <- rowSums(!is.na(values))
    required <- .pm_required_valid_items(length(scale_items), min_valid)
    eligible <- valid_count >= required

    scored_values <- values
    if (impute == "person_mean") {
      row_mean <- rowMeans(scored_values, na.rm = TRUE)
      for (i in seq_len(nrow(scored_values))) {
        if (eligible[i] && anyNA(scored_values[i, ])) {
          scored_values[i, is.na(scored_values[i, ])] <- row_mean[i]
        }
      }
    }

    score <- if (method == "mean") {
      rowMeans(scored_values, na.rm = FALSE)
    } else {
      rowSums(scored_values, na.rm = FALSE)
    }
    score[!eligible] <- NA_real_
    scores[[scale]] <<- score

    alpha_value <- .pm_alpha(item_data[scale_items])
    data.frame(
      scale = scale,
      n_items = length(scale_items),
      min_valid = required,
      alpha = alpha_value,
      mean_score = mean(score, na.rm = TRUE),
      sd_score = stats::sd(score, na.rm = TRUE),
      n_scored = sum(!is.na(score)),
      stringsAsFactors = FALSE
    )
  })

  scores$row <- NULL
  reliability <- do.call(rbind, reliability)
  rownames(reliability) <- NULL

  result <- list(
    scores = scores,
    reliability = reliability,
    keys = keys,
    method = method,
    impute = impute,
    min_valid = min_valid,
    reversed_items = reverse,
    call = match.call()
  )
  class(result) <- "psychomatic_scores"
  result
}

.pm_prepare_numeric_items <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }
  out <- as.data.frame(data)
  numeric_cols <- vapply(out, is.numeric, logical(1))
  if (!any(numeric_cols)) {
    stop("`data` must contain at least one numeric item.", call. = FALSE)
  }
  out[numeric_cols]
}

.pm_has_variation <- function(x) {
  observed <- x[!is.na(x)]
  length(unique(observed)) > 1L
}

.pm_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || stats::sd(x) == 0) return(NA_real_)
  z <- (x - mean(x)) / stats::sd(x)
  mean(z^3)
}

.pm_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L || stats::sd(x) == 0) return(NA_real_)
  z <- (x - mean(x)) / stats::sd(x)
  mean(z^4) - 3
}

.pm_required_valid_items <- function(n_items, min_valid) {
  if (!is.numeric(min_valid) || length(min_valid) != 1L || is.na(min_valid) ||
      min_valid <= 0) {
    stop("`min_valid` must be a positive number.", call. = FALSE)
  }
  if (min_valid <= 1) {
    max(1L, ceiling(n_items * min_valid))
  } else {
    if (min_valid > n_items) {
      stop("`min_valid` cannot exceed the number of items in a scale.",
           call. = FALSE)
    }
    as.integer(ceiling(min_valid))
  }
}

.pm_alpha <- function(data) {
  if (ncol(data) < 2L || any(!vapply(data, .pm_has_variation, logical(1)))) {
    return(NA_real_)
  }
  tryCatch(
    suppressWarnings(psych::alpha(data, check.keys = FALSE)$total$raw_alpha),
    error = function(e) NA_real_
  )
}
