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
#' @param language Output language: `"eng"` (default) or `"esp"`.
#' @param report Generate a report, or supply an output path.
#' @param format Report format: `"excel"`, `"word"`, or both.
#' @param file_name Output file name, without extension when using both formats.
#' @param overwrite Allow replacing an existing report (default `FALSE`).
#' @param verbose Show results directly in the console (default `TRUE`).
#' @details Displayed percentages are on a 0--100 scale. The returned data frame
#'   retains its original English column names and 0--1 proportions for compatibility.
#'
#' @return A data frame of class `psychomatic_item_screen`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' screen_items(MY_DATA[, 3:14])
#' }
#'
#' @export
screen_items <- function(data,
                         digits = 3,
                         floor_cutoff = 0.80,
                         ceiling_cutoff = 0.80,
                         min_item_rest = -0.10,
                         language = "eng", report = FALSE, format = "excel",
                         file_name = NULL, overwrite = FALSE, verbose = TRUE) {
  language <- .pm_normalize_language(language)
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
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
  .pm_finish_analysis(out, output, language)
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
#' @details Repeated selections are applied once. Non-missing responses must be
#'   finite and fall within the effective scoring range. When reversed columns
#'   are appended, an existing target name is treated as a collision and is
#'   never overwritten silently.
#'
#' @return A data frame with attributes `reversed_items` and `reverse_range`.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' reverse_items(MY_DATA, items = c("R2", "R4"), min_score = 1, max_score = 5)
#' }
#'
#' @export
reverse_items <- function(data,
                          items = NULL,
                          min_score = NULL,
                          max_score = NULL,
                          suffix = "_r",
                          overwrite = FALSE,
                          detect = FALSE) {
  out <- .pm_as_data_frame(data)
  numeric_items <- .pm_prepare_numeric_items(out)

  selected <- character(0)
  if (!is.null(items)) {
    selected <- as.character(items)
  }
  if (isTRUE(detect)) {
    screened <- screen_items(numeric_items, verbose = FALSE)
    selected <- unique(c(selected, screened$item[screened$probable_reverse]))
  }
  selected <- unique(selected)

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
  .pm_assert_finite_item_values(out, selected, "reverse scoring")

  if (!is.character(suffix) || length(suffix) != 1L || is.na(suffix)) {
    stop("`suffix` must be a character scalar.", call. = FALSE)
  }
  if (!isTRUE(overwrite)) {
    output_names <- paste0(selected, suffix)
    collisions <- intersect(output_names, names(out))
    if (length(collisions) > 0L) {
      stop(
        "Reverse-scored column already exists: ",
        paste(collisions, collapse = ", "),
        call. = FALSE
      )
    }
  }

  observed <- unlist(out[selected], use.names = FALSE)
  observed <- observed[!is.na(observed)]
  if (length(observed) == 0L && (is.null(min_score) || is.null(max_score))) {
    stop("`min_score` and `max_score` are required when selected items are all missing.",
         call. = FALSE)
  }
  if (is.null(min_score)) min_score <- min(observed)
  if (is.null(max_score)) max_score <- max(observed)
  min_score <- .pm_validate_scalar_number(min_score, "min_score")
  max_score <- .pm_validate_scalar_number(max_score, "max_score")
  if (min_score >= max_score) {
    stop("`min_score` and `max_score` must define an increasing numeric range.",
         call. = FALSE)
  }
  outside_range <- vapply(
    out[selected],
    function(x) any(!is.na(x) & (x < min_score | x > max_score)),
    logical(1)
  )
  if (any(outside_range)) {
    stop(
      sprintf(
        "Observed values outside the reverse-scoring range [%s, %s] for items: %s.",
        min_score,
        max_score,
        paste(selected[outside_range], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (item in selected) {
    reversed <- min_score + max_score - out[[item]]
    if (is.integer(out[[item]]) &&
        all(is.na(reversed) | reversed == trunc(reversed))) {
      reversed <- as.integer(reversed)
    }
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
#' @details Eligible incomplete rows are prorated: mean scores use the available
#'   item mean and sum scores multiply that mean by the scale length. Missing
#'   values may remain `NA`, but non-finite responses are rejected before any
#'   scale is scored. Duplicate display labels receive stable unique keys
#'   recorded in `scale_map`.
#'
#' @return A list of class `psychomatic_scores` with score and reliability
#'   tables and a `scale_map` linking unique keys to display labels.
#'
#' @examples
#' # Replace MY_DATA with your data frame and adapt the item names or columns.
#' \dontrun{
#' score_scale(MY_DATA, items = c("R1", "R2", "R3", "R4"))
#' }
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
  if (inherits(data, "psychomatic_spec")) {
    item_data <- data$data[data$item_index]
    names(item_data) <- data$items
  } else {
    item_data <- .pm_as_data_frame(data)
  }

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
      anyNA(names(keys)) || any(!nzchar(names(keys)))) {
    stop("`keys` must be a named list of item names.", call. = FALSE)
  }

  scale_labels <- names(keys)
  scale_keys <- make.unique(scale_labels, sep = "_")
  normalized_keys <- unname(keys)
  names(normalized_keys) <- scale_keys
  scale_item_sets <- vector("list", length(keys))
  for (i in seq_along(keys)) {
    scale_items <- as.character(keys[[i]])
    missing_items <- setdiff(scale_items, names(item_data))
    if (length(missing_items) > 0L) {
      stop("Items not found for scale `", scale_labels[i], "`: ",
           paste(missing_items, collapse = ", "), call. = FALSE)
    }
    if (!all(vapply(item_data[scale_items], is.numeric, logical(1)))) {
      stop("All scoring items must be numeric.", call. = FALSE)
    }
    scale_item_sets[[i]] <- scale_items
  }
  .pm_assert_finite_item_values(
    item_data,
    unique(unlist(scale_item_sets, use.names = FALSE)),
    "scale scoring"
  )

  score_columns <- vector("list", length(keys))
  reliability <- vector("list", length(keys))

  for (i in seq_along(keys)) {
    scale_key <- scale_keys[i]
    scale_label <- scale_labels[i]
    scale_items <- scale_item_sets[[i]]

    values <- as.matrix(item_data[scale_items])
    valid_count <- rowSums(!is.na(values))
    required <- .pm_required_valid_items(length(scale_items), min_valid)
    eligible <- valid_count >= required

    person_mean <- rowMeans(values, na.rm = TRUE)
    person_mean[valid_count == 0L] <- NA_real_
    score <- if (method == "mean") person_mean else person_mean * length(scale_items)
    score[!eligible] <- NA_real_
    score_columns[[i]] <- score

    alpha_value <- .pm_alpha(item_data[scale_items])
    reliability[[i]] <- data.frame(
      scale = scale_key,
      label = scale_label,
      n_items = length(scale_items),
      min_valid = required,
      alpha = alpha_value,
      mean_score = mean(score, na.rm = TRUE),
      sd_score = stats::sd(score, na.rm = TRUE),
      n_scored = sum(!is.na(score)),
      stringsAsFactors = FALSE
    )
  }

  names(score_columns) <- scale_keys
  scores <- as.data.frame(score_columns, check.names = FALSE, optional = TRUE)
  reliability <- do.call(rbind, reliability)
  rownames(reliability) <- NULL
  scale_map <- data.frame(
    key = scale_keys,
    label = scale_labels,
    stringsAsFactors = FALSE
  )

  result <- list(
    scores = scores,
    reliability = reliability,
    scale_map = scale_map,
    keys = normalized_keys,
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
  if (inherits(data, "psychomatic_spec")) {
    specification <- data
    data <- specification$data[specification$item_index]
    names(data) <- specification$items
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }
  out <- as.data.frame(data, check.names = FALSE, optional = TRUE)
  ordered_cols <- vapply(out, is.ordered, logical(1))
  if (any(ordered_cols)) {
    out[ordered_cols] <- lapply(out[ordered_cols], as.integer)
  }
  numeric_cols <- vapply(out, is.numeric, logical(1))
  if (!any(numeric_cols)) {
    stop("`data` must contain at least one numeric item.", call. = FALSE)
  }
  out[numeric_cols]
}

.pm_assert_finite_item_values <- function(data, items, context) {
  invalid <- vapply(
    data[items],
    function(x) any(is.nan(x) | (!is.na(x) & !.pm_available_values(x))),
    logical(1)
  )
  if (any(invalid)) {
    stop(
      sprintf(
        "All non-missing values used for %s must be finite; check items: %s.",
        context,
        paste(items[invalid], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(data)
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
  tryCatch({
    invisible(utils::capture.output(
      result <- suppressWarnings(suppressMessages(psych::alpha(
        data,
        check.keys = FALSE
      )))
    ))
    result$total$raw_alpha
  }, error = function(e) NA_real_)
}
