# PsychoMatic 0.3.0

## CRAN-Readiness

- Updated `DESCRIPTION` with a more informative title, description, DOI-backed
  methodological references, maintainer metadata, `VignetteBuilder`, and
  `LazyData`.
- Reduced strong dependencies by moving report-only and workflow-specific
  packages to `Suggests` and loading them conditionally.
- Regenerated English-facing manual pages for exported functions.
- Replaced the old `NEWS` file with `NEWS.md`.
- Added lightweight vignettes for EFA, CFA, measurement invariance, and
  alignment invariance.
- Added small bundled datasets for reproducible examples:
  `psychomatic_continuous`, `psychomatic_ordinal`, and
  `psychomatic_multigroup`.
- Ensured examples, tests, and report exports use temporary files where
  appropriate.

## New Functions

- Added `screen_items()` for item-quality screening.
- Added `reverse_items()` for reverse scoring with explicit or detected items.
- Added `score_scale()` for mean/sum scale scores with minimum valid-response
  rules and optional person-mean imputation.
- Added `psychomatic_report()` for integrated Excel or Word reports.
- Added `irt_auto()` for basic Rasch, 2PL, and graded response IRT workflows,
  including item-level interpretation, an overall model summary, and APA 7
  methodological references.
- Added `dif_auto()` for lightweight logistic DIF screening by group, with
  item-level interpretation, an overall DIF summary, and APA 7 methodological
  references.
- Added `compare_models_auto()` for comparing alternative CFA models with
  fit-index-based recommendations and support for WLSMV-style scaled indices.

## Object Structure And Summaries

- Added `summary()` methods for `efa_auto`, `cfa_auto`, and
  `factorial_invariance_auto`.
- Added more predictable structured tables to EFA, CFA, and invariance result
  objects.
- Preserved bilingual output support through `language = "esp"` and
  `language = "eng"` while keeping the main documentation in English.
- Aligned `export_efa()` with `psychomatic_report()` so EFA exports use the
  same structured reporting path as the rest of the package.

## Refactoring

- Split large EFA and CFA files into focused helpers, methods, exports, and
  summary utilities.
- Kept legacy Spanish aliases for backward compatibility, but documented the
  English argument names as the recommended interface.

## Testing

- Expanded tests for item screening, reverse scoring, scale scoring, report
  export, IRT, DIF, CFA model comparison, datasets, summaries, and error paths.
- Verified the release candidate with `R CMD check --as-cran`.

# PsychoMatic 0.2.0

## New Features

- Added inter-item correlation matrices with Pearson, polychoric, and
  tetrachoric options.
- Added Excel export for correlation matrices.

## Fixes

- Fixed the label for confidence intervals in corrected factor loadings.
