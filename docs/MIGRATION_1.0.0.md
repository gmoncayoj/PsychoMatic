# Migrating from PsychoMatic 0.3.0 to 1.0.0

The 0.4.0 development series is incorporated into 1.0.0. This version identifies
the project's first stable API; the version number alone does not establish
that a release has been published on CRAN.

## Existing code

Previous public functions and compatible positional arguments are retained.
Review results before reusing an archived analysis: estimator selection,
ordinal-data preparation, diagnostics and presentation have received
substantial corrections documented in NEWS.

- In CFA and complete invariance, `result$warnings` is a structured table.
  Legacy code expecting text can use `result$warnings_list`, a transitional
  alias retained during the deprecation cycle.
- New results generally inherit from `psychomatic_result`, separating
  requested/effective settings, results, diagnostics and decisions.
  Display precision does not alter stored values.
- Factorial invariance and the ordinal procedure with Wu–Estabrook
  identification have different APIs. Choose the procedure appropriate for
  your model and explain any parameter releases in partial invariance.
- Demonstration datasets remain development fixtures, outside the installed
  API. Help examples are templates for researchers' own data; vignettes
  provide reproducible demonstrations.
- Automatic reports are created in `getwd()` with unique names unless an
  explicit path is supplied. `overwrite = TRUE` permits replacing existing files.

## Added functions

`psychomatic_spec()` declares structure and coding; `reliability_auto()`
centralizes reliability analysis; `ordinal_invariance_auto()` and
`invariance_partial()` extend invariance workflows; `network_auto()` and
`network_plots()` integrate network analysis and figures.

`esem_auto()` fits an exploratory measurement block with lavaan. For geomin,
select columns using `items = 2:13` or column names; do not pass a data frame
as `items`. Target rotation requires a theoretically specified list of primary
item memberships in `factors`. Both rotations are oblique; global fit does
not determine which rotation to choose. `seed = NULL` preserves the engine's
random-number behavior; an explicit integer seed enables reproducibility and
restores the previous RNG state. The final interpretation section depends on
the results and their admissibility and requires substantive review before
inclusion in a manuscript.

## Verification before updating a project

Save your code, R/dependency versions and previous results. First install the
candidate in an isolated library, run your usual analyses and compare data
preparation, samples, estimators and the fit indices actually used.
Complete EN/ES reference manuals are available in [manuals/](manuals/README.md),
and summarized candidate evidence is in [VALIDATION.md](VALIDATION.md).
