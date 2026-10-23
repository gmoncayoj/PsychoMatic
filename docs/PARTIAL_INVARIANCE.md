# Partial measurement invariance

`invariance_partial()` defaults to `workflow = "factorial"`, matching the
current `factorial_invariance_auto()` preparation, estimator selection, native
lavaan identification, missing-data policy, constraint sequence and fit-index
selection. The shared preparation helper prevents categorical-score conversion
or WLSMV ordered defaults from drifting between functions.

## Native factorial workflow

The complete reference model uses `lavaan::cfa()` with the requested target's
`group.equal`: loadings (metric), loadings plus intercepts (scalar), or loadings,
intercepts and residuals (strict). Its partial counterpart adds the explicitly
requested `group.partial`. WLSMV/ULSMV default to `ordered = TRUE`, marker
identification and native delta parameterization. `ordered = FALSE` requests
numeric-score estimation. ML/MLR/ULS use numeric scores, including explicit
conversion of declared categorical responses, as in complete invariance.

`fit_variant = "scaled"` is the default. For corrected estimators, CFI/RMSEA,
changes, admissibility and decisions use scaled indices. Explicit `"robust"`
uses robust CFI/RMSEA throughout, without replacing unavailable values.
Corrected chi-square, df and p stay scaled; SRMR stays base. Uncorrected
estimators use base indices with the default selection.

```r
partial <- invariance_partial(
  DATA_OK, group = "Sexo", model = RESIL,
  stage = "metric", estimator = "WLSMV",
  free = "INGE =~ R3",
  rationale = "Replace with an independently justified substantive rationale.",
  fit_variant = "scaled", language = "esp"
)
```

This example specifies a release; it is not a recommendation to free R3 in the
user's instrument. Rationale is mandatory. No modification-index search or
automatic release is performed.

## Preserved safeguards

Requested parameters must belong to the model and an active equality family.
The fitted parameter table must confirm their equality in the complete model
and actual release in the partial model. Marker loadings or fixed parameters
cannot be described as released merely because they appear in `free`.

Releases from earlier stages propagate to the matched prerequisite, preserving
exact nesting. At least two fully constrained anchor items per factor remain
required by default. Case indices, group coding, estimator, identification,
ordered indicators and missing-data policy are frozen across comparisons.
Unsupported configural fit, inadmissibility, unavailable indices, ineffective
releases or lack of additional degrees of freedom prevent a supported decision.

Native WLSMV scalar/strict calls can lack independent restrictions. Alignment
with lavaan does not turn those calls into an identified ordinal threshold
invariance test. For that procedure, use the explicit ordinal workflow.

## Identified ordinal option

`workflow = "ordinal"` retains the prior semTools procedure with theta
parameterization, `std.lv` and Wu-Estabrook categorical identification. It
supports threshold releases, loading releases, and the previously supported
latent-response scalar/strict extensions. For binary indicators, identified
joint scalar/strict stages are required. semTools is required for this option;
native partial estimation only requires lavaan. The independent
`ordinal_invariance_auto()` API retains its three-stage complete sequence.

```r
partial_ordinal <- invariance_partial(
  DATA_OK, "Sexo", RESIL, free = "R3 | t2",
  rationale = "Replace with a defensible threshold-specific rationale.",
  workflow = "ordinal", stage = "thresholds", estimator = "WLSMV",
  ordered = TRUE, fit_variant = "scaled", language = "esp"
)
```

The previous public positional arguments are preserved; `fit_variant`,
`workflow`, `missing` and `min_group_n` are appended. Existing threshold calls
must now specify `workflow = "ordinal"` explicitly. Console and readable
Word/Excel reports identify the selected workflow and variant. Audit tables
retain all numerical variants and release diagnostics.

## Sources

Native equality and partial-release syntax follows the official
[lavaan multigroup tutorial](https://lavaan.ugent.be/tutorial/groups.html).
Ordered estimation follows the
[lavaan categorical-data tutorial](https://lavaan.ugent.be/tutorial/cat.html).
The identified ordinal method and its references are documented in
`ORDINAL_INVARIANCE.md`.

## Verification

The regression coverage includes native configural, full metric and partial
metric fit measures, scaled/robust reports and bilingual Word/Excel exports.
See [VALIDATION.md](VALIDATION.md) for the exact GPL archive's complete Windows
release and R-devel checks, including PDF/HTML manuals. [TESTING.md](TESTING.md)
explains the distribution suite and optional intensive regressions.
