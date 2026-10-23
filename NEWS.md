# PsychoMatic 1.0.0

## Release scope

- Change the package license from MIT to GPL version 3 or later
  (`GPL (>= 3)`), copyright 2026 José Gamarra-Moncayo.

- Consolidates the unreleased 0.4.0 development series into the first stable
  1.0.0 API. The preceding published release is 0.3.0.
- Complete English and Spanish reference PDFs are generated reproducibly from
  every documented topic; editorial sources stay outside the source tarball.
- See `docs/MIGRATION_1.0.0.md` for warning-table compatibility and workflows.

- `esem_auto()` now defaults to `seed = NULL`, using the session RNG unless
  the investigator supplies a reproducible seed. Its console omits preparation,
  internal maps and versions while preserving them in objects and exports.
  References precede the final Interpretation section in both languages;
  the generic closing manuscript sentence was removed.
  Negative degrees of freedom independently block interpretation, even without
  a native warning; zero degrees of freedom remain a separate saturation case.

- Add `esem_auto()` for a single measurement block with target or geomin oblique
  rotation, native rotated inference, input mappings and admissibility gates.
  Raise minimum lavaan to 0.7-3. Provide descriptive CFA comparison, bilingual
  console/Word/Excel reports and a final numerical manuscript draft. Documentation
  explains rotation choice, soft targets, epsilon and inference limitations.


- `invariance_partial()` now defaults to the same native lavaan settings,
  item preparation and metric/scalar/strict constraints as
  `factorial_invariance_auto()`. Explicit WLSMV/ULSMV default to ordered
  indicators, marker identification and delta parameterization. Existing
  positional arguments are preserved.
- Added `fit_variant` (scaled by default, robust on request), `missing` and
  `min_group_n`. Variant selection propagates to tables, comparisons,
  admissibility gates, interpretations and bilingual console/Word/Excel reports;
  unavailable robust indices have no fallback.
- Identified threshold and latent-response partial invariance remains available
  through explicit `workflow = "ordinal"`. semTools is required only for this
  workflow. Researcher rationale, retained anchors, effective release and exact
  matched-prerequisite nesting checks remain in place.

- Simplify traditional invariance console and reports: display observed response
  options instead of the measurement-level label and omit the item coding table.
  Raw coding metadata and all model estimation settings remain available.

- Reproduce native lavaan calls in the traditional four-stage workflow.
  `factorial_invariance_auto()` now always fits only configural, metric, scalar
  and strict models, constraining loadings, intercepts and residuals. Explicit
  WLSMV/ULSMV default to `ordered = TRUE`, marker identification and delta
  parameterization, matching native CFA defaults. Automatic selection uses
  WLSMV for detected ordinal responses and ML/MLR otherwise. Explicit ML/MLR/ULS
  analyze numeric scores. Ordinal comparisons without positive additional df
  remain unevaluable. Use `ordinal_invariance_auto()` for identified ordinal stages.
- Add `fit_variant = c("scaled", "robust")` to both invariance APIs. Scaled is
  the default; explicit robust selection governs CFI/RMSEA, changes, decisions
  and reports. Corrected chi-square/df/p remain scaled; SRMR is base. Requested
  unavailable robust indices are not silently substituted.

- Add `ordinal_invariance_auto()` with exactly configural, equal thresholds,
  and equal thresholds and loadings stages, using Wu-Estabrook identification.
  Restrict estimators to WLSMV/ULSMV, mark untestable threshold stages explicitly,
  and provide bilingual reader tables and Word/Excel reports with APA references
  to Wu and Estabrook (2016), Svetina et al. (2020), and Chen (2007).
- Preserve signed fit changes and use the loading-specific SRMR guideline when
  adding loading equality to the ordinal sequence. Chen cutoffs remain
  descriptive extrapolations for ordinal outcomes.

This is a local development candidate until the release checks and publication
process are completed.

## Default report destination

- Save reports in the current R working directory (`getwd()`) when no output
  path is supplied, including per-analysis reports, integrated reports and
  EFA/CFA exports. Generate unique timestamped names and share one base name
  for Word/Excel pairs. Preserve explicit paths and overwrite protection.

## Network manuscript presentation and APA typography

- Enlarge fixed node labels and stability values; add `label_cex` for plots and
  exports. Add explicit `model`, `algorithm`, and `bootstrap_type` choices to
  the network workflow, with accurate methods and TMFG edge captions.
  Defaults remain GLASSO, consensus Louvain and resampling. BGGM availability
  depends on EGAnet (disabled in 2.4.1); never silently substitute an estimator.
- Harmonize community colors with distinct blue, coral, green, purple and gold
  fills; use opaque node/point backgrounds and matching legends. Keep positive
  edges neutral and negative edges dashed. Loading maps add a signed color key
  and contrasting black/white cell labels without changing estimates.
- Include signed cross-loadings beneath primary loadings and a complete
  item stability table beneath community structural consistency.
- Add combined EGA/item-stability panels, colored replication dot plots and
  soft community fills. Two-group networks mark complete BH loading tests,
  retaining unavailable tests and exploratory configural interpretation.
- Classify loading magnitude using unrounded values and explain rounded
  boundary values in reader output.
- Default network output emphasizes community membership, bootstrap structural
  consistency, centrality stability and separate NCT/loading interpretations.
- Retain full technical conditions in the object and `print(..., detail="audit")`;
  translate unavailable assessments and inadmissible CFA into scientific caveats.
- Display only BH-significant loading differences in the compact item table;
  missing tests and unsupported configural structure remain explicit.
- Share partial APA reference italics across analyses in compatible consoles,
  Word and Excel, while preserving plain reference strings and numeric results.
- Select relevant publication figures by default; retain all diagnostic types
  on request. Add community legends, consistent colors, signed edges and
  BH-aware loading difference plots. Separate TEFI comparison sets into panels.

## Network runtime and correlation corrections

- Add ordinal `corr="polychoric"` to network estimation and resampling.
- Show stage/substage progress, elapsed seconds and optional-stage interruptions.
- Avoid unused consensus Louvain in fixed-partition loading permutations and
  respect requested consensus iterations in hierarchical bootstrap.
- Serialize a self-contained bootnet estimator for Windows workers, independent
  of other installed PsychoMatic versions.
- Distinguish network hierarchy/TEFI from supplementary CFA and explain negative
  variances, WLSMV AIC/BIC and second-order components with three factors.

## Psychometric networks (2026-10-03)

- Make UVA explicitly optional: NULL (default) skips it, `"diagnostic"` reports
  candidates and reasons without removal, and `"removed"` authorizes justified
  native reduction. Remove interactive prompts in every execution environment.
  Record proposed/actual removals, representatives, scores, ties and the decision
  in console and bilingual Word/Excel reports.
- Show regularized partial correlations on EGA edges and item labels in nodes.
  Plot primary revised loadings on item-dimension edges; preserve all cross-loadings
  in the separate loading matrix figure.
- Accept integer-valued seeds in EGAnet resampling and permutation calls; use
  session RNG by default (`seed=NULL`) and preserve caller RNG with explicit seeds.
- Show only UVA pairs meeting the declared cutoff; keep full pairs and package
  versions as internal diagnostics. Explain unavailable RI partitions and skip
  their structural bootstrap without blocking main stability or invariance.
- Attribute warnings/errors to their analysis stages and omit empty reasons.
- Add `network_auto()` with optional pooled UVA, Spearman/EBICglasso consensus Louvain EGA,
  revised network loadings, centrality and entropy diagnostics.
- Preserve original unrecoded responses for riEGA; evaluate identifiable correlated,
  second-order and bifactor candidates with explicit exploratory recommendations.
- Add bootEGA structural consistency and bootnet edge precision/case-drop stability.
- Compare all observed group pairs using NCT and EGA loading permutations, with
  full-precision raw p values and four declared BH families.
- Add cached publication figures, bilingual Word/Excel reports and an APA 7
  methodological reference registry. No participant rows appear in reader reports.

## EFA omega consistency (2026-09-26)

- Use the shared standalone reliability engine for EFA omega per dimension,
  including measurement-level detection and inadmissible-model checks.
- Explain the independent one-factor calculation in console and reports;
  retain CFA reliability from the user-specified joint model unchanged.

## Model-comparison estimator and index correction (2026-09-25)

- Default model comparisons to automatic estimator selection using the same
  selection helper as CFA, with an explicit estimator override preserved.
- Display estimator, selection source and rationale. Use scaled CFI/TLI/RMSEA
  for corrected estimators in both the comparison and reader output; retain
  base SRMR and do not silently replace unavailable scaled values.
- Allow explicit `fit_variant = "scaled"`, `"robust"`, or `"base"` for both
  the displayed comparison and recommendation inputs.
- Suppress dependency startup messages while loading MVN, without suppressing
  model warnings or errors.

## Coefficient interpretation and help revision (2026-09-23)

- Show descriptive adequacy and explicit cutoffs for reliability (.70) and AVE
  (.50), preserving computational failures and the original numeric results.
- Simplify item-analysis sample-size display; retain per-fit sample sizes internally.
- Link numbered CFA suggestions to a dedicated notes section.
- Complete Spanish invariance labels and missing-data warnings; include the
  Hirschfeld and von Brachel DOI.
- Replace programming guards and simulations in public help examples with
  adaptable user-data syntax.

## Reader-output corrections

- Avoid duplicate console titles; add bilingual direct console and Excel/Word
  reporting to `screen_items`, preserving its raw data-frame columns.
- Shorten measurement-invariance output without changing decisions or diagnostics.
- Present alignment tables with compact headers, separate flagged groups and
  explicit group names, latent means and latent standard deviations.

## Language and statistical reporting revision (2026-09-22)

- Public analysis functions now default to English; Spanish is explicit.
- `cormat()` defaults to automatic selection using category count and observed
  distribution shape, explains the decision, names the method, and cites sources.
- CFA computes Mardia regardless of estimator, distinguishing observed ordinal
  codes from latent-response assumptions. Alpha is labelled Alpha (Cronbach).
- Invariance uses scaled indices for corrected estimators (base SRMR), with no
  silent robust fallback, and consistently computes current-minus-previous deltas.
- Configural gates and change limits are visible. Removed automatic classification
  of unknown Chen scenarios; ordinal extrapolation and threshold identification
  are explained with the supplied Hirschfeld/von Brachel and Bowen/Masa sources.

## Minimal-interface revision (2026-09-17)

- Removed public `check_items()`, `dimensionality_sensitivity()` and
  `audit_report()` entry points; validators and audit metadata remain internal.
- Moved demo datasets to test fixtures, removed their data catalog entries,
  and made installed examples self-contained. S3 methods remain internal and
  are hidden from the regular help index.
- Added automatic congeneric omega, corrected item-rest correlations, and
  reliability after item deletion. Ordinal omega targets observed scores;
  ordinal alpha targets latent responses. Non-evaluable results retain reasons.
- Added observed and ordinal alpha per eligible CFA factor/group. Clarified
  standardized loading intervals, omitted redundant hierarchy columns, explained
  why Mardia is not applied to ordinal models, localized inadmissible conclusions,
  and retained original lavaan fit warnings and vcov diagnostics.
- Simplified EFA variance percentages and correlation output. Invariance labels
  use delta symbols with explicit subtraction direction; variant priority is
  unchanged (robust, scaled, base) and is disclosed in the report.
- The earlier development notes below describe the history, not the current
  public interface. No publication or remote Git update is implied.

## New workflows

- Added `psychomatic_spec()` to declare item columns, measurement levels,
  ordered-category mappings, reverse keys, ranges, grouping, metadata, and
  missing-data intent in one reusable object.
- Added `check_items()` with stable issue codes for missing/non-finite values,
  variation, sparse or empty categories, names, ranges, pairwise overlap, and
  correlation-matrix validity.
- Added `reliability_auto()` with explicitly labelled observed-score,
  ordinal, two-item, CFA-based, second-order, and bifactor estimands where the
  supplied data/model support them.
- Added `dimensionality_sensitivity()` to compare parallel-analysis retention
  results across correlation methods, quantiles, and seeds without declaring an
  automatic winner.
- Added researcher-directed `invariance_partial()`. It requires exact parameter
  releases and a substantive rationale, protects anchor items, and never runs
  an automatic modification-index search.
- Added `audit_report()` to extract serializable settings, provenance,
  conditions, diagnostics, decisions, references, result tables, and report
  artifacts without refitting models.

## Statistical and implementation corrections

- Unified item classification across specifications, checks, correlations,
  EFA, reliability, and ordinal model workflows. Ordered-factor levels are
  authoritative; ambiguous nominal encodings fail closed.
- Corrected reverse-key deduplication, append-name collisions, non-finite input
  validation, and incomplete-row proration for mean and sum scale scores.
- Correlation matrices now retain requested/effective methods, pairwise sample
  sizes, and structured diagnostics. No positive-definite repair occurs unless
  explicitly requested where supported.
- Descriptive and analytical fields retain full-precision raw values; `digits`
  is restricted to display and export formatting.
- Reworked EFA iteration so every item-set change is followed by a complete
  refit and reclassification. Parallel analysis can retain zero factors, uses
  the intended factor-eigenvalue reference, restores caller RNG state, and
  avoids categorical “optimal solution” language.
- Updated CFA normality integration for `MVN >= 6.3`; retained base, scaled, and
  robust fit evidence; and separated requested from effective estimators.
- Corrected two-item reliability to preserve the sign of the correlation.
  CFA reliability now uses current public `semTools` APIs and full Lambda, Psi,
  and Theta matrices, with group/factor-specific non-evaluable states.
- Strengthened CFA admissibility checks, bifactor-structure validation,
  second-order/bifactor reliability mapping, modification-syntax preservation,
  and structured condition capture.
- Rebuilt complete invariance around one frozen analytic sample and explicit
  fit/admissibility/comparability gates. Corrected fit-change directions and
  records the effective estimator, missing policy, parameterization,
  identification, and equality families at every stage.
- Ordinal complete and partial invariance use identification-aware semTools
  syntax. Chen thresholds are exposed as scenario-dependent heuristics, never
  as formal proof of invariance.
- Alignment now matches matrices and weights by semantic identities, rejects
  incomplete/ambiguous names and protected argument overrides, preserves
  missingness diagnostics, captures engine output/conditions, and fails closed
  for malformed, non-finite, or nonconverged results.
- Automatic binary IRT now selects 2PL; Rasch is explicit. Mixed item banks use
  item-specific types, raw slopes are preserved, and only applicable parameter
  cells enter admissibility checks.
- DIF matching scores exclude the focal item and undeclared metadata, freeze
  one analysis sample per comparison, and fail closed on separation,
  nonconvergence, non-nesting, or degenerate probabilities. Polytomous items now
  default to `not_evaluable` rather than silent median dichotomization.
- CFA model comparison preserves positional candidate identity, effective
  estimator/estimand metadata, common fit-index variants, and comparability
  gates. AIC/BIC are suppressed outside the unique comparable ranking set.
- Report and export adapters now preserve diagnostic data frames, sanitize
  unsafe live model objects recursively, avoid duplicate table expansion,
  validate all formats/targets before writing, retain worksheet-name mappings,
  and apply Excel number formats without changing numeric cell values.

## API, objects, and compatibility

- Introduced the common `psychomatic_result` envelope: call, requested and
  effective settings, data summary, results, structured warnings/errors,
  diagnostics, decision log, and conclusion.
- Added or completed concise `print()`, structured `summary()`, and selected
  base `plot()` methods for the new workflows and invariance results.
- Existing public function names and valid 0.3.x argument positions are
  retained. New controls were added as trailing arguments where applicable.
- CFA and complete-invariance `$warnings` are now structured tables. Use the
  transitional `$warnings_list` character alias for 0.3.x code during the
  deprecation cycle.
- Single-group second-order reliability retains its historical factor-indexed
  view. Multigroup output uses group-then-factor identity; canonical long rows
  are available in `$reliability$rows`.
- Spanish argument aliases in EFA remain available for compatibility; the
  documented English names are preferred for new code.

## Data, examples, and release infrastructure

- Moved demonstration datasets to development test fixtures; they are not
  exported datasets. Preserved reproducible generators and simulation metadata
  for numerical regression checks, including `likert_demo` and `invariance_demo`.
- Added twelve installed example scripts covering the principal 1.0.0
  workflows without personal paths, attached packages, network access, or
  repository writes.
- Expanded tests with numerical contrasts against public `psych`, `lavaan`,
  `semTools`, `MVN`, `mirt`, and `sirt` APIs where applicable.
- Added a cross-platform R CMD check workflow for Ubuntu, Windows, and macOS.

## Reader-facing usability


* Analysis functions display reader-facing results directly, including assigned calls.
  `verbose = FALSE` supports silent workflows; `digits` controls presentation only.
* Per-analysis `report`, `format`, `file_name`, and `overwrite` arguments provide
  Excel/Word reports without a separate exporter. Integrated reports default to
  analysis sections; `detail = "audit"` preserves technical evidence tables.
* EFA restores comprehensive diagnostics, fit, variance, loadings, item allocation,
  reliability, and references. CFA defaults to explicitly labelled scaled indices;
  robust and base estimates remain available.
* Reliability explains score scope and calculation, separating AVE from reliability.
  Invariance retains non-evaluable stages and actual engine errors visibly.
* Corrected the Wu and Estabrook (2016) reference and formatted references as prose.

# PsychoMatic 0.3.0

## CRAN readiness

- Updated `DESCRIPTION`, maintainer metadata, vignette configuration, and
  DOI-backed methodological references.
- Moved report-only and workflow-specific packages to conditional suggested
  dependencies.
- Added English manual pages, lightweight vignettes, bundled example datasets,
  and temporary-path report examples.

## New functions

- Added `screen_items()`, `reverse_items()`, `score_scale()`,
  `psychomatic_report()`, `irt_auto()`, `dif_auto()`, and
  `compare_models_auto()`.
- Added structured summaries for EFA, CFA, and complete invariance.
- Aligned EFA export with the integrated report writer and retained bilingual
  output through `language = "esp"` and `language = "eng"`.

# PsychoMatic 0.2.0

- Added Pearson, polychoric, and tetrachoric inter-item correlation matrices.
- Added Excel export for correlation matrices.
- Corrected the confidence-interval label for factor loadings.
