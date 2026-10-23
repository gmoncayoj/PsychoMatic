# PsychoMatic

<img width="1536" height="1024" alt="Package image" src="docs/assets/Package_image.png" />

## Automated psychometric workflows with minimal syntax

`PsychoMatic` helps applied researchers move from item checks to scoring,
factor analysis, reliability, measurement invariance, IRT/DIF screening, and
reproducible reports without assembling long procedural scripts. Automation is
used to reduce syntax, not to replace measurement theory or substantive
judgment.

The development tree is currently preparing version 1.0.0. Until that version
is released, `install.packages("PsychoMatic")` installs the public CRAN version;
a locally built `PsychoMatic_1.0.0.tar.gz` can be installed explicitly for
candidate testing.

The complete English and Spanish reference manuals, including every public
function and documented method, are in the repository at `docs/manuals/README.md`.
See `docs/MIGRATION_1.0.0.md` in the repository before updating code
written for 0.3.x.

## One call: analysis, console results, optional report

With `report=TRUE` and no `file_name`, reports are saved in the current R
working directory (`getwd()`), using unique timestamped names. Word and Excel
from the same call share a base name. Explicit paths in `file_name` (or in
`report` where that argument accepts a path) are respected. Set the R working directory to the desired project folder
before analysis; PsychoMatic does not change it or redirect reports to a
temporary directory.

Psychometric networks are available through `network_auto()`: UVA, consensus
Louvain EGA with Spearman/EBICglasso (or Pearson/polychoric correlations), revised loadings, TEFI, riEGA, hierEGA,
bootEGA/bootnet, and group comparisons with NCT and EGA loading permutations.

```r
net <- network_auto(MY_DATA, items = 1:12, group = "sex", seed = 2026)
plot(net, type = "network")
plot(net, type = "loadings") # primary loadings on item-dimension edges
network_plots(net, "network_figures", formats = c("pdf", "svg", "png"))
psychomatic_report(net, format = c("word", "excel"), file_name = "Network")
```

The default network output emphasizes community membership, bootstrap structural
consistency and separate group comparisons of structure, global strength and
item loadings. Signed cross-loadings follow the primary loading table, and all
item assignment stabilities follow community structural consistency. Only
BH-significant loading differences enter the compact item
table; unavailable tests and scientific limitations remain explicit. Full
diagnostics are available with `print(net, detail="audit")` or
`psychomatic_report(net, detail="audit")`.

Reference typography is shared across analyses: journal names and volumes, and
the book titles used by the package, appear in italics in Word, Excel and
compatible interactive consoles. Logs retain plain text. Use
`options(PsychoMatic.ansi=FALSE)` to disable console styling, or `TRUE` to force
ANSI italics when supported by the terminal and font.

Network reports and `network_plots(net, directory)` select relevant scientific
figures by default. Unavailable diagnostics remain accessible with explicit
types, e.g. `types=network_plots(net)`. Figures use community legends, shared
group positions and weight scales, and BH-aware item loading comparisons.
With bootstrap results, `ega_stability` combines the EGA network and a colored
replication dot plot. Standalone `network` and `item_stability` remain available:
`network_plots(net, "network_figures", types=c("ega_stability", "item_stability"))`.
Community fills are opaque and use the same distinct palette in networks,
replication plots and legends. Positive edges are neutral; negative edges are
dashed. Loading maps retain a signed scale and choose contrasting cell text.
Node labels use a fixed larger font, without automatic shrinking. Increase it
for a reduced publication width with `network_plots(net, "figures", label_cex=1.5)`.
Stability values and markers are also larger.

`network_auto()` accepts `model="glasso"` (default), `"TMFG"`, or `"BGGM"`,
`algorithm="louvain"` (default), `"leiden"`, or `"walktrap"`, and
`bootstrap_type="resampling"` (default) or `"parametric"`.
The selected estimator is used in bootstraps and comparisons; community
detection is used in empirical/group EGA and both hierarchy orders. Loading
permutations keep the original partition fixed and avoid unused consensus.
EGAnet 2.4.1 disables BGGM internally: requesting it produces the native error,
never a silent fallback. TMFG uses filtered correlations, not GLASSO partial
correlations; gamma applies only to GLASSO. UVA uses its separate native
redundancy estimator. Parametric bootEGA simulates continuous normal responses
and does not preserve ordinal categories. No automatic best-bootstrap selection
is claimed; bootnet nonparametric edge and case-drop bootstraps are separate.
Native limitations remain explicit: TMFG hierarchy orders need at least four
nodes, and EGAnet can fail in riEGA bootstraps for some selected options. These
optional stages retain explanatory errors without replacing the user's methods.
For a single pair of groups, thick orange node borders flag complete loading
tests with BH p < .05; gray borders indicate no difference detected and squares
mark unavailable tests. These marks do not establish equivalence or configural
invariance, and remain exploratory when configural support is absent.
PDF/SVG are vector formats; PNG defaults to 300 dpi. Check the target journal's
required dimensions and resolution before submission.

Optional dependencies must be installed before analysis. Use `language="esp"`
for Spanish. UVA is optional: `uva=NULL` (default) skips redundancy analysis,
`uva="diagnostic"` reports candidate pairs, wTO and removal reasons while keeping
all items, and `uva="removed"` authorizes justified native item removal.
There are no interactive questions.
Explicit reverse keys require theoretical `ranges`; riEGA uses
original responses. Defaults use 1000 bootstrap replicates and 5000 permutations.
See `vignette("network-workflow", package="PsychoMatic")` for assumptions, BH
families, exploratory hierarchy recommendations and interpretation limits.

```r
library(PsychoMatic)
continuous_items <- lavaan::simulateData(
  "general =~ 0.7*item1 + 0.7*item2 + 0.7*item3 + 0.7*item4 + 0.7*item5 + 0.7*item6",
  sample.nobs = 200, seed = 2401
)
desc_auto(continuous_items, digits = 3)
efa_auto(continuous_items, seed = 123)

# The same call also writes a reader-facing workbook or Word report.
efa_auto(continuous_items, seed = 123, report = TRUE,
         format = "excel", file_name = "Exploratory-analysis.xlsx")
```

Analyses print their results directly, even when assigned to an object. No
translation argument is needed for English; use `language = "esp"` for Spanish.
`cormat()` defaults to `type = "auto"`, using categories and distribution-shape
flags and reporting the selected method, rationale, and references. No
`print()` or `summary()` is required. Use `verbose = FALSE` for silent scripted
work. `digits` controls presentation only; stored estimates retain precision.
Use `format = "word"` or `c("excel", "word")` for other output formats, and
`overwrite = TRUE` only when replacement is intended. EFA keeps its historical
item-removal default; `refine = "suggest"` retains all items for researcher review.

Reader reports group results, decisions, non-evaluable outcomes, and references
by analysis. `psychomatic_report(..., detail = "audit")` is optional for full
technical metadata. Specifications and exporters remain available for advanced
workflows, but are not prerequisites for an analysis.

CFA displays scaled fit indices by default, with explicit variant labels.
`fit_variant = "robust"` requests robust indices. Both variants and the base
estimates remain in the returned object; they are not interchangeable.

Measurement invariance defaults to scaled indices for corrected estimators and
base indices for uncorrected estimators (SRMR is base). `fit_variant = "robust"`
selects robust CFI/RMSEA for tables, changes and decisions, without a silent fallback.
Every delta is current minus previous. Ordinal threshold stages retain the
Wu--Estabrook identification sequence; the report explains why some other
protocols incorporate thresholds at the scalar stage instead.

## Optional explicit data specification

```r
set.seed(2401)
latent <- matrix(stats::rnorm(300 * 3), ncol = 3)
ordinal_items <- as.data.frame(sapply(seq_len(12), function(j) {
  as.integer(cut(0.8 * latent[, ceiling(j / 4)] + stats::rnorm(300, sd = 0.6),
                 breaks = c(-Inf, -1, -0.3, 0.3, 1, Inf), labels = FALSE))
}))
names(ordinal_items) <- sprintf("item%02d", 1:12)
attr(ordinal_items, "reverse_keys") <- c("item03", "item07", "item11")
ordinal_items[attr(ordinal_items, "reverse_keys")] <-
  lapply(ordinal_items[attr(ordinal_items, "reverse_keys")], function(x) 6L - x)
attr(ordinal_items, "factor_structure") <- list(
  factor1 = names(ordinal_items)[1:4], factor2 = names(ordinal_items)[5:8],
  factor3 = names(ordinal_items)[9:12]
)

specification <- PsychoMatic::psychomatic_spec(
  ordinal_items,
  level = "ordinal",
  reverse = attr(ordinal_items, "reverse_keys"),
  ranges = stats::setNames(rep(list(c(1, 5)), 12), names(ordinal_items)),
  language = "eng"
)

scores <- PsychoMatic::score_scale(
  specification,
  keys = attr(ordinal_items, "factor_structure"),
  reverse = attr(ordinal_items, "reverse_keys"),
  min_score = 1,
  max_score = 5
)

utils::head(scores$scores)
```

`psychomatic_spec()` makes item selection, measurement level, reverse keys,
ranges, grouping, and missing-data intent explicit once. Downstream functions
can consume the specification directly.

## Public workflows

| Task | Main function |
| --- | --- |
| Data contract | `psychomatic_spec()` |
| Descriptives and correlations | `desc_auto()`, `cormat()` |
| Screening, reverse scoring, scale scores | `screen_items()`, `reverse_items()`, `score_scale()` |
| Reliability | `reliability_auto()` |
| EFA | `efa_auto()`, `export_efa()` |
| ESEM measurement model | `esem_auto()` (target/geomin, native lavaan, manuscript draft) |
| CFA and model comparison | `cfa_auto()`, `compare_models_auto()`, `export_cfa()` |
| Complete measurement invariance | `factorial_invariance_auto()` |
| Ordinal threshold and loading invariance | `ordinal_invariance_auto()` |
| Researcher-directed partial invariance | `invariance_partial()` |
| Alignment optimization | `inv_align_auto()` |
| IRT and DIF screening | `irt_auto()`, `dif_auto()` |
| Reports | `psychomatic_report()` |

## Reproducible examples

`esem_auto()` estimates one exploratory measurement block with freely estimated
cross-loadings (lavaan >= 0.7-3). Use target rotation when primary item dimensions
are theoretically specified; use geomin when that organization is uncertain.
Both are oblique; rotation choice cannot be justified by global fit alone.
See `vignette("esem-workflow")` and `inst/examples/EXAMPLE_ESEM_AUTO.R` for
executable ordinal target and continuous geomin examples. Console, Word and
Excel end with a data-based manuscript draft, also available as
`result$manuscript_text`; failed/inadmissible fits receive limitation statements.

Examples simulate their inputs with explicit seeds. The package does not
install demonstration datasets; use your own study data for substantive work.

## Reliability and dimensionality

```r
items <- PsychoMatic::reverse_items(
  ordinal_items,
  items = attr(ordinal_items, "reverse_keys"),
  min_score = 1,
  max_score = 5,
  overwrite = TRUE
)

reliability <- PsychoMatic::reliability_auto(
  items,
  items = attr(ordinal_items, "factor_structure")$factor1,
  level = "ordinal",
  language = "eng"
)

```

Here reliability is estimated for the first simulated domain rather than for
the known three-domain item set as a whole. `reliability_auto()` labels the
estimand and method for every coefficient. An alpha or omega value is not
evidence of unidimensionality, tau-equivalence, or validity.
Inspect parallel-analysis stability in `efa_auto()` and increase simulation
iterations for final analyses.

## EFA and CFA

```r
efa <- PsychoMatic::efa_auto(
  items,
  refine = "suggest",
  parallel_iter = 100,
  seed = 2401,
  verbose = FALSE,
  language = "eng"
)

continuous_items <- lavaan::simulateData(
  "general =~ 0.7*item1 + 0.7*item2 + 0.7*item3 + 0.7*item4 + 0.7*item5 + 0.7*item6",
  sample.nobs = 200, seed = 2401
)
model <- "general =~ item1 + item2 + item3 + item4 + item5 + item6"
cfa <- PsychoMatic::cfa_auto(
  continuous_items,
  model = model,
  estimator = "ML",
  language = "eng"
)
```

Use `refine = "suggest"` when item removal has not been preregistered or
substantively justified. Automated loading and fit cutoffs are screening
heuristics, not formal decision rules. CFA output records the requested and
effective estimator where the engine distinguishes them.

## Three distinct invariance workflows

Complete, partial, and alignment invariance answer related but different
questions and should not be treated as interchangeable stages of one automatic
procedure.

### Complete invariance

```r
grouped_items <- lavaan::simulateData(
  "general =~ 0.7*mg1 + 0.7*mg2 + 0.7*mg3 + 0.7*mg4 + 0.7*mg5 + 0.7*mg6",
  sample.nobs = 400, seed = 2402
)
grouped_items$group <- factor(rep(c("A", "B"), each = 200))
model <- "general =~ mg1 + mg2 + mg3 + mg4 + mg5 + mg6"

complete <- PsychoMatic::factorial_invariance_auto(
  grouped_items,
  group = "group",
  model = model,
  estimator = "ML",
  language = "eng"
)
```

The function freezes one analytic sample and fits native lavaan configural,
metric, scalar and strict calls. Restrictions accumulate on loadings, intercepts
and residual variances. Explicit WLSMV/ULSMV default to `ordered = TRUE`,
marker identification, delta parameterization and listwise handling. Automatic
selection uses WLSMV for detected ordinal responses and ML/MLR for continuous
responses according to group normality. Explicit ML/MLR/ULS analyze numeric
scores; `ordered = FALSE` explicitly requests numeric-score estimation.
No separate threshold-equality stage is fitted here. Native ordinal scalar/strict
calls may not add independent restrictions; those comparisons are marked
unevaluable rather than declared invariant. Use `ordinal_invariance_auto()` for
the identified three-stage threshold/loading workflow. Chen-style changes are
descriptive heuristics; a supported stage is not proof of invariance.

The reader summary shows the number of observed response options instead of
the measurement-level label. Different counts are listed when they vary across
items. The item coding table is omitted from console and reports; coding remains
available in `effective_settings$score_coding`. An observed count does not infer
response categories that were never present in the analytic data. To reproduce
WLSMV results in native lavaan, use `ordered = TRUE`, lavaan's default marker
identification and delta parameterization, and the same missing-data policy,
sample and equality constraints. `fit_variant = "scaled"` is the default;
`fit_variant = "robust"` selects robust CFI/RMSEA, including their changes and
decisions. Chi-square/df/p remain scaled for corrected estimators; SRMR is base.
Unavailable requested robust indices remain missing.

```r
fit <- factorial_invariance_auto(MY_DATA, group = "Sexo", model = RESIL,
  estimator = "WLSMV", fit_variant = "scaled", language = "esp")
# To select robust CFI and RMSEA instead:
fit_robust <- factorial_invariance_auto(MY_DATA, group = "Sexo", model = RESIL,
  estimator = "WLSMV", fit_variant = "robust", language = "esp")
```

### Ordinal threshold and loading invariance

For a three-stage ordinal analysis, use `ordinal_invariance_auto()`:

```r
ordinal_fit <- PsychoMatic::ordinal_invariance_auto(
  MY_DATA, group = "Sexo", model = "F =~ R1 + R2 + R3 + R4",
  ordered = TRUE, estimator = "WLSMV", language = "eng",
  report = TRUE, format = c("word", "excel")
)
```

Its stages are configural, equal thresholds, and equal thresholds and loadings.
`"auto"` selects WLSMV; ULSMV is also accepted. ML/MLR and nonordinal
indicators are rejected. The identification follows Wu and Estabrook (2016)
and the tutorial of Svetina, Rutkowski, and Rutkowski (2020), DOI
`10.1080/10705511.2019.1602776` (published online in 2019, volume dated 2020).
The table includes chi-square, df, CFI, RMSEA, SRMR and signed changes
(current minus preceding evaluable model). Chen (2007) cutoffs are descriptive
extrapolations to ordinal outcomes, not proof of invariance. Binary and
three-category indicators may make the separate threshold test unidentifiable;
such stages are explicitly marked unevaluable. This three-stage sequence does
not establish equal latent-response intercepts or residuals.

### Researcher-directed partial invariance

```r
set.seed(2402)
latent <- matrix(stats::rnorm(600 * 3), ncol = 3)
group <- factor(rep(c("A", "B"), each = 300))
ordinal_groups <- as.data.frame(sapply(seq_len(9), function(j) {
  response <- 0.8 * latent[, ceiling(j / 3)] + stats::rnorm(600, sd = 0.6)
  values <- integer(600)
  for (g in levels(group)) {
    thresholds <- c(-Inf, -1, -0.3, 0.3, 1, Inf)
    if (j == 3 && g == "B") thresholds[4] <- 0.8
    values[group == g] <- cut(response[group == g], thresholds, labels = FALSE)
  }
  values
}))
names(ordinal_groups) <- paste0("i", 1:9)
ordinal_groups[] <- lapply(ordinal_groups, ordered, levels = 1:5)
ordinal_groups$group <- group
ordinal_model <- paste(
  "factor1 =~ i1 + i2 + i3",
  "factor2 =~ i4 + i5 + i6",
  "factor3 =~ i7 + i8 + i9",
  sep = "\n"
)

if (requireNamespace("semTools", quietly = TRUE)) {
  partial <- PsychoMatic::invariance_partial(
    ordinal_groups,
    group = "group",
    model = ordinal_model,
    free = "i3 | t3",
    rationale = "The simulation preregisters this threshold difference.",
    stage = "thresholds",
    workflow = "ordinal",
    measurement = "ordinal",
    estimator = "WLSMV",
    ordered = paste0("i", 1:9),
    language = "eng"
  )
}
```

`invariance_partial()` never searches modification indices or frees parameters
automatically. The researcher must name every release and provide a substantive
rationale. The function checks equality families, effective releases, nesting,
model admissibility, and at least two remaining anchor items per affected
factor by default.

Its default `workflow = "factorial"` uses the same native lavaan preparation,
marker identification, constraints and fit-index selection as
`factorial_invariance_auto()`. WLSMV/ULSMV default to ordered indicators;
`fit_variant = "scaled"` is the default and `"robust"` can be requested.
Use `stage = "metric"`, `"scalar"`, or `"strict"` and specify releases in
`free` (for example, `"factor1 =~ i2"`). Releases from lower stages also apply
to the matched prerequisite. The threshold example above explicitly chooses
`workflow = "ordinal"` to retain identified semTools estimation. Unavailable
indices, ineffective releases and nonnested comparisons cannot support partial
invariance. Both workflows accept `missing` and `min_group_n`.

### Alignment

`inv_align_auto()` accepts item-level data or parameter matrices and requires
the optional `sirt` package. Alignment is especially relevant when conventional
exact-invariance sequences become impractical across many groups. Its R-squared
and percentage summaries remain approximate diagnostics, not authorization to
compare latent means without further scrutiny.

## Ordered items, missing data, and non-evaluable results

- Supply ordered factors, `ordered =`, or an explicit `psychomatic_spec()`;
  numeric category codes are not always self-explanatory.
- Inspect requested and effective missing-data settings and analytic case
  counts. Different estimators support different missing-data policies.
- Ordinal models record parameterization and identification choices. Binary
  and polytomous invariance stages are not mechanically identical to continuous
  metric/scalar stages.
- A status such as `not_evaluable` is an analytical result. It must not be
  recoded as “no problem,” “no DIF,” or support for invariance.
- `dif_auto()` is a lightweight binary logistic screen. Polytomous items default
  to `not_evaluable`; explicit median dichotomization is a compatibility option,
  not a recommended substitute for ordinal or IRT DIF methods.
- `irt_auto()` automates a bounded set of `mirt` models and retains convergence,
  admissibility, category maps, and raw coefficients. It does not automate item
  calibration design, local-dependence studies, or linking.

## Structured results, audit, and reporting

Core 1.0.0 workflows return `psychomatic_result` objects with requested and
effective settings, data summaries, raw results, structured warnings/errors,
diagnostics, decision logs, and cautious conclusions. Analytical values retain
full precision (raw values); `digits` affects display or exported formatting.

```r

if (requireNamespace("openxlsx", quietly = TRUE)) {
  path <- tempfile("psychomatic-report-", fileext = ".xlsx")
  PsychoMatic::psychomatic_report(
    list(reliability = reliability, efa = efa),
    format = "excel",
    file_name = path,
    overwrite = FALSE
  )
}
```

Report writers preflight all target paths and refuse to overwrite by
default. Use `overwrite = TRUE` only after deliberately selecting an existing
file.

## Migration from 0.3.x

- CFA and complete-invariance `$warnings` are structured condition tables.
  Transitional `$warnings_list` aliases preserve the former character view.
- EFA, CFA, invariance, IRT, DIF, and reporting retain compatible legacy fields
  where practical. Console results require no extra calls; `summary()` and
  `$tables` remain optional for programmatic access.
- Raw numerical fields are no longer rounded for presentation. Apply explicit
  rounding only when displaying a result.
- `dif_auto()` no longer silently median-dichotomizes polytomous items.
- Automatic binary IRT selects 2PL; Rasch must be requested explicitly.
- Report paths no longer overwrite existing files unless explicitly allowed.

Ten complete, noninteractive workflows are installed under
`inst/examples/` in the source package.

## Development tests

The source package runs a bounded distribution suite. Exhaustive native failure
probes, missing-data combinations, and large report/permutation regressions are
retained in an explicitly enabled extended suite:

```r
Sys.setenv(PSYCHOMATIC_EXTENDED_TESTS = "true")
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```

Unset `PSYCHOMATIC_EXTENDED_TESTS` to run the distribution suite. The CI workflow
enables the complete suite on every configured platform. See
`docs/TESTING.md` in the repository for the scope and evidence.

## License and support

PsychoMatic is licensed under the GNU General Public License, version 3 or later
(GPL-3.0-or-later). Copyright (c) 2026 José Gamarra-Moncayo.
See https://www.gnu.org/licenses/gpl-3.0.html and `LICENSE.note`.
Please report reproducible problems through
the [GitHub issue tracker](https://github.com/gmoncayoj/PsychoMatic/issues).
