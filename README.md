# PsychoMatic

<img width="1536" height="1024" alt="Package image" src="https://github.com/user-attachments/assets/77d9c2d8-31b8-4680-961b-5f260a71e85f" />

## Automated psychometric analysis with minimal syntax

`PsychoMatic` is an R package for automated psychometric workflows. Its main
goal is to let applied researchers run common psychometric analyses with the
least possible syntax, while still returning transparent results, structured
tables, bilingual reports, and literature-informed decision support.

The package is designed for workflows where the researcher wants to move from
item screening to factor analysis, reliability, scoring, measurement
invariance, IRT/DIF screening with interpretation, and reporting without
writing long procedural code for each step.

## Design Philosophy

PsychoMatic prioritizes:

- Minimal syntax for common psychometric analyses.
- Conservative automated defaults based on psychometric practice.
- English-facing API and documentation for CRAN-readiness.
- Bilingual output through `language = "esp"` or `language = "eng"`.
- Structured return objects that can be inspected, tested, summarized, and
  exported.
- Optional advanced arguments for researchers who need more control.

## Installation

```r
install.packages("remotes")
remotes::install_github("gmoncayoj/PsychoMatic")
```

## Minimal Workflow

```r
library(PsychoMatic)

data(psychomatic_ordinal)
data(psychomatic_continuous)

# Screen item quality
screen_items(psychomatic_ordinal)

# Compute a total score with minimal syntax
scores <- score_scale(psychomatic_ordinal)

# Run automated EFA
efa_result <- efa_auto(psychomatic_ordinal, verbose = FALSE, language = "eng")
summary(efa_result)

# Compare CFA models
models <- list(
  one_factor = "general =~ item1 + item2 + item3 + item4 + item5 + item6",
  two_factor = "
    factor1 =~ item1 + item2 + item3
    factor2 =~ item4 + item5 + item6
  "
)
compare_models_auto(psychomatic_continuous, models = models)
```

## Core Functions

| Task | Function |
| --- | --- |
| Item descriptives | `desc_auto()` |
| Item screening | `screen_items()` |
| Reverse scoring | `reverse_items()` |
| Scale scoring | `score_scale()` |
| Inter-item correlations | `cormat()` |
| Exploratory factor analysis | `efa_auto()`, `export_efa()` |
| Confirmatory factor analysis | `cfa_auto()`, `export_cfa()` |
| Measurement invariance | `factorial_invariance_auto()` |
| Alignment invariance | `inv_align_auto()` |
| IRT modeling and item interpretation | `irt_auto()` |
| DIF screening and interpretation | `dif_auto()` |
| CFA model comparison | `compare_models_auto()` |
| Integrated reports | `psychomatic_report()` |

## Example Datasets

`PsychoMatic` includes small simulated datasets for reproducible examples:

| Dataset | Purpose |
| --- | --- |
| `psychomatic_continuous` | Continuous item responses for CFA examples. |
| `psychomatic_ordinal` | Five-category ordinal item responses for EFA, scoring, and IRT examples. |
| `psychomatic_multigroup` | Multigroup item responses for invariance workflows. |

## Descriptives And Item Screening

```r
data(psychomatic_ordinal)

desc_auto(psychomatic_ordinal, language = "eng")
screen_items(psychomatic_ordinal)
```

## Scoring

```r
keys <- list(
  total = names(psychomatic_ordinal)
)

score_scale(
  psychomatic_ordinal,
  keys = keys,
  method = "mean",
  min_valid = 0.80
)
```

## Exploratory Factor Analysis

```r
efa_result <- efa_auto(
  psychomatic_ordinal,
  rotation = "oblique",
  verbose = FALSE,
  language = "eng"
)

summary(efa_result)
```

## Confirmatory Factor Analysis

```r
data(psychomatic_continuous)

model <- "
factor1 =~ item1 + item2 + item3
factor2 =~ item4 + item5 + item6
"

cfa_result <- cfa_auto(
  psychomatic_continuous,
  model = model,
  estimator = "ML",
  language = "eng"
)

summary(cfa_result)
```

## Measurement Invariance

```r
data(psychomatic_multigroup)

model <- "
factor1 =~ mg1 + mg2 + mg3
factor2 =~ mg4 + mg5 + mg6
"

invariance <- factorial_invariance_auto(
  psychomatic_multigroup,
  group = "group",
  model = model,
  estimator = "ML",
  language = "eng",
  report = FALSE
)

summary(invariance)
```

## IRT And DIF Screening

```r
data(psychomatic_ordinal)
data(psychomatic_multigroup)

irt_result <- irt_auto(psychomatic_ordinal[, 1:4], model = "grm", language = "eng")
irt_result$model_summary
irt_result$interpretation
irt_result$references

dif_result <- dif_auto(
  psychomatic_multigroup[paste0("mg", 1:6)],
  group = psychomatic_multigroup$group,
  language = "eng"
)
dif_result$summary
dif_result$interpretation
dif_result$references
```

## Reporting

Reports are explicit actions. In examples and tests, use temporary files so the
package does not write uncontrolled files to the working directory.

```r
report_path <- tempfile("psychomatic-report-", fileext = ".xlsx")

psychomatic_report(
  screening = screen_items(psychomatic_ordinal),
  scores = score_scale(psychomatic_ordinal),
  format = "excel",
  file_name = report_path,
  language = "eng"
)
```

## CRAN-Readiness

Version 0.3.0 focuses on CRAN-readiness:

- English documentation and manual pages.
- Smaller, focused R files for EFA/CFA helpers, methods, and exports.
- Optional dependencies loaded only when needed.
- Executable tests for success and error paths.
- Minimal vignettes using bundled datasets.
- Clean `R CMD check --as-cran` on the release candidate tarball.

## Development Status

`PsychoMatic` is in active development. The package aims to support applied
psychometric research while keeping the interface simple enough for teaching,
graduate research, and reproducible analytical reports.

## Author

Jose Gamarra-Moncayo  
Psychology professor, Faculty of Medicine, Universidad Catolica Santo Toribio de Mogrovejo  
Email: gamarramoncayoj@gmail.com

Suggestions and bug reports are welcome through the GitHub issue tracker.

## License

MIT
