# PsychoMatic

<img width="1536" height="1024" alt="Package image" src="https://github.com/user-attachments/assets/77d9c2d8-31b8-4680-961b-5f260a71e85f" />

### Automated psychometric analysis in R

`PsychoMatic` is an R package that automates psychometric workflows using SEM-based procedures and literature-informed heuristics. It covers descriptives, EFA, CFA, reliability, measurement invariance, and alignment, with flexible settings, bilingual reporting, and APA 7-style references to support analytical decisions.

## What PsychoMatic can do?

PsychoMatic makes the best decisions for you so you can run your analyses with confidence!  
However, you can also customize it according to your own specifications; in any case, it will provide you with the final report in the language of your choice. It was initially designed to promote its use in Spanish, but since most of the literature worldwide is in English, it is also capable of providing the report for each analysis in that language.

| Analysis | Function |
|--------|-------------|
| Descriptives | `desc_auto()` |
| Inter-item correlations | `cormat()` |
| EFA + reliability | `efa_auto()`, `export_efa()` |
| CFA + reliability | `cfa_auto()`, `export_cfa()` |
| Factorial invariance | `factorial_invariance_auto()` |
| Alignment invariance | `inv_align_auto()` |

## Development status

PsychoMatic is currently in active development and is open to suggestions for improvement.

## Installation

```r
install.packages("remotes")
remotes::install_github("gmoncayoj/PsychoMatic")
```

## Descriptive analysis of items

### Main arguments

| Argument   | Description | Default  |
|------------|-------------|----------|
| `data`     | A data frame containing the variables to analyze. Only numeric columns are processed. | n/a |
| `digits`   | Number of decimal places for rounding descriptive statistics (mean, SD, skewness, kurtosis). Percentages always use at least 2 decimals. | `2` |
| `language` | Output language for column labels and messages. Options: `"esp"` (Spanish) or `"eng"` (English). | `"esp"` |
| `report`   | If `TRUE`, exports the results to a timestamped `.xlsx` file with two sheets: descriptives and a notes sheet with missing value count. Requires the `writexl` package. | `FALSE` |

## Inter-item correlation matrix

### Main arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Data frame or matrix containing the items | n/a |
| `type` | Correlation type: `"pearson"`, `"poly"` for polychoric, or `"tetra"` for tetrachoric | `"pearson"` |
| `report` | If `TRUE`, exports the matrix to Excel | `FALSE` |
| `file_name` | Optional Excel file name | `NULL` |
| `digits` | Number of decimal places used in the Excel report | `3` |

## Exploratory factor analysis + reliability

### Main arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Data frame or matrix containing the items (numeric variables) | n/a |
| `rotation` | Rotation type: `"oblique"` (oblimin) or `"orthogonal"` (varimax) | `"oblique"` |
| `min_loading` | Minimum acceptable factor loading; items below this threshold are flagged or removed | `0.30` |
| `min_communality` | Minimum acceptable communality (`h2`); items below this threshold are flagged or removed | `0.30` |
| `min_cross_loading_diff` | Minimum difference between the two highest loadings to avoid cross-loading flags | `0.15` |
| `max_iter` | Maximum number of iterative refinement cycles | `20` |
| `verbose` | If `TRUE`, prints the full analysis report to the console | `TRUE` |
| `export_format` | Optional export format for results (for example, `"excel"` or `"word"`) | `NULL` |
| `file_name` | Optional base name for exported file(s) | `NULL` |
| `language` | Report language: `"esp"` (Spanish) or `"eng"` (English) | `"esp"` |

Previous Spanish argument aliases are still accepted temporarily for backward compatibility, but the English names above are now the recommended interface.

## Confirmatory factor analysis + reliability

### Main arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Data frame containing item responses | n/a |
| `model` | Model specified using `lavaan` syntax | n/a |
| `ordered` | Vector of ordinal items (or `TRUE` for all items) | `NULL` |
| `estimator` | Estimator (`"ML"`, `"MLR"`, `"WLSMV"`, etc.). If `NULL`, it is selected automatically | `NULL` |
| `mi_threshold` | Minimum threshold for reporting modification indices | `10` |
| `language` | Report language: `"esp"` or `"eng"` | `"esp"` |

Use `export_cfa()` to export the resulting CFA object to Excel or Word with an English API aligned with the rest of the package.

## Factorial invariance

### Main arguments

| Argument | Description | Default value |
|----------|-------------|---------------|
| `data` | Data frame containing the items and the grouping variable | n/a |
| `group` | Name of the grouping variable (string or unquoted symbol) | n/a |
| `model` | CFA model specified in `lavaan` syntax using `=~` | n/a |
| `language` | Report language: `"esp"` (Spanish) or `"eng"` (English) | `"esp"` |
| `estimator` | Estimator to use: `"auto"`, `"ML"`, `"MLR"`, `"ULS"`, or `"WLSMV"` | `"auto"` |
| `ordered` | Override ordinal treatment: `TRUE`, `FALSE`, or `NULL` (auto-detected from data) | `NULL` |
| `report` | Export to Excel: `FALSE` (none), `TRUE` (default filename), or a custom file path string | `FALSE` |

## Alignment invariance

### Main arguments

| Argument | Description | Default value |
|----------|-------------|---------------|
| `data` | Data frame containing only numeric item columns, one row per respondent. Use together with `group` | `NULL` |
| `group` | Vector of group labels with one value per row in `data` | `NULL` |
| `lambda` | Groups x items matrix of pre-estimated factor loadings. Use together with `nu` | `NULL` |
| `nu` | Groups x items matrix of pre-estimated intercepts. Use together with `lambda` | `NULL` |
| `language` | Report language: `"esp"` (Spanish) or `"eng"` (English) | `"esp"` |
| `model` | Configural CFA model when `data`/`group` are supplied: `"2PM"` (2-parameter) or `"1PM"` (1-parameter) | `"2PM"` |
| `align.scale` | Named numeric vector `c(lambda, nu)` controlling scale penalty in the alignment loss function | `c(lambda = 0.40, nu = 0.20)` |
| `align.pow` | Named numeric vector `c(lambda, nu)` controlling power penalty in the alignment loss function | `c(lambda = 0.25, nu = 0.25)` |
| `parm_tol` | Tolerance used by `sirt::invariance_alignment_constraints()` to classify parameters as noninvariant; defaults to `align.scale` | Same as `align.scale` |
| `noninvariance_cutoff` | Percentage threshold above which parameters (loadings or intercepts) are flagged as noninvariant | `25` |
| `sampling_weights` | Optional numeric vector of sampling weights, one per row in `data` | `NULL` |
| `wgt` | Optional weight matrix (groups x items) or per-group vector passed to `sirt::invariance.alignment()` | `NULL` |
| `digits` | Number of decimal places used in the printed report | `3` |
| `config_args` | Named list of additional arguments passed to `sirt::invariance_alignment_cfa_config()` | `list()` |
| `alignment_args` | Named list of additional arguments passed to `sirt::invariance.alignment()` | `list()` |

## Author

Jose Gamarra-Moncayo  
Psychology professor, Faculty of Medicine, Universidad Catolica Santo Toribio de Mogrovejo  
Email: gamarramoncayoj@gmail.com  
Please feel free to contact the author to offer suggestions and/or report any bugs in the package.

## License

MIT
