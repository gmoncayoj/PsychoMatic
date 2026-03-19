# ⚡PsychoMatic

<img width="1536" height="1024" alt="Package image" src="https://github.com/user-attachments/assets/77d9c2d8-31b8-4680-961b-5f260a71e85f" />

### 🚀 Automated psychometric analysis in R

`PsychoMatic` is an R package that automates psychometric workflows using SEM-based procedures and literature-informed heuristics. It covers descriptives, EFA, CFA, reliability, measurement invariance, and alignment, with flexible settings, bilingual reporting, and APA 7–style references to support analytical decisions.

## 🧠 What PsychoMatic can do?

PsychoMatic makes the best decisions for you so you can run your analyses with confidence!  
However, you can also customize it according to your own specifications; in any case, it will provide you with the final report in the language of your choice. It was initially designed to promote its use in Spanish, but since most of the literature worldwide is in English, it is also capable of providing the report for each analysis in that language.

| Analysis | Function |
|--------|-------------|
| Descriptives | `desc_auto()` |
| EFA + reliability | `efa_auto()` |
| CFA + reliability | `cfa_auto` |
| Invariance | `factorial_invariance_auto()` |
| Alignment | `inv_align_auto()` |

## 🚧 Development status

PsychoMatic is currently in active development and is open to suggestions for improvement.

## ⚙️ Installation

```r
install.packages("remotes")
remotes::install_github("gmoncayoj/PsychoMatic")
```

## 📊 Example of univariate descriptive analysis of items

```r

library(PsychoMatic)

### Basic use (Spanish, no report)
desc_auto(base)

### With English labels
desc_auto(base, language = "eng") # If you want Spanish, enter “esp”

### Adjust decimal places
desc_auto(base, digits = 2)

### Export report to Excel
desc_auto(base, report = TRUE, language = "eng")
```

Using `desc_auto()`, you can calculate the mean, standard deviation, skewness, kurtosis, and response rate per item.

## 🔍 Example of exploratory factor analysis + reliability

```r
library(PsychoMatic)

result <- efa_auto(data,
                  rotacion = "oblicua", #You can change the rotation to “orthogonal” if needed
                  language = "eng") #For Spanish, enter “esp”
result
```

---

### Export results

```r
# Export to Excel
exportar_efa(result, formato = "excel", archivo = "results_afe")

# Export to Word
exportar_efa(result, formato = "word", archivo = "results_afe")
```

### Main arguments

| Argument | Description | Default |
|-----------|-------------|---------------|
| `datos` | Data frame or matrix containing the items (numeric variables) | — |
| `rotacion` | Rotation type: `"oblicua"` (oblimin) or `"ortogonal"` (varimax) | `"oblicua"` |
| `carga_min` | Minimum acceptable factor loading; items below this threshold are flagged or removed | `0.30` |
| `comunalidad_min` | Minimum acceptable communality (h²); items below this threshold are flagged or removed | `0.30` |
| `dif_cargas_cruzadas` | Minimum difference between the two highest loadings to avoid cross-loading flags | `0.15` |
| `max_iter` | Maximum number of iterative refinement cycles | `20` |
| `verbose` | If `TRUE`, prints the full analysis report to the console | `TRUE` |
| `exportar` | Optional export format for results (e.g., `"excel"`, `"word"`) | `NULL` |
| `archivo` | Optional base name for exported file(s) | `NULL` |
| `language` | Report language: `"esp"` (Spanish) or `"eng"` (English) | `"esp"` |

## 📐 Example of confirmatory factor analysis + reliability

```r

library(PsychoMatic)

# Define the model in Lavaan syntax
# Unidimensional
model <- 'F =~ i1 + i2 + i3 + i4 + i5'

# Correlated factors
model <- '
  F1 =~ i1 + i2 + i3
  F2 =~ i4 + i5 + i6
'

# Second order
model <- '
  F1 =~ i1 + i2 + i3
  F2 =~ i4 + i5 + i6
  SO  =~ F1 + F2
'

# Bifactor
model <- '
  GEN  =~ i1 + i2 + i3 + i4 + i5 + i6
  F1 =~ i1 + i2 + i3
  F2 =~ i4 + i5 + i6
'

# Run the automated AFC
result <- cfa_auto(data = data, model = model)

# View the complete psychometric report
print(result)
```

💡 PsychoMatic can automatically detect which model you are about to run. If it detects a second-order model, it will report the higher-order omega reliability. If it detects a two-factor model, it will also calculate the hierarchical omega, ECV, and PUC indices, and provide you with a final interpretation of the strength of the general factor!

### Main arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `data` | Data frame containing item responses | — |
| `model` | Model specified using `lavaan` syntax | — |
| `ordered` | Vector of ordinal items (or `TRUE` for all items) | `NULL` |
| `estimator` | Estimator (`"ML"`, `"MLR"`, `"WLSMV"`, etc.). If `NULL`, it is selected automatically | `NULL` |
| `mi_threshold` | Minimum threshold for reporting modification indices | `10` |
| `language` | Report language: `"esp"` or `"eng"` | `"esp"` |

## ⚖️ Example of factorial invariance

```r

library(PsychoMatic)

# Define your CFA model in lavaan syntax
model <- "
  F1 =~ item1 + item2 + item3
  F2 =~ item4 + item5 + item6
"

# Run — automatic estimator selection, output in Spanish (default)
result <- factorial_invariance_auto(
  data  = my_data,
  group = "group_var",
  model = model
)

# Common variants

# English/spanish output + save Excel report to working directory
result <- factorial_invariance_auto(
  data     = my_data,
  group    = "group_var",
  model    = model,
  language = "eng", #If you need the Spanish language, change “esp”
  report   = TRUE
)

# Force a specific estimator
result <- factorial_invariance_auto(
  data      = my_data,
  group     = "group_var",
  model     = model,
  estimator = "MLR"
)

# Save report to a custom path
result <- factorial_invariance_auto(
  data   = my_data,
  group  = "group_var",
  model  = model,
  report = "C:/results/invariance_report.xlsx"
)

```

💡 When `estimator` is set to “auto”, the function automatically selects between WLSMV, ML, and MLR based on the number of response categories and multivariate normality (Mardia's test).

### Main arguments

| Argument | Description | Default value |
|-----------|-------------|---------------|
| `data` | Data frame containing the items and the grouping variable | — |
| `group` | Name of the grouping variable (string or unquoted symbol) | — |
| `model` | CFA model specified in `lavaan` syntax using `=~` | — |
| `language` | Report language: `"esp"` (Spanish) or `"eng"` (English) | `"esp"` |
| `estimator` | Estimator to use: `"auto"`, `"ML"`, `"MLR"`, `"ULS"`, or `"WLSMV"` | `"auto"` |
| `ordered` | Override ordinal treatment: `TRUE`, `FALSE`, or `NULL` (auto-detected from data) | `NULL` |
| `report` | Export to Excel: `FALSE` (none), `TRUE` (default filename), or a custom file path string | `FALSE` |

## 🎯 Example of alignment invariance

```r

library(PsychoMatic)

## Option 1 — From raw item data and a grouping variable

# items_df: data frame with numeric item columns only
# group_var: vector indicating group membership (one value per row)

res <- inv_align_auto(
  data     = data,
  group    = group_var,
  language = "eng"
)

print(res)       # Full alignment report
summary(res)     # Summary table (R2, % noninvariant, invariance verdict)

## Option 2 — From pre-estimated λ and ν matrices

# lambda and nu: groups × items numeric matrices
# (rows = groups, columns = items)

res <- inv_align_auto(
  lambda   = lambda_mat,
  nu       = nu_mat,
  language = "eng"
)

print(res)

```

### Main arguments

| Argument | Description | Default value |
|-----------|-------------|---------------|
| `data` | Data frame containing only numeric item columns, one row per respondent. Use together with `group` | `NULL` |
| `group` | Vector of group labels with one value per row in `data` | `NULL` |
| `lambda` | Groups × items matrix of pre-estimated factor loadings. Use together with `nu` | `NULL` |
| `nu` | Groups × items matrix of pre-estimated intercepts. Use together with `lambda` | `NULL` |
| `language` | Report language: `"esp"` (Spanish) or `"eng"` (English) | `"esp"` |
| `model` | Configural CFA model when `data`/`group` are supplied: `"2PM"` (2-parameter) or `"1PM"` (1-parameter) | `"2PM"` |
| `align.scale` | Named numeric vector `c(lambda, nu)` controlling scale penalty in the alignment loss function | `c(lambda = 0.40, nu = 0.20)` |
| `align.pow` | Named numeric vector `c(lambda, nu)` controlling power penalty in the alignment loss function | `c(lambda = 0.25, nu = 0.25)` |
| `parm_tol` | Tolerance used by `sirt::invariance_alignment_constraints()` to classify parameters as noninvariant; defaults to `align.scale` | Same as `align.scale` |
| `noninvariance_cutoff` | Percentage threshold above which parameters (loadings or intercepts) are flagged as noninvariant | `25` |
| `sampling_weights` | Optional numeric vector of sampling weights, one per row in `data` | `NULL` |
| `wgt` | Optional weight matrix (groups × items) or per-group vector passed to `sirt::invariance.alignment()` | `NULL` |
| `digits` | Number of decimal places used in the printed report | `3` |
| `config_args` | Named list of additional arguments passed to `sirt::invariance_alignment_cfa_config()` | `list()` |
| `alignment_args` | Named list of additional arguments passed to `sirt::invariance.alignment()` | `list()` |

## 👤 Author

Jose Gamarra Moncayo  
Email: gamarramoncayoj@gmail.com  
Please feel free to contact the author to offer suggestions and/or report any bugs in the package.

## 📄 License

MIT