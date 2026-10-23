# PsychoMatic

<img src="docs/assets/Package_image.png" alt="PsychoMatic" width="480" />

**Psychometric analysis and bilingual reports in R.**

PsychoMatic 1.0.1 simplifies descriptives, correlations, factor analysis
(EFA, CFA and ESEM), reliability, measurement invariance, IRT/DIF and network
analysis. Results are available in the console and as Word or Excel reports.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes") # Only needed once.
remotes::install_github("gmoncayoj/PsychoMatic")
```

For the published CRAN version, use `install.packages("PsychoMatic")`.

## Quick start

An ESEM analysis with three factors and geomin rotation:

```r
library(PsychoMatic)

result <- esem_auto(
  data = MY_DATA,
  items = 2:13,
  n_factors = 3,
  rotation = "geomin",
  ordered = FALSE
)
```

Replace `MY_DATA` with your data frame and `2:13` with your item columns.
This example treats items as continuous; use `ordered = TRUE` for ordinal items.
Set `language = "esp"` for Spanish output.

Psychometric networks with Word and Excel reports:

```r
network <- network_auto(
  data = MY_DATA,
  items = 2:13,
  seed = 2026,
  report = TRUE,
  format = c("word", "excel"),
  file_name = "Network"
)
```

Reports are saved in your working directory. Optional dependencies are listed
in the [network guide](https://github.com/gmoncayoj/PsychoMatic/blob/main/vignettes/network-workflow.Rmd).

## Documentation

- [Reference manuals · English / Español](https://github.com/gmoncayoj/PsychoMatic/blob/main/docs/manuals/README.md)
- [Guides and examples](https://github.com/gmoncayoj/PsychoMatic/tree/main/docs)
  · [Migration from 0.3.x](https://github.com/gmoncayoj/PsychoMatic/blob/main/docs/MIGRATION_1.0.0.md)

## License and support

GPL version 3 or later. Copyright © 2026 José Gamarra-Moncayo.
[License notice](https://github.com/gmoncayoj/PsychoMatic/blob/main/LICENSE.note)
· [Report an issue](https://github.com/gmoncayoj/PsychoMatic/issues).
