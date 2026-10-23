# Reproducible tests

The distribution suite retains input validation, decisions, native model
fitting, result equivalence and reporting checks. Extensive regressions remain
in the same test tree and are enabled explicitly. No statistical engine or
its defaults are changed to speed up checks.

The extended selection includes exhaustive ESEM missing-data matrices, a large
matrix of partial-invariance variants, native inadmissible/nonconvergent
solutions, repeated complete exports, and network integrations with more
replications, seed repetition and all figures. Their test bodies are retained.
The ordinary selection retains continuous and ordinal ESEM target/geomin
native equivalence, invariance, unevaluable states, bilingual reports and a
native network smoke test. Known failures are not skipped.

## Distribution suite

```r
Sys.unsetenv("PSYCHOMATIC_EXTENDED_TESTS")
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```

The tarball runs this same selection through `R CMD check --as-cran`. Skips
identified as `Extended regression` are deliberate and are distinguished from
skips caused by missing dependencies.

## Complete suite

```r
Sys.setenv(PSYCHOMATIC_EXTENDED_TESTS = "true")
pkgload::load_all(".")
testthat::test_dir("tests/testthat")
```

The value must be exactly `true`. It does not depend on `interactive()`,
`NOT_CRAN`, the R version or computer speed. The workflow
`.github/workflows/R-CMD-check.yaml` sets it explicitly to retain complete
native regressions on Linux, Windows and macOS. Configuring a workflow does
not demonstrate that it has run on GitHub.

## Evidence and limitations

The [validation summary](VALIDATION.md) identifies the checked GPL candidate.
The maintainer retains per-case profiles, test-selection policy checks,
reverse-dependency reviews, HTTP verification and check results locally.
Each platform's archive check must identify the tarball SHA256; a pending
submission does not count as a successful check. Reviewing the CRAN index
does not detect private dependencies or external repositories.

Compact smoke tests reopen Word and Excel reports in both languages and
check numerical precision. They also run network bootstrapping and
permutations with native engines and few replications. This small budget
verifies software integration; it is not a recommended configuration for
analyzing research data.
