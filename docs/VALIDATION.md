# PsychoMatic validation

## Maintenance 1.0.1: CRAN M1 correction

Local preparation began on October 23, 2026. CRAN already publishes 1.0.0;
its additional M1mac report identifies four failures in textual expectations
for equivalent `/var` and `/private/var` paths. The canonical-path corrections
already on main were incorporated into the new 1.0.1 tarball, with additional
regressions for repeated separators and report creation/overwriting.

The 44 statistical-engine files and 47 Rd topics remain identical to CRAN
1.0.0. The API is unchanged. The 1.0.0 artifacts and PDFs are retained as history.
Checked archive: `PsychoMatic_1.0.1.tar.gz`, 452715 bytes; SHA256
`30890BE9873CBBB2878FDA79B98C89D53AFF4C3D221D5339C32F08204D9535AD`.
The complete `R CMD check --as-cran` on Windows 11 / R 4.6.1 finished with
**0 errors, 0 warnings and 1 NOTE**: `Days since last update: 0`.
This administrative NOTE is explained by CRAN requesting a correction to
1.0.0 on its publication day; incoming checks were not disabled.

Tests passed **3415 assertions**, with no failures or warnings and 13 deliberate
intensive-test skips in the distribution profile. None of the four reported
tests was skipped. Examples, vignettes and PDF/HTML manuals completed
successfully. Independent review found no blocking issues.

Complete 1.0.1 manuals: EN 49 pages, ES 51 pages; 47 topics and 22 functions.
SHA256 EN: `229C7E32AEC444DAB818E89A851B8547D1BF87B78C738F91C35C09AD3827EC57`.
SHA256 ES: `E44F2CD9B446CC6242D3970EA5DDF11715592E30AB44853624BAA62650DC9DCA`.

MacBuilder checked the same archive on Apple M1 / macOS Tahoe 26.6 / R 4.6.1
Patched: **Status OK, 0 errors, 0 warnings and 0 NOTEs**, 3416 passed assertions
and the same 13 expected skips. Test sources in the remote report match the
checked source. See the
[MacBuilder report](https://mac.R-project.org/macbuilder/results/1792787197-53b23dc7692f1423/).
The service ran its R-release profile despite R-devel being requested through
the official form value; this result does not establish R-devel coverage.
Its standard check also does not replace the complete Windows `--as-cran`
check described above. Multiplatform CI results are available through
[GitHub Actions](https://github.com/gmoncayoj/PsychoMatic/actions/workflows/R-CMD-check.yaml)
and apply to the source commit identified in each run. Checks of the exact
tarball are documented separately above. Editorial resubmission of 1.0.1
to CRAN has not been performed as part of this preparation.

## History: GPL 1.0.0 candidate

Checks completed on October 23, 2026. This historical section describes the
candidate licensed under GPL version 3 or later, copyright 2026
José Gamarra-Moncayo. It does not constitute confirmation of CRAN publication.

## Checked archive

`PsychoMatic_1.0.0.tar.gz`, 459506 bytes.
SHA256: `A062B6F2C4107BA1D1207907A2DBE61C890E21DB259AB10C0C10E4C68BA5B9BC`.
The hash identifies the checked archive; rebuilding can produce a different
hash because of metadata, even when the code is unchanged.

| Environment | Result for the same GPL tarball |
| --- | --- |
| Windows 11 x64, R 4.6.1 | Complete `R CMD check --as-cran`, Status OK |
| Windows x64, R-devel r90748 (2026-10-22) | Win-builder, Status OK |

Both checks completed without errors, warnings or NOTEs and included tests,
examples, vignettes and PDF/HTML manuals. The
[Win-builder report](https://win-builder.r-project.org/AyEWPJy5LMt4/00check.log)
documents the remote result; a copy and the local report are retained in the
maintainer's evidence because service links are temporary.

The locally installed suite passed **3402 assertions**, with zero failures
or warnings and **13 deliberate intensive-test skips**. Those regressions
remain available with `PSYCHOMATIC_EXTENDED_TESTS=true`; see
[TESTING.md](TESTING.md). Tests took 583.53 seconds in that run. Duration varies
with load and environment; this is neither a benchmark nor a duration guarantee.

The GPL candidate's 44 R files, 47 help topics and 60 test files are identical
to those of the preceding MIT candidate. The license change was verified inside
the new archive. Linux/macOS checks of the previous candidate are not attributed
to the new GPL tarball. The GitHub Actions matrix was configured for Linux,
Windows, macOS and multiple R versions, but had not run on GitHub when these
historical checks were closed.

## Manuals and reverse dependencies

The [complete manuals](manuals/README.md) cover all 22 public functions and
47 help topics; matching signatures and examples were verified across languages.

| Manual | Pages | SHA256 |
| --- | --- | --- |
| English | 49 | `53B5915329A98DF9464F123C1FC288CD4DA30AA078F393E8AFA6D6B566AE5148` |
| Spanish | 51 | `76B10102D01418C8369978D183FD35AF30413F36C5CBA85794FA7DAAEE197C5A` |

The CRAN index reviewed on October 22, 2026 contained 25528 packages and showed
no reverse dependencies for PsychoMatic in Depends, Imports, LinkingTo,
Suggests or Enhances, including recursive dependencies. This does not cover
private code or external repositories. The lavaan tutorial and its official
sources also returned HTTP200 with the expected content on that date.
