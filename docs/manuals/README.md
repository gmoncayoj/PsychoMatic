# Complete reference manuals

The 1.0.1 maintenance edition retains all 22 functions and 47 topics from
1.0.0. The correction affects path tests, with no API or method changes.
The updated manuals were built and verified in both languages.

- [PsychoMatic_1.0.1_EN.pdf](PsychoMatic_1.0.1_EN.pdf): updated English reference.
- [PsychoMatic_1.0.1_ES.pdf](PsychoMatic_1.0.1_ES.pdf): updated Spanish translation.

The 1.0.0 PDFs are retained to preserve README links already included in the
package published on CRAN.

The October 23, 2026 edition of both manuals includes GPL version 3 or later
and copyright 2026 José Gamarra-Moncayo. The maintainer retains earlier MIT
PDFs as local historical evidence. Hashes for the corresponding 1.0.0 and
1.0.1 editions are recorded in [VALIDATION.md](../VALIDATION.md).

- `PsychoMatic_1.0.0_EN.pdf`: canonical English reference for 1.0.0.
- `PsychoMatic_1.0.0_ES.pdf`: complete Spanish translation for 1.0.0.

Both languages include every topic in `man/`, including those marked internal
by R: 47 topics covering all 22 exported functions and documented methods.
Signatures, options and R syntax are identical across languages. Help examples
use templates that researchers adapt to their data frames.

`man/` remains the installed English help and the source CRAN uses to generate
its official manual. The Spanish PDF is supplementary project documentation;
it does not replace that help or imply CRAN editorial approval. Both PDFs
use the native Rd/LaTeX reference-manual layout.

## Rebuilding

From the repository root, with R and a LaTeX distribution available:

```sh
Rscript docs/manuals/check_manuals.R
Rscript docs/manuals/build_manuals.R
```

When the API changes, update roxygen and the canonical help first, then
synchronize the Spanish prose in `es/man/`. `check_manuals.R` requires matching
topics, names/aliases, signatures, arguments and executable examples, as well
as coverage of every export. `coverage.csv` records hashes for both sources.
The builder preserves syntax and technical labels, translates Spanish-manual
headings and uses an explicit edition date.

Temporary files are created in `work/`, which Git ignores. The entire `docs/`
directory is excluded from the source package through `.Rbuildignore`,
avoiding duplicate help files or redundant PDFs in the CRAN submission.
