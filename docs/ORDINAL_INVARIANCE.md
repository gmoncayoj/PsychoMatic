# Ordinal measurement invariance

Implemented the independent three-stage ordinal API in
`R/ordinal_invariance_auto.R`, using the shared native fitting and diagnostic
engine in `R/INV_FAC_AUTO.R`. It fits exactly configural, equal thresholds,
and equal thresholds and loadings, with theta parameterization and
`ID.cat = "Wu.Estabrook.2016"`. Estimators are WLSMV (automatic default) and
ULSMV; ML/MLR and nonordinal model indicators are rejected.

The third stage is a test of thresholds and loadings, not a declaration of
equal latent-response intercepts or residuals. Separate threshold equality
with zero additional df is unevaluable; the next identified comparison uses
the last evaluable prerequisite. Failed configural fit or inadmissible models
prevent favorable downstream conclusions. Fit changes are signed differences;
fit improvements are not counted as deterioration. Chen's .010 CFI / .015
RMSEA / .030 loading-SRMR / .010 other-SRMR guidelines are descriptive
extrapolations to ordinal estimation. Custom criteria remain available.

Reader tables and Word/Excel exports show the three levels and fit changes.
References follow existing APA typography. The Svetina et al. article is cited
as 2020 (volume 27, issue 1), although it appeared online in 2019.

## Verification

The regression coverage includes ordinal semTools agreement, ML/MLR rejection,
ordinal-input checks, two/three-category identification, ULSMV, group
symbols/vectors, signed changes and bilingual Word/Excel reporting. See
[VALIDATION.md](VALIDATION.md) for the exact GPL archive's complete Windows
release and R-devel checks and [TESTING.md](TESTING.md) for the distinction
between the distribution suite and optional intensive regressions.

## Separation from native factorial invariance

`factorial_invariance_auto()` fits native lavaan configural, metric,
scalar and strict calls. WLSMV/ULSMV default to ordered responses, marker
identification and delta parameterization. Thresholds are estimated internally
but no threshold-equality stage is added. Native ordinal scalar/strict calls
may not add independent restrictions and do not establish ordinal invariance.
The ordinal API retains its identified Wu-Estabrook three-stage sequence.
Both APIs default to scaled CFI/RMSEA for corrected estimators, and accept
`fit_variant = "robust"`; this controls fit changes and interpretations too.
SRMR is base and corrected chi-square/df/p stay scaled. Unavailable requested
robust indices remain missing, with no substitution.
`invariance_partial()` defaults to the native factorial
workflow too. Identified threshold and latent-response releases require its
explicit `workflow = "ordinal"` option. The categorical identification
primitives remain tested separately from the default native workflow.

## References

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement
invariance. *Structural Equation Modeling: A Multidisciplinary Journal, 14*(3),
464-504. https://doi.org/10.1080/10705510701301834

Svetina, D., Rutkowski, L., & Rutkowski, D. (2020). Multiple-group invariance with
categorical outcomes using updated guidelines: An illustration using Mplus and
the lavaan/semTools packages. *Structural Equation Modeling: A Multidisciplinary
Journal, 27*(1), 111-130. https://doi.org/10.1080/10705511.2019.1602776

Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor analysis
models of different levels of invariance for ordered categorical outcomes.
*Psychometrika, 81*(4), 1014-1045. https://doi.org/10.1007/s11336-016-9506-0
