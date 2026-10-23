.pm_network_versions <- function() {
  packages <- c("PsychoMatic", "EGAnet", "bootnet", "NetworkComparisonTest", "lavaan", "igraph", "qgraph")
  data.frame(package = packages, version = vapply(packages, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed", character(1)))
}

.pm_network_limitations <- function() {
  c("EGA partitions and hierarchical candidates are exploratory. Fit on the same sample is not independent confirmation.",
    "UVA reduction is conditional item selection; substantive content and retained representatives require review.",
    "riEGA requires original unrecoded responses; a fitted random intercept does not prove wording effects or identify inverse items.",
    "Revised loading magnitudes: |loading| <.20 below small, .20 small, .35 moderate, .50 large; heuristic interpretation, not cutoffs for item deletion.",
    "TEFI is relative fit on the same correlation matrix and item set. genTEFI and residual TEFI have separate comparison sets.",
    "Strength uses absolute edge weights. Shortest paths use inverse absolute weights; closeness is undefined for disconnected networks.",
    "Centrality interpretation requires sampling stability: CS <.25 is inadequate; >=.50 is preferable. Betweenness and closeness are often unstable.",
    "NCT compares network structure and global strength. EGA loading permutation tests concern metric network invariance; neither establishes CFA scalar invariance.",
    "BH controls false discovery rate within each explicitly declared family across all tested group pairs. Nonsignificant results do not establish equivalence.",
    "Small bootstrap/permutation counts are for diagnostics. Publication defaults are 1000 bootstraps and 5000 permutations; report Monte Carlo uncertainty.",
    "Listwise deletion is explicit; pairwise correlations can be non-positive-definite. Missing-data mechanisms require substantive assessment.")
}

.pm_network_references <- function() {
  c(
    "Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B (Methodological), 57(1), 289\u2013300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
    "Bulut, H. C., Bulut, O., & Clelland, A. (2025). A psychometric network analysis approach for detecting item wording effects in self-report measures across subgroups. Field Methods, 37(1), 20\u201336. https://doi.org/10.1177/1525822X241247444",
    "Christensen, A. P., & Golino, H. (2021). Estimating the stability of psychological dimensions via bootstrap exploratory graph analysis: A Monte Carlo simulation and tutorial. Psych, 3(3), 479\u2013500. https://doi.org/10.3390/psych3030032",
    "Christensen, A. P., Garrido, L. E., & Golino, H. (2023). Unique variable analysis: A network psychometrics method to detect local dependence. Multivariate Behavioral Research, 58(6), 1165\u20131182. https://doi.org/10.1080/00273171.2023.2194606",
    "Christensen, A. P., Garrido, L. E., Guerra-Pe\u00f1a, K., & Golino, H. (2024). Comparing community detection algorithms in psychometric networks: A Monte Carlo simulation. Behavior Research Methods, 56(3), 1485\u20131505. https://doi.org/10.3758/s13428-023-02106-4",
    "Christensen, A. P., Golino, H., Abad, F. J., & Garrido, L. E. (2025). Revised network loadings. Behavior Research Methods, 57, Article 114. https://doi.org/10.3758/s13428-025-02640-3",
    "Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks and their accuracy: A tutorial paper. Behavior Research Methods, 50(1), 195\u2013212. https://doi.org/10.3758/s13428-017-0862-1",
    "Faraci, P., & Nasonte, G. (2026). Disentangling substantive and method variance in mixed-worded scales: An empirical application of the random intercept item factor analysis (RIIFA). Behavior Research Methods, 58, Article 299. https://doi.org/10.3758/s13428-026-03163-1",
    "Fonseca-Pedrero, E., & Mu\u00f1iz, J. (2025). An\u00e1lisis de redes en la medici\u00f3n psicol\u00f3gica: Fundamentos. Acci\u00f3n Psicol\u00f3gica, 22(1), 87\u2013100. https://doi.org/10.5944/ap.22.1.43296",
    "Garcia-Pardina, A., Abad, F. J., Christensen, A. P., Golino, H., & Garrido, L. E. (2024). Dimensionality assessment in the presence of wording effects: A network psychometric and factorial approach. Behavior Research Methods, 56(6), 6179\u20136197. https://doi.org/10.3758/s13428-024-02348-w",
    "Golino, H., Moulder, R., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2021). Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables. Multivariate Behavioral Research, 56(6), 874\u2013902. https://doi.org/10.1080/00273171.2020.1779642",
    "Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Mart\u00ednez-Molina, A. (2020). Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial. Psychological Methods, 25(3), 292\u2013320. https://doi.org/10.1037/met0000255",
    "Isvoranu, A. M., & Epskamp, S. (2023). Which estimation method to choose in network psychometrics? Deriving guidelines for applied researchers. Psychological Methods, 28(4), 925\u2013946. https://doi.org/10.1037/met0000439",
    "Jamison, L., Christensen, A. P., & Golino, H. F. (2024). Metric invariance in exploratory graph analysis via permutation testing. Methodology, 20(2), Article e12877. https://doi.org/10.5964/meth.12877",
    "Jim\u00e9nez, M., Abad, F. J., Garcia-Garzon, E., Golino, H., Christensen, A. P., & Garrido, L. E. (2025). Dimensionality assessment in bifactor structures with multiple general factors: A network psychometrics approach. Psychological Methods, 30(4), 770\u2013792. https://doi.org/10.1037/met0000590",
    "Shi, D., Christensen, A. P., Day, E. A., Golino, H. F., & Garrido, L. E. (2025). Exploring estimation procedures for reducing dimensionality in psychological network modeling. Multivariate Behavioral Research, 60(2), 184\u2013210. https://doi.org/10.1080/00273171.2024.2395941",
    "van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2023). Comparing network structures on three aspects: A permutation test. Psychological Methods, 28(6), 1273\u20131285. https://doi.org/10.1037/met0000476"
  )
}
