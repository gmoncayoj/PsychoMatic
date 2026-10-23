# PsychoMatic 1.0.0: alignment from group-level parameter matrices


item_names <- paste0("i", 1:9)
groups <- c("A", "B")
# These are illustrative configural-model parameter matrices, not raw item means.
lambda <- rbind(
  A = c(0.78, 0.73, 0.76, 0.72, 0.80, 0.69, 0.77, 0.74, 0.71),
  B = c(0.76, 0.75, 0.72, 0.74, 0.78, 0.71, 0.75, 0.76, 0.70)
)
nu <- rbind(
  A = rep(0, length(item_names)),
  B = c(0.03, -0.02, 0.45, 0.01, -0.03, 0.02, 0.00, 0.02, -0.01)
)
dimnames(lambda) <- dimnames(nu) <- list(groups, item_names)

if (!requireNamespace("sirt", quietly = TRUE)) {
  message("Skipped alignment: optional package 'sirt' is not installed.")
} else {
  alignment <- PsychoMatic::inv_align_auto(
    lambda = lambda,
    nu = nu,
    fixed = TRUE,
    meth = 2,
    language = "eng"
  )
}
