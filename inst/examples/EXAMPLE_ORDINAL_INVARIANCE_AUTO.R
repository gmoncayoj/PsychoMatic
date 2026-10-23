# PsychoMatic: three-stage ordinal invariance with reproducible simulated data.
set.seed(7114)
n <- 600L
latent <- stats::rnorm(n)
items <- as.data.frame(sapply(seq_len(6L), function(j) {
  response <- .75 * latent + stats::rnorm(n, sd = sqrt(1 - .75^2))
  ordered(cut(response, breaks = c(-Inf, -1, -.3, .3, 1, Inf), labels = FALSE),
    levels = 1:5)
}))
names(items) <- paste0("R", seq_len(6L))
items$Sexo <- factor(rep(c("Femenino", "Masculino"), each = n / 2L))
model <- "F =~ R1 + R2 + R3 + R4 + R5 + R6"
ordinal_fit <- PsychoMatic::ordinal_invariance_auto(
  items, group = "Sexo", model = model, estimator = "WLSMV", language = "eng"
)
summary(ordinal_fit)$fit
# Reports are saved in getwd() unless an explicit file_name is supplied:
# PsychoMatic::psychomatic_report(ordinal_fit, format = c("word", "excel"))
