# Regenerate the two internal regression fixtures from the
# canonical package generators. Run this script from the package root.

source(file.path("R", "demo_data_helpers.R"), local = TRUE)

if (!dir.exists(file.path("tests", "fixtures"))) {
  stop("Run data-raw/generate_demo_data.R from the PsychoMatic package root.")
}

likert_demo <- .pm_generate_likert_demo()
invariance_demo <- .pm_generate_invariance_demo()

save(likert_demo, file = file.path("tests", "fixtures", "likert_demo.rda"),
     compress = "xz", version = 3)
save(invariance_demo, file = file.path("tests", "fixtures", "invariance_demo.rda"),
     compress = "xz", version = 3)
