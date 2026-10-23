test_that("English reader headers do not expose internal deletion names", {
  raw <- data.frame(item_rest_r = .6, alpha_if_deleted = .8,
    alpha_ordinal_if_deleted = .9, omega_if_deleted = .85,
    p_value = .01, CI_lower = .7, CI_upper = .9)
  shown <- PsychoMatic:::.pm_label_table(raw, "eng")
  expect_identical(names(shown), c("Item-rest r", "Alpha if deleted",
    "Ordinal alpha if deleted", "Omega if deleted", "p", "Lower CI", "Upper CI"))
  expect_identical(unname(shown), unname(raw))
  expect_identical(names(raw)[2], "alpha_if_deleted")
})
