.pm_select_cfa_estimator <- function(max_cat, min_cat, n_complete,
                                     mvn_diagnostic, estimator = NULL,
                                     ordered = NULL, language = "eng") {
  is_english <- identical(language, "eng")
  if (is.null(estimator)) {
    if (max_cat <= 5) {
      # ------ Items ordinales ------
      if (is.null(ordered)) ordered <- TRUE

      if (n_complete < 200) {
        estimator <- "ULSMV"
        estimator_reason <- if (is_english) {
          sprintf(
            "Ordinal items (%d-%d cat.) + n=%d (<200) -> ULSMV (recommended for small ordinal samples)",
            min_cat, max_cat, n_complete
          )
        } else {
          sprintf(
            "Items ordinales (%d-%d cat.) + n=%d (<200) -> ULSMV (recomendado para muestras pequenas ordinales)",
            min_cat, max_cat, n_complete
          )
        }
      } else {
        estimator <- "WLSMV"
        estimator_reason <- if (is_english) {
          sprintf(
            "Ordinal items (%d-%d cat.) + n=%d -> WLSMV (robust WLS for ordinal data)",
            min_cat, max_cat, n_complete
          )
        } else {
          sprintf(
            "Items ordinales (%d-%d cat.) + n=%d -> WLSMV (WLS robusto para datos ordinales)",
            min_cat, max_cat, n_complete
          )
        }
      }

    } else if (max_cat <= 7 && min_cat <= 5) {
      # ------ Mezcla ------
      if (is.null(ordered)) ordered <- TRUE
      estimator <- "WLSMV"
      estimator_reason <- if (is_english) {
        sprintf(
          "Mixed item categories (%d-%d) -> WLSMV (ordinal treatment recommended)",
          min_cat, max_cat
        )
      } else {
        sprintf(
          "Mezcla de items con %d-%d categorias -> WLSMV (tratamiento ordinal recomendado)",
          min_cat, max_cat
        )
      }

    } else {
      # ------ Items continuos (>= 6 categorias) ------
      if (is.null(ordered)) ordered <- FALSE

      if (identical(mvn_diagnostic$status, "ok")) {
        mardia_p_skew <- mvn_diagnostic$p_skew
        mardia_p_kurt <- mvn_diagnostic$p_kurt
        mardia_skew_stat <- mvn_diagnostic$skew_stat
        mardia_kurt_stat <- mvn_diagnostic$kurt_stat
        is_mvn <- isTRUE(mvn_diagnostic$normal)

        if (is_mvn) {
          estimator <- "ML"
          estimator_reason <- if (is_english) {
            sprintf(
              "Continuous items (%d-%d cat.), multivariate normality NOT rejected (Mardia: p_skew=%.4f, p_kurt=%.4f) -> ML",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          } else {
            sprintf(
              "Items continuos (%d-%d cat.), normalidad multivariada NO rechazada (Mardia: p_asim=%.4f, p_curt=%.4f) -> ML",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          }
        } else {
          estimator <- "MLR"
          estimator_reason <- if (is_english) {
            sprintf(
              "Continuous items (%d-%d cat.), multivariate normality REJECTED (Mardia: p_skew=%.4f, p_kurt=%.4f) -> MLR",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          } else {
            sprintf(
              "Items continuos (%d-%d cat.), normalidad multivariada RECHAZADA (Mardia: p_asim=%.4f, p_curt=%.4f) -> MLR",
              min_cat, max_cat, mardia_p_skew, mardia_p_kurt
            )
          }
        }
      } else {
        estimator <- "MLR"
        estimator_reason <- if (is_english) {
          sprintf(
            "Continuous items (%d-%d cat.), normality test unavailable -> MLR as a precaution",
            min_cat, max_cat
          )
        } else {
          sprintf(
            "Items continuos (%d-%d cat.), prueba de normalidad no disponible -> MLR como precaucion",
            min_cat, max_cat
          )
        }
      }
    }
  } else {
    estimator_reason <- if (is_english) {
      "Estimator specified by the user"
    } else {
      "Estimador especificado por el usuario"
    }
    if (is.null(ordered)) {
      ordered <- if (max_cat <= 5) TRUE else FALSE
    }
  }
  list(estimator = estimator, ordered = ordered, reason = estimator_reason)
}
