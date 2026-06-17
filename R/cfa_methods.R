# METODO PRINT
# =============================================================================

.print_cfa_auto_legacy <- function(x, digits = 3, ...) {

  sep  <- paste(rep("=", 70), collapse = "")
  sep2 <- paste(rep("-", 70), collapse = "")

  cat("\n", sep, "\n", sep = "")
  cat("       CONFIRMATORY FACTOR ANALYSIS - PSYCHOMETRIC REPORT\n")
  cat(sep, "\n\n", sep = "")

  # --- WARNINGS ---
  if (length(x$warnings) > 0) {
    cat("*** WARNINGS ***\n")
    cat(sep2, "\n", sep = "")
    for (w in x$warnings) {
      cat("  [!] ", w, "\n", sep = "")
    }
    cat("\n")
  }

  # --- MODEL INFORMATION ---
  cat("MODEL INFORMATION\n")
  cat(sep2, "\n", sep = "")
  cat("  Model type           :", x$model_type, "\n")
  cat("  Factors              :", paste(x$factors, collapse = ", "), "\n")
  cat("  Observed items       :", x$n_items, "\n")
  cat("  N (total)            :", x$n, "\n")
  if (x$pct_missing > 0) {
    cat("  N (completos)        :", x$n_complete,
        sprintf("(%.1f%% with missing data)\n", x$pct_missing))
  }
  cat("  Response categories  :", paste0(min(x$n_categories), "-",
                                          max(x$n_categories)), "\n")
  cat("  Estimator            :", x$estimator, "\n")
  cat("  Rationale            :", x$estimator_reason, "\n")
  cat("  Converged            :",
      ifelse(x$converged, "Yes", "*** NO ***"), "\n")
  cat("\n")

  # --- MULTIVARIATE NORMALITY ---
  if (!is.na(x$mardia$p_skew)) {
    cat("MULTIVARIATE NORMALITY (Mardia)\n")
    cat(sep2, "\n", sep = "")
    cat(sprintf("  Skewness  : Statistic = %.2f, p = %.4f  %s\n",
                x$mardia$skew_stat, x$mardia$p_skew,
                ifelse(x$mardia$p_skew < 0.05, "[Rejected]", "[Not rejected]")))
    cat(sprintf("  Kurtosis  : Statistic = %.2f, p = %.4f  %s\n",
                x$mardia$kurt_stat, x$mardia$p_kurt,
                ifelse(x$mardia$p_kurt < 0.05, "[Rejected]", "[Not rejected]")))
    cat("\n")
  }

  # --- FIT INDICES ---
  cat("FIT INDICES\n")
  cat(sep2, "\n", sep = "")
  cat("  Nota:", x$fi_variant_note, "\n\n")
  fi <- x$fit_indices
  for (i in seq_len(nrow(fi))) {
    val <- fi$Valor[i]
    ev  <- fi$Evaluacion[i]
    if (is.na(val)) next
    eval_str <- if (nchar(ev) > 0) paste0("  [", ev, "]") else ""
    if (fi$Indice[i] %in% c("gl")) {
      cat(sprintf("  %-22s: %d%s\n", fi$Indice[i], as.integer(val), eval_str))
    } else {
      cat(sprintf("  %-22s: %.3f%s\n", fi$Indice[i], val, eval_str))
    }
  }
  cat("\n")

  # --- FACTOR LOADINGS ---
  cat("STANDARDIZED FACTOR LOADINGS\n")
  cat(sep2, "\n", sep = "")
  ld <- x$loadings
  for (f in unique(ld$Factor)) {
    cat(sprintf("\n  Factor: %s\n", f))
    sub <- ld[ld$Factor == f, ]
    for (j in seq_len(nrow(sub))) {
      alerta <- if (nchar(sub$Alerta[j]) > 0) paste0("  ", sub$Alerta[j]) else ""
      cat(sprintf("    %-18s %6.3f  (SE=%.3f, z=%.2f, p=%.4f)  95%% CI [%.3f, %.3f]%s\n",
                  sub$Item[j], sub$Lambda[j], sub$SE[j], sub$z[j], sub$p[j],
                  sub$IC_inf[j], sub$IC_sup[j], alerta))
    }
  }
  cat("\n")

  # --- Second-order loadings, if applicable ---
  if (!is.null(x$so_loadings) && nrow(x$so_loadings) > 0) {
    cat("SECOND-ORDER LOADINGS\n")
    cat(sep2, "\n", sep = "")
    so <- x$so_loadings
    for (j in seq_len(nrow(so))) {
      cat(sprintf("  %s -> %s : %.3f  (SE=%.3f, z=%.2f, p=%.4f)\n",
                  so$Factor_SO[j], so$Factor_PO[j], so$Lambda[j],
                  so$SE[j], so$z[j], so$p[j]))
    }
    cat("\n")
  }

  # --- R-SQUARED ---
  if (!is.null(x$r_squared)) {
    cat("R-SQUARED BY ITEM\n")
    cat(sep2, "\n", sep = "")
    r2 <- x$r_squared
    if (is.list(r2)) {
      for (nm in names(r2)) {
        r2_vec <- r2[[nm]]
        for (j in seq_along(r2_vec)) {
          cat(sprintf("    %-18s %.3f\n", names(r2_vec)[j], r2_vec[j]))
        }
      }
    } else {
      for (j in seq_along(r2)) {
        cat(sprintf("    %-18s %.3f\n", names(r2)[j], r2[j]))
      }
    }
    cat("\n")
  }

  # --- FACTOR CORRELATIONS ---
  if (!is.null(x$factor_correlations)) {
    cat("FACTOR CORRELATIONS\n")
    cat(sep2, "\n", sep = "")
    fc <- x$factor_correlations
    for (j in seq_len(nrow(fc))) {
      alerta <- if (nchar(fc$Alerta[j]) > 0) paste0("  *** ", fc$Alerta[j]) else ""
      cat(sprintf("  %-12s <-> %-12s : %6.3f  (SE=%.3f, p=%.4f)%s\n",
                  fc$Factor1[j], fc$Factor2[j], fc$r[j], fc$SE[j], fc$p[j], alerta))
    }
    cat("\n")
  }

  # --- RELIABILITY ---
  cat("RELIABILITY\n")
  cat(sep2, "\n", sep = "")

  rel <- x$reliability

  if (rel$type == "error") {
    cat("  Calculation error: ", rel$message, "\n")

  } else if (rel$type == "standard") {
    cat("  Method: McDonald's omega and Cronbach's alpha (semTools::reliability)\n\n")
    rv <- rel$values
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    etiquetas <- c(omega = "Omega", alpha = "Alpha", avevar = "AVE")
    for (f in colnames(rv)) {
      if (f == "total" && ncol(rv) > 2) next
      cat(sprintf("  Factor: %s\n", f))
      for (m in show_rows) {
        cat(sprintf("    %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
      cat("\n")
    }

  } else if (rel$type == "second_order") {
    cat("  Second-order model\n\n")
    cat("  -- First order, by dimension --\n")
    rv <- rel$first_order
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    etiquetas <- c(omega = "Omega", alpha = "Alpha", avevar = "AVE")
    for (f in colnames(rv)) {
      if (f == "total") next
      cat(sprintf("    Factor: %s\n", f))
      for (m in show_rows) {
        cat(sprintf("      %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
    }
    cat("\n  -- Second order --\n")
    if (!is.null(rel$second_order)) {
      for (sof in rel$so_factors) {
        l2 <- rel$second_order[[sof]]
        if (!is.null(l2)) {
          cat(sprintf("    Factor: %s\n", sof))
          for (nm in names(l2)) {
            cat(sprintf("      %-25s: %.3f\n", nm, l2[nm]))
          }
          cat("\n")
        } else {
          cat(sprintf("    Factor: %s - Could not be computed\n", sof))
        }
      }
    } else {
      cat("    Could not be computed; inspect the model specification.\n")
    }
    cat("\n")

  } else if (rel$type == "bifactor") {
    cat("  Bifactor model - specialized indices\n\n")
    cat(sprintf("  Omega total             : %.3f\n", rel$omega_total))
    cat(sprintf("  Hierarchical omega (w_h): %.3f  [%s]\n",
                rel$omega_h, rel$general_factor))
    cat(sprintf("  ECV                     : %.3f\n", rel$ecv))
    cat(sprintf("  PUC                     : %.3f\n", rel$puc))
    cat(sprintf("  H (general factor)      : %.3f\n", rel$h_general))

    cat("\n  By specific factor:\n")
    cat(sprintf("  %-20s  %8s  %8s\n", "Factor", "Omega_hs", "H"))
    for (i in seq_along(rel$omega_hs)) {
      cat(sprintf("  %-20s  %8.3f  %8.3f\n",
                  names(rel$omega_hs)[i], rel$omega_hs[i], rel$h_specific[i]))
    }

    cat("\n  Interpretation:\n")
    if (rel$omega_h >= 0.80) {
      cat("    Omega_h >= .80: strong general factor.\n")
      cat("    -> The total score mainly reflects the general factor.\n")
    } else if (rel$omega_h >= 0.50) {
      cat("    Omega_h between .50 and .80: moderate general factor.\n")
      cat("    -> Evaluate whether subscales add incremental information.\n")
    } else {
      cat("    Omega_h < .50: weak general factor.\n")
      cat("    -> Subscale scores may be more informative.\n")
    }
    # Interpretación conjunta ECV/PUC basada en Rodriguez et al. (2016) y Reise et al. (2013)
    if (rel$puc > 0.80 && rel$ecv > 0.60 && rel$omega_h >= 0.70) {
      cat("    PUC > .80, ECV > .60, and Omega_h >= .70: evidence supports\n")
      cat("    essential unidimensionality (Reise et al., 2013).\n")
    } else if (rel$ecv > 0.70 && rel$puc > 0.70) {
      cat("    ECV > .70 and PUC > .70: evidence may support\n")
      cat("    unidimensionality (Rodriguez et al., 2016).\n")
    } else if (rel$ecv < 0.60) {
      cat("    ECV < .60: substantial multidimensionality is present.\n")
    } else {
      cat(sprintf("    ECV = %.3f, PUC = %.3f: No se cumplen los umbrales para\n", rel$ecv, rel$puc))
      cat("    support unidimensionality. Interpret subscales cautiously.\n")
    }
    for (i in seq_along(rel$h_specific)) {
      if (rel$h_specific[i] < 0.50) {
        cat(sprintf("    H(%s) < .50: specific factor has limited replicability.\n",
                    names(rel$h_specific)[i]))
      }
    }
    cat("\n")
  }

  # --- RESPECIFICATIONS ---
  if (length(x$respecifications) > 0) {
    cat("RESPECIFICATION SUGGESTIONS\n")
    cat(sep2, "\n", sep = "")
    cat("  (Consider only with theoretical or substantive justification)\n\n")
    for (i in seq_along(x$respecifications)) {
      r <- x$respecifications[[i]]
      respec_type <- if (!is.null(r$type)) r$type else r$tipo
      respec_syntax <- if (!is.null(r$syntax)) r$syntax else r$sintaxis
      respec_note <- if (!is.null(r$note)) r$note else r$nota
      cat(sprintf("  %d. [%s] %s\n", i, respec_type, respec_syntax))
      cat(sprintf("     IM = %.2f, EPC = %.3f\n", r$mi, r$epc))
      cat(sprintf("     %s\n\n", respec_note))
    }
  } else {
    cat("RESPECIFICATIONS\n")
    cat(sep2, "\n", sep = "")
    cat("  No modification indices above the threshold were found.\n\n")
  }

  # --- REFERENCES ---
  cat("REFERENCES (APA 7)\n")
  cat(sep2, "\n", sep = "")
  cat("  (Only references relevant to the analysis performed are included)\n\n")
  refs <- .build_references(x)
  for (i in seq_along(refs)) {
    cat(sprintf("  [%d] %s\n\n", i, refs[i]))
  }

  cat(sep, "\n", sep = "")
  cat("lavaan object available in $fit | All indices in $all_fit_measures\n")
  cat(sep, "\n\n", sep = "")

  invisible(x)
}

#' Print a `cfa_auto` object
#'
#' @param x Object returned by `cfa_auto()`.
#' @param digits Number of decimal places shown in the printed report.
#' @param ... Additional unused arguments.
#'
#' @rdname cfa_auto
#' @method print cfa_auto
#' @export
print.cfa_auto <- function(x, digits = 3, ...) {
  language <- if (!is.null(x$language)) .normalize_language(x$language) else "esp"
  is_english <- identical(language, "eng")
  txt <- .cfa_auto_labels(language)

  sep  <- paste(rep("=", 70), collapse = "")
  sep2 <- paste(rep("-", 70), collapse = "")

  cat("\n", sep, "\n", sep = "")
  cat(txt$report_title, "\n")
  cat(sep, "\n\n", sep = "")

  if (length(x$warnings) > 0) {
    cat(txt$warnings_title, "\n")
    cat(sep2, "\n", sep = "")
    for (w in x$warnings) {
      cat("  [!] ", w, "\n", sep = "")
    }
    cat("\n")
  }

  cat(txt$model_info_title, "\n")
  cat(sep2, "\n", sep = "")
  cat(txt$model_type_label, x$model_type, "\n")
  cat(txt$factors_label, paste(x$factors, collapse = ", "), "\n")
  cat(txt$observed_items_label, x$n_items, "\n")
  cat(txt$n_total_label, x$n, "\n")
  if (x$pct_missing > 0) {
    cat(txt$n_complete_label, x$n_complete, txt$missing_suffix(x$pct_missing), "\n")
  }
  cat(txt$categories_label, paste0(min(x$n_categories), "-", max(x$n_categories)), "\n")
  cat(txt$estimator_label, x$estimator, "\n")
  cat(txt$justification_label, x$estimator_reason, "\n")
  cat(txt$converged_label, ifelse(x$converged, txt$yes, txt$no), "\n")
  cat("\n")

  if (!is.na(x$mardia$p_skew)) {
    cat(txt$mardia_title, "\n")
    cat(sep2, "\n", sep = "")
    cat(sprintf("%s = %.2f, p = %.4f  %s\n",
                txt$skew_label,
                x$mardia$skew_stat, x$mardia$p_skew,
                ifelse(x$mardia$p_skew < 0.05, txt$rejected, txt$not_rejected)))
    cat(sprintf("%s = %.2f, p = %.4f  %s\n",
                txt$kurt_label,
                x$mardia$kurt_stat, x$mardia$p_kurt,
                ifelse(x$mardia$p_kurt < 0.05, txt$rejected, txt$not_rejected)))
    cat("\n")
  }

  cat(txt$fit_title, "\n")
  cat(sep2, "\n", sep = "")
  cat(txt$note_label, x$fi_variant_note, "\n\n")
  fi <- x$fit_indices
  metric_col <- if ("Metric" %in% names(fi)) fi$Metric else rep(NA_character_, nrow(fi))
  for (i in seq_len(nrow(fi))) {
    val <- fi$Valor[i]
    ev  <- fi$Evaluacion[i]
    if (is.na(val)) next
    eval_str <- if (nchar(ev) > 0) paste0("  [", ev, "]") else ""
    if ((!is.na(metric_col[i]) && metric_col[i] == "df") || fi$Indice[i] %in% c("gl", "df")) {
      cat(sprintf("  %-22s: %d%s\n", fi$Indice[i], as.integer(val), eval_str))
    } else {
      cat(sprintf("  %-22s: %.3f%s\n", fi$Indice[i], val, eval_str))
    }
  }
  cat("\n")

  cat(txt$loadings_title, "\n")
  cat(sep2, "\n", sep = "")
  ld <- x$loadings
  for (f in unique(ld$Factor)) {
    cat(sprintf("\n%s %s\n", txt$factor_label, f))
    sub <- ld[ld$Factor == f, ]
    for (j in seq_len(nrow(sub))) {
      alerta <- if (nchar(sub$Alerta[j]) > 0) paste0("  ", sub$Alerta[j]) else ""
      cat(sprintf("    %-18s %6.3f  (SE=%.3f, z=%.2f, p=%.4f)  %s [%.3f, %.3f]%s\n",
                  sub$Item[j], sub$Lambda[j], sub$SE[j], sub$z[j], sub$p[j],
                  txt$ci_label, sub$IC_inf[j], sub$IC_sup[j], alerta))
    }
  }
  cat("\n")

  if (!is.null(x$so_loadings) && nrow(x$so_loadings) > 0) {
    cat(txt$second_order_loadings_title, "\n")
    cat(sep2, "\n", sep = "")
    so <- x$so_loadings
    for (j in seq_len(nrow(so))) {
      cat(sprintf("  %s -> %s : %.3f  (SE=%.3f, z=%.2f, p=%.4f)\n",
                  so$Factor_SO[j], so$Factor_PO[j], so$Lambda[j],
                  so$SE[j], so$z[j], so$p[j]))
    }
    cat("\n")
  }

  if (!is.null(x$r_squared)) {
    cat(txt$r2_title, "\n")
    cat(sep2, "\n", sep = "")
    r2 <- x$r_squared
    if (is.list(r2)) {
      for (nm in names(r2)) {
        r2_vec <- r2[[nm]]
        for (j in seq_along(r2_vec)) {
          cat(sprintf("    %-18s %.3f\n", names(r2_vec)[j], r2_vec[j]))
        }
      }
    } else {
      for (j in seq_along(r2)) {
        cat(sprintf("    %-18s %.3f\n", names(r2)[j], r2[j]))
      }
    }
    cat("\n")
  }

  if (!is.null(x$factor_correlations)) {
    cat(txt$factor_corr_title, "\n")
    cat(sep2, "\n", sep = "")
    fc <- x$factor_correlations
    for (j in seq_len(nrow(fc))) {
      alerta <- if (nchar(fc$Alerta[j]) > 0) paste0("  *** ", fc$Alerta[j]) else ""
      cat(sprintf("  %-12s <-> %-12s : %6.3f  (SE=%.3f, p=%.4f)%s\n",
                  fc$Factor1[j], fc$Factor2[j], fc$r[j], fc$SE[j], fc$p[j], alerta))
    }
    cat("\n")
  }

  cat(txt$reliability_title, "\n")
  cat(sep2, "\n", sep = "")

  rel <- x$reliability

  if (rel$type == "error") {
    if (is_english) {
      cat("  Error while computing: ", rel$message, "\n")
    } else {
      cat("  Error al calcular: ", rel$message, "\n")
    }

  } else if (rel$type == "standard") {
    if (is_english) {
      cat("  Method: McDonald's Omega and Cronbach's Alpha (semTools::reliability)\n\n")
      etiquetas <- c(omega = "Omega", alpha = "Alpha", avevar = "AVE")
    } else {
      cat("  Metodo: Omega de McDonald y Alfa de Cronbach (semTools::reliability)\n\n")
      etiquetas <- c(omega = "Omega", alpha = "Alfa", avevar = "AVE")
    }
    rv <- rel$values
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    for (f in colnames(rv)) {
      if (f == "total" && ncol(rv) > 2) next
      cat(sprintf("%s %s\n", txt$factor_label, f))
      for (m in show_rows) {
        cat(sprintf("    %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
      cat("\n")
    }

  } else if (rel$type == "second_order") {
    if (is_english) {
      cat("  Second-order model\n\n")
      cat("  -- First-order (by dimension) --\n")
      etiquetas <- c(omega = "Omega", alpha = "Alpha", avevar = "AVE")
    } else {
      cat("  Modelo de Segundo Orden\n\n")
      cat("  -- Primer Orden (por dimension) --\n")
      etiquetas <- c(omega = "Omega", alpha = "Alfa", avevar = "AVE")
    }
    rv <- rel$first_order
    show_rows <- intersect(rownames(rv), c("omega", "alpha", "avevar"))
    for (f in colnames(rv)) {
      if (f == "total") next
      cat(sprintf("    Factor: %s\n", f))
      for (m in show_rows) {
        cat(sprintf("      %-12s: %.3f\n", etiquetas[m], rv[m, f]))
      }
    }
    if (is_english) {
      cat("\n  -- Second-order --\n")
    } else {
      cat("\n  -- Segundo Orden --\n")
    }
    if (!is.null(rel$second_order)) {
      for (sof in rel$so_factors) {
        l2 <- rel$second_order[[sof]]
        if (!is.null(l2)) {
          cat(sprintf("    Factor: %s\n", sof))
          for (nm in names(l2)) {
            cat(sprintf("      %-25s: %.3f\n", nm, l2[nm]))
          }
          cat("\n")
        } else if (is_english) {
          cat(sprintf("    Factor: %s - Could not be computed\n", sof))
        } else {
          cat(sprintf("    Factor: %s - No se pudo calcular\n", sof))
        }
      }
    } else if (is_english) {
      cat("    Could not be computed (check model specification)\n")
    } else {
      cat("    No se pudo calcular (revisar especificacion del modelo)\n")
    }
    cat("\n")

  } else if (rel$type == "bifactor") {
    if (is_english) {
      cat("  Bifactor model - Specialized indices\n\n")
      cat(sprintf("  Total omega             : %.3f\n", rel$omega_total))
      cat(sprintf("  Hierarchical omega (w_h): %.3f  [%s]\n",
                  rel$omega_h, rel$general_factor))
      cat(sprintf("  ECV                     : %.3f\n", rel$ecv))
      cat(sprintf("  PUC                     : %.3f\n", rel$puc))
      cat(sprintf("  H (general factor)      : %.3f\n", rel$h_general))
      cat("\n  By specific factor:\n")
    } else {
      cat("  Modelo Bifactor - Indices Especializados\n\n")
      cat(sprintf("  Omega total             : %.3f\n", rel$omega_total))
      cat(sprintf("  Omega jerarquico (w_h)  : %.3f  [%s]\n",
                  rel$omega_h, rel$general_factor))
      cat(sprintf("  ECV                     : %.3f\n", rel$ecv))
      cat(sprintf("  PUC                     : %.3f\n", rel$puc))
      cat(sprintf("  H (factor general)      : %.3f\n", rel$h_general))
      cat("\n  Por factor especifico:\n")
    }
    cat(sprintf("  %-20s  %8s  %8s\n", "Factor", "Omega_hs", "H"))
    for (i in seq_along(rel$omega_hs)) {
      cat(sprintf("  %-20s  %8.3f  %8.3f\n",
                  names(rel$omega_hs)[i], rel$omega_hs[i], rel$h_specific[i]))
    }

    if (is_english) {
      cat("\n  Interpretation:\n")
    } else {
      cat("\n  Interpretacion:\n")
    }
    if (rel$omega_h >= 0.80) {
      if (is_english) {
        cat("    Omega_h >= .80: STRONG general factor.\n")
        cat("    -> The total score mainly reflects the general factor.\n")
      } else {
        cat("    Omega_h >= .80: Factor general FUERTE.\n")
        cat("    -> La puntuacion total refleja predominantemente el factor general.\n")
      }
    } else if (rel$omega_h >= 0.50) {
      if (is_english) {
        cat("    Omega_h between .50-.80: MODERATE general factor.\n")
        cat("    -> Evaluate whether the subscales add incremental information.\n")
      } else {
        cat("    Omega_h entre .50-.80: Factor general MODERADO.\n")
        cat("    -> Evaluar si las subescalas aportan informacion incremental.\n")
      }
    } else if (is_english) {
      cat("    Omega_h < .50: WEAK general factor.\n")
      cat("    -> Subscale scores may be more informative.\n")
    } else {
      cat("    Omega_h < .50: Factor general DEBIL.\n")
      cat("    -> Las puntuaciones por subescala pueden ser mas informativas.\n")
    }

    if (rel$puc > 0.80 && rel$ecv > 0.60 && rel$omega_h >= 0.70) {
      if (is_english) {
        cat("    PUC > .80, ECV > .60, and Omega_h >= .70: Evidence favorable to\n")
        cat("    unidimensionality (Reise et al., 2013).\n")
      } else {
        cat("    PUC > .80, ECV > .60 y Omega_h >= .70: Evidencia favorable a la\n")
        cat("    unidimensionalidad (Reise et al., 2013).\n")
      }
    } else if (rel$ecv > 0.70 && rel$puc > 0.70) {
      if (is_english) {
        cat("    ECV > .70 and PUC > .70: One could conclude in favor of\n")
        cat("    unidimensionality (Rodriguez et al., 2016).\n")
      } else {
        cat("    ECV > .70 y PUC > .70: Se podria concluir a favor de la\n")
        cat("    unidimensionalidad (Rodriguez et al., 2016).\n")
      }
    } else if (rel$ecv < 0.60) {
      if (is_english) {
        cat("    ECV < .60: Substantial multidimensionality is present.\n")
      } else {
        cat("    ECV < .60: Multidimensionalidad sustancial presente.\n")
      }
    } else if (is_english) {
      cat(sprintf("    ECV = %.3f, PUC = %.3f: The thresholds for concluding\n", rel$ecv, rel$puc))
      cat("    in favor of unidimensionality are not met. Interpret subscales\n")
      cat("    with caution.\n")
    } else {
      cat(sprintf("    ECV = %.3f, PUC = %.3f: No se cumplen los umbrales para\n", rel$ecv, rel$puc))
      cat("    concluir a favor de la unidimensionalidad. Interpretar subescalas\n")
      cat("    con cautela.\n")
    }
    for (i in seq_along(rel$h_specific)) {
      if (rel$h_specific[i] < 0.50) {
        if (is_english) {
          cat(sprintf("    H(%s) < .50: Specific factor shows low replicability.\n",
                      names(rel$h_specific)[i]))
        } else {
          cat(sprintf("    H(%s) < .50: Factor especifico poco replicable.\n",
                      names(rel$h_specific)[i]))
        }
      }
    }
    cat("\n")
  }

  if (length(x$respecifications) > 0) {
    cat(txt$respec_title, "\n")
    cat(sep2, "\n", sep = "")
    cat(txt$respec_note, "\n\n")
    for (i in seq_along(x$respecifications)) {
      r <- x$respecifications[[i]]
      cat(sprintf("  %d. [%s] %s\n", i, r$tipo, r$sintaxis))
      cat(sprintf("     %s = %.2f, EPC = %.3f\n",
                  if (is_english) "MI" else "IM", r$mi, r$epc))
      cat(sprintf("     %s\n\n", r$nota))
    }
  } else {
    cat(txt$respec_none_title, "\n")
    cat(sep2, "\n", sep = "")
    cat(txt$respec_none, "\n\n")
  }

  cat(txt$references_title, "\n")
  cat(sep2, "\n", sep = "")
  cat(txt$references_note, "\n\n")
  refs <- .build_references(x)
  for (i in seq_along(refs)) {
    cat(sprintf("  [%d] %s\n\n", i, refs[i]))
  }

  cat(sep, "\n", sep = "")
  cat(txt$footer, "\n")
  cat(sep, "\n\n", sep = "")

  invisible(x)
}

.cfa_respecifications_to_df <- function(x) {
  respec <- if (!is.null(x$respecification_suggestions)) {
    x$respecification_suggestions
  } else {
    x$respecifications
  }

  if (is.null(respec) || length(respec) == 0) {
    return(data.frame())
  }

  do.call(
    rbind,
    lapply(respec, function(item) {
      data.frame(
        Type = if (!is.null(item$type)) item$type else item$tipo,
        Syntax = if (!is.null(item$syntax)) item$syntax else item$sintaxis,
        MI = item$mi,
        EPC = item$epc,
        SameFactor = if (!is.null(item$same_factor)) item$same_factor else item$mismo_factor,
        Note = if (!is.null(item$note)) item$note else item$nota,
        stringsAsFactors = FALSE
      )
    })
  )
}

.cfa_export_excel <- function(result, file_name) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required to export Excel reports.", call. = FALSE)
  }

  createWorkbook <- openxlsx::createWorkbook
  addWorksheet <- openxlsx::addWorksheet
  writeData <- openxlsx::writeData
  saveWorkbook <- openxlsx::saveWorkbook

  wb <- createWorkbook()
  fit_df <- if (!is.null(result$fit_summary)) result$fit_summary else result$fit_indices
  if (!"Index" %in% names(fit_df) && "Indice" %in% names(fit_df)) fit_df$Index <- fit_df$Indice
  if (!"Value" %in% names(fit_df) && "Valor" %in% names(fit_df)) fit_df$Value <- fit_df$Valor
  if (!"Evaluation" %in% names(fit_df) && "Evaluacion" %in% names(fit_df)) fit_df$Evaluation <- fit_df$Evaluacion

  loadings_df <- if (!is.null(result$factor_loadings)) result$factor_loadings else result$loadings
  if (!"Loading" %in% names(loadings_df) && "Lambda" %in% names(loadings_df)) loadings_df$Loading <- loadings_df$Lambda
  if (!"p_value" %in% names(loadings_df) && "p" %in% names(loadings_df)) loadings_df$p_value <- loadings_df$p
  if (!"CI_lower" %in% names(loadings_df) && "IC_inf" %in% names(loadings_df)) loadings_df$CI_lower <- loadings_df$IC_inf
  if (!"CI_upper" %in% names(loadings_df) && "IC_sup" %in% names(loadings_df)) loadings_df$CI_upper <- loadings_df$IC_sup
  if (!"Alert" %in% names(loadings_df) && "Alerta" %in% names(loadings_df)) loadings_df$Alert <- loadings_df$Alerta

  higher_order_df <- if (!is.null(result$higher_order_loadings)) result$higher_order_loadings else result$loadings_ho
  if (!is.null(higher_order_df) && nrow(higher_order_df) > 0) {
    if (!"HigherOrderFactor" %in% names(higher_order_df) && "Factor_SO" %in% names(higher_order_df)) higher_order_df$HigherOrderFactor <- higher_order_df$Factor_SO
    if (!"LowerOrderFactor" %in% names(higher_order_df) && "Factor_PO" %in% names(higher_order_df)) higher_order_df$LowerOrderFactor <- higher_order_df$Factor_PO
    if (!"Loading" %in% names(higher_order_df) && "Lambda" %in% names(higher_order_df)) higher_order_df$Loading <- higher_order_df$Lambda
    if (!"p_value" %in% names(higher_order_df) && "p" %in% names(higher_order_df)) higher_order_df$p_value <- higher_order_df$p
  }

  corr_df <- if (!is.null(result$factor_correlation_table)) result$factor_correlation_table else result$factor_correlations
  if (!is.null(corr_df) && nrow(corr_df) > 0) {
    if (!"Correlation" %in% names(corr_df) && "r" %in% names(corr_df)) corr_df$Correlation <- corr_df$r
    if (!"p_value" %in% names(corr_df) && "p" %in% names(corr_df)) corr_df$p_value <- corr_df$p
    if (!"Alert" %in% names(corr_df) && "Alerta" %in% names(corr_df)) corr_df$Alert <- corr_df$Alerta
  }

  warnings_vec <- if (!is.null(result$warnings_list)) result$warnings_list else result$warnings

  summary_df <- data.frame(
    Element = c(
      "Model type", "Factors", "Observed items", "Sample size",
      "Complete cases", "Missing rows (%)", "Estimator",
      "Estimator rationale", "Converged", "Language"
    ),
    Value = c(
      result$model_type,
      paste(result$factors, collapse = ", "),
      result$n_items,
      result$n,
      result$n_complete,
      result$pct_missing,
      result$estimator,
      result$estimator_reason,
      result$converged,
      result$language
    ),
    stringsAsFactors = FALSE
  )

  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", summary_df)

  addWorksheet(wb, "Fit Indices")
  writeData(wb, "Fit Indices", fit_df[, c("Metric", "Index", "Value", "Evaluation"), drop = FALSE])

  addWorksheet(wb, "Loadings")
  writeData(
    wb,
    "Loadings",
    loadings_df[, c("Factor", "Item", "Loading", "SE", "z", "p_value", "CI_lower", "CI_upper", "Alert"), drop = FALSE]
  )

  if (!is.null(higher_order_df) && nrow(higher_order_df) > 0) {
    addWorksheet(wb, "Higher Order")
    writeData(
      wb,
      "Higher Order",
      higher_order_df[, c("HigherOrderFactor", "LowerOrderFactor", "Loading", "SE", "z", "p_value"), drop = FALSE]
    )
  }

  if (!is.null(corr_df) && nrow(corr_df) > 0) {
    addWorksheet(wb, "Factor Corr")
    writeData(
      wb,
      "Factor Corr",
      corr_df[, c("Factor1", "Factor2", "Correlation", "SE", "z", "p_value", "Alert"), drop = FALSE]
    )
  }

  if (!is.null(result$modification_indices) && nrow(result$modification_indices) > 0) {
    addWorksheet(wb, "Mod Indices")
    writeData(wb, "Mod Indices", result$modification_indices)
  }

  respec_df <- .cfa_respecifications_to_df(result)
  if (nrow(respec_df) > 0) {
    addWorksheet(wb, "Respec")
    writeData(wb, "Respec", respec_df)
  }

  if (length(warnings_vec) > 0) {
    addWorksheet(wb, "Warnings")
    writeData(wb, "Warnings", data.frame(Warning = warnings_vec, stringsAsFactors = FALSE))
  }

  if (!grepl("\\.xlsx$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".xlsx")
  }

  saveWorkbook(wb, file_name, overwrite = TRUE)
  invisible(normalizePath(file_name))
}

.cfa_export_word <- function(result, file_name, digits = 3) {
  if (!requireNamespace("officer", quietly = TRUE)) {
    stop("Package 'officer' is required to export Word reports.", call. = FALSE)
  }

  read_docx <- officer::read_docx
  body_add_par <- officer::body_add_par

  doc <- read_docx()
  report_lines <- capture.output(print(result, digits = digits))

  for (line in report_lines) {
    doc <- body_add_par(doc, if (nzchar(line)) line else " ", style = "Normal")
  }

  if (!grepl("\\.docx$", file_name, ignore.case = TRUE)) {
    file_name <- paste0(file_name, ".docx")
  }

  print(doc, target = file_name)
  invisible(normalizePath(file_name))
}

#' Export CFA results to Excel or Word
#'
#' @param result Object of class `"cfa_auto"` returned by `cfa_auto()`.
#' @param format `"excel"` or `"word"` (you can provide both).
#' @param file_name Optional base file name without extension. If `NULL`, a
#'   timestamped name is generated automatically.
#' @param digits Number of decimal places used when exporting the Word report.
#'
#' @return The input `result`, invisibly.
#'
#' @examples
#' if (identical(Sys.getenv("PSYCHOMATIC_RUN_SLOW_EXAMPLES"), "true") &&
#'     requireNamespace("semTools", quietly = TRUE) &&
#'     requireNamespace("MVN", quietly = TRUE) &&
#'     requireNamespace("openxlsx", quietly = TRUE)) {
#'   model <- "visual =~ x1 + x2 + x3"
#'   items <- lavaan::HolzingerSwineford1939[1:80, c("x1", "x2", "x3")]
#'   result <- cfa_auto(items, model = model, estimator = "ML", language = "eng")
#'   export_cfa(result, format = "excel", file_name = tempfile("cfa-", fileext = ".xlsx"))
#' }
#'
#' @export
export_cfa <- function(result, format = c("excel", "word"), file_name = NULL, digits = 3) {
  if (!inherits(result, "cfa_auto")) {
    stop("'result' must be an object of class 'cfa_auto'.", call. = FALSE)
  }

  format <- match.arg(format, c("excel", "word"), several.ok = TRUE)

  if (is.null(file_name)) {
    file_name <- paste0("CFA_Results_", base::format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  for (fmt in format) {
    if (fmt == "excel") {
      .cfa_export_excel(result, file_name)
    } else {
      .cfa_export_word(result, file_name, digits = digits)
    }
  }

  invisible(result)
}

# Internal legacy alias for `export_cfa()`
#'
#' @param resultado Deprecated. Use `result`.
#' @param formato Deprecated. Use `format`.
#' @param archivo Deprecated. Use `file_name`.
#' @param digits Number of decimal places used when exporting the Word report.
#'
#' @return The input `resultado`, invisibly.
#' @noRd
exportar_cfa <- function(resultado, formato = c("excel", "word"), archivo = NULL, digits = 3) {
  warning("`exportar_cfa()` is deprecated; use `export_cfa()` instead.", call. = FALSE)
  export_cfa(result = resultado, format = formato, file_name = archivo, digits = digits)
}
