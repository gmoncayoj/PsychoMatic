# ============================================================
#                   FUNCION PRINCIPAL
# ============================================================

#' Automated exploratory factor analysis
#'
#' Runs an automated exploratory factor analysis workflow with estimator
#' selection, parallel analysis, rotation, iterative item screening, reliability,
#' and structured output tables.
#'
#' @param data Data frame or matrix containing numeric item responses.
#' @param rotation Rotation family: `"oblique"` (oblimin) or `"orthogonal"`
#'   (varimax).
#' @param min_loading Minimum acceptable factor loading.
#' @param min_communality Minimum acceptable communality.
#' @param min_cross_loading_diff Minimum difference between the two largest
#'   loadings required to avoid a cross-loading flag.
#' @param max_iter Maximum number of iterative item-removal steps.
#' @param verbose If `TRUE`, prints a console report.
#' @param export_format Optional export format passed to `export_efa()`.
#' @param file_name Optional base file name for exported files.
#' @param language Report language: `"esp"` or `"eng"`.
#' @param seed Optional single integer seed for reproducible random operations.
#'   Defaults to `NULL`, so no seed is set unless explicitly requested.
#' @param ... Deprecated legacy aliases kept temporarily for backward
#'   compatibility. New code should use the English argument names documented
#'   above.
#'
#' @return Object of class `efa_auto` with model results and tidy tables.
#'
#' @examples
#' set.seed(1)
#' factor_score <- rnorm(80)
#' items <- data.frame(
#'   item1 = factor_score + rnorm(80, sd = 0.4),
#'   item2 = factor_score + rnorm(80, sd = 0.4),
#'   item3 = factor_score + rnorm(80, sd = 0.4),
#'   item4 = factor_score + rnorm(80, sd = 0.4)
#' )
#' result <- efa_auto(items, max_iter = 1, verbose = FALSE, language = "eng")
#'
#' @export

efa_auto <- function(data,
                     rotation = c("oblique", "orthogonal"),
                     min_loading = 0.30,
                     min_communality = 0.30,
                     min_cross_loading_diff = 0.15,
                     max_iter = 20,
                     verbose = TRUE,
                     export_format = NULL,
                     file_name = NULL,
                     language = "esp",
                     seed = NULL,
                     ...) {

  data_missing <- missing(data)
  rotation_missing <- missing(rotation)
  min_loading_missing <- missing(min_loading)
  min_communality_missing <- missing(min_communality)
  min_cross_loading_diff_missing <- missing(min_cross_loading_diff)
  export_format_missing <- missing(export_format)
  file_name_missing <- missing(file_name)

  legacy_args <- list(...)
  legacy_names <- names(legacy_args)

  if (length(legacy_args) > 0 && (is.null(legacy_names) || any(!nzchar(legacy_names)))) {
    stop("All additional arguments passed to `efa_auto()` must be named.", call. = FALSE)
  }

  if ("datos" %in% legacy_names) {
    if (!data_missing) {
      stop("Use only one of `data` or deprecated `datos`.", call. = FALSE)
    }
    data <- legacy_args$datos
    legacy_args$datos <- NULL
    .efa_warn_deprecated_arg("datos", "data")
  }

  if ("rotacion" %in% legacy_names) {
    if (!rotation_missing) {
      stop("Use only one of `rotation` or deprecated `rotacion`.", call. = FALSE)
    }
    rotation <- legacy_args$rotacion
    legacy_args$rotacion <- NULL
    .efa_warn_deprecated_arg("rotacion", "rotation")
  }

  if ("carga_min" %in% legacy_names) {
    if (!min_loading_missing) {
      stop("Use only one of `min_loading` or deprecated `carga_min`.", call. = FALSE)
    }
    min_loading <- legacy_args$carga_min
    legacy_args$carga_min <- NULL
    .efa_warn_deprecated_arg("carga_min", "min_loading")
  }

  if ("comunalidad_min" %in% legacy_names) {
    if (!min_communality_missing) {
      stop("Use only one of `min_communality` or deprecated `comunalidad_min`.", call. = FALSE)
    }
    min_communality <- legacy_args$comunalidad_min
    legacy_args$comunalidad_min <- NULL
    .efa_warn_deprecated_arg("comunalidad_min", "min_communality")
  }

  if ("dif_cargas_cruzadas" %in% legacy_names) {
    if (!min_cross_loading_diff_missing) {
      stop("Use only one of `min_cross_loading_diff` or deprecated `dif_cargas_cruzadas`.", call. = FALSE)
    }
    min_cross_loading_diff <- legacy_args$dif_cargas_cruzadas
    legacy_args$dif_cargas_cruzadas <- NULL
    .efa_warn_deprecated_arg("dif_cargas_cruzadas", "min_cross_loading_diff")
  }

  if ("exportar" %in% legacy_names) {
    if (!export_format_missing) {
      stop("Use only one of `export_format` or deprecated `exportar`.", call. = FALSE)
    }
    export_format <- legacy_args$exportar
    legacy_args$exportar <- NULL
    .efa_warn_deprecated_arg("exportar", "export_format")
  }

  if ("archivo" %in% legacy_names) {
    if (!file_name_missing) {
      stop("Use only one of `file_name` or deprecated `archivo`.", call. = FALSE)
    }
    file_name <- legacy_args$archivo
    legacy_args$archivo <- NULL
    .efa_warn_deprecated_arg("archivo", "file_name")
  }

  if (data_missing && is.null(data)) {
    stop("Argument `data` is required.", call. = FALSE)
  }

  if (length(legacy_args) > 0) {
    stop(
      sprintf("Unused argument(s): %s", paste(names(legacy_args), collapse = ", ")),
      call. = FALSE
    )
  }

  language <- .efa_normalize_language(language)
  rotation <- .efa_normalize_rotation(rotation)
  txt <- .efa_labels(language)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
        seed != floor(seed)) {
      stop("Argument `seed` must be NULL or a single non-negative integer.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  if (!requireNamespace("psych", quietly = TRUE))
    stop(txt$required_psych)
  if (!requireNamespace("GPArotation", quietly = TRUE))
    stop(txt$required_gpa)

  sep <- paste0(rep("=", 65), collapse = "")
  rotacion_label <- ifelse(rotation == "oblique", txt$oblique, txt$orthogonal)

  # ================================================================
  # 1. VALIDACION Y PREPARACION DE DATOS
  # ================================================================
  if (!is.data.frame(data) && !is.matrix(data))
    stop(txt$invalid_data)

  datos <- as.data.frame(data)
  cols_num <- sapply(datos, is.numeric)

  if (sum(cols_num) < 3)
    stop(txt$need_numeric)

  if (sum(!cols_num) > 0 && verbose)
    cat(txt$excluded_non_numeric, sum(!cols_num), txt$non_numeric_cols, "\n")

  datos <- datos[, cols_num, drop = FALSE]

  n_original <- nrow(datos)
  datos <- na.omit(datos)
  n <- nrow(datos)
  p <- ncol(datos)

  if (verbose && n_original != n)
    cat(txt$removed_rows, n_original - n, txt$rows_with_na,
        txt$final_sample, n, "\n\n")

  # ================================================================
  # 2. CARACTERISTICAS DE LOS DATOS
  # ================================================================
  n_categorias <- sapply(datos, function(x) length(unique(x)))
  mediana_cat <- median(n_categorias)
  es_ordinal <- mediana_cat <= 5

  if (verbose) {
    cat(sep, "\n")
    cat(txt$analysis_title, "\n")
    cat(sep, "\n\n")
    cat(txt$data_title, "\n")
    cat(txt$sample_size, n, "\n")
    cat(txt$n_items, p, "\n")
    cat(txt$median_categories, mediana_cat, "\n")
    cat(txt$data_type,
        ifelse(es_ordinal, txt$data_type_ordinal, txt$data_type_continuous), "\n")
    cat(txt$categories_by_item, "\n")
    cat_tbl <- data.frame(
      Item = names(n_categorias),
      Categorias = as.integer(n_categorias)
    )
    for (i in seq_len(nrow(cat_tbl)))
      cat(sprintf("    %-25s %d\n", cat_tbl$Item[i], cat_tbl$Categorias[i]))
  }

  # ================================================================
  # 3. NORMALIDAD MULTIVARIADA (Test de Mardia)
  # ================================================================
  normalidad <- .efa_mardia(datos)

  if (verbose) {
    cat("\n", txt$normality_title, "\n", sep = "")
    if (language == "eng") {
      cat(sprintf("  Skewness: chi2 = %.3f, df = %.0f, p = %.4f\n",
                  normalidad$skewness_chi, normalidad$skewness_df,
                  normalidad$skewness_p))
      cat(sprintf("  Kurtosis:  z = %.3f, p = %.4f\n",
                  normalidad$kurtosis_z, normalidad$kurtosis_p))
    } else {
      cat(sprintf("  Asimetria: chi2 = %.3f, gl = %.0f, p = %.4f\n",
                  normalidad$skewness_chi, normalidad$skewness_df,
                  normalidad$skewness_p))
      cat(sprintf("  Curtosis:  z = %.3f, p = %.4f\n",
                  normalidad$kurtosis_z, normalidad$kurtosis_p))
    }
    cat(txt$normality_conclusion,
        ifelse(normalidad$normal, txt$normality_ok, txt$normality_bad))
  }

  # ================================================================
  # 4. SELECCION DEL ESTIMADOR
  # ================================================================
  est <- .efa_seleccionar_estimador(es_ordinal, normalidad$normal, n, p, language)

  if (verbose) {
    cat("\n", txt$estimator_title, "\n", sep = "")
    cat(txt$estimator_label, est$nombre, "\n")
    cat(txt$justification_label, est$justificacion, "\n")
    cat(txt$correlation_label, est$tipo_cor, "\n")
  }

  # ================================================================
  # 5. MATRIZ DE CORRELACIONES
  # ================================================================
  if (est$usar_policorica) {
    if (verbose) cat(txt$calculating_polychoric)
    cor_result <- tryCatch(
      polychoric(datos),
      error = function(e) {
        warning(txt$polychoric_warning, e$message)
        list(rho = cor(datos))
      }
    )
    R <- cor_result$rho
  } else {
    R <- cor(datos)
  }

  # Verificar que R es definida positiva
  eigvals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0)) {
    if (verbose) cat(txt$non_pd_matrix)
    R <- as.matrix(nearPD(R, corr = TRUE)$mat)
  }

  # ================================================================
  # 6. KMO Y TEST DE BARTLETT (sobre items originales)
  # ================================================================
  kmo_original <- KMO(R)
  bartlett_original <- cortest.bartlett(R, n = n)

  if (verbose) {
    cat("\n", txt$adequacy_title, "\n", sep = "")
    if (language == "eng") {
      cat(sprintf("  Overall KMO: %.3f - %s\n",
                  kmo_original$MSA, .efa_interpretar_kmo(kmo_original$MSA, language)))
      cat(sprintf("  Bartlett: chi2 = %.3f, df = %.0f, p = %s\n",
                  bartlett_original$chisq, bartlett_original$df,
                  formatC(bartlett_original$p.value, format = "e", digits = 4)))
    } else {
      cat(sprintf("  KMO global: %.3f - %s\n",
                  kmo_original$MSA, .efa_interpretar_kmo(kmo_original$MSA, language)))
      cat(sprintf("  Bartlett: chi2 = %.3f, gl = %.0f, p = %s\n",
                  bartlett_original$chisq, bartlett_original$df,
                  formatC(bartlett_original$p.value, format = "e", digits = 4)))
    }

    kmo_bajo <- kmo_original$MSAi[kmo_original$MSAi < 0.50]
    if (length(kmo_bajo) > 0) {
      cat(txt$kmo_warning_items)
      for (i in seq_along(kmo_bajo))
        cat(sprintf("    %-25s %.3f\n", names(kmo_bajo)[i], kmo_bajo[i]))
    }
  }

  if (kmo_original$MSA < 0.50)
    warning(txt$kmo_warning_global)

  # ================================================================
  # 7. ANALISIS PARALELO
  # ================================================================
  if (verbose) cat("\n", txt$parallel_title, "\n", sep = "")

  ap <- suppressWarnings(
    fa.parallel(R, n.obs = n, fa = "fa", fm = est$fm, plot = FALSE)
  )
  n_factores <- max(1, ap$nfact)

  if (verbose) cat(txt$optimal_factors, n_factores, "\n")

  # ================================================================
  # 8. SELECCION DE ROTACION
  # ================================================================
  if (n_factores == 1) {
    tipo_rot <- "none"
    if (verbose) cat(txt$no_rotation)
  } else {
    tipo_rot <- ifelse(rotation == "oblique", "oblimin", "varimax")
    if (verbose)
      cat(sprintf(txt$rotation_label, tipo_rot, rotacion_label))
  }

  # ================================================================
  # 9. PROCESO ITERATIVO DE REFINAMIENTO
  # ================================================================
  if (verbose) {
    cat("\n", sep, "\n", sep = "")
    cat(txt$refinement_title, "\n")
    cat(sep, "\n")
  }

  datos_actual <- datos
  R_actual <- R
  items_eliminados <- data.frame(
    Item = character(), Razon = character(),
    Iteracion = integer(), stringsAsFactors = FALSE
  )
  iteracion <- 0
  efa <- NULL

  while (iteracion < max_iter) {
    iteracion <- iteracion + 1
    p_actual <- ncol(datos_actual)

    # Recalcular correlaciones y numero de factores tras eliminacion
    if (iteracion > 1) {
      if (est$usar_policorica) {
        R_actual <- tryCatch(
          polychoric(datos_actual)$rho,
          error = function(e) cor(datos_actual)
        )
      } else {
        R_actual <- cor(datos_actual)
      }

      # Asegurar definida positiva
      eigvals_act <- eigen(R_actual, symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals_act <= 0))
        R_actual <- as.matrix(nearPD(R_actual, corr = TRUE)$mat)

      ap_new <- suppressWarnings(
        fa.parallel(R_actual, n.obs = n, fa = "fa", fm = est$fm, plot = FALSE)
      )
      nf_new <- max(1, ap_new$nfact)

      if (nf_new != n_factores) {
        if (verbose)
          cat(sprintf(txt$factors_adjusted, n_factores, nf_new))
        n_factores <- nf_new
        tipo_rot <- if (n_factores == 1) "none" else
          ifelse(rotation == "oblique", "oblimin", "varimax")
      }
    }

    # Verificar items suficientes
    if (p_actual < n_factores * 3) {
      if (n_factores > 1) {
        n_factores <- max(1, floor(p_actual / 3))
        tipo_rot <- if (n_factores == 1) "none" else
          ifelse(rotation == "oblique", "oblimin", "varimax")
        if (verbose)
          cat(txt$insufficient_items, n_factores, txt$factor_suffix)
      } else {
        break
      }
    }

    # Ejecutar AFE
    efa <- tryCatch(
      suppressWarnings(
        fa(R_actual, nfactors = n_factores, n.obs = n,
           rotate = tipo_rot, fm = est$fm)
      ),
      error = function(e) {
        if (est$fm != "minres") {
          if (verbose)
            cat(txt$estimator_error_note, est$fm, txt$switching_minres)
          suppressWarnings(
            fa(R_actual, nfactors = n_factores, n.obs = n,
               rotate = tipo_rot, fm = "minres")
          )
        } else {
          stop(if (language == "eng") "Error in EFA: " else "Error en AFE: ",
               e$message)
        }
      }
    )

    cargas <- as.data.frame(unclass(efa$loadings))
    comunalidades <- efa$communality

    if (verbose)
      cat(sprintf(txt$iteration_title, iteracion, p_actual, n_factores))

    # Evaluar items
    problemas <- .efa_evaluar_items(
      cargas, comunalidades, min_loading,
      min_communality, min_cross_loading_diff, n_factores, language
    )

    if (nrow(problemas) == 0) {
      if (verbose) cat(txt$no_problem_items)
      break
    }

    if (verbose) {
      cat(txt$problematic_items)
      for (i in 1:nrow(problemas))
        cat(sprintf("    - %-20s %s\n",
                    problemas$Item[i], problemas$Razon[i]))
    }

    # Eliminar el item mas problematico
    item_out <- problemas$Item[1]
    razon_out <- problemas$Razon[1]

    items_eliminados <- rbind(items_eliminados, data.frame(
      Item = item_out, Razon = razon_out,
      Iteracion = iteracion, stringsAsFactors = FALSE
    ))

    datos_actual <- datos_actual[, colnames(datos_actual) != item_out,
                                 drop = FALSE]

    if (verbose)
      cat(sprintf(txt$removed_item, item_out, razon_out))

    if (ncol(datos_actual) < 3) {
      warning(txt$fewer_than_three)
      break
    }
  }

  if (is.null(efa))
    stop(txt$no_model)

  # ================================================================
  # 10. RESULTADOS FINALES
  # ================================================================
  cargas_final <- as.data.frame(unclass(efa$loadings))
  asignacion <- apply(abs(cargas_final), 1, which.max)

  distribucion <- data.frame(
    Item = rownames(cargas_final),
    Dimension = paste0("F", asignacion),
    Carga_Principal = sapply(seq_along(asignacion), function(i)
      cargas_final[i, asignacion[i]]),
    Comunalidad = efa$communality[rownames(cargas_final)],
    stringsAsFactors = FALSE
  )
  distribucion <- distribucion[order(distribucion$Dimension,
                                     -abs(distribucion$Carga_Principal)), ]
  rownames(distribucion) <- NULL

  var_explicada <- efa$Vaccounted
  chi2 <- if (!is.null(efa$STATISTIC)) efa$STATISTIC else NA
  gl   <- if (!is.null(efa$dof)) efa$dof else NA
  p_chi2 <- if (!is.null(efa$PVAL)) efa$PVAL else NA

  # Recalcular KMO y Bartlett sobre items finales
  R_final <- R_actual
  kmo_final <- KMO(R_final)
  bartlett_final <- cortest.bartlett(R_final, n = n)

  if (verbose) {
    cat("\n", sep, "\n", sep = "")
    cat(txt$final_results_title, "\n")
    cat(sep, "\n\n")

    cat(txt$fit_indices_title, "\n")
    if (!is.na(chi2))
      cat(sprintf(txt$chi_square, chi2))
    if (!is.na(gl))
      cat(sprintf(txt$degrees_freedom, as.integer(gl)))
    if (!is.na(p_chi2))
      cat(sprintf(txt$p_value, formatC(p_chi2, format = "e", digits = 4)))
    if (!is.null(efa$TLI))
      cat(sprintf(txt$tli, efa$TLI))
    if (!is.null(efa$RMSEA))
      cat(sprintf(txt$rmsea, efa$RMSEA[1], efa$RMSEA[2], efa$RMSEA[3]))
    if (!is.null(efa$BIC))
      cat(sprintf(txt$bic, efa$BIC))

    cat(sprintf(txt$final_kmo,
                kmo_final$MSA, .efa_interpretar_kmo(kmo_final$MSA, language)))
    cat(sprintf(txt$final_bartlett,
                bartlett_final$chisq, bartlett_final$df,
                formatC(bartlett_final$p.value, format = "e", digits = 4)))

    # --- Varianza explicada ---
    cat("\n", txt$variance_title, "\n", sep = "")
    if (is.matrix(var_explicada)) {
      for (j in 1:ncol(var_explicada)) {
        if ("Proportion Var" %in% rownames(var_explicada))
          cat(sprintf(txt$factor_variance, j,
                      var_explicada["Proportion Var", j] * 100))
      }
      if ("Cumulative Var" %in% rownames(var_explicada))
        cat(sprintf(txt$total_variance,
                    var_explicada["Cumulative Var",
                                  ncol(var_explicada)] * 100))
    }

    # --- Matriz de cargas factoriales ---
    cat("\n", txt$loadings_title, "\n", sep = "")
    cargas_print <- cargas_final
    for (i in 1:nrow(cargas_print))
      for (j in 1:ncol(cargas_print))
        cargas_print[i, j] <- ifelse(
          abs(cargas_final[i, j]) < min_loading,
          NA_real_,
          round(cargas_final[i, j], 3)
        )
    cargas_print$h2 <- round(efa$communality[rownames(cargas_final)], 3)
    print(cargas_print, na.print = "   ")
    cat(txt$loadings_note_left, min_loading, txt$loadings_note_right)

    # --- Distribucion de items ---
    cat("\n", txt$distribution_title, "\n", sep = "")
    for (f in sort(unique(distribucion$Dimension))) {
      items_f <- distribucion[distribucion$Dimension == f, ]
      cat(sprintf(txt$distribution_factor, f, nrow(items_f)))
      for (i in 1:nrow(items_f))
        cat(sprintf(txt$distribution_item,
                    items_f$Item[i], items_f$Carga_Principal[i],
                    items_f$Comunalidad[i]))
    }

    # --- Correlaciones entre factores ---
    if (rotation == "oblique" && n_factores > 1 && !is.null(efa$Phi)) {
      cat("\n", txt$factor_corr_title, "\n", sep = "")
      print(round(efa$Phi, 3))
    }

    # --- Items eliminados ---
    if (nrow(items_eliminados) > 0) {
      cat("\n", txt$removed_items_title, "\n", sep = "")
      for (i in 1:nrow(items_eliminados))
        cat(sprintf(txt$removed_items_line,
                    i, items_eliminados$Item[i],
                    items_eliminados$Razon[i],
                    items_eliminados$Iteracion[i]))
      cat(sprintf(txt$removed_items_total,
                  nrow(items_eliminados), p))
    } else {
      cat(txt$no_items_removed)
    }
  }

  # ================================================================
  # 11. CONFIABILIDAD POR DIMENSION
  # ================================================================
  if (verbose) cat("\n", txt$reliability_title, "\n", sep = "")

  fiabilidad <- list()
  for (f in sort(unique(distribucion$Dimension))) {
    items_f <- distribucion$Item[distribucion$Dimension == f]

    if (length(items_f) >= 3) {
      alpha_f <- tryCatch(
        alpha(datos_actual[, items_f, drop = FALSE], check.keys = TRUE),
        error = function(e) NULL
      )

      omega_val <- tryCatch({
        invisible(capture.output(
          om <- omega(datos_actual[, items_f, drop = FALSE],
                      nfactors = 1, plot = FALSE)
        ))
        om$omega.tot
      }, error = function(e) NA_real_)

      if (!is.null(alpha_f)) {
        fiabilidad[[f]] <- list(
          alpha = alpha_f$total$raw_alpha,
          omega = omega_val,
          items = items_f,
          n_items = length(items_f),
          alpha_drop = alpha_f$alpha.drop
        )

        if (verbose) {
          if (language == "eng") {
            cat(sprintf("\n  %s (%d items): alpha = %.3f",
                        f, length(items_f), alpha_f$total$raw_alpha))
          } else {
            cat(sprintf("\n  %s (%d items): alpha = %.3f",
                        f, length(items_f), alpha_f$total$raw_alpha))
          }
          if (!is.na(omega_val))
            cat(sprintf(", omega = %.3f", omega_val))
          cat("\n")

          # Items cuya eliminacion mejoraria alpha
          if (length(items_f) > 3 && !is.null(alpha_f$alpha.drop)) {
            mejora <- alpha_f$alpha.drop[
              alpha_f$alpha.drop$raw_alpha > alpha_f$total$raw_alpha,
              , drop = FALSE
            ]
            if (nrow(mejora) > 0) {
              cat(txt$alpha_improve)
              for (r in rownames(mejora))
                cat(sprintf(txt$alpha_if_removed,
                            r, mejora[r, "raw_alpha"],
                            mejora[r, "raw_alpha"] - alpha_f$total$raw_alpha))
            }
          }
        }
      } else {
        fiabilidad[[f]] <- list(
          alpha = NA, omega = NA,
          items = items_f, n_items = length(items_f)
        )
        if (verbose)
          cat(sprintf(txt$reliability_error, f))
      }

    } else if (length(items_f) == 2) {
      # Con 2 items, solo Spearman-Brown
      cor_2 <- cor(datos_actual[, items_f[1]], datos_actual[, items_f[2]],
                   use = "complete.obs")
      sb <- (2 * abs(cor_2)) / (1 + abs(cor_2))
      fiabilidad[[f]] <- list(
        alpha = NA, omega = NA, spearman_brown = sb,
        items = items_f, n_items = 2
      )
      if (verbose)
        cat(sprintf(txt$spearman_brown, f, sb, cor_2))
    } else {
      fiabilidad[[f]] <- list(
        alpha = NA, omega = NA,
        items = items_f, n_items = length(items_f)
      )
      if (verbose)
        cat(sprintf(txt$reliability_na, f, length(items_f)))
    }
  }

  # ================================================================
  # 12. VEREDICTO FINAL
  # ================================================================
  sugerencias <- .efa_sugerencias(
    efa, distribucion, fiabilidad, kmo_final, bartlett_final,
    n_factores, n, items_eliminados, min_loading, var_explicada, p, language
  )

  if (verbose) {
    cat("\n", sep, "\n", sep = "")
    if (length(sugerencias) == 0) {
      cat(txt$optimal_solution_title)
      cat(txt$optimal_solution_1)
      cat(txt$optimal_solution_2)
    } else {
      cat(txt$respec_title)
      cat(sep, "\n\n")
      for (i in seq_along(sugerencias))
        cat(sprintf("  %d. %s\n", i, sugerencias[[i]]))
      cat("\n")
    }
    cat(sep, "\n")

    # --- Referencias bibliograficas ---
    cat(txt$refs_title)
    cat(txt$refs_intro)

    cat(txt$refs_estimator)
    cat("    Lloret-Segura, S., Ferreres-Traver, A., Hernandez-Baeza, A., &\n")
    cat("      Tomas-Marco, I. (2014). El analisis factorial exploratorio de los\n")
    cat("      items: una guia practica, revisada y actualizada. Anales de\n")
    cat("      Psicologia, 30(3), 1151-1169.\n")
    cat("      https://doi.org/10.6018/analesps.30.3.199361\n\n")

    cat("    Ferrando, P. J., & Lorenzo-Seva, U. (2017). Program FACTOR at 10:\n")
    cat("      Origins, development and future directions. Psicothema, 29(2),\n")
    cat("      236-240. https://doi.org/10.7334/psicothema2016.304\n\n")

    cat(txt$refs_polychoric)
    cat("    Holgado-Tello, F. P., Chacon-Moscoso, S., Barbero-Garcia, I., &\n")
    cat("      Vila-Abad, E. (2010). Polychoric versus Pearson correlations in\n")
    cat("      exploratory and confirmatory factor analysis of ordinal variables.\n")
    cat("      Quality & Quantity, 44(1), 153-166.\n")
    cat("      https://doi.org/10.1007/s11135-008-9190-y\n\n")

    cat(txt$refs_parallel)
    cat("    Horn, J. L. (1965). A rationale and test for the number of factors\n")
    cat("      in factor analysis. Psychometrika, 30(2), 179-185.\n")
    cat("      https://doi.org/10.1007/BF02289447\n\n")

    cat("    Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality\n")
    cat("      assessment of ordered polytomous items with parallel analysis.\n")
    cat("      Psychological Methods, 16(2), 209-220.\n")
    cat("      https://doi.org/10.1037/a0023353\n\n")

    cat(txt$refs_normality)
    cat("    Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis\n")
    cat("      with applications. Biometrika, 57(3), 519-530.\n")
    cat("      https://doi.org/10.1093/biomet/57.3.519\n\n")

    cat(txt$refs_kmo)
    cat("    Kaiser, H. F. (1974). An index of factorial simplicity. Psychometrika,\n")
    cat("      39(1), 31-36. https://doi.org/10.1007/BF02291575\n\n")

    cat(txt$refs_elimination)
    cat("    Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019).\n")
    cat("      Multivariate data analysis (8th ed.). Cengage.\n\n")

    cat("    Watkins, M. W. (2018). Exploratory factor analysis: A guide to best\n")
    cat("      practice. Journal of Black Psychology, 44(3), 219-246.\n")
    cat("      https://doi.org/10.1177/0095798418771807\n\n")

    cat("    Howard, M. C. (2016). A review of exploratory factor analysis\n")
    cat("      decisions and overview of current practices. International Journal\n")
    cat("      of Human-Computer Interaction, 32(1), 51-62.\n")
    cat("      https://doi.org/10.1080/10447318.2015.1087664\n\n")

    cat("    Costello, A. B., & Osborne, J. W. (2005). Best practices in\n")
    cat("      exploratory factor analysis: Four recommendations for getting the\n")
    cat("      most from your analysis. Practical Assessment, Research, and\n")
    cat("      Evaluation, 10(7), 1-9. https://doi.org/10.7275/jyj1-4868\n\n")

    cat(txt$refs_reliability)
    cat("    McDonald, R. P. (1999). Test theory: A unified treatment. Lawrence\n")
    cat("      Erlbaum Associates.\n\n")

    cat("    Dunn, T. J., Baguley, T., & Brunsden, V. (2014). From alpha to\n")
    cat("      omega: A practical solution to the pervasive problem of internal\n")
    cat("      consistency estimation. British Journal of Psychology, 105(3),\n")
    cat("      399-412. https://doi.org/10.1111/bjop.12046\n\n")

    cat(txt$refs_software)
    cat("    Revelle, W. (2024). psych: Procedures for psychological,\n")
    cat("      psychometric, and personality research. R package.\n")
    cat("      https://CRAN.R-project.org/package=psych\n\n")

    cat(sep, "\n")
  }

  # ================================================================
  # OBJETO DE RETORNO
  # ================================================================
  resultado <- list(
    efa                  = efa,
    cargas               = cargas_final,
    comunalidades        = efa$communality,
    distribucion         = distribucion,
    language             = language,
    seed                 = seed,
    n_factores           = n_factores,
    estimador            = est,
    rotation             = tipo_rot,
    rotacion             = tipo_rot,
    kmo_original         = kmo_original,
    kmo_final            = kmo_final,
    bartlett_original    = bartlett_original,
    bartlett_final       = bartlett_final,
    chi_cuadrado         = chi2,
    grados_libertad      = gl,
    p_valor              = p_chi2,
    varianza_explicada   = var_explicada,
    normalidad           = normalidad,
    fiabilidad           = fiabilidad,
    items_eliminados     = items_eliminados,
    correlacion_factores = if (!is.null(efa$Phi)) efa$Phi else NULL,
    sugerencias          = sugerencias,
    solucion_optima      = length(sugerencias) == 0,
    datos_finales        = datos_actual,
    n                    = n,
    p_original           = p,
    p_final              = ncol(datos_actual)
  )

  resultado <- .efa_add_tidy_fields(resultado, call = match.call())
  class(resultado) <- "efa_auto"

  # Exportar si se solicito
  if (!is.null(export_format)) {
    export_efa(resultado, format = export_format, file_name = file_name)
  }

  invisible(resultado)
}
