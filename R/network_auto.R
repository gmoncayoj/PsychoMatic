#' Automated exploratory psychometric network analysis
#'
#' Optionally runs unique variable analysis before EGA (by default Spearman,
#' EBICglasso and consensus Louvain), revised network loadings, TEFI, wording and hierarchy diagnostics,
#' structural bootstrap and group comparisons. Optional packages are required
#' only for requested stages. See the network workflow vignette for assumptions.
#' @param data Numeric/ordered item data or a `psychomatic_spec`.
#' @param items Item columns; defaults to eligible columns, excluding `group`.
#' @param group Group column (e.g., sex), or declared in the specification.
#' @param reverse Known reverse-keyed item columns. Never inferred from signs.
#' @param ranges Theoretical response ranges for reverse items.
#' @param uva NULL (default) skips UVA entirely. `"diagnostic"` reports redundant
#'   pairs, wTO and justified removal proposals without removing items.
#'   `"removed"` authorizes native representative selection and removal, with
#'   recorded criteria, scores and reasons. No mode requests interactive input.
#' @param uva_cutoff Weighted topological overlap threshold, default .25.
#' @param corr `"spearman"` (default), `"pearson"`, or `"polychoric"` for ordinal
#'   items (including binary items), up to 11 categories per item in EGAnet.
#'   Correlations are re-estimated in resamples.
#' @param gamma EBIC penalty, default .5.
#' @param model Network estimator: `"glasso"` (default), `"TMFG"`, or `"BGGM"`.
#'   BGGM requires an EGAnet version that supports it; EGAnet 2.4.1 disables it.
#'   No unavailable estimator is silently replaced. Gamma applies only to GLASSO.
#' @param algorithm Community detection: `"louvain"` (default), `"leiden"`,
#'   or `"walktrap"`. Used in empirical/group EGA and both hierEGA orders.
#' @param bootstrap_type bootEGA bootstrap: `"resampling"` (default) or
#'   `"parametric"`. Applies also to configural, riEGA and hierEGA bootstraps.
#'   Parametric generates continuous multivariate-normal data; it does not
#'   preserve ordinal response categories. No automatic superiority selection.
#'   bootnet edge precision and case-drop procedures remain separate.
#' @param missing `"listwise"` (default) or `"pairwise"` for correlations.
#' @param wording Run riEGA on original, unrecoded responses.
#' @param hierarchy Run hierEGA and evaluate identifiable CFA candidates.
#' @param bootstrap Run bootEGA and bootnet accuracy/case-drop analyses.
#' @param invariance Run NCT and EGA loading permutation tests if a group exists.
#' @param n_boot Replications per bootstrap analysis, default 1000.
#' @param n_perm Permutations per comparison, default 5000.
#' @param consensus_iter Louvain consensus iterations, default 1000.
#' @param n_cores Workers for bootEGA, bootnet and loading permutations, default 1.
#'   NCT remains sequential in the supported version.
#' @param seed Optional integer-valued seed. NULL (default) draws an analysis seed
#'   from the session RNG and records it in `effective_settings$seed`. An explicit
#'   seed is reproducible and restores the caller's RNG state.
#' @param language `"eng"` (default) or `"esp"`.
#' @param report Export a report, logical or output path.
#' @param format Report format(s), `"excel"` and/or `"word"`.
#' @param file_name Report output basename/path.
#' @param overwrite Allow overwriting existing report files, default FALSE.
#' @param digits Display precision; numerical results retain full precision.
#' @param verbose Print stage progress, elapsed time and the analysis summary.
#' @return A `psychomatic_network` and `psychomatic_result` object. Failed optional
#'   stages remain explicit in `results$stages` and `errors`; partial results are
#'   never described as a completed analysis. Original responses remain available
#'   in `diagnostics$input` but are excluded from reader-facing reports.
#' @details The default console and reports emphasize communities, structural
#'   consistency and group comparisons. Non-significance does not establish
#'   equivalence. Complete diagnostic tables remain in the object and can be
#'   displayed with `print(fit, detail="audit")` or exported with
#'   `psychomatic_report(fit, detail="audit")`. APA reference italics are shared
#'   across analyses in Word, Excel and compatible interactive consoles; plain
#'   logs retain unformatted text. `options(PsychoMatic.ansi=FALSE)` disables
#'   console styling; `TRUE` forces it for terminals supporting ANSI italics.
#' @examples
#' \dontrun{
#' fit <- network_auto(MY_DATA, items = 1:12, group = "sex")
#' plot(fit, type = "network")
#' psychomatic_report(fit, format = c("word", "excel"), file_name = "network")
#' }
#' @export
network_auto <- function(data, items = NULL, group = NULL, reverse = NULL, ranges = NULL,
    uva = NULL, uva_cutoff = .25,
    corr = c("spearman", "pearson", "polychoric"), gamma = .5, missing = c("listwise", "pairwise"),
    wording = TRUE, hierarchy = TRUE, bootstrap = TRUE, invariance = TRUE,
    n_boot = 1000L, n_perm = 5000L, consensus_iter = 1000L, n_cores = 1L,
    seed = NULL, language = "eng", report = FALSE, format = "excel", file_name = NULL,
    overwrite = FALSE, digits = 3L, verbose = TRUE,
    model = c("glasso", "TMFG", "BGGM"),
    algorithm = c("louvain", "leiden", "walktrap"),
    bootstrap_type = c("resampling", "parametric")) {
  call <- match.call()
  if (!is.null(uva) && (!is.character(uva) || length(uva)!=1L ||
      is.na(uva) || !uva %in% c("diagnostic","removed")))
    stop('`uva` must be NULL, "diagnostic", or "removed".',call.=FALSE)
  corr <- match.arg(corr); missing <- match.arg(missing)
  model <- match.arg(model); algorithm <- match.arg(algorithm)
  bootstrap_type <- match.arg(bootstrap_type)
  for (arg in c("n_boot", "n_perm", "consensus_iter", "n_cores"))
    .pm_validate_scalar_number(get(arg), arg, 1, .Machine$integer.max, integer = TRUE)
  if (!is.null(seed)) {
    .pm_validate_scalar_number(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
    # EGAnet's compiled reproducible_seeds() accepts double, not integer storage.
    seed <- as.double(seed)
  }
  .pm_validate_scalar_number(gamma, "gamma", 0, 1)
  .pm_validate_scalar_number(uva_cutoff, "uva_cutoff", 0, 1)
  for (arg in c("wording", "hierarchy", "bootstrap", "invariance"))
    if (!is.logical(get(arg)) || length(get(arg)) != 1L || is.na(get(arg)))
      stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  output <- .pm_output_options(report, format, file_name, overwrite, digits, verbose)
  language <- .pm_normalize_language(language)
  .pm_network_require("EGAnet")
  if (bootstrap) .pm_network_require("bootnet")
  if (invariance && (!is.null(group) || (inherits(data,"psychomatic_spec") && !is.null(data$group))))
    .pm_network_require("NetworkComparisonTest")
  input <- .pm_network_input(data, items, group, reverse, ranges, missing)
  if(corr=="polychoric" && any(input$spec$measurement_level[input$spec$items] != "ordinal"))
    stop('`corr="polychoric"` requires ordinal items. Declare their levels in psychomatic_spec; continuous items cannot be discretized automatically.',call.=FALSE)
  # Native EGAnet NULL seeds can use a generator independent of R's RNG. Derive
  # one seed from R, then pass an explicit double to every native stage instead.
  analysis_seed <- if (is.null(seed)) as.double(sample.int(.Machine$integer.max, 1L)) else seed
  settings <- list(corr = corr, gamma = gamma, missing = missing, model = model,
    algorithm = algorithm, bootstrap_type = bootstrap_type,
    consensus_iter = consensus_iter, seed = analysis_seed,
    seed_source = if (is.null(seed)) "session_generated" else "explicit",
    n_boot = n_boot, n_perm = n_perm, n_cores = n_cores, language = language,
    uva = uva, uva_cutoff = uva_cutoff, reverse = input$spec$reverse,
    verbose=verbose,ordinal_categories=if(corr=="polychoric")
      .pm_network_corr_args(corr,input$keyed)$ordinal.categories else NULL)
  requested_settings <- settings
  requested_settings["seed"] <- list(seed)
  .pm_network_seed(analysis_seed, {
    stages <- data.frame(stage = character(), status = character(), reason = character(),elapsed_seconds=numeric())
    interrupted <- FALSE
    warns <- errs <- list()
    warning_stages <- error_stages <- character()
    run <- function(name, expr, required = FALSE) {
      if(interrupted) {
        stages <<- rbind(stages,data.frame(stage=name,status="skipped",
          reason="Skipped after user interruption.",elapsed_seconds=0))
        return(NULL)
      }
      .pm_network_progress(settings,paste(if(language=="esp")"Iniciando"else"Starting",name))
      started <- proc.time()[["elapsed"]]
      capture <- tryCatch(.pm_capture_conditions(expr),interrupt=function(e)
        list(value=NULL,warnings=list(),messages=list(),error=e))
      elapsed <- proc.time()[["elapsed"]]-started
      warns <<- c(warns, capture$warnings)
      warning_stages <<- c(warning_stages, rep(name, length(capture$warnings)))
      if (!is.null(capture$error)) {
        if (required) stop(capture$error)
        interrupted <<- .pm_network_interrupted(capture$error)
        errs <<- c(errs, list(capture$error))
        error_stages <<- c(error_stages, name)
      }
      not_evaluable <- is.list(capture$value) && identical(capture$value$status, "not_evaluable")
      stages <<- rbind(stages, data.frame(stage = name,
        status = if (interrupted) "interrupted" else if (!is.null(capture$error)) "failed" else if (not_evaluable) "not_evaluable" else "completed",
        reason = if (!is.null(capture$error)) conditionMessage(capture$error) else
          if (not_evaluable) capture$value$reason else "",elapsed_seconds=elapsed))
      .pm_network_progress(settings,sprintf("%s %s (%.1f s; %s)",
        if(language=="esp")"Finaliz\u00f3"else"Finished",name,elapsed,tail(stages$status,1)))
      capture$value
    }
    uv <- NULL
    if(is.null(uva)) {
      stages <- rbind(stages,data.frame(stage="UVA",status="skipped",
        reason="UVA was not requested (uva=NULL); redundancy was not assessed.",elapsed_seconds=0))
    } else uv <- run("UVA", do.call(EGAnet::UVA,c(list(data = as.matrix(input$keyed),
      gamma = gamma, na.data = missing, cut.off = uva_cutoff,
      reduce = TRUE, reduce.method = "remove", auto = TRUE),
      .pm_network_corr_args(corr,input$keyed,settings$ordinal_categories))), TRUE)
    selection <- .pm_network_uva_select(uv,names(input$keyed),uva,uva_cutoff)
    retained <- selection$retained
    settings["uva"] <- list(selection$effective_mode)
    settings$uva_decision <- selection$decision
    if (length(retained) < 3L) stop("UVA retained fewer than three items; review redundancy or use uva='diagnostic'.", call. = FALSE)
    keyed <- input$keyed[retained]; raw <- input$raw[retained]
    ega <- run("EGA", .pm_network_ega(keyed, settings), TRUE)
    if (all(is.na(ega$wc))) stop("EGA found no evaluable community structure. Review sample size, correlations and item selection.", call. = FALSE)
    loads <- run("loadings", .pm_network_loadings(ega), TRUE)
    centrality <- run("centrality", .pm_network_centrality(ega$network), TRUE)
    tefi <- data.frame(model = "EGA", index = "TEFI", value = as.numeric(ega$TEFI),
      comparison_set = "keyed_retained", comparable = TRUE)
    word <- if (wording) run("riEGA", .pm_network_wording(raw, keyed, ega, settings)) else NULL
    hier <- if (hierarchy) run("hierarchy", .pm_network_hierarchy(keyed, ega, input$spec, settings)) else NULL
    boots <- if (bootstrap) run("bootstrap", .pm_network_bootstrap(keyed, ega, settings)) else NULL
    if (bootstrap && !is.null(word)) {
      if (identical(word$status, "not_evaluable")) {
        run("riEGA stability", .pm_not_evaluable(
          "riEGA stability was not computed because the empirical residual partition is not evaluable."))
      } else word$stability <- run("riEGA stability",
        .pm_network_boot_summary(.pm_network_boot_ega(raw,settings,"riEGA"), word$adjusted$EGA$wc))
    }
    if (bootstrap && !is.null(hier)) {
      hier$stability <- run("hierEGA stability", .pm_network_hier_bootstrap(keyed,settings))
      if (hier$recommendation != "indeterminate") {
        required_orders <- if(hier$recommendation=="correlated")"lower_order"else c("lower_order","higher_order")
        stability <- hier$stability$dimensions
        supported <- !is.null(stability) && nrow(stability)>0 &&
          all(is.finite(stability$structural_consistency[stability$order %in% required_orders]) &
            stability$structural_consistency[stability$order %in% required_orders] >= .70)
        if (!supported) {
          hier$fit_candidate <- hier$recommendation
          hier$recommendation <- "indeterminate"
          hier$reason <- "Fit candidate lacked >=.70 structural consistency in required bootstrap orders."
        }
      }
    }
    inv <- if (invariance && !is.null(input$groups))
      run("invariance", .pm_network_invariance(keyed, input$groups, ega, settings)) else NULL
    if (!is.null(word)) tefi <- rbind(tefi, word$tefi)
    if (!is.null(hier) && !is.null(hier$tefi)) tefi <- rbind(tefi, hier$tefi)
    decisions <- rbind(
      .pm_decision("uva", if(is.null(uva))"NULL"else uva, paste(retained, collapse = ", "),
        paste(.pm_network_uva_decision_text(selection$decision),
          "All groups share the same retained items. Bootstrap inference uses this fixed item set.")),
      .pm_decision("estimator", paste(corr,model,algorithm), .pm_network_method(settings),
        "Prespecified estimator; no claim of universal superiority. Correlations are recomputed inside bootstrap and permutations."),
      .pm_decision("bootstrap_type", bootstrap_type, bootstrap_type,
        "User-selected bootEGA generator; no automatic superiority selection. Parametric simulates continuous normal responses; resampling draws observed rows with replacement. bootnet uses nonparametric and case-drop bootstraps separately."),
      .pm_decision("wording", paste(input$spec$reverse, collapse = ", "), "original raw responses for riEGA",
        "Only declared keys and theoretical ranges reverse items; signs are diagnostic, not semantic evidence."))
    result <- .pm_result("psychomatic_network", call, requested_settings, settings,
      list(n = nrow(keyed), p_original = ncol(input$keyed), p_retained = ncol(keyed),
        omitted = input$omitted, group = input$spec$group),
      list(uva = uv, uva_selection = selection, retained_items = retained,
        removed_items = selection$removed,
        ega = ega, loadings = loads, centrality = centrality, tefi = tefi,
        wording = word, hierarchy = hier, bootstrap = boots, invariance = inv, stages = stages),
      warns, errs, list(input = input, versions = .pm_network_versions(),
        limitations = .pm_network_limitations()), decisions,
      list(status = if (any(!stages$status %in% c("completed","skipped"))) "partial" else "completed",
        dimensions = ega$n.dim))
    result$warnings$stage <- warning_stages
    result$errors$stage <- error_stages
    result$language <- language
    result$references <- .pm_network_references()
    .pm_finish_analysis(result, output, language)
  })
}

#' Print a psychometric network result
#' @param x A `psychomatic_network` object.
#' @param detail `"report"` (default) displays a manuscript-oriented summary;
#'   `"audit"` displays complete diagnostic tables and conditions.
#' @param ... Additional arguments, currently unused.
#' @return `x`, invisibly.
#' @method print psychomatic_network
#' @export
print.psychomatic_network <- function(x, detail = c("report", "audit"), ...) {
  .pm_console_report(x, .pm_null_default(attr(x, "display_digits"), 3L), x$language,
    detail = match.arg(detail))
  invisible(x)
}
