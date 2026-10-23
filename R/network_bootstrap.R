.pm_network_boot_ega <- function(data, settings, variant = "EGA") {
  args <- c(list(data = as.matrix(data), na.data = settings$missing,
    seed = settings$seed,
    iter = settings$n_boot, ncores = settings$n_cores,
    type = .pm_null_default(settings$bootstrap_type,"resampling"),
    EGA.type = variant, plot.itemStability = FALSE, typicalStructure = FALSE, verbose = FALSE),
    .pm_network_native_options(settings,hierarchy=variant=="hierEGA"),
    .pm_network_corr_args(settings$corr,data,settings$ordinal_categories))
  if (variant == "hierEGA") {
    # hierEGA resets consensus.iter to 1000 when `algorithm` is supplied.
    args$algorithm <- NULL
    args <- c(args,list(loading.method="revised",scores="network",impute="none"))
  }
  tryCatch(do.call(EGAnet::bootEGA,args),error=function(e) {
    message <- conditionMessage(e)
    if(variant=="riEGA" && grepl("invalid 'use' argument",message,fixed=TRUE))
      stop("Installed EGAnet could not evaluate riEGA bootstrap with the selected options: invalid native correlation missing-data argument. The main EGA bootstrap is separate; no estimator or bootstrap type was substituted.",call.=FALSE)
    if(variant=="hierEGA" && identical(settings$model,"TMFG") && grepl("invalid 'nrow' value",message,fixed=TRUE))
      stop("TMFG requires at least four nodes in an estimated network order. A hierarchical bootstrap replicate could not be estimated; consider GLASSO or hierarchy=FALSE. No estimator was substituted.",call.=FALSE)
    stop(e)
  })
}

.pm_network_hier_bootstrap <- function(data, settings) {
  native <- .pm_network_boot_ega(data,settings,"hierEGA")
  tables <- lapply(c("lower_order","higher_order"),function(order) {
    s <- native$stability[[order]]
    consistency <- s$dimension.stability$structural.consistency
    item <- s$item.stability$item.stability$empirical.dimensions
    list(dimensions=data.frame(order=order,dimension=names(consistency),structural_consistency=as.numeric(consistency)),
      items=data.frame(order=order,item=names(item),stability=as.numeric(item)),
      requested=settings$n_boot,valid=sum(is.finite(native[[order]]$boot.ndim)))
  })
  names(tables) <- c("lower_order","higher_order")
  list(native=native,dimensions=do.call(rbind,lapply(tables,`[[`,"dimensions")),
    items=do.call(rbind,lapply(tables,`[[`,"items")),
    counts=data.frame(order=names(tables),requested=settings$n_boot,valid=vapply(tables,`[[`,numeric(1),"valid")))
}

.pm_network_boot_summary <- function(native, structure = NULL) {
  valid <- is.finite(native$boot.ndim) & native$boot.ndim > 0 &
    apply(native$boot.wc, 1L, function(wc) any(!is.na(wc)))
  requested <- native$iter
  if (!any(valid)) stop("No usable bootstrap community partitions; structural stability is not evaluable.", call. = FALSE)
  usable <- native
  for (field in c("bootGraphs", "bootCorrs", "boot.ndim", "TEFI"))
    usable[[field]] <- native[[field]][valid]
  usable$boot.wc <- native$boot.wc[valid, , drop = FALSE]
  usable$iter <- sum(valid)
  dim <- EGAnet::dimensionStability(usable, IS.plot = FALSE, structure = structure)
  item <- dim$item.stability$item.stability$empirical.dimensions
  consistency <- dim$dimension.stability$structural.consistency
  average <- dim$dimension.stability$average.item.stability
  frequency <- as.data.frame(table(native$boot.ndim[valid]) / sum(valid))
  names(frequency) <- c("dimensions", "frequency")
  list(native = native,
    items = data.frame(item = names(item), stability = as.numeric(item)),
    dimensions = data.frame(dimension = names(consistency), structural_consistency = as.numeric(consistency),
      average_item_stability = as.numeric(average)),
    frequency = frequency, requested = requested,
    valid = sum(valid), failed = sum(!valid))
}

.pm_network_bootstrap <- function(data, ega, settings) {
  .pm_network_progress(settings,sprintf("bootEGA (%s): %d replicates; %d workers",
    .pm_null_default(settings$bootstrap_type,"resampling"),settings$n_boot,settings$n_cores))
  structural <- .pm_network_boot_summary(.pm_network_boot_ega(data, settings), ega$wc)
  # The function is serialized with the estimator settings; ranks/correlations
  # must be recomputed on every bootstrap sample, not bootstrapped from a fixed S.
  fun <- .pm_network_estimator_function(settings)
  estimated <- bootnet::estimateNetwork(as.matrix(data), fun = fun, verbose = FALSE,
    labels = colnames(data), weighted = TRUE, signed = TRUE)
  .pm_network_progress(settings,sprintf("bootnet: %d edge resamples",settings$n_boot))
  accuracy <- bootnet::bootnet(estimated, nBoots = settings$n_boot,
    nCores = settings$n_cores, type = "nonparametric", statistics = "edge", verbose = FALSE)
  .pm_network_progress(settings,sprintf("bootnet: %d case-drop resamples",settings$n_boot))
  case <- bootnet::bootnet(estimated, nBoots = settings$n_boot,
    nCores = settings$n_cores, type = "case", statistics = c("strength", "betweenness", "closeness"), verbose = FALSE)
  cs_capture <- .pm_capture_conditions(bootnet::corStability(case, statistics = c("strength", "betweenness", "closeness"), verbose = FALSE))
  cs <- cs_capture$value
  cs_table <- if (is.null(cs)) data.frame(metric = c("strength", "betweenness", "closeness"), CS = NA_real_) else
    data.frame(metric = names(cs), CS = as.numeric(cs))
  cs_table$interpretation <- ifelse(is.na(cs_table$CS), "not_evaluable",
    ifelse(cs_table$CS < .25, "unstable", ifelse(cs_table$CS < .50, "caution", "preferable")))
  if (!all(.pm_network_centrality(estimated$graph)$connected)) {
    cs_table$CS[cs_table$metric == "closeness"] <- NA_real_
    cs_table$interpretation[cs_table$metric == "closeness"] <- "not_evaluable"
  }
  edges <- accuracy$bootTable[accuracy$bootTable$type == "edge", , drop = FALSE]
  edge_table <- if (!nrow(edges)) data.frame() else {
    ids <- unique(edges$id)
    do.call(rbind, lapply(ids, function(id) {
      values <- edges$value[edges$id == id]
      data.frame(edge = id, mean = mean(values, na.rm = TRUE),
        lower = unname(stats::quantile(values, .025, na.rm = TRUE)),
        upper = unname(stats::quantile(values, .975, na.rm = TRUE)), valid = sum(is.finite(values)))
    }))
  }
  list(structural = structural, estimated = estimated, accuracy = accuracy, case = case,
    cs = cs_table, cs_error = if (is.null(cs_capture$error)) NULL else conditionMessage(cs_capture$error),
    edge_intervals = edge_table, counts = data.frame(stage = c("bootEGA", "edge", "case"),
      requested = settings$n_boot,
      valid = c(structural$valid, length(accuracy$boots), length(case$boots))),
    conditional_on_uva = identical(settings$uva,"removed"))
}
